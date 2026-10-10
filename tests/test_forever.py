import copy
import unittest
from fractions import Fraction as Q
from dynacomp.uniform import SparseLinearModel, uniform_enrichment, sharp_box_bound, _project
from dynacomp.forever import all_time_bound, replay_all_time, moment_rows, all_time_enrichment


def rod(n=8, k="1/4", loss="1/20"):
    k, loss = Q(k), Q(loss)
    entries = []
    for i in range(n):
        near = [j for j in (i-1, i+1) if 0 <= j < n]
        entries.append([i, i, str(1 - k*len(near) - loss)])
        entries += [[i, j, str(k)] for j in near]
    return SparseLinearModel.from_spec({"n": n, "entries": entries, "B": [["1"]] + [["0"]]*(n-1),
                                        "C": [["1/2", "1/2"] + ["0"]*(n-2)]})


class AllTimeTests(unittest.TestCase):
    def setUp(self):
        self.m = rod()
        self.n = self.m.operator.n
        self.box = ([1]*self.n, ["1/2"]*self.n, [1], [1])

    def candidate(self, rows):
        v, e, g, h, d = _project(self.m, rows)
        return v, e, g, h

    def exact_sharp(self, cand, t):
        c, r, a, b = self.box
        return Q(sharp_box_bound(self.m, *cand, c, r, [b]*t, [a]*t)["bound"][t])

    def test_bracket_contains_exact_finite_suprema(self):
        cand = self.candidate(list(self.m.C) + [self.m.operator.pullback(self.m.C)[0]])
        p = all_time_bound(self.m, *cand, *self.box)
        self.assertEqual(p["status"], "converged")
        lower, upper = Q(p["lower"]), Q(p["upper"])
        self.assertLess(upper - lower, Q(1, 10**5))
        for t in range(0, 40, 3):
            self.assertLessEqual(self.exact_sharp(cand, t), upper)
        self.assertLessEqual(lower, self.exact_sharp(cand, p["worst_time_lower_bound"]))
        self.assertTrue(replay_all_time(p))

    def test_refutation_is_an_attained_violation(self):
        cand = self.candidate(list(self.m.C))
        p = all_time_bound(self.m, *cand, *self.box, tolerance="1/100")
        self.assertEqual(p["status"], "refuted")
        t = p["refutation"]["time"]
        self.assertGreater(self.exact_sharp(cand, t), Q(1, 100))
        self.assertTrue(replay_all_time(p))

    def test_moments_match_steady_state(self):
        rows = list(self.m.C) + moment_rows(self.m, 1)
        p = all_time_bound(self.m, *self.candidate(rows), *self.box)
        self.assertEqual(p["steady_state_error"], ["0"])
        q = all_time_bound(self.m, *self.candidate(list(self.m.C)), *self.box)
        self.assertNotEqual(q["steady_state_error"], ["0"])

    def test_enrichment_certifies_and_finite_horizon_does_not_suffice(self):
        tol = "1/20"
        p = all_time_enrichment(self.m, *self.box, tol)
        self.assertEqual(p["status"], "certified")
        self.assertLessEqual(Q(p["upper"]), Q(tol))
        self.assertTrue(replay_all_time(p))
        short = uniform_enrichment(self.m, self.box[0], self.box[1], [self.box[3]]*3, tol,
                                   input_centers=[self.box[2]]*3, strategy="prefix")
        q = all_time_bound(self.m, short["V"], short["E"], short["G"], short["H"], *self.box, tolerance=tol)
        if short["rank"] < p["reduced_order"]:
            self.assertEqual(q["status"], "refuted")

    def test_tampering_rejected(self):
        p = all_time_bound(self.m, *self.candidate(list(self.m.C)), *self.box, tolerance="1/100")
        bad = copy.deepcopy(p)
        bad["lower"] = str(Q(bad["lower"]) * 2)
        with self.assertRaises(ValueError):
            replay_all_time(bad)
        bad = copy.deepcopy(p)
        bad["gamma"] = "1/1000"
        with self.assertRaises(ValueError):
            replay_all_time(bad)
