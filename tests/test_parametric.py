import random
import unittest
from fractions import Fraction as Q
from dynacomp import upoly as P
from dynacomp.parametric import parametric_analysis, replay_parametric, classify, poly_matrix, rank_at, _krylov_stack


class PolynomialTests(unittest.TestCase):
    def test_roots_rational_and_irrational(self):
        p = P.mul(P.mul((-2, 0, 1), (-1, 3)), P.mul((1, 1), (1, 1)))
        rational, irrational = P.real_roots(p)
        self.assertEqual(rational, [Q(-1), Q(1, 3)])
        self.assertEqual(len(irrational), 2)
        for lo, hi in irrational:
            self.assertLess(lo, hi)
            self.assertTrue(lo*lo < 2 < hi*hi or hi*hi < 2 < lo*lo)

    def test_random_root_counts(self):
        rng = random.Random(4)
        for _ in range(30):
            roots = sorted({Q(rng.randrange(-20, 21), rng.randrange(1, 6)) for _ in range(rng.randrange(1, 6))})
            p = (1,)
            for r in roots:
                p = P.mul(p, (-r, 1))
            p = P.mul(p, (1, 0, 1))  # no real roots added
            self.assertEqual(P.real_roots(p), (roots, []))


class ParametricTests(unittest.TestCase):
    def test_accidental_closure_detected(self):
        p = parametric_analysis([[1, [0, 1]], [0, 2]], [[0], [1]], [[1, 0]])
        self.assertEqual(p["observable"]["generic_rank"], 2)
        self.assertEqual([x["value"] for x in p["observable"]["rational_exceptions"]], ["0"])
        self.assertTrue(classify(p, 0)["accidental"])
        self.assertFalse(classify(p, "1/7")["accidental"])
        self.assertEqual(p["minimal"]["rational_exceptions"][0]["rank"], 0)
        self.assertTrue(replay_parametric(p))

    def test_irrational_exceptions(self):
        # O = [[1, 1], [t^2-2, 1]] has determinant 3 - t^2: observability drops exactly at +-sqrt(3).
        p = parametric_analysis([[0, 1], [[-2, 0, 1], 0]], [[1], [0]], [[1, 1]])
        prof = p["observable"]
        self.assertEqual(prof["generic_rank"], 2)
        self.assertEqual(prof["rational_exceptions"], [])
        self.assertEqual(prof["real_exception_count"], 2)
        drop = tuple(Q(x) for x in prof["drop_polynomial"])
        for item in prof["irrational_exceptions"]:
            lo, hi = (Q(x) for x in item["interval"])
            self.assertTrue(lo*lo < 3 < hi*hi or hi*hi < 3 < lo*lo)
            self.assertLess(P.evaluate(drop, lo)*P.evaluate(drop, hi), 0)
        self.assertTrue(replay_parametric(p))

    def test_random_drop_points_match_pointwise_ranks(self):
        rng = random.Random(77)
        for _ in range(12):
            n = rng.randrange(2, 4)
            A = [[rng.choice([0, 1, -1, [0, 1], [1, -1], [0, 0, 1]]) for _ in range(n)] for _ in range(n)]
            B = [[rng.choice([0, 1])] for _ in range(n)]
            C = [[rng.choice([0, 1, [0, 1]]) for _ in range(n)]]
            p = parametric_analysis(A, B, C)
            self.assertTrue(replay_parametric(p))
            a, c = poly_matrix(A), poly_matrix(C)
            o = _krylov_stack(a, c, n)
            exceptions = {Q(x["value"]) for x in p["observable"]["rational_exceptions"]}
            for t in (Q(k, d) for k in range(-6, 7) for d in (1, 2, 3)):
                expected = p["observable"]["generic_rank"] - (t in exceptions)
                if t in exceptions:
                    self.assertLess(rank_at(o, t), p["observable"]["generic_rank"])
                else:
                    self.assertEqual(rank_at(o, t), expected)
