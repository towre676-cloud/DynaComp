import copy
import itertools
import random
import unittest
from fractions import Fraction as Q
from dynacomp.uniform import SparseLinearModel, sharp_box_bound, _project
from dynacomp.zonotope import sharp_zonotope_bound, box_generators, replay_zonotope
from dynacomp.rational import mv, addv, mul, matrix


def random_model(rng, n, m):
    entries = [[i, j, str(Q(rng.randrange(-3, 4), 6))] for i in range(n) for j in range(n) if rng.random() < 0.6]
    return SparseLinearModel.from_spec({"n": n, "entries": entries,
                                        "B": [[str(rng.randrange(-1, 2)) for _ in range(m)] for _ in range(n)],
                                        "C": [[str(rng.randrange(0, 2)) for _ in range(n)]]})


class ZonotopeTests(unittest.TestCase):
    def test_boxes_reproduce_sharp_box_bound(self):
        rng = random.Random(21)
        for _ in range(10):
            n, m, T = rng.randrange(2, 5), rng.randrange(1, 3), rng.randrange(1, 6)
            model = random_model(rng, n, m)
            if not any(model.C[0]):
                continue
            v, e, g, h, _ = _project(model, [model.C[0]])
            c = [Q(rng.randrange(-2, 3)) for _ in range(n)]
            r = [Q(rng.randrange(0, 3), 2) for _ in range(n)]
            a = [[Q(rng.randrange(-1, 2)) for _ in range(m)] for _ in range(T)]
            b = [[Q(rng.randrange(0, 2)) for _ in range(m)] for _ in range(T)]
            box = sharp_box_bound(model, v, e, g, h, c, r, b, a)
            zono = sharp_zonotope_bound(model, v, e, g, h, c, a, box_generators(n, r, b))
            self.assertEqual(box["bound"], zono["bound"])
            self.assertTrue(replay_zonotope(zono))

    def test_vertex_enumeration_matches(self):
        rng = random.Random(5)
        for _ in range(8):
            n, m, T = 3, 2, 4
            model = random_model(rng, n, m)
            if not any(model.C[0]):
                continue
            v, e, g, h, d = _project(model, [model.C[0]])
            gens = [{"x0": {k: Q(rng.randrange(-2, 3), 2) for k in range(n)}},
                    {"inputs": {t: [1, -1] for t in range(T)}},
                    {"x0": {0: 1}, "inputs": {1: [Q(1, 2), 0]}}]
            c, a = [1, 0, -1], [[0, 0]]*T
            z = sharp_zonotope_bound(model, v, e, g, h, c, a, gens)
            best = Q(0)
            for signs in itertools.product((-1, 1), repeat=len(gens)):
                x = [Q(x) for x in c]
                word = [list(map(Q, row)) for row in a]
                for s, gen in zip(signs, gens):
                    for k, w in gen.get("x0", {}).items():
                        x[k] += s*Q(w)
                    for t, vec in gen.get("inputs", {}).items():
                        word[t] = [p + s*Q(w) for p, w in zip(word[t], vec)]
                xs, zs = tuple(x), mv(e, tuple(x))
                for t in range(T+1):
                    best = max(best, abs(mv(model.C, xs)[0] - mv(d, zs)[0]))
                    if t < T:
                        xs = addv(model.operator.step(xs), mv(model.B, word[t]))
                        zs = addv(mv(g, zs), mv(h, word[t]))
            self.assertEqual(Q(z["worst_bound"]), best)

    def test_correlation_never_loosens_and_tamper_rejected(self):
        rng = random.Random(9)
        model = random_model(rng, 4, 1)
        model = SparseLinearModel(model.operator, model.B, ((1, 1, 0, 0),))
        v, e, g, h, _ = _project(model, [model.C[0]])
        box = sharp_zonotope_bound(model, v, e, g, h, [0]*4, [[0]]*5, box_generators(4, [1]*4, [[1]]*5))
        shared = sharp_zonotope_bound(model, v, e, g, h, [0]*4, [[0]]*5,
                                      [{"x0": {k: 1 for k in range(4)}}, {"inputs": {t: [1] for t in range(5)}}])
        self.assertLessEqual(Q(shared["worst_bound"]), Q(box["worst_bound"]))
        bad = copy.deepcopy(shared)
        bad["worst_bound"] = str(Q(bad["worst_bound"]) + 1)
        with self.assertRaises(ValueError):
            replay_zonotope(bad)
