import copy
import random
import unittest
from fractions import Fraction as Q
from dynacomp import LinearModel, compile_linear
from dynacomp.realization import minimal_realization, replay_realization, markov_parameter, simulate_zero_state
from dynacomp.rational import matrix, mul, mv, addv, rref


def hankel_rank(m, size):
    rows = []
    for i in range(size):
        blocks = [markov_parameter(m.A, m.B, m.C, i+j) for j in range(size)]
        for r in range(len(m.C)):
            rows.append(tuple(x for b in blocks for x in b[r]))
    return len(rref(rows)[1])


class RealizationTests(unittest.TestCase):
    def test_unreachable_mode_removed_beyond_observability(self):
        # Observable from C but never excited by B: observability keeps 2, realization keeps 1.
        m = LinearModel.from_rows([[1, 0], [0, 2]], [[1], [0]], [[1, 1]])
        self.assertEqual(len(compile_linear(m).E), 2)
        p = minimal_realization(m)
        self.assertEqual(p["minimal_order"], 1)
        self.assertTrue(replay_realization(p))

    def test_random_markov_parameters_and_hankel_rank(self):
        rng = random.Random(1729)
        for n in range(1, 6):
            for _ in range(8):
                A = [[Q(rng.randrange(-2, 3), rng.choice([1, 2, 3])) for _ in range(n)] for _ in range(n)]
                # Sparse B and C make cancellations common.
                B = [[rng.choice([0, 0, 1, -1]) for _ in range(2)] for _ in range(n)]
                C = [[rng.choice([0, 0, 1, 2]) for _ in range(n)] for _ in range(2)]
                if not any(map(any, B)) or not any(map(any, C)):
                    continue
                m = LinearModel.from_rows(A, B, C)
                try:
                    p = minimal_realization(m)
                except ValueError:
                    self.assertEqual(hankel_rank(m, n), 0)
                    continue
                self.assertTrue(replay_realization(p))
                self.assertEqual(p["minimal_order"], hankel_rank(m, n))
                G, H, D = (matrix(p[k]) for k in ("G", "H", "D"))
                for k in range(2*n+1):
                    self.assertEqual(markov_parameter(m.A, m.B, m.C, k), markov_parameter(G, H, D, k))
                u = [[rng.randrange(-3, 4), rng.randrange(-3, 4)] for _ in range(6)]
                x, ys = (Q(0),)*n, [mv(m.C, (Q(0),)*n)]
                for step in u:
                    x = addv(mv(m.A, x), mv(m.B, step))
                    ys.append(mv(m.C, x))
                self.assertEqual(simulate_zero_state(G, H, D, u), tuple(ys))

    def test_mutation_and_false_minimality_rejected(self):
        m = LinearModel.from_rows([[0, 1, 0], [0, 0, 1], [1, 0, 0]], [[1], [0], [0]], [[1, 0, 0]])
        good = minimal_realization(m)
        self.assertTrue(replay_realization(good))
        for key in ("P", "F", "K", "E", "G", "H", "D"):
            bad = copy.deepcopy(good)
            bad[key][0][0] = str(Q(bad[key][0][0]) + 1)
            with self.assertRaises(ValueError):
                replay_realization(bad)
        bad = copy.deepcopy(good)
        bad["hankel_rows"][1] = bad["hankel_rows"][0]
        with self.assertRaises(ValueError):
            replay_realization(bad)

    def test_mode_excitation_orders(self):
        # A driven cascade needs both modes; a diagonal model whose actuator excites one observed mode needs one.
        m = LinearModel.from_rows([["1/2", 0], [1, "1/3"]], [[1], [0]], [[0, 1]])
        self.assertEqual(minimal_realization(m)["minimal_order"], 2)
        m = LinearModel.from_rows([["1/2", 0], [0, "1/3"]], [[0], [1]], [[1, 1]])
        self.assertEqual(minimal_realization(m)["minimal_order"], 1)
