import copy
import random
import unittest
from fractions import Fraction as Q
from dynacomp import LinearModel, compile_linear, replay_linear, closure_witness, distinguishing_horizon
from dynacomp.rational import matrix, mv, addv, mul, row_basis, nullspace


class LinearTests(unittest.TestCase):
    def test_hidden_driver_repair(self):
        m = LinearModel.from_rows([[1, 1, 0], [0, 1, 0], [0, 0, 2]], [[0], [1], [0]], [[1, 0, 0]])
        r = compile_linear(m)
        self.assertEqual(len(r.E), 2)
        self.assertEqual(len(r.rounds), 1)
        self.assertTrue(replay_linear(r.packet()))
        w = r.rounds[0]["witness"]
        self.assertEqual(w["Ev"], ["0"])
        self.assertTrue(any(Q(v) for v in w["EAv"]))

    def test_forced_random_trajectories_and_minimality(self):
        rng = random.Random(8128)
        for n in range(1, 6):
            for _ in range(12):
                m = LinearModel.from_rows([[Q(rng.randrange(-2, 3), 3) for _ in range(n)] for _ in range(n)],
                    [[rng.randrange(-2, 3)] for _ in range(n)], [[1]+[0]*(n-1)])
                r = compile_linear(m)
                self.assertTrue(replay_linear(r.packet()))
                rows, ca = [], m.C
                for _ in range(n):
                    rows.extend(ca)
                    ca = mul(ca, m.A)
                self.assertEqual(len(r.E), len(row_basis(rows)))
                x = tuple(Q(rng.randrange(-4, 5)) for _ in range(n))
                controls = [[rng.randrange(-2, 3)] for _ in range(8)]
                outputs = [mv(m.C, x)]
                x0 = x
                for u in controls:
                    x = addv(mv(m.A, x), mv(m.B, u))
                    outputs.append(mv(m.C, x))
                self.assertEqual(r.simulate(x0, controls), tuple(outputs))

    def test_unobservable_modes_are_removed_even_if_unstable(self):
        m = LinearModel.from_rows([[1, 0], [0, 9]], [[1], [3]], [[1, 0]])
        r = compile_linear(m)
        self.assertEqual(len(r.E), 1)
        self.assertIsNone(distinguishing_horizon(m, [0, 1]))

    def test_earliest_horizon(self):
        m = LinearModel.from_rows([[0, 1, 0], [0, 0, 1], [0, 0, 0]], [[0], [0], [0]], [[1, 0, 0]])
        self.assertEqual(distinguishing_horizon(m, [0, 0, 1])["horizon"], 2)

    def test_identity_periodic_and_oblique_modes(self):
        for A in ([[1, 1], [0, 0]], [[0, -1], [1, 0]], [[1, 1], [0, 1]]):
            m = LinearModel.from_rows(A, [[0], [0]], [[1, 0]])
            r = compile_linear(m)
            self.assertTrue(replay_linear(r.packet()))
            self.assertEqual(r.simulate([2, 3], [[0]]*10)[-1], self.full(m, [2, 3], 10))

    @staticmethod
    def full(m, x, steps):
        for _ in range(steps):
            x = mv(m.A, x)
        return mv(m.C, x)

    def test_certificate_mutation_rejected(self):
        m = LinearModel.from_rows([[1, 1], [0, 1]], [[0], [1]], [[1, 0]])
        for key in ("E", "G", "H", "D", "C", "A", "B"):
            p = copy.deepcopy(compile_linear(m).packet())
            p[key][0][0] = str(Q(p[key][0][0])+1)
            with self.assertRaises(ValueError):
                replay_linear(p)

    def test_bad_dimensions_and_float_inputs(self):
        with self.assertRaises(TypeError):
            LinearModel.from_rows([[0.5]], [[0]], [[1]])
        with self.assertRaises(ValueError):
            LinearModel.from_rows([[1, 2]], [[0]], [[1]])
        with self.assertRaises(ValueError):
            compile_linear(LinearModel.from_rows([[1]], [[0]], [[0]]))

    def test_no_closure_witness_for_closed_encoder(self):
        self.assertIsNone(closure_witness([[1, 0], [0, 2]], [[1, 0]]))

    def test_kernel_matches_all_future_indistinguishability(self):
        m = LinearModel.from_rows([[1, 0, 1], [0, 2, 0], [0, 0, 1]], [[0], [0], [0]], [[1, 0, 0]])
        r = compile_linear(m)
        for v in nullspace(r.E):
            self.assertIsNone(distinguishing_horizon(m, v))
