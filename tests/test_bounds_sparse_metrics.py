import random
import unittest
from fractions import Fraction as Q
from dynacomp import LinearModel, finite_horizon_bound, usable_horizon
from dynacomp.rational import mv, addv, matrix
from dynacomp.metrics import metric_margin, discover_metric
from dynacomp.sparse import SparseOperator, compile_sparse, replay_sparse


class BoundsTests(unittest.TestCase):
    def test_random_forcing_initial_error_and_transition_defects(self):
        rng = random.Random(29)
        for _ in range(30):
            A = [[Q(rng.randrange(-3, 4), 4) for _ in range(3)] for _ in range(3)]
            m = LinearModel.from_rows(A, [[1], [2], [-1]], [[1, 0, 1]])
            V, G, H = [[1], [0], [0]], [[Q(1, 2)]], [[1]]
            x, z = (Q(1), Q(-2), Q(3)), (Q(1),)
            inputs = [[rng.randrange(-2, 3)] for _ in range(6)]
            packet = finite_horizon_bound(m, V, G, H, {"full": x, "reduced": z}, inputs)
            for t in range(7):
                difference = abs(mv(m.C, x)[0]-mv(matrix([[1]]), z)[0])
                self.assertLessEqual(difference, Q(packet["bound"][t]))
                if t < 6:
                    x = addv(mv(m.A, x), mv(m.B, inputs[t]))
                    z = addv(mv(matrix(G), z), mv(matrix(H), inputs[t]))

    def test_nonnormal_transient_is_not_hidden(self):
        m = LinearModel.from_rows([["1/2", 8], [0, "1/2"]], [[0], [0]], [[1, 0]])
        p = finite_horizon_bound(m, [[1], [0]], [["1/2"]], [[0]],
                                 {"full": [0, 1], "reduced": [0]}, [[0]]*8)
        self.assertEqual(Q(p["bound"][0]), 0)
        self.assertEqual(Q(p["bound"][1]), 8)
        self.assertEqual(usable_horizon(p["bound"], 1), 0)

    def test_exact_reduction_zero_bound(self):
        m = LinearModel.from_rows([["1/2", 0], [0, 2]], [[1], [0]], [[1, 0]])
        p = finite_horizon_bound(m, [[1], [0]], [["1/2"]], [[1]],
                                 {"full": [2, 5], "reduced": [2]}, [[1]]*5)
        self.assertEqual(p["bound"], ["0"]*6)

    def test_prefix_horizon_and_negative_tolerance(self):
        self.assertEqual(usable_horizon([0, 3, 0], 1), 0)
        self.assertEqual(usable_horizon([2], 1), -1)
        with self.assertRaises(ValueError):
            usable_horizon([0], -1)


class SparseMetricTests(unittest.TestCase):
    def test_ring_total_closure_and_trajectory(self):
        n = 17
        entries = [(i, i, "1/2") for i in range(n)] + [(i, (i+1)%n, "1/4") for i in range(n)] + [(i, (i-1)%n, "1/4") for i in range(n)]
        op = SparseOperator.from_entries(n, entries)
        p = compile_sparse(op, [[1]]*n, [[1]*n])
        self.assertEqual(p["G"], [["1"]])
        self.assertTrue(replay_sparse(p))
        x = list(range(n))
        for _ in range(10):
            self.assertEqual(sum(op.step(x)), sum(x))
            x = op.step(x)

    def test_defect_destroys_aggregate_closure(self):
        op = SparseOperator.from_entries(5, [(i, i, 1) for i in range(5)]+[(0, 0, 1)])
        p = compile_sparse(op, [[0]]*5, [[1]*5])
        self.assertEqual(len(p["E"]), 2)
        self.assertTrue(replay_sparse(p))

    def test_sparse_duplicate_entries_and_rank_limit(self):
        op = SparseOperator.from_entries(2, [(0, 0, 1), (0, 0, -1), (0, 1, 1)])
        self.assertEqual(len(op.entries), 1)
        with self.assertRaises(OverflowError):
            compile_sparse(op, [[0], [0]], [[1, 0]], max_rank=1)

    def test_metric_can_certify_nonnormal_stability(self):
        A = [["1/2", 8], [0, "1/2"]]
        self.assertFalse(metric_margin(A, [1, 1])["contractive"])
        self.assertTrue(metric_margin(A, [32, 1])["contractive"])
        self.assertTrue(discover_metric(A)["contractive"])

    def test_metric_search_failure_not_instability_proof(self):
        self.assertFalse(discover_metric([[0, -1], [1, 0]])["contractive"])
        with self.assertRaises(ValueError):
            metric_margin([[1]], [0])
