import unittest
from fractions import Fraction as Q
from dynacomp import LinearModel, compile_linear
from dynacomp.fastforward import repeat_input_block, repeat_affine
from dynacomp.structure import aggregate_partition
from dynacomp.sparse import SparseOperator, compile_sparse
from dynacomp.rational import matrix, mul, coordinates


class StructureTests(unittest.TestCase):
    def test_regular_ring_one_aggregate(self):
        op = SparseOperator.from_entries(12, [(i, (i+1)%12, 1) for i in range(12)])
        self.assertEqual(len(aggregate_partition(op, [0]*12)["E"]), 1)

    def test_defect_discovers_reflection_classes(self):
        n = 16
        entries = [(i, i, "1/2") for i in range(n)]
        entries += [(i, (i+1)%n, "1/4") for i in range(n)]
        entries += [(i, (i-1)%n, "1/4") for i in range(n)]
        entries += [(0, 0, "1/100")]
        op = SparseOperator.from_entries(n, entries)
        discovered = aggregate_partition(op, [0]*n)
        self.assertEqual(len(discovered["E"]), 9)
        e = matrix(discovered["E"])
        g = matrix([coordinates(row, e) for row in op.pullback(e)])
        self.assertEqual(op.pullback(e), mul(g, e))
        for i in range(1, n):
            self.assertEqual(discovered["classes"][i], discovered["classes"][n-i])

    def test_initial_target_groups_are_preserved(self):
        op = SparseOperator.from_entries(4, [(i, i, 1) for i in range(4)])
        self.assertEqual(aggregate_partition(op, [0, 0, 1, 1])["classes"], [0, 0, 1, 1])


class FastForwardTests(unittest.TestCase):
    def test_trillion_steps_exact(self):
        m = LinearModel.from_rows([[1, 1], [0, 1]], [[0], [1]], [[1, 0]])
        r = compile_linear(m)
        t = 10**12
        self.assertEqual(repeat_input_block(r, [2, 3], [[1]], t), (Q(2+3*t+t*(t-1)//2),))

    def test_control_word_matches_replay(self):
        m = LinearModel.from_rows([["1/2", 1], [0, "1/3"]], [[1], [2]], [[1, 0]])
        r = compile_linear(m)
        inputs = [[2], [-1], [3]]
        for k in (0, 1, 2, 9):
            self.assertEqual(repeat_input_block(r, [1, 2], inputs, k), r.simulate([1, 2], inputs*k)[-1])

    def test_negative_repetition_rejected(self):
        with self.assertRaises(ValueError):
            repeat_affine([[1]], [0], -1)
