import unittest
from dynacomp import LinearModel
from dynacomp.proposals import reduced_candidate
from dynacomp.rational import matrix, mul, identity


class ProposalTests(unittest.TestCase):
    def test_oblique_left_inverse(self):
        m = LinearModel.from_rows([[1, 1], [0, 1]], [[0], [1]], [[1, 0]])
        p = reduced_candidate(m, [[1], [2]])
        self.assertEqual(mul(matrix(p["E"]), matrix(p["V"])), identity(1))
        self.assertEqual(p["status"], "candidate_requires_error_bound")

    def test_dependent_lift_rejected(self):
        m = LinearModel.from_rows([[1, 0], [0, 1]], [[0], [0]], [[1, 0]])
        with self.assertRaises(ValueError):
            reduced_candidate(m, [[1, 1], [2, 2]])
