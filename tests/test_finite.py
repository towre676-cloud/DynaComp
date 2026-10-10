import copy
import itertools
import unittest
from dynacomp.finite import FiniteModel, compile_finite, replay_finite, distinguish, outcome


class FiniteTests(unittest.TestCase):
    def test_all_324_two_state_partial_machines(self):
        count = 0
        words = [word for length in range(4) for word in itertools.product(range(2), repeat=length)]
        for observations in itertools.product(range(2), repeat=2):
            for values in itertools.product((None, 0, 1), repeat=4):
                model = FiniteModel(observations, (values[:2], values[2:]))
                packet = compile_finite(model)
                self.assertTrue(replay_finite(packet))
                same = all(outcome(model, 0, w) == outcome(model, 1, w) for w in words)
                self.assertEqual(packet["classes"][0] == packet["classes"][1], same)
                count += 1
        self.assertEqual(count, 324)

    def test_larger_machine_requires_pair_specific_experiments(self):
        model = FiniteModel((0, 0, 0, 1, 0, 0),
                            ((1, 2), (3, 1), (2, 3), (3, 3), (5, 2), (3, 5)))
        p = compile_finite(model)
        self.assertTrue(replay_finite(p))
        self.assertEqual(p["classes"][0], p["classes"][4])
        self.assertEqual(p["classes"][1], p["classes"][5])
        self.assertNotEqual(distinguish(model, 0, 2), distinguish(model, 0, 1))

    def test_missing_witness_and_forged_quotient_rejected(self):
        m = FiniteModel((0, 1), ((0,), (1,)))
        p = compile_finite(m)
        p["distinguishing_witnesses"] = []
        with self.assertRaises(ValueError):
            replay_finite(p)
        p = compile_finite(m)
        p["quotient_transitions"][0][0] = 1
        with self.assertRaises(ValueError):
            replay_finite(p)

    def test_definedness_distinguishes_equal_labels(self):
        m = FiniteModel((0, 0), ((None,), (1,)))
        self.assertEqual(distinguish(m, 0, 1), [0])
        self.assertTrue(replay_finite(compile_finite(m)))

    def test_invalid_successor(self):
        with self.assertRaises(ValueError):
            FiniteModel((0,), ((2,),))
