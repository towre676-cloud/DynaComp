import copy
import itertools
import random
import unittest
from fractions import Fraction as Q
from dynacomp.switched import (compile_switched, replay_switched, WeightedAutomaton, minimize_automaton,
                               replay_automaton, check_equivalence, replay_equivalence)
from dynacomp.rational import matrix, mv, mul


def words(alphabet, length):
    for k in range(length+1):
        yield from itertools.product(alphabet, repeat=k)


def random_automaton(rng, n, letters=("a", "b")):
    modes = {s: [[Q(rng.choice([0, 0, 1, -1, 2]), rng.choice([1, 2])) for _ in range(n)] for _ in range(n)]
             for s in letters}
    return WeightedAutomaton([rng.choice([0, 1]) for _ in range(n)], modes, [rng.choice([0, 1, -1]) for _ in range(n)])


class SwitchedTests(unittest.TestCase):
    def test_closure_under_every_mode_sequence(self):
        modes = {"heat": [[1, 1, 0], [0, 1, 0], [0, 0, 3]], "cool": [["1/2", 0, 0], [0, 1, 0], [0, 0, 5]]}
        p = compile_switched(modes, [[1, 0, 0]], {"heat": [[0], [1], [1]], "cool": [[0], [0], [1]]})
        self.assertEqual(len(p["E"]), 2)
        self.assertTrue(replay_switched(p))
        e, d = matrix(p["E"]), matrix(p["D"])
        rng = random.Random(5)
        for _ in range(20):
            x = tuple(Q(rng.randrange(-3, 4)) for _ in range(3))
            z = mv(e, x)
            for _ in range(6):
                s = rng.choice(["heat", "cool"])
                u = (Q(rng.randrange(-2, 3)),)
                x = tuple(a+b for a, b in zip(mv(matrix(modes[s]), x), mv(matrix(p["inputs"][s]), u)))
                z = tuple(a+b for a, b in zip(mv(matrix(p["G"][s]), z), mv(matrix(p["H"][s]), u)))
                self.assertEqual(mv(d, z), x[:1])

    def test_single_mode_misses_what_switching_needs(self):
        modes = {"a": [[1, 0], [0, 1]], "b": [[0, 1], [1, 0]]}
        self.assertEqual(len(compile_switched({"a": modes["a"]}, [[1, 0]])["E"]), 1)
        self.assertEqual(len(compile_switched(modes, [[1, 0]])["E"]), 2)

    def test_minimization_matches_brute_force_weights_and_hankel(self):
        rng = random.Random(2718)
        for n in range(1, 5):
            for _ in range(10):
                aut = random_automaton(rng, n)
                p = minimize_automaton(aut)
                if p["status"] == "zero_series":
                    self.assertTrue(all(aut.weight(w) == 0 for w in words("ab", n)))
                    continue
                self.assertTrue(replay_automaton(p))
                small = WeightedAutomaton.from_spec(p["minimal"])
                self.assertLessEqual(small.n, n)
                for w in words("ab", n+1):
                    self.assertEqual(aut.weight(w), small.weight(w))

    def test_redundant_hmm_collapses(self):
        # Two copies of the same biased coin, split evenly: a 2-state HMM with a 1-state series.
        half, coin = Q(1, 2), {"H": [["3/8", "3/8"], ["3/8", "3/8"]], "T": [["1/8", "1/8"], ["1/8", "1/8"]]}
        aut = WeightedAutomaton([half, half], coin, [1, 1])
        p = minimize_automaton(aut)
        self.assertEqual(p["order"], 1)
        self.assertEqual(aut.weight("HHT"), Q(3, 4)**2*Q(1, 4))

    def test_equivalence_certificate_and_shortest_counterexample(self):
        rng = random.Random(31)
        for _ in range(40):
            first, second = random_automaton(rng, rng.randrange(1, 4)), random_automaton(rng, rng.randrange(1, 4))
            p = check_equivalence(first, second)
            self.assertTrue(replay_equivalence(p))
            differing = [w for w in words("ab", first.n+second.n) if first.weight(w) != second.weight(w)]
            if p["equivalent"]:
                self.assertEqual(differing, [])
            else:
                self.assertEqual(len(p["word"]), len(differing[0]))
        rng = random.Random(9)
        aut = random_automaton(rng, 3)
        while minimize_automaton(aut)["status"] != "minimized":
            aut = random_automaton(rng, 3)
        p = check_equivalence(aut, WeightedAutomaton.from_spec(minimize_automaton(aut)["minimal"]))
        self.assertTrue(p["equivalent"])
        self.assertTrue(replay_equivalence(p))

    def test_forged_certificates_rejected(self):
        aut = WeightedAutomaton([1, 0], {"a": [[0, 1], [1, 0]]}, [1, 0])
        p = minimize_automaton(aut)
        bad = copy.deepcopy(p)
        bad["minimal"]["eta"][0] = str(Q(bad["minimal"]["eta"][0]) + 1)
        with self.assertRaises(ValueError):
            replay_automaton(bad)
        other = WeightedAutomaton([1, 0], {"a": [[0, 1], [1, 0]]}, [1, 1])
        q = check_equivalence(aut, other)
        self.assertFalse(q["equivalent"])
        q["word"] = []  # both automata weigh the empty word 1
        with self.assertRaises(ValueError):
            replay_equivalence(q)
        r = check_equivalence(aut, aut)
        r["G"]["a"][0][0] = str(Q(r["G"]["a"][0][0]) + 1)
        with self.assertRaises(ValueError):
            replay_equivalence(r)
