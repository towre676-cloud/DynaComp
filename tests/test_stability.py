import copy
import random
import unittest
from fractions import Fraction as Q
from dynacomp.stability import (schur_stable, hurwitz_stable, positive_definite, quadratic_certificate,
                                replay_quadratic, common_quadratic_certificate, replay_common_quadratic,
                                sqrt_up, charpoly)
from dynacomp.rational import matrix, mul, transpose, rref, identity


def lyapunov_solution(a, continuous=False):
    """Exact P solving A^T P A - P = -I (or A^T P + P A = -I) via the n^2 linear system; None if singular."""
    n = len(a)
    rows = []
    for i in range(n):
        for j in range(n):
            row = [Q(0)]*(n*n)
            for k in range(n):
                for l in range(n):
                    if continuous:
                        coef = (a[k][i] if l == j else 0) + (a[l][j] if k == i else 0)
                    else:
                        coef = a[k][i]*a[l][j] - (1 if (k == i and l == j) else 0)
                    row[k*n+l] += coef
            rows.append(row + [Q(-1 if i == j else 0)])
    reduced, pivots = rref(rows)
    if len(pivots) < n*n or n*n in pivots:
        return None
    sol = [Q(0)]*(n*n)
    for r, p in zip(reduced, pivots):
        sol[p] = r[-1]
    return tuple(tuple(sol[i*n+j] for j in range(n)) for i in range(n))


def random_matrix(rng, n):
    return matrix([[Q(rng.randrange(-6, 7), rng.choice([2, 3, 4, 5, 8])) for _ in range(n)] for _ in range(n)])


class StabilityTests(unittest.TestCase):
    def test_charpoly_cayley_hamilton(self):
        rng = random.Random(1)
        for n in range(1, 5):
            a = random_matrix(rng, n)
            p = charpoly(a)
            total, power = [[Q(0)]*n for _ in range(n)], identity(n)
            for c in p:
                total = [[x + c*y for x, y in zip(r, s)] for r, s in zip(total, power)]
                power = mul(power, a)
            self.assertTrue(all(x == 0 for row in total for x in row))

    def test_schur_cohn_matches_exact_lyapunov_theorem(self):
        rng = random.Random(7)
        seen = set()
        for n in range(1, 5):
            for _ in range(40):
                a = random_matrix(rng, n)
                p = lyapunov_solution(a)
                lyapunov = p is not None and positive_definite(p)
                self.assertEqual(schur_stable(a)["stable"], lyapunov)
                seen.add(lyapunov)
        self.assertEqual(seen, {True, False})

    def test_hurwitz_matches_continuous_lyapunov_theorem(self):
        rng = random.Random(11)
        seen = set()
        for n in range(1, 4):
            for _ in range(40):
                a = random_matrix(rng, n)
                p = lyapunov_solution(a, continuous=True)
                lyapunov = p is not None and positive_definite(p)
                self.assertEqual(hurwitz_stable(a)["stable"], lyapunov)
                seen.add(lyapunov)
        self.assertEqual(seen, {True, False})

    def test_boundary_cases(self):
        self.assertFalse(schur_stable([[1]])["stable"])
        self.assertFalse(schur_stable([[0, 1], [-1, 0]])["stable"])  # eigenvalues on the circle
        self.assertTrue(schur_stable([["-99/100"]])["stable"])
        self.assertFalse(hurwitz_stable([[0]])["stable"])
        self.assertFalse(hurwitz_stable([[1, 0], [0, -1]])["stable"])

    def test_quadratic_certificates_for_stable_matrices(self):
        rng = random.Random(3)
        certified, stable = 0, 0
        for n in range(1, 5):
            for _ in range(25):
                a = random_matrix(rng, n)
                c = quadratic_certificate(a)
                if schur_stable(a)["stable"]:
                    stable += 1
                    self.assertIn(c["status"], ("certified", "proposal_failed"))
                    if c["status"] == "certified":
                        certified += 1
                        self.assertTrue(replay_quadratic(c))
                else:
                    self.assertEqual(c["status"], "not_schur_stable")
        self.assertGreater(stable, 15)
        self.assertEqual(certified, stable)  # every stable matrix in this corpus is certified

    def test_nonnormal_metric_is_near_spectral(self):
        c = quadratic_certificate([["1/2", 8], [0, "1/2"]])
        self.assertEqual(c["status"], "certified")
        self.assertLess(Q(c["gamma"]), Q(3, 10))  # spectral radius squared is 1/4
        bad = copy.deepcopy(c)
        bad["gamma"] = "1/5"
        with self.assertRaises(ValueError):
            replay_quadratic(bad)

    def test_common_quadratic_switching(self):
        good = {"a": [["1/2", "1/4"], [0, "1/2"]], "b": [["1/2", 0], ["1/4", "1/2"]]}
        c = common_quadratic_certificate(good)
        self.assertEqual(c["status"], "certified")
        self.assertTrue(replay_common_quadratic(c))
        # Each mode is stable, but alternating them is not: no common P may exist.
        bad = {"a": [["1/2", 1], [0, "1/2"]], "b": [["1/2", 0], [1, "1/2"]]}
        self.assertNotEqual(common_quadratic_certificate(bad)["status"], "certified")
        self.assertFalse(schur_stable(mul(matrix(bad["a"]), matrix(bad["b"])))["stable"])

    def test_sqrt_up(self):
        for x in (Q(2), Q(1, 3), Q(10**30, 7), Q(0), Q(49, 4)):
            r = sqrt_up(x)
            self.assertGreaterEqual(r*r, x)
            self.assertLess(r*r - x, max(x, 1)*Q(1, 2**40))


class ContinuousCertificateTests(unittest.TestCase):
    def test_continuous_certificates(self):
        from dynacomp.stability import continuous_certificate, replay_continuous
        rng = random.Random(17)
        stable = 0
        for n in range(1, 4):
            for _ in range(20):
                a = random_matrix(rng, n)
                c = continuous_certificate(a)
                if hurwitz_stable(a)["stable"]:
                    stable += 1
                    self.assertEqual(c["status"], "certified")
                    self.assertTrue(replay_continuous(c))
                else:
                    self.assertEqual(c["status"], "not_hurwitz_stable")
        self.assertGreater(stable, 5)
