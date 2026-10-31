import math
import unittest

from matrix_trace_condition_utils import (
    trace,
    determinant,
    first_invariant,
    second_invariant,
    third_invariant,
    principal_invariants,
)


class TraceTests(unittest.TestCase):
    def test_identity(self):
        self.assertEqual(trace([[1, 0, 0], [0, 1, 0], [0, 0, 1]]), 3)

    def test_general(self):
        self.assertEqual(trace([[1, 2, 3], [4, 5, 6], [7, 8, 9]]), 15)

    def test_negative_entries(self):
        self.assertEqual(trace([[-1, 0, 0], [0, -2, 0], [0, 0, -3]]), -6)


class DeterminantTests(unittest.TestCase):
    def test_identity(self):
        self.assertEqual(determinant([[1, 0, 0], [0, 1, 0], [0, 0, 1]]), 1)

    def test_singular(self):
        # Rows are linearly dependent, so the determinant must be zero.
        self.assertEqual(determinant([[1, 2, 3], [4, 5, 6], [7, 8, 9]]), 0)

    def test_known_value(self):
        # det([[2, 0, 0], [0, 3, 0], [0, 0, 4]]) == 24
        self.assertEqual(determinant([[2, 0, 0], [0, 3, 0], [0, 0, 4]]), 24)

    def test_negative_determinant(self):
        # Swapping two rows flips the sign of the determinant.
        self.assertEqual(determinant([[0, 1, 0], [1, 0, 0], [0, 0, 1]]), -1)


class InvariantTests(unittest.TestCase):
    def test_first_invariant_is_trace(self):
        m = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
        self.assertEqual(first_invariant(m), trace(m))

    def test_third_invariant_is_determinant(self):
        m = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
        self.assertEqual(third_invariant(m), determinant(m))

    def test_second_invariant_identity(self):
        # For the identity, each principal 2x2 minor is 1, so I2 == 3.
        self.assertEqual(second_invariant([[1, 0, 0], [0, 1, 0], [0, 0, 1]]), 3)

    def test_second_invariant_diagonal(self):
        # diag(2, 3, 4): minors are 6, 8, 12 -> I2 = 26.
        self.assertEqual(second_invariant([[2, 0, 0], [0, 3, 0], [0, 0, 4]]), 26)

    def test_second_invariant_general(self):
        m = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
        # (1*5 - 2*4) + (1*9 - 3*7) + (5*9 - 6*8)
        # = -3 + -12 + -3 = -18
        self.assertEqual(second_invariant(m), -18)

    def test_principal_invariants_tuple(self):
        m = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
        self.assertEqual(principal_invariants(m), (15, -18, 0))

    def test_principal_invariants_match_individuals(self):
        m = [[2, -1, 0], [3, 4, 1], [0, 5, -2]]
        i1, i2, i3 = principal_invariants(m)
        self.assertEqual(i1, first_invariant(m))
        self.assertEqual(i2, second_invariant(m))
        self.assertEqual(i3, third_invariant(m))

    def test_characteristic_polynomial_relation(self):
        # The invariants are exactly the coefficients of the characteristic
        # polynomial, so evaluating it at a known eigenvalue must give zero.
        m = [[2, 0, 0], [0, 3, 0], [0, 0, 4]]
        i1, i2, i3 = principal_invariants(m)
        lam = 3
        value = lam**3 - i1 * lam**2 + i2 * lam - i3
        self.assertEqual(value, 0)

    def test_float_tolerance(self):
        m = [[0.1, 0.2, 0.3], [0.4, 0.5, 0.6], [0.7, 0.8, 0.9]]
        i1, i2, i3 = principal_invariants(m)
        self.assertAlmostEqual(i1, 1.5)
        self.assertAlmostEqual(i2, -0.18)
        self.assertAlmostEqual(i3, 0.0)


if __name__ == "__main__":
    unittest.main()
