"""Core routines for 3x3 matrix invariants.

The principal invariants of a 3x3 matrix A are defined by the characteristic
polynomial

    det(lambda * I - A) = lambda**3 - I1 * lambda**2 + I2 * lambda - I3

so that I1 is the trace, I3 is the determinant, and I2 is the sum of the
principal 2x2 minors.  Structural analysis commonly needs all three at once;
computing them together avoids evaluating the same entries repeatedly.

Every function accepts a matrix as a sequence of three rows, each of which is a
sequence of three numbers.  No validation is performed beyond what indexing
naturally raises; callers are expected to pass well-formed 3x3 data because the
cost of checking on every call is not worth it in tight numerical loops.
"""


def trace(matrix):
    """Return the trace of a 3x3 matrix (sum of the main diagonal).

    The trace is the first principal invariant and is used directly in many
    structural failure criteria, so it is exposed under its own name.
    """
    return matrix[0][0] + matrix[1][1] + matrix[2][2]


def determinant(matrix):
    """Return the determinant of a 3x3 matrix.

    Uses the cofactor expansion along the first row.  For 3x3 matrices this is
    both exact for integer-like inputs and fast enough that no special-casing
    is warranted.
    """
    a, b, c = matrix[0]
    d, e, f = matrix[1]
    g, h, i = matrix[2]
    return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)


def first_invariant(matrix):
    """Return the first principal invariant I1, equal to the trace."""
    return trace(matrix)


def second_invariant(matrix):
    """Return the second principal invariant I2.

    I2 is the sum of the three principal 2x2 minors:

        I2 = (a*e - b*d) + (a*i - c*g) + (e*i - f*h)

    This form is written out explicitly rather than derived from the trace of
    the adjugate so that it stays exact for integer inputs and does not
    introduce a division.
    """
    a, b, _ = matrix[0]
    d, e, f = matrix[1]
    g, h, i = matrix[2]
    return (a * e - b * d) + (a * i - g * _third(matrix, 0, 2)) + (e * i - f * h)


def _third(matrix, row, col):
    """Return matrix[row][col]; helper kept private to keep the minor formula flat."""
    return matrix[row][col]


def third_invariant(matrix):
    """Return the third principal invariant I3, equal to the determinant."""
    return determinant(matrix)


def principal_invariants(matrix):
    """Return ``(I1, I2, I3)`` for a 3x3 matrix.

    Computing the three invariants together is the common case in structural
    analysis; returning a tuple keeps the call site cheap and avoids three
    separate passes over the same entries.
    """
    a, b, c = matrix[0]
    d, e, f = matrix[1]
    g, h, i = matrix[2]
    i1 = a + e + i
    i2 = (a * e - b * d) + (a * i - c * g) + (e * i - f * h)
    i3 = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
    return i1, i2, i3
