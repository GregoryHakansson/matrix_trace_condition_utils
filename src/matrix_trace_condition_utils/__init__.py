"""Utilities for the trace, determinant, and principal invariants of 3x3 matrices.

All functions operate on nested sequences of three rows of three numeric values.
The module is deliberately dependency-free and works with any objects that
support the arithmetic operations used (int, float, Fraction, Decimal, ...).
"""

from .core import (
    trace,
    determinant,
    first_invariant,
    second_invariant,
    third_invariant,
    principal_invariants,
)

__all__ = [
    "trace",
    "determinant",
    "first_invariant",
    "second_invariant",
    "third_invariant",
    "principal_invariants",
]
