# matrix_trace_condition_utils

Computes the trace, determinant, and principal invariants of a 3x3 matrix for structural analysis.

```python
from matrix_trace_condition_utils import (
    trace,
    determinant,
    first_invariant,
    second_invariant,
    third_invariant,
    principal_invariants,
)

m = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]

trace(m)                 # 15
determinant(m)           # 0
first_invariant(m)       # 15
second_invariant(m)      # -18
third_invariant(m)       # 0
principal_invariants(m)  # (15, -18, 0)
```

## Why this exists

Structural checks often need all three principal invariants of a 3x3 stress or strain tensor at once, and the standard formulas are easy to get wrong when written out by hand. The library keeps one clear definition of each quantity and computes them with exact integer arithmetic where the inputs allow it. The trade-off is that it does not validate shapes or types: it assumes a well-formed 3x3 sequence of three rows of three numbers, because checking on every call would slow down the numerical loops where these functions are used.

## Edge cases

The invariants are defined by the characteristic polynomial

    det(lambda * I - A) = lambda**3 - I1 * lambda**2 + I2 * lambda - I3

so `first_invariant` is exactly `trace`, `third_invariant` is exactly `determinant`, and `second_invariant` is the sum of the three principal 2x2 minors. Inputs are not copied or coerced, so passing a matrix of `float` values returns `float` results and small rounding errors are possible; compare floating results with a tolerance rather than equality.
