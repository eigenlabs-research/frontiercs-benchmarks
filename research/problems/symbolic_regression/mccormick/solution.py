import numpy as np


class Solution:
    """Baseline: predict a constant (the mean of y). Replace with a real symbolic regressor.

    solve(X, y) returns a dict with:
      - "expression": a sympy-parseable string in variables x1..xn (here just the mean constant), and
      - "predictions": one value per input row. We supply predictions explicitly because the evaluator
        uses them directly; a bare constant expression alone would evaluate to a scalar and fail the
        per-row shape check. A constant baseline scores low, leaving the benchmark honestly unclaimed.
    """

    def solve(self, X, y):
        X = np.asarray(X, dtype=float)
        y = np.asarray(y, dtype=float)
        mean = float(np.mean(y)) if y.size else 0.0
        preds = np.full(X.shape[0], mean, dtype=float)
        return {"expression": repr(mean), "predictions": preds}
