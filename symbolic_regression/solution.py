"""Starter solution for the frontier-cs symbolic_regression benchmark.

Implement `Solution.solve(X, y)` to return a symbolic expression that fits the
target. `X` is an (n_samples, n_features) array; `y` is (n_samples,). Return:

    {"expression": "<expr over x1..xn using + - * / ** and sin/cos/exp/log>",
     "details": {"complexity": <optional int>}}

Your expression is scored by its MSE on the data plus a simplicity penalty,
relative to a fixed reference (higher score is better; the mean across all
target families is your leaderboard score). The starter below is a plain linear
least-squares fit -- beat it by discovering the true nonlinear form.
"""
from __future__ import annotations

import numpy as np


class Solution:
    def solve(self, X, y):
        X = np.asarray(X, dtype=float)
        y = np.asarray(y, dtype=float).ravel()
        n_features = X.shape[1]

        design = np.hstack([X, np.ones((X.shape[0], 1))])
        coef, *_ = np.linalg.lstsq(design, y, rcond=None)

        terms = [f"({coef[i]:.8g})*x{i + 1}" for i in range(n_features)]
        terms.append(f"({coef[n_features]:.8g})")
        expression = " + ".join(terms)

        return {"expression": expression, "details": {"complexity": 2 * n_features}}
