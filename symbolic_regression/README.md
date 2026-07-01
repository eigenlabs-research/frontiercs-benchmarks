# symbolic_regression

A [frontier-cs](https://openfrontiercs.com) benchmark, adapted from
[FrontierCS/Frontier-CS](https://github.com/FrontierCS/Frontier-CS) `research/problems/symbolic_regression` (MIT).

## Problem

Discover a closed-form symbolic expression `f(x1..xn)` that fits each target
function with **low error** and **low complexity**. Five target families are
scored: `mccormick`, `mixed_polyexp_4d`, `peaks`, `ripple`, `sincos`.

## What you edit

Only `solution.py` — a `Solution` class whose `solve(X, y)` returns:

```python
{"expression": "sin(x1 + x2) + (x1 - x2)**2 - 1.5*x1 + 2.5*x2 + 1",
 "details": {"complexity": 8}}   # complexity optional; computed if omitted
```

Allowed operators: `+ - * / **` and functions `sin, cos, exp, log` over
`x1..xn`.

## How it's graded

Per target: `score = 100 · clamp((m_base − MSE)/(m_base − m_ref), 0, 1) · 0.99^max(complexity − C_ref, 0)`.
Your leaderboard score is the **mean** across the five targets; the per-target
breakdown is shown as metrics. Higher is better; an invalid/unparseable
expression scores 0.

## Run it

```bash
frontier-cs clone symbolic_regression
cd frontiercs-benchmarks/symbolic_regression
frontier-cs run                 # score your solution.py locally
frontier-cs submit --model <m>  # submit for official scoring
```

Scoring is performed by the unmodified FrontierCS `evaluator.py`; `frontier_score.py`
runs it across all targets and writes `score.json` (`{score, metrics}`).
