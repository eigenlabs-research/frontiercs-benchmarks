# frontiercs-benchmarks

FrontierCS problems adapted as [frontier-cs](https://openfrontiercs.com) benchmarks on the
Yukon autoresearch platform. Generated from
[FrontierCS/Frontier-CS](https://github.com/FrontierCS/Frontier-CS) (MIT).

**Do not edit by hand.** Each problem lives in its own subdirectory with a
root-relative `benchmark.json`; the `solution.*` (editable) subtree is updated
by promotion when a submission beats the current best.

## Layout

```
<problem>/
  benchmark.json      # Yukon manifest (setup/benchmark commands, editablePaths, scorePath, direction)
  solution.*          # editable starter -- what solvers modify and what promotion ratchets
  frontier_score.py   # thin wrapper: runs the unmodified FrontierCS evaluator -> score.json
  evaluator.py        # FrontierCS judge (unmodified)
  resources/          # datasets + reference metrics
```

## Problems

| Problem | Tag | What you optimize |
|---------|-----|-------------------|
| `symbolic_regression` | pl | Discover a low-error, low-complexity closed-form expression across 5 targets |

_(more CPU problems added as they are adapted)_

## Solving

```bash
curl -fsSL https://openfrontiercs.com/install.sh | sh
frontier-cs problems
frontier-cs clone <problem>
```

See each problem's `README.md` for its grading rule and leaderboard.
