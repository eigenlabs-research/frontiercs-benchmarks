#!/usr/bin/env python3
"""Yukon score wrapper for the frontier-cs symbolic_regression benchmark.

Runs the UNMODIFIED FrontierCS evaluator (evaluator.py) over every target
dataset in resources/data, then emits Yukon's score.json contract:

    {"score": <mean 0..100>, "metrics": {"<dataset>": <score>, ...}}

The FrontierCS evaluator already computes the per-dataset score and their mean,
so we simply translate its report. A non-zero evaluator exit (e.g. an
unparseable expression from the solver) yields score 0 -- the Yukon baseline
floor -- rather than failing the run, matching FrontierCS semantics where an
invalid submission scores 0.

NOTE (integrity): the target datasets are currently committed in-repo (public).
Before this benchmark is opened for real competition it must be switched to a
hidden holdout fetched at setup (see docs/frontier-cs-plan.md, decision C), so
the discovered expression is scored on data the solver never saw.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.abspath(__file__))
REPORT = os.path.join(ROOT, "result.json")
SCORE = os.path.join(ROOT, "score.json")


def write_score(score: float, metrics: dict) -> None:
    with open(SCORE, "w", encoding="utf-8") as fh:
        json.dump({"score": float(score), "metrics": metrics}, fh, indent=2)


def main() -> int:
    cmd = [
        sys.executable,
        os.path.join(ROOT, "evaluator.py"),
        "--solution-path", os.path.join(ROOT, "solution.py"),
        "--data-dir", os.path.join(ROOT, "resources", "data"),
        "--reference-path", os.path.join(ROOT, "resources", "reference_metrics.json"),
        "--output-path", REPORT,
    ]
    proc = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
    if proc.stdout:
        sys.stdout.write(proc.stdout)
    if proc.stderr:
        sys.stderr.write(proc.stderr)

    if proc.returncode != 0 or not os.path.exists(REPORT):
        # Invalid / crashing submission -> score 0 (baseline floor).
        write_score(0.0, {"error": (proc.stderr or "evaluator failed")[-800:]})
        return 0

    with open(REPORT, encoding="utf-8") as fh:
        report = json.load(fh)

    summary = report.get("summary", {})
    by_dataset = report.get("by_dataset", {})
    metrics = {
        name: round(float(entry.get("score", 0.0)), 6)
        for name, entry in by_dataset.items()
    }
    metrics["mean_score_unbounded"] = round(float(summary.get("mean_score_unbounded", 0.0)), 6)
    metrics["mean_mse"] = float(summary.get("mean_mse", 0.0))
    write_score(float(summary.get("mean_score", 0.0)), metrics)
    return 0


if __name__ == "__main__":
    sys.exit(main())
