#!/usr/bin/env bash
# Yukon glue for the FrontierCS **symbolic_regression** research family ONLY. Scoring stays FrontierCS:
# this runs FrontierCS's own evaluator.py and translates its report into Yukon's score.json.
#
# NOTE: research families do NOT share an evaluator CLI. symbolic_regression uses
# --solution-path/--data-dir/--reference-path/--output-path; most other families use --spec-path (and no
# --data-dir/--reference-path), and a few have no argparse at all. This runner is therefore
# family-specific by design — render-glue dispatches per family; a new family needs its own runner.
#
# The evaluator's real CLI (from upstream) is:
#   evaluator.py --solution-path <solution.py> --data-dir <csv dir> \
#                --reference-path <reference_metrics.json> --output-path <report.json>
# It writes a JSON report to --output-path whose `summary.mean_score` is the aggregate on a 0..100
# scale (compute_score clamps to [0,100]), and also prints "<mean_score> <mean_score_unbounded>" to
# stdout. We read the report (authoritative), normalize /100 -> [0,1] (matching algorithmic), and write
# score.json. Scoring stays FrontierCS — we never recompute it.
set -euo pipefail

LEAF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SOLUTION="$LEAF_DIR/solution.py"
DATA_DIR="$LEAF_DIR/resources/data"
REFERENCE="$LEAF_DIR/resources/reference_metrics.json"
REPORT="$LEAF_DIR/.yukon/eval-report.json"
SCORE_OUT="$LEAF_DIR/.yukon/score.json"

[[ -f "$SOLUTION" ]]              || { echo "run.sh: missing $SOLUTION" >&2; exit 1; }
[[ -f "$LEAF_DIR/evaluator.py" ]] || { echo "run.sh: missing evaluator.py in leaf" >&2; exit 1; }
[[ -d "$DATA_DIR" ]]             || { echo "run.sh: missing data dir $DATA_DIR" >&2; exit 1; }
[[ -f "$REFERENCE" ]]            || { echo "run.sh: missing reference metrics $REFERENCE" >&2; exit 1; }

echo "run.sh: evaluating $SOLUTION against $(ls "$DATA_DIR"/*.csv 2>/dev/null | wc -l | tr -d ' ') dataset(s)" >&2

# Run FrontierCS's evaluator. It writes the report to --output-path; its stdout/stderr stream to our log.
CBL_LOG_LEVEL=WARNING python3 "$LEAF_DIR/evaluator.py" \
  --solution-path "$SOLUTION" \
  --data-dir "$DATA_DIR" \
  --reference-path "$REFERENCE" \
  --output-path "$REPORT" >&2

# Translate the evaluator's report -> Yukon score.json. summary.mean_score is already in [0,1];
# carry the full report as metrics. Fail CLOSED with context (never a silent 0).
python3 - "$REPORT" "$SCORE_OUT" <<'PY'
import json, sys
report_path, out = sys.argv[1], sys.argv[2]
try:
    report = json.load(open(report_path))
    summary = report.get("summary") or {}
    if "mean_score" not in summary:
        raise KeyError("report missing summary.mean_score")
    raw = float(summary["mean_score"])   # FrontierCS 0..100 (compute_score clamps to [0,100])
    score = raw / 100.0                    # -> [0,1], matching algorithmic
    # The evaluator also reports an UNbounded aggregate (see header: it prints "<mean_score>
    # <mean_score_unbounded>"). Surface it as scoreRatioUnbounded so a solver beating the reference can
    # exceed 100%; the backend ranks/displays it, falling back to `score` when absent.
    raw_unbounded = summary.get("mean_score_unbounded")
    score_unbounded = float(raw_unbounded) / 100.0 if raw_unbounded is not None else score
    json.dump({"score": score, "scoreRatioUnbounded": score_unbounded, "metrics": report}, open(out, "w"))
    by = report.get("by_dataset") or {}
    print(f"run.sh: ===== mean_score={raw:.4f}/100 -> {score:.6f} (unbounded {score_unbounded:.6f}) over {len(by)} dataset(s) =====", file=sys.stderr)
    for name, entry in by.items():
        print(f"run.sh:   {name}: score={entry.get('score')} mse={entry.get('mse')} complexity={entry.get('complexity')}", file=sys.stderr)
except Exception as e:  # noqa: BLE001 — fail closed with context, never a silent 0
    print(f"run.sh: could not extract score from evaluator report: {e}", file=sys.stderr)
    sys.exit(1)
PY
echo "run.sh: wrote $SCORE_OUT" >&2
