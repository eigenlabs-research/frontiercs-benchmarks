#!/usr/bin/env bash
# Yukon glue for a FrontierCS research leaf. Invokes FrontierCS's own evaluator and
# translates its output into Yukon's score.json. Scoring stays FrontierCS.
set -euo pipefail
LEAF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SOLUTION="$LEAF_DIR/solution.py"
SPEC="$LEAF_DIR/resources/submission_spec.json"
[[ -f "$SOLUTION" ]] || { echo "run.sh: missing $SOLUTION" >&2; exit 1; }
EVAL_JSON="$(CBL_LOG_LEVEL=WARNING python3 "$LEAF_DIR/evaluator.py" --solution "$SOLUTION" --spec "$SPEC")"
python3 - "$EVAL_JSON" "$LEAF_DIR/.yukon/score.json" <<'PY'
import json, sys
try:
    raw = sys.argv[1]
    p = json.loads(raw)
    if "score" not in p:
        raise KeyError('"score" field missing')
    score_val = float(p["score"])
    json.dump({"score": score_val, "metrics": p}, open(sys.argv[2], "w"))
except (json.JSONDecodeError, KeyError, ValueError, TypeError) as e:
    print(f"run.sh: could not extract score from evaluator output: {e}", file=sys.stderr)
    print(f"Raw output (first 500 chars): {raw[:500]}", file=sys.stderr)
    sys.exit(1)
PY
echo "run.sh: wrote $LEAF_DIR/.yukon/score.json" >&2
