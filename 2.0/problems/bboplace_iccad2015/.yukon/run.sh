#!/usr/bin/env bash
# Yukon glue for a FrontierCS 2.0 leaf. Scoring stays FrontierCS: this runs the leaf's own
# evaluator.py (which scores the solver's solution.py) and writes score.json.
#
# 2.0's evaluator needs the benchmark data the upstream "judge_image" ships, exposed at
# $BBOPLACE_ROOT (default /opt/bboplace-bench). It can be provided three ways:
#   1. BBOPLACE_ROOT already points at the data (platform-mounted the judge image) — used as-is.
#   2. BBOPLACE_ROOT unset but docker + the judge image are available here — we extract it (in-script).
#   3. Neither — fail CLOSED with a clear message (never a silent 0 score).
set -euo pipefail

LEAF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SOLUTION="$LEAF_DIR/solution.py"
SCORE_OUT="$LEAF_DIR/.yukon/score.json"
JUDGE_IMAGE="${JUDGE_IMAGE:-ghcr.io/frontiercs/frontiercs-bboplace-data:2026-06-ispd-iccad}"        # rendered per leaf by the generator (from config.yaml)

[[ -s "$SOLUTION" ]] || { echo "run.sh: missing or empty $SOLUTION" >&2; exit 1; }

CONTAINER=""
cleanup() { [[ -n "$CONTAINER" ]] && docker rm -f "$CONTAINER" >/dev/null 2>&1 || true; }
trap cleanup EXIT

# (2) self-provide the data from the judge image if no BBOPLACE_ROOT was mounted.
if [[ -z "${BBOPLACE_ROOT:-}" ]]; then
  if command -v docker >/dev/null 2>&1 && [[ -n "$JUDGE_IMAGE" && "$JUDGE_IMAGE" != "ghcr.io/frontiercs/frontiercs-bboplace-data:2026-06-ispd-iccad" ]]; then
    EXTRACT_DIR="$LEAF_DIR/.yukon/bboplace-data"
    mkdir -p "$EXTRACT_DIR"
    # Create (not run) a container from the data image and copy its /opt/bboplace-bench out.
    CONTAINER="$(docker create "$JUDGE_IMAGE" 2>/dev/null || true)"
    if [[ -n "$CONTAINER" ]] && docker cp "$CONTAINER:/opt/bboplace-bench/." "$EXTRACT_DIR" >/dev/null 2>&1; then
      BBOPLACE_ROOT="$EXTRACT_DIR"
    fi
  fi
fi

if [[ -z "${BBOPLACE_ROOT:-}" ]]; then
  echo "run.sh: 2.0 benchmark data unavailable (set BBOPLACE_ROOT to the judge-image data, or provide docker + the judge image for in-script extraction). Failing closed." >&2
  exit 1
fi

# Run FrontierCS's own evaluator; it prints a JSON payload with a score. Point it at the solver's
# solution (upstream evaluate.sh hardcodes a Harbor mount path; we pass the leaf-local file instead).
EVAL_JSON="$(BBOPLACE_ROOT="$BBOPLACE_ROOT" python3 "$LEAF_DIR/evaluator.py" "$SOLUTION")"

python3 - "$EVAL_JSON" "$SCORE_OUT" <<'PY'
import json, sys
try:
    raw = sys.argv[1]
    p = json.loads(raw)
    if "score" not in p:
        raise KeyError('"score" field missing')
    json.dump({"score": float(p["score"]), "metrics": p}, open(sys.argv[2], "w"))
except (json.JSONDecodeError, KeyError, ValueError, TypeError) as e:
    print(f"run.sh: could not extract score from evaluator output: {e}", file=sys.stderr)
    print(f"Raw output (first 500 chars): {raw[:500]}", file=sys.stderr)
    sys.exit(1)
PY
echo "run.sh: wrote $SCORE_OUT" >&2
