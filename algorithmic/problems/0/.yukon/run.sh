#!/usr/bin/env bash
# Yukon glue for a FrontierCS algorithmic leaf. Scoring stays FrontierCS: this submits the solver's
# C++ to FrontierCS's own judge (a Node server backed by a go-judge sandbox) and writes score.json.
#
# The judge is a SERVICE, not a one-shot CLI, so it must be reachable at $JUDGE_URL. Three ways it can be:
#   1. JUDGE_URL already set (a platform-provided sidecar) — used as-is.
#   2. JUDGE_URL unset but the judge tooling is present here — we boot it (in-script), set JUDGE_URL.
#   3. Neither — fail CLOSED with a clear message (never a silent 0 score).
#
# Contract for a sidecar (platform) or self-boot: a judge exposing GET /health, POST /submit
# {pid,lang,code} -> {sid}, GET /result/:sid -> {status,score}. score is 0..100; we normalize to [0,1].
set -euo pipefail

LEAF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SOLUTION="$LEAF_DIR/solution.cpp"
SCORE_OUT="$LEAF_DIR/.yukon/score.json"
PROBLEM_ID="${PROBLEM_ID:-0}"          # rendered per leaf by the generator
MAX_POLL_TIME="${MAX_POLL_TIME:-600}"

[[ -s "$SOLUTION" ]] || { echo "run.sh: missing or empty $SOLUTION" >&2; exit 1; }

JUDGE_STARTED_BY_US=""
cleanup() { [[ -n "$JUDGE_STARTED_BY_US" ]] && kill "$JUDGE_STARTED_BY_US" 2>/dev/null || true; }
trap cleanup EXIT

# (2) self-boot if no sidecar was provided and the judge is present in the repo.
if [[ -z "${JUDGE_URL:-}" ]]; then
  JUDGE_DIR="$LEAF_DIR/../../judge"    # algorithmic/judge, shared across problems
  if [[ -d "$JUDGE_DIR" ]] && command -v node >/dev/null 2>&1; then
    ( cd "$JUDGE_DIR" && npm install --silent >/dev/null 2>&1 || true )
    # The judge's Node server also needs a go-judge sandbox at 127.0.0.1:5050; if go-judge is
    # available, start it, else the judge boot will fail its /health and we exit closed below.
    command -v go-judge >/dev/null 2>&1 && ( go-judge >/dev/null 2>&1 & )
    ( cd "$JUDGE_DIR" && node src/index.js >/dev/null 2>&1 & ) && JUDGE_STARTED_BY_US=$!
    JUDGE_URL="http://127.0.0.1:${JUDGE_PORT:-8082}"
  fi
fi

if [[ -z "${JUDGE_URL:-}" ]]; then
  echo "run.sh: no judge service available (set JUDGE_URL to a running FrontierCS judge, or provide the judge tooling for in-script boot). Failing closed." >&2
  exit 1
fi

# Wait for the judge to be healthy.
deadline=$(( $(date +%s) + 60 ))
until curl -fsS "$JUDGE_URL/health" >/dev/null 2>&1; do
  [[ $(date +%s) -ge $deadline ]] && { echo "run.sh: judge at $JUDGE_URL not healthy within 60s" >&2; exit 1; }
  sleep 2
done

# Submit the solution, then poll /result until terminal, then normalize score to [0,1].
python3 - "$JUDGE_URL" "$PROBLEM_ID" "$SOLUTION" "$MAX_POLL_TIME" "$SCORE_OUT" <<'PY'
import json, sys, time, urllib.request, urllib.parse
judge, pid, sol_path, timeout, out = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), sys.argv[5]
code = open(sol_path).read()

def post_submit():
    body = urllib.parse.urlencode({"pid": pid, "lang": "cpp", "code": code}).encode()
    req = urllib.request.Request(f"{judge}/submit", data=body,
                                 headers={"Content-Type": "application/x-www-form-urlencoded"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read())["sid"]

def get_result(sid):
    with urllib.request.urlopen(f"{judge}/result/{sid}", timeout=60) as r:
        return json.loads(r.read())

try:
    sid = post_submit()
    deadline = time.time() + timeout
    result = None
    while time.time() < deadline:
        result = get_result(sid)
        if result.get("status") in ("done", "finished", "error", "completed"):
            break
        time.sleep(3)
    if result is None:
        raise TimeoutError("no result before deadline")
    raw = result.get("score")
    if raw is None:
        raise KeyError('judge result missing "score"')
    score = float(raw) / 100.0   # FrontierCS scores 0..100 -> [0,1]
    json.dump({"score": score, "metrics": result}, open(out, "w"))
except Exception as e:  # noqa: BLE001 — fail closed with context, never a silent 0
    print(f"run.sh: judge scoring failed: {e}", file=sys.stderr)
    sys.exit(1)
PY
echo "run.sh: wrote $SCORE_OUT" >&2
