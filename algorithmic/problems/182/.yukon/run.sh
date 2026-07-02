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
PROBLEM_ID="${PROBLEM_ID:-182}"          # rendered per leaf by the generator
MAX_POLL_TIME="${MAX_POLL_TIME:-600}"

[[ -s "$SOLUTION" ]] || { echo "run.sh: missing or empty $SOLUTION" >&2; exit 1; }

# Track every process we self-boot so cleanup can reap them all. We background at PARENT level
# (not inside a subshell) so `$!` is actually set — a subshell-backgrounded `&` leaves `$!` unset,
# which aborts under `set -u`, and its PID would be unreapable (leaking go-judge on :5050).
JUDGE_PIDS=()
cleanup() { for pid in "${JUDGE_PIDS[@]:-}"; do [[ -n "$pid" ]] && kill "$pid" 2>/dev/null || true; done; }
trap cleanup EXIT

# (2) self-boot if no sidecar was provided and the judge is present in the repo.
if [[ -z "${JUDGE_URL:-}" ]]; then
  JUDGE_DIR="$LEAF_DIR/../../judge"    # algorithmic/judge, shared across problems
  if [[ -d "$JUDGE_DIR" ]] && command -v node >/dev/null 2>&1; then
    ( cd "$JUDGE_DIR" && npm install --silent >/dev/null 2>&1 || true )
    # The judge's Node server also needs a go-judge sandbox at 127.0.0.1:5050; if go-judge is
    # available, start it, else the judge boot will fail its /health and we exit closed below.
    if command -v go-judge >/dev/null 2>&1; then
      go-judge >/dev/null 2>&1 &
      JUDGE_PIDS+=("$!")
    fi
    ( cd "$JUDGE_DIR" && exec node src/index.js >/dev/null 2>&1 ) &
    JUDGE_PIDS+=("$!")
    JUDGE_URL="http://127.0.0.1:${JUDGE_PORT:-8081}"
  fi
fi

if [[ -z "${JUDGE_URL:-}" ]]; then
  echo "run.sh: no judge service available (set JUDGE_URL to a running FrontierCS judge, or provide the judge tooling for in-script boot). Failing closed." >&2
  exit 1
fi

# Wait for the judge to be healthy. Use python3 (guaranteed present in the sandbox; curl may not be).
deadline=$(( $(date +%s) + 60 ))
until python3 -c "import sys,urllib.request; urllib.request.urlopen(sys.argv[1].rstrip('/')+'/health',timeout=5)" "$JUDGE_URL" >/dev/null 2>&1; do
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
    print(f"run.sh: submitting solution ({len(code)} bytes) to {judge} for problem {pid}", file=sys.stderr)
    sid = post_submit()
    print(f"run.sh: submission id={sid}; polling /result (timeout {timeout}s)", file=sys.stderr)
    deadline = time.time() + timeout
    result = None
    last_status = None
    polls = 0
    while time.time() < deadline:
        result = get_result(sid)
        polls += 1
        st = result.get("status")
        if st != last_status:
            print(f"run.sh: [poll {polls}] status={st}", file=sys.stderr)
            last_status = st
        if st in ("done", "finished", "error", "completed"):
            break
        time.sleep(3)
    if result is None:
        raise TimeoutError("no result before deadline")
    raw = result.get("score")
    if raw is None:
        raise KeyError('judge result missing "score"')
    score = float(raw) / 100.0   # FrontierCS scores 0..100 -> [0,1]
    # CI-style per-case breakdown from the judge result, so the run log shows every test case.
    cases = result.get("cases") or (result.get("metrics") or {}).get("cases") or []
    print(f"run.sh: ===== judge result: {raw}/100 -> {score:.4f} over {len(cases)} case(s) =====", file=sys.stderr)
    for i, c in enumerate(cases):
        ok = "ok" if c.get("ok") else "FAIL"
        msg = str(c.get("msg", "")).strip().replace("\n", " ")[:80]
        print(f"run.sh:   case {i:>3} [{ok}] time={c.get('time','?')} mem={c.get('memory','?')} | {msg}", file=sys.stderr)
    json.dump({"score": score, "metrics": result}, open(out, "w"))
except Exception as e:  # noqa: BLE001 — fail closed with context, never a silent 0
    print(f"run.sh: judge scoring failed: {e}", file=sys.stderr)
    sys.exit(1)
PY
echo "run.sh: wrote $SCORE_OUT" >&2
