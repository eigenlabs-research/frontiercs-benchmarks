#!/usr/bin/env bash
# Yukon glue for a FrontierCS algorithmic leaf. Scoring stays FrontierCS: this submits the solver's
# C++ to FrontierCS's own judge (a Node server backed by a go-judge sandbox) and writes score.json.
#
# The judge is a SERVICE, not a one-shot CLI, so it must be reachable at $JUDGE_URL. Three ways it can be:
#   1. JUDGE_URL already set (a platform-provided sidecar) — used as-is.
#   2. JUDGE_URL unset but the FULL judge tooling is present here — we boot it (in-script), set JUDGE_URL.
#   3. Neither — fail CLOSED with an actionable message (never a silent 0 score). This is the normal
#      outcome on a solver's own machine: the judge is Linux-only tooling that isn't shipped in the
#      cloned leaf, so local scoring isn't possible — solvers score by submitting instead.
#
# Contract for a sidecar (platform) or self-boot: a judge exposing GET /health, POST /submit
# {pid,lang,code} -> {sid}, GET /result/:sid -> {status,score}. score is 0..100; we normalize to [0,1].
set -euo pipefail

LEAF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SOLUTION="$LEAF_DIR/solution.cpp"
SCORE_OUT="$LEAF_DIR/.yukon/score.json"
PROBLEM_ID="${PROBLEM_ID:-2}"          # rendered per leaf by the generator
MAX_POLL_TIME="${MAX_POLL_TIME:-600}"

[[ -s "$SOLUTION" ]] || { echo "run.sh: missing or empty $SOLUTION" >&2; exit 1; }

# Track every process we self-boot so cleanup can reap them all. We background at PARENT level
# (not inside a subshell) so `$!` is actually set — a subshell-backgrounded `&` leaves `$!` unset,
# which aborts under `set -u`, and its PID would be unreapable (leaking go-judge on :5050).
JUDGE_PIDS=()
cleanup() { for pid in "${JUDGE_PIDS[@]:-}"; do [[ -n "$pid" ]] && kill "$pid" 2>/dev/null || true; done; }
trap cleanup EXIT

# (2) self-boot if no sidecar was provided AND the complete judge tooling is present here. We only
# start it when every piece the boot needs is actually available — the Node judge server entrypoint,
# node, and the go-judge sandbox binary. Booting on a partial toolchain (e.g. the leaf ships the
# judge sources but not the go-judge binary, or the entrypoint isn't materialized) can never pass
# /health, so it would only surface later as a misleading "judge not healthy within 60s". Requiring
# the full set here lets us instead fail fast below with an actionable message.
if [[ -z "${JUDGE_URL:-}" ]]; then
  JUDGE_DIR="$LEAF_DIR/../../judge"    # algorithmic/judge, shared across problems
  JUDGE_ENTRYPOINT="$JUDGE_DIR/src/index.js"
  if [[ -f "$JUDGE_ENTRYPOINT" ]] && command -v node >/dev/null 2>&1 && command -v go-judge >/dev/null 2>&1; then
    ( cd "$JUDGE_DIR" && npm install --silent >/dev/null 2>&1 || true )
    # The judge's Node server needs a go-judge sandbox at 127.0.0.1:5050; start it (presence checked above).
    go-judge >/dev/null 2>&1 &
    JUDGE_PIDS+=("$!")
    ( cd "$JUDGE_DIR" && exec node "$JUDGE_ENTRYPOINT" >/dev/null 2>&1 ) &
    JUDGE_PIDS+=("$!")
    JUDGE_URL="http://127.0.0.1:${JUDGE_PORT:-8081}"
  fi
fi

if [[ -z "${JUDGE_URL:-}" ]]; then
  # No sidecar and the judge can't be booted here. On a solver's own machine this is expected:
  # the FrontierCS judge is Linux-only tooling (go-judge sandbox) that isn't shipped in the cloned
  # leaf, so algorithmic problems can't be scored locally. Say so plainly and point at `submit`,
  # rather than proceeding into a misleading health-check timeout. Still fails CLOSED (never a 0).
  echo "run.sh: local scoring isn't available for algorithmic problems on this machine — the FrontierCS judge (a Linux-only go-judge sandbox) isn't part of the cloned repo. Score by submitting instead: \`frontiercs submit --model <name> --note <text>\` (the platform runs the official judge in its scoring sandbox). To score locally you'd need a running FrontierCS judge reachable via JUDGE_URL. Failing closed." >&2
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
import json, os, sys, time, urllib.request, urllib.parse, urllib.error
judge, pid, sol_path, timeout, out = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), sys.argv[5]
code = open(sol_path).read()

# In prod the platform injects `Authorization: Bearer <token>` at the network layer, so run.sh sends
# no auth. For a LOCAL judge there's no such proxy, so if JUDGE_TOKEN is present in the environment we
# attach it ourselves. Harmless in prod (unset in the scoring sandbox → no header, proxy still injects).
_token = (os.environ.get("JUDGE_TOKEN") or "").strip()
def _auth_headers(extra=None):
    h = dict(extra or {})
    if _token:
        h["Authorization"] = f"Bearer {_token}"
    return h

def judge_open(req):
    # Single choke point for the judge's authenticated calls. /health is unauthenticated, so a
    # JUDGE_URL pointed at the hosted judge passes the health probe but rejects /submit (and can gate
    # /result too, or a token can expire mid-poll) without the solver token — which the platform
    # injects only inside its scoring sandboxes. Translate that (401/403) into an actionable message
    # instead of surfacing a bare "HTTP Error 401" from whichever call hit the auth boundary first.
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            return json.loads(r.read())
    except urllib.error.HTTPError as e:
        if e.code in (401, 403):
            raise RuntimeError(
                f"judge rejected the request ({e.code}). The FrontierCS judge at {judge} requires a "
                f"solver token that is only injected inside the platform's scoring sandboxes, so it can't "
                f"be scored from here. Score by submitting instead: `frontiercs submit`."
            ) from e
        raise

def post_submit():
    body = urllib.parse.urlencode({"pid": pid, "lang": "cpp", "code": code}).encode()
    req = urllib.request.Request(f"{judge}/submit", data=body,
                                 headers=_auth_headers({"Content-Type": "application/x-www-form-urlencoded"}))
    return judge_open(req)["sid"]

def get_result(sid):
    return judge_open(urllib.request.Request(f"{judge}/result/{sid}", headers=_auth_headers()))

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
    # The judge also reports an UNbounded aggregate (per-case ratios not clamped to <=1), so a solver
    # that beats the reference scores >1.0. Surface it as scoreRatioUnbounded; the backend ranks/displays
    # it (falling back to `score` when absent). Older judges omit it -> fall back to the bounded score.
    raw_unbounded = result.get("scoreUnbounded")
    score_unbounded = float(raw_unbounded) / 100.0 if raw_unbounded is not None else score
    # CI-style per-case breakdown from the judge result, so the run log shows every test case.
    cases = result.get("cases") or (result.get("metrics") or {}).get("cases") or []
    print(f"run.sh: ===== judge result: {raw}/100 -> {score:.4f} (unbounded {score_unbounded:.4f}) over {len(cases)} case(s) =====", file=sys.stderr)
    for i, c in enumerate(cases):
        ok = "ok" if c.get("ok") else "FAIL"
        msg = str(c.get("msg", "")).strip().replace("\n", " ")[:80]
        print(f"run.sh:   case {i:>3} [{ok}] time={c.get('time','?')} mem={c.get('memory','?')} | {msg}", file=sys.stderr)
    json.dump({"score": score, "scoreRatioUnbounded": score_unbounded, "metrics": result}, open(out, "w"))
except Exception as e:  # noqa: BLE001 — fail closed with context, never a silent 0
    print(f"run.sh: judge scoring failed: {e}", file=sys.stderr)
    sys.exit(1)
PY
echo "run.sh: wrote $SCORE_OUT" >&2
