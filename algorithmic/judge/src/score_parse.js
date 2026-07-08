// Per-case score derivation for the FrontierCS algorithmic judge.
//
// SINGLE SOURCE OF THE TRUST RULE. This file exists in three places and score_parse.js is the
// authority for all of them:
//   1. scripts/frontiercs/judge-image/judge/src/  (this file) — the yukon PROD judge image.
//   2. algorithmic/judge/src/ in eigenlabs-research/frontiercs-benchmarks — the LOCAL judge solvers
//      run. It is GENERATED, not hand-maintained: scripts/frontiercs/generate.ts copies the pristine
//      UPSTREAM judge, then patch-judge.ts writes THIS file in and re-applies the fix to
//      judge_engine.js. So the benchmarks copy is byte-identical to this one by construction.
//   3. UPSTREAM FrontierCS/Frontier-CS (algorithmic/judge/src/) — still VULNERABLE; we don't control
//      it, hence the deterministic post-copy patch in (2). The bug should be reported upstream.
// judge_engine.js has legitimate per-copy drift (prod-only cleanup/zero-case guards), so it is NOT
// byte-identical across copies; the trust INVARIANT is what must hold, and it is gated by
// tests/unit/frontiercs-judge-trust-invariant.test.ts. Keep this file identical across (1) and (2).
//
// SECURITY-CRITICAL: optimization problems report their score via "Ratio:" / "RatioUnbounded:"
// markers, but those markers must be trusted ONLY when they come from the judge's own checker or
// interactor — never from the contestant program's own stdout/stderr. A contestant can print
// anything to its own streams; if we parsed markers from there, a program that never produces a
// valid answer could exit nonzero (so the checker never runs), print "RatioUnbounded: 9.0000" to
// its stderr, and be credited a 9x score. That exact exploit topped the prod traveling-santa (pid
// 44) leaderboard. When the case result carries `msgFromContestant: true`, `msg` is the contestant's
// own stream and MUST NOT be mined for score markers — such a run scores 0 (it failed to pass the
// checker, by definition).
//
// The judge sets `r.msgFromContestant = true` on every early-return where the contestant program
// failed to run to a clean exit (so the checker/interactor never scored it). On the normal path,
// `msg` is the checker/interactor output and `msgFromContestant` is falsy → markers are trusted.

// Returns { scoreRatio, scoreRatioUnbounded } for a judged case result `r`:
//   r.ok               — did the checker/interactor accept the run (exit 0)
//   r.msg              — checker/interactor output (trusted) OR contestant stream (untrusted)
//   r.msgFromContestant— when true, `msg` is the contestant's own stream; ignore markers, score 0
export function scoreCaseResult(r) {
    // Untrusted message (contestant stream): the checker did not score this run. It failed — score 0.
    if (r.msgFromContestant) {
        return { scoreRatio: 0, scoreRatioUnbounded: 0 };
    }
    // CAVEAT: on the trusted path we mine the first marker anywhere in `msg`; this assumes the
    // checker/interactor never echoes untrusted contestant text (which could contain a marker) back
    // into its own output — true for the current testlib checkers, but a constraint on future ones.
    const matchBounded = r.msg.match(/Ratio: ([\d.]+)/);
    const matchUnbounded = r.msg.match(/RatioUnbounded: ([\d.]+)/);
    const scoreRatio = matchBounded ? parseFloat(matchBounded[1]) : (r.ok ? 1.0 : 0);
    const scoreRatioUnbounded = matchUnbounded ? parseFloat(matchUnbounded[1]) : scoreRatio;
    return { scoreRatio, scoreRatioUnbounded };
}
