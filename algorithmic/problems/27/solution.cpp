#include <bits/stdc++.h>
using namespace std;

/* -------------------------------------------------------------
   1.  Finite field GF(q)  –  q = p^e   (prime power)
   ------------------------------------------------------------- */
struct Field {
    int p = 0, e = 0, q = 0;               // q = p^e
    vector<vector<int>> add, mul;          // tables
    vector<int> neg, inv;                  // additive / multiplicative inverse
};

/* ---------- helpers ------------------------------------------------ */
static bool isPrime(int x) {
    if (x < 2) return false;
    for (int d = 2; d * d <= x; ++d) if (x % d == 0) return false;
    return true;
}
static pair<int,int> factorPrimePower(int q) {    // q = p^e
    for (int p = 2; p * p <= q; ++p) if (q % p == 0 && isPrime(p)) {
        int e = 0, t = q;
        while (t % p == 0) { t /= p; ++e; }
        if (t == 1) return {p, e};
    }
    return {q, 1};
}

/* ---------- polynomial utilities (coefficients in Zp) ------------- */
using Poly = vector<int>;

static void trim(Poly &a) {
    while (!a.empty() && a.back() == 0) a.pop_back();
}
static int degree(const Poly &a) { return (int)a.size() - 1; }

static Poly poly_mod(const Poly &a, const Poly &b, int p) {
    Poly r = a;
    int db = degree(b);
    while (degree(r) >= db && db >= 0) {
        int shift = degree(r) - db;
        int coeff = r.back();                // leading coeff (b is monic)
        for (int i = 0; i <= db; ++i) {
            int idx = i + shift;
            r[idx] = (r[idx] - coeff * b[i]) % p;
            if (r[idx] < 0) r[idx] += p;
        }
        trim(r);
    }
    return r;
}

/* ---------- irreducibility test (deterministic, tiny fields) ----- */
static bool isIrreducible(const Poly &f, int p) {
    int e = degree(f);
    for (int d = 1; d <= e / 2; ++d) {
        long long total = 1;
        for (int i = 0; i < d; ++i) total *= p;
        for (long long mask = 0; mask < total; ++mask) {
            Poly g(d + 1, 0);
            long long tmp = mask;
            for (int i = 0; i < d; ++i) {
                g[i] = int(tmp % p);
                tmp /= p;
            }
            g[d] = 1;                     // monic
            if (poly_mod(f, g, p).empty()) return false;
        }
    }
    return true;
}
static Poly findIrreduciblePoly(int p, int e) {
    long long total = 1;
    for (int i = 0; i < e; ++i) total *= p;
    for (long long mask = 0; mask < total; ++mask) {
        Poly f(e + 1, 0);
        long long tmp = mask;
        for (int i = 0; i < e; ++i) {
            f[i] = int(tmp % p);
            tmp /= p;
        }
        f[e] = 1;
        if (isIrreducible(f, p)) return f;
    }
    assert(false);
    return {};
}

/* ---------- build field tables for a given prime power q ----------- */
static Field buildField(int q) {
    Field F;  F.q = q;
    auto [p, e] = factorPrimePower(q);
    F.p = p;  F.e = e;
    F.add.assign(q, vector<int>(q));
    F.mul.assign(q, vector<int>(q));
    F.neg.assign(q, 0);
    F.inv.assign(q, 0);

    if (e == 1) {                 // prime field – trivial tables
        for (int a = 0; a < q; ++a)
            for (int b = 0; b < q; ++b) {
                int s = a + b; if (s >= p) s -= p;
                F.add[a][b] = s;
                F.mul[a][b] = int(1LL * a * b % p);
            }
        for (int a = 0; a < q; ++a)
            for (int b = 0; b < q; ++b)
                if (F.add[a][b] == 0) { F.neg[a] = b; break; }
        F.inv[1] = 1;
        for (int a = 2; a < q; ++a)
            for (int b = 1; b < q; ++b)
                if (F.mul[a][b] == 1) { F.inv[a] = b; break; }
        return F;
    }

    /* ---------- extension field ------------------------------------- */
    Poly irr = findIrreduciblePoly(p, e);          // monic irreducible of degree e

    // powers of p for digit encoding
    vector<int> powp(e);
    powp[0] = 1;
    for (int i = 1; i < e; ++i) powp[i] = powp[i-1] * p;

    auto decode = [&](int x, array<int,8> &out) {
        for (int i = 0; i < e; ++i) {
            out[i] = x % p;
            x /= p;
        }
    };

    // addition – componentwise
    for (int a = 0; a < q; ++a)
        for (int b = 0; b < q; ++b) {
            array<int,8> A{}, B{}, C{};
            decode(a, A); decode(b, B);
            int v = 0;
            for (int i = 0; i < e; ++i) {
                int t = A[i] + B[i];
                if (t >= p) t -= p;
                C[i] = t;
                v += C[i] * powp[i];
            }
            F.add[a][b] = v;
        }

    // additive inverses
    for (int a = 0; a < q; ++a)
        for (int b = 0; b < q; ++b)
            if (F.add[a][b] == 0) { F.neg[a] = b; break; }

    // multiplication – polynomial multiplication modulo irr
    for (int a = 0; a < q; ++a)
        for (int b = 0; b < q; ++b) {
            array<int,8> A{}, B{};
            decode(a, A); decode(b, B);
            int maxDeg = 2 * e - 2;
            vector<int> prod(maxDeg + 1, 0);
            for (int i = 0; i < e; ++i)
                for (int j = 0; j < e; ++j) {
                    int cur = prod[i + j] + A[i] * B[j];
                    prod[i + j] = cur % p;
                }
            // reduction modulo irr (irr is monic)
            for (int d = maxDeg; d >= e; --d) {
                int coeff = prod[d];
                if (!coeff) continue;
                int shift = d - e;
                for (int i = 0; i <= e; ++i) {
                    int idx = i + shift;
                    int cur = prod[idx] - coeff * irr[i];
                    cur %= p; if (cur < 0) cur += p;
                    prod[idx] = cur;
                }
            }
            int v = 0;
            for (int i = 0; i < e; ++i) v += prod[i] * powp[i];
            F.mul[a][b] = v;
        }

    // multiplicative inverses (0 unused)
    F.inv[1] = 1;
    for (int a = 2; a < q; ++a)
        for (int b = 1; b < q; ++b)
            if (F.mul[a][b] == 1) { F.inv[a] = b; break; }

    return F;
}

/* -------------------------------------------------------------
   2.  Counting / emitting helpers – affine design
   ------------------------------------------------------------- */
static long long countAffine(const Field &F, bool rowsArePoints,
                             int R, int C) {
    int q = F.q;
    long long cnt = 0;
    if (rowsArePoints) {                 // rows = points (x,y)
        for (int r = 0; r < R; ++r) {
            int x = r / q, y = r % q;
            // non‑vertical lines
            for (int a = 0; a < q; ++a) {
                int ax = F.mul[a][x];
                int b  = F.add[y][F.neg[ax]];
                int col = a * q + b;                     // 0 … q²-1
                if (col < C) ++cnt;
            }
            // vertical line x = const
            int col = q * q + x;                         // q² … q²+q-1
            if (col < C) ++cnt;
        }
    } else {                              // rows = lines, columns = points
        for (int r = 0; r < R; ++r) {
            if (r < q * q) {                // non‑vertical line y = a·x + b
                int a = r / q;
                int b = r % q;
                for (int x = 0; x < q; ++x) {
                    int ax = F.mul[a][x];
                    int y  = F.add[ax][b];
                    int col = x * q + y;               // point index
                    if (col < C) ++cnt;
                }
            } else {                         // vertical line x = const
                int x = r - q * q;
                for (int y = 0; y < q; ++y) {
                    int col = x * q + y;
                    if (col < C) ++cnt;
                }
            }
        }
    }
    return cnt;
}
static void emitAffine(const Field &F, bool rowsArePoints,
                       int R, int C,
                       vector<pair<int,int>> &out) {
    int q = F.q;
    if (rowsArePoints) {
        for (int r = 0; r < R; ++r) {
            int x = r / q, y = r % q;
            for (int a = 0; a < q; ++a) {
                int ax = F.mul[a][x];
                int b  = F.add[y][F.neg[ax]];
                int col = a * q + b;
                if (col < C) out.emplace_back(r + 1, col + 1);
            }
            int col = q * q + x;
            if (col < C) out.emplace_back(r + 1, col + 1);
        }
    } else {
        for (int r = 0; r < R; ++r) {
            if (r < q * q) {
                int a = r / q, b = r % q;
                for (int x = 0; x < q; ++x) {
                    int ax = F.mul[a][x];
                    int y  = F.add[ax][b];
                    int col = x * q + y;
                    if (col < C) out.emplace_back(r + 1, col + 1);
                }
            } else {
                int x = r - q * q;
                for (int y = 0; y < q; ++y) {
                    int col = x * q + y;
                    if (col < C) out.emplace_back(r + 1, col + 1);
                }
            }
        }
    }
}

/* -------------------------------------------------------------
   3.  Counting / emitting helpers – projective plane
   ------------------------------------------------------------- */
static long long countProjective(const Field &F, int R, int C) {
    int q = F.q;
    long long cnt = 0;
    for (int r = 0; r < R; ++r) {
        if (r < q * q) {                     // ordinary point (x,y)
            int x = r / q, y = r % q;
            for (int a = 0; a < q; ++a) {
                int ax = F.mul[a][x];
                int b  = F.add[y][F.neg[ax]];
                int col = a * q + b;                     // non‑vertical line
                if (col < C) ++cnt;
            }
            int col = q * q + x;                         // vertical line
            if (col < C) ++cnt;
        } else {                              // point at infinity
            int d = r - q * q;
            if (d < q) {                     // direction ∞ (slope‑∞)
                for (int b = 0; b < q; ++b) {
                    int col = d * q + b;
                    if (col < C) ++cnt;
                }
            } else {                         // vertical‑∞
                for (int c = 0; c < q; ++c) {
                    int col = q * q + c;
                    if (col < C) ++cnt;
                }
            }
            int col = q * q + q;               // line at infinity
            if (col < C) ++cnt;
        }
    }
    return cnt;
}
static void emitProjective(const Field &F,
                           int R, int C,
                           vector<pair<int,int>> &out) {
    int q = F.q;
    for (int r = 0; r < R; ++r) {
        if (r < q * q) {                     // ordinary point (x,y)
            int x = r / q, y = r % q;
            for (int a = 0; a < q; ++a) {
                int ax = F.mul[a][x];
                int b  = F.add[y][F.neg[ax]];
                int col = a * q + b;                     // non‑vertical line
                if (col < C) out.emplace_back(r + 1, col + 1);
            }
            int col = q * q + x;                         // vertical line
            if (col < C) out.emplace_back(r + 1, col + 1);
        } else {                              // point at infinity
            int d = r - q * q;
            if (d < q) {                     // direction ∞
                for (int b = 0; b < q; ++b) {
                    int col = d * q + b;
                    if (col < C) out.emplace_back(r + 1, col + 1);
                }
            } else {                         // vertical‑∞
                for (int c = 0; c < q; ++c) {
                    int col = q * q + c;
                    if (col < C) out.emplace_back(r + 1, col + 1);
                }
            }
            int col = q * q + q;               // line at infinity
            if (col < C) out.emplace_back(r + 1, col + 1);
        }
    }
}

/* -------------------------------------------------------------
   4.  Counting / emitting helpers – polynomial design
   ------------------------------------------------------------- */
static long long countPolynomial(const Field &F, int p, int R, int C) {
    // rows = (i,x), 0≤i<p, x∈GF(q)
    // columns = (j,y), 0≤j<p, y∈GF(q)
    // edge (i,x)-(j,y) iff y = j·x + i   (mod q)
    long long cnt = 0;
    int q = F.q;
    for (int i = 0; i < p; ++i) {
        for (int x = 0; x < q; ++x) {
            int r = i * q + x;               // row id (0‑based)
            if (r >= R) continue;
            for (int j = 0; j < p; ++j) {
                int y = F.add[F.mul[j][x]][i];
                int c = j * q + y;           // column id
                if (c < C) ++cnt;
            }
        }
    }
    return cnt;
}
static void emitPolynomial(const Field &F, int p,
                           int R, int C,
                           vector<pair<int,int>> &out) {
    int q = F.q;
    for (int i = 0; i < p; ++i) {
        for (int x = 0; x < q; ++x) {
            int r = i * q + x;               // row id
            if (r >= R) continue;
            for (int j = 0; j < p; ++j) {
                int y = F.add[F.mul[j][x]][i];
                int c = j * q + y;           // column id
                if (c < C) out.emplace_back(r + 1, c + 1);
            }
        }
    }
}

/* -------------------------------------------------------------
   5.   Candidate description
   ------------------------------------------------------------- */
struct Candidate {
    long long cells = -1;               // number of ones after augmentation
    enum Kind { STAR, AFFINE, PROJECTIVE, POLYNOMIAL } kind = STAR;
    bool rowsArePoints = true;          // used by AFFINE / PROJECTIVE
    bool transpose = false;             // PROJECTIVE only
    int q = 0;                          // order of the geometry (0 for star)
    int p = 0;                          // order for polynomial design (p ≤ q)
    int R = 0, C = 0;                   // rows / columns actually used
    const Field *F = nullptr;
};

/* -------------------------------------------------------------
   6.   Main
   ------------------------------------------------------------- */
int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n, m;
    if (!(cin >> n >> m)) return 0;

    /* ---------- trivial star (one full row + one full column) ----- */
    Candidate best;
    best.cells = (long long)n + m - 1;
    best.kind = Candidate::STAR;

    /* ---------- enumerate all relevant prime powers ---------------- */
    const int LIM = 320;               // because q² ≤ n·m ≤ 1e5
    vector<int> primePowers;
    vector<bool> isPrime(LIM + 1, true);
    isPrime[0] = isPrime[1] = false;
    for (int i = 2; i * i <= LIM; ++i)
        if (isPrime[i])
            for (int j = i * i; j <= LIM; j += i) isPrime[j] = false;
    for (int p = 2; p <= LIM; ++p) if (isPrime[p]) {
        long long cur = p;
        while (cur <= LIM) {
            primePowers.push_back((int)cur);
            cur *= p;
        }
    }
    sort(primePowers.begin(), primePowers.end());
    primePowers.erase(unique(primePowers.begin(), primePowers.end()),
                      primePowers.end());

    vector<Field> fields; fields.reserve(primePowers.size());

    for (int q : primePowers) {
        fields.emplace_back(buildField(q));
        const Field *F = &fields.back();

        /* ---------- affine (rows = points) ---------- */
        {
            int maxRows = q * q;
            int maxCols = q * q + q;
            int R = min(n, maxRows);
            int C = min(m, maxCols);
            long long cur = countAffine(*F, true, R, C);
            long long final = cur + (m - C);
            if (final > best.cells) {
                best = {final, Candidate::AFFINE, true, false,
                        q, 0, R, C, F};
            }
        }
        /* ---------- affine (rows = lines) ----------- */
        {
            int maxRows = q * q + q;
            int maxCols = q * q;
            int R = min(n, maxRows);
            int C = min(m, maxCols);
            long long cur = countAffine(*F, false, R, C);
            long long final = cur + (m - C);
            if (final > best.cells) {
                best = {final, Candidate::AFFINE, false, false,
                        q, 0, R, C, F};
            }
        }
        /* ---------- projective (rows = points) ----- */
        {
            int total = q * q + q + 1;
            int R = min(n, total);
            int C = min(m, total);
            long long cur = countProjective(*F, R, C);
            long long final = cur + (m - C);
            if (final > best.cells) {
                best = {final, Candidate::PROJECTIVE, true, false,
                        q, 0, R, C, F};
            }
        }
        /* ---------- projective (transpose) ---------- */
        {
            int total = q * q + q + 1;
            int R = min(n, total);
            int C = min(m, total);
            long long cur = countProjective(*F, R, C);
            long long final = cur + (m - C);
            if (final > best.cells) {
                best = {final, Candidate::PROJECTIVE, true, true,
                        q, 0, R, C, F};
            }
        }

        /* ---------- polynomial design (p = ceil(sqrt(m))) ---------- */
        {
            // we only need a p such that p ≤ q and p·q ≥ m
            int p = (int)ceil(sqrt((double)m));
            if (p <= q) {
                int maxRows = p * q;
                int maxCols = p * q;
                int R = min(n, maxRows);
                int C = min(m, maxCols);
                long long cur = countPolynomial(*F, p, R, C);
                long long final = cur + (m - C);
                if (final > best.cells) {
                    best = {final, Candidate::POLYNOMIAL, false, false,
                            q, p, R, C, F};
                }
            }
        }
    }

    /* ---------- build the selected construction ------------------- */
    vector<pair<int,int>> answer;
    answer.reserve((size_t)best.cells + 10);

    // bitset per row, list of rows per column (for fast checking)
    int WORDS = (m + 63) >> 6;
    vector<vector<uint64_t>> bits(n, vector<uint64_t>(WORDS, 0));
    vector<vector<int>> colRows(m);

    auto setCell = [&](int r0, int c0) {
        bits[r0][c0 >> 6] |= 1ULL << (c0 & 63);
        colRows[c0].push_back(r0);
        answer.emplace_back(r0 + 1, c0 + 1);
    };

    if (best.kind == Candidate::STAR) {
        // whole first row
        for (int c = 0; c < m; ++c) setCell(0, c);
        // remaining rows only column 0
        for (int r = 1; r < n; ++r) setCell(r, 0);
    } else if (best.kind == Candidate::AFFINE) {
        emitAffine(*best.F, best.rowsArePoints, best.R, best.C, answer);
        for (auto &pr : answer) {
            int r = pr.first - 1, c = pr.second - 1;
            bits[r][c >> 6] |= 1ULL << (c & 63);
            colRows[c].push_back(r);
        }
        // add star row for all remaining columns
        for (int c = best.C; c < m; ++c) setCell(0, c);
    } else if (best.kind == Candidate::PROJECTIVE) {
        vector<pair<int,int>> temp;
        emitProjective(*best.F, best.R, best.C, temp);
        if (best.transpose)
            for (auto &pr : temp) swap(pr.first, pr.second);
        for (auto &pr : temp) {
            int r = pr.first - 1, c = pr.second - 1;
            if (r < n && c < m) {
                bits[r][c >> 6] |= 1ULL << (c & 63);
                colRows[c].push_back(r);
                answer.emplace_back(r + 1, c + 1);
            }
        }
        // star row (row 0) for the rest
        for (int c = best.C; c < m; ++c) setCell(0, c);
    } else { // POLYNOMIAL
        vector<pair<int,int>> temp;
        emitPolynomial(*best.F, best.p, best.R, best.C, temp);
        for (auto &pr : temp) {
            int r = pr.first - 1, c = pr.second - 1;
            if (r < n && c < m) {
                bits[r][c >> 6] |= 1ULL << (c & 63);
                colRows[c].push_back(r);
                answer.emplace_back(r + 1, c + 1);
            }
        }
        // star row for leftovers
        for (int c = best.C; c < m; ++c) setCell(0, c);
    }

    /* ---------- greedy augmentation (maximal C4‑free) ---------- */
    vector<pair<int,int>> empty;
    empty.reserve((size_t)n * m - answer.size());
    for (int r = 0; r < n; ++r)
        for (int c = 0; c < m; ++c)
            if ( (bits[r][c >> 6] & (1ULL << (c & 63))) == 0 )
                empty.emplace_back(r, c);

    mt19937 rng( (uint64_t)n * 1234567ULL + (uint64_t)m );
    shuffle(empty.begin(), empty.end(), rng);

    auto canAdd = [&](int r, int c) -> bool {
        const auto &rowsHere = colRows[c];
        for (int r2 : rowsHere) {
            const auto &A = bits[r];
            const auto &B = bits[r2];
            for (int w = 0; w < WORDS; ++w)
                if (A[w] & B[w]) return false;   // they already share a column
        }
        return true;
    };

    for (auto &pr : empty) {
        int r = pr.first, c = pr.second;
        if (canAdd(r, c)) setCell(r, c);
    }

    /* ---------- output ------------------------------------------- */
    cout << answer.size() << '\n';
    for (auto &p : answer)
        cout << p.first << ' ' << p.second << '\n';
    return 0;
}