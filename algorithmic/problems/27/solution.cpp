#include <bits/stdc++.h>
using namespace std;

/*-------------------------------------------------------------
   1.  Finite field  GF(q)   (prime power q ≤ 317)
   ------------------------------------------------------------*/
struct Field {
    int p = 0, e = 0, q = 0;               // q = p^e
    vector<vector<int>> add, mul;          // tables
    vector<int> neg, inv;                  // additive / multiplicative inverses
};

/* primality test */
static bool isPrime(int x) {
    if (x < 2) return false;
    for (int d = 2; d * d <= x; ++d) if (x % d == 0) return false;
    return true;
}

/* factor a prime power q = p^e (e ≥ 1) */
static pair<int,int> factorPrimePower(int q) {
    for (int p = 2; p * p <= q; ++p) if (q % p == 0 && isPrime(p)) {
        int e = 0, t = q;
        while (t % p == 0) { t /= p; ++e; }
        return {p, e};
    }
    return {q, 1};                     // q itself is prime
}

/* ----------  irreducible polynomial for extensions (e == 2) ---------- */
using Poly = vector<int>;

static void trim(Poly &a) {
    while (!a.empty() && a.back() == 0) a.pop_back();
}
static int deg(const Poly &a) { return (int)a.size() - 1; }

static Poly poly_mod(Poly a, const Poly &b, int p) {
    int db = deg(b);
    while (deg(a) >= db && db >= 0) {
        int shift = deg(a) - db;
        int coeff = a.back();
        for (int i = 0; i <= db; ++i) {
            int idx = i + shift;
            a[idx] = (a[idx] - coeff * b[i]) % p;
            if (a[idx] < 0) a[idx] += p;
        }
        trim(a);
    }
    return a;
}
static bool isIrreducible(const Poly &f, int p) {
    int d = deg(f);
    for (int k = 1; k <= d / 2; ++k) {
        long long tot = 1;
        for (int i = 0; i < k; ++i) tot *= p;
        for (long long mask = 0; mask < tot; ++mask) {
            Poly g(k + 1, 0);
            long long tmp = mask;
            for (int i = 0; i < k; ++i) { g[i] = int(tmp % p); tmp /= p; }
            g[k] = 1;                                   // monic
            if (poly_mod(f, g, p).empty()) return false;
        }
    }
    return true;
}
static Poly findIrredPoly(int p, int e) {          // only e==2 is needed here
    long long tot = 1;
    for (int i = 0; i < e; ++i) tot *= p;
    for (long long mask = 0; mask < tot; ++mask) {
        Poly f(e + 1, 0);
        long long tmp = mask;
        for (int i = 0; i < e; ++i) { f[i] = int(tmp % p); tmp /= p; }
        f[e] = 1;                                   // monic
        if (isIrreducible(f, p)) return f;
    }
    assert(false);
    return {};
}

/* ----------- build field tables for a given prime power q ----------- */
static Field buildField(int q) {
    Field F;  F.q = q;
    auto [p, e] = factorPrimePower(q);
    F.p = p;  F.e = e;
    F.add.assign(q, vector<int>(q));
    F.mul.assign(q, vector<int>(q));
    F.neg.assign(q, 0);
    F.inv.assign(q, 0);

    if (e == 1) {                                 // prime field
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

    /* ----------- extension field (e == 2) ----------------- */
    Poly irr = findIrredPoly(p, e);               // monic irreducible degree e
    vector<int> powp(e);
    powp[0] = 1;
    for (int i = 1; i < e; ++i) powp[i] = powp[i-1] * p;

    auto decode = [&](int x, array<int,8> &out) {
        for (int i = 0; i < e; ++i) {
            out[i] = x % p; x /= p;
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

    for (int a = 0; a < q; ++a)
        for (int b = 0; b < q; ++b)
            if (F.add[a][b] == 0) { F.neg[a] = b; break; }

    // multiplication – polynomial multiplication mod irr
    for (int a = 0; a < q; ++a)
        for (int b = 0; b < q; ++b) {
            array<int,8> A{}, B{};
            decode(a, A); decode(b, B);
            vector<int> prod(2 * e - 1, 0);
            for (int i = 0; i < e; ++i)
                for (int j = 0; j < e; ++j) {
                    prod[i + j] = (prod[i + j] + A[i] * B[j]) % p;
                }
            // reduction
            for (int d = 2 * e - 2; d >= e; --d) {
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

    F.inv[1] = 1;
    for (int a = 2; a < q; ++a)
        for (int b = 1; b < q; ++b)
            if (F.mul[a][b] == 1) { F.inv[a] = b; break; }

    return F;
}

/*-------------------------------------------------------------
   2.   Counting / emitting helpers for the families
   -------------------------------------------------------------*/

static long long countAffineOptim(const Field &F, bool rowsArePoints,
                                  int R, int C) {
    int q = F.q;
    int maxCols = q * q + q;                     // total columns in affine plane
    if (rowsArePoints) {
        vector<int> cnt(maxCols, 0);
        for (int r = 0; r < R; ++r) {
            int x = r / q, y = r % q;
            for (int a = 0; a < q; ++a) {
                int ax = F.mul[a][x];
                int b  = F.add[y][F.neg[ax]];
                int col = a * q + b;                     // non‑vertical line
                ++cnt[col];
            }
            ++cnt[q * q + x];           // vertical line x = const
        }
        nth_element(cnt.begin(), cnt.begin() + C, cnt.end(),
                    greater<int>());
        long long sum = 0;
        for (int i = 0; i < C; ++i) sum += cnt[i];
        return sum;
    } else {
        // rows are lines, columns are points
        int maxPoints = q * q;
        vector<int> cnt(maxPoints, 0);
        for (int r = 0; r < R; ++r) {
            if (r < q * q) {                         // ordinary line
                int a = r / q, b = r % q;
                for (int x = 0; x < q; ++x) {
                    int ax = F.mul[a][x];
                    int y  = F.add[ax][b];
                    int pt = x * q + y;
                    ++cnt[pt];
                }
            } else {                                 // vertical line
                int x = r - q * q;
                for (int y = 0; y < q; ++y) {
                    int pt = x * q + y;
                    ++cnt[pt];
                }
            }
        }
        nth_element(cnt.begin(), cnt.begin() + C, cnt.end(),
                    greater<int>());
        long long sum = 0;
        for (int i = 0; i < C; ++i) sum += cnt[i];
        return sum;
    }
}

/* emission (the same selection of columns) */
static void emitAffineOptim(const Field &F, bool rowsArePoints,
                            int R, int C,
                            vector<pair<int,int>> &out) {
    int q = F.q;
    int maxCols = q * q + q;
    if (rowsArePoints) {
        vector<int> cnt(maxCols, 0);
        for (int r = 0; r < R; ++r) {
            int x = r / q, y = r % q;
            for (int a = 0; a < q; ++a) {
                int ax = F.mul[a][x];
                int b  = F.add[y][F.neg[ax]];
                int col = a * q + b;
                ++cnt[col];
            }
            ++cnt[q * q + x];
        }
        vector<int> idx(maxCols);
        iota(idx.begin(), idx.end(), 0);
        nth_element(idx.begin(), idx.begin() + C, idx.end(),
                    [&](int i, int j){ return cnt[i] > cnt[j]; });
        vector<char> keep(maxCols, 0);
        for (int i = 0; i < C; ++i) keep[idx[i]] = 1;

        for (int r = 0; r < R; ++r) {
            int x = r / q, y = r % q;
            for (int a = 0; a < q; ++a) {
                int ax = F.mul[a][x];
                int b  = F.add[y][F.neg[ax]];
                int col = a * q + b;
                if (keep[col])
                    out.emplace_back(r + 1, col + 1);
            }
            int col = q * q + x;
            if (keep[col]) out.emplace_back(r + 1, col + 1);
        }
    } else {
        // rows = lines, columns = points
        int maxPts = q * q;
        vector<int> cnt(maxPts, 0);
        for (int r = 0; r < R; ++r) {
            if (r < q * q) {
                int a = r / q, b = r % q;
                for (int x = 0; x < q; ++x) {
                    int ax = F.mul[a][x];
                    int y  = F.add[ax][b];
                    int pt = x * q + y;
                    ++cnt[pt];
                }
            } else {
                int x = r - q * q;
                for (int y = 0; y < q; ++y) {
                    int pt = x * q + y;
                    ++cnt[pt];
                }
            }
        }
        vector<int> idx(maxPts);
        iota(idx.begin(), idx.end(), 0);
        nth_element(idx.begin(), idx.begin() + C, idx.end(),
                    [&](int i, int j){ return cnt[i] > cnt[j]; });
        vector<char> keep(maxPts, 0);
        for (int i = 0; i < C; ++i) keep[idx[i]] = 1;

        for (int r = 0; r < R; ++r) {
            if (r < q * q) {
                int a = r / q, b = r % q;
                for (int x = 0; x < q; ++x) {
                    int ax = F.mul[a][x];
                    int y  = F.add[ax][b];
                    int pt = x * q + y;
                    if (keep[pt])
                        out.emplace_back(r + 1, pt + 1);
                }
            } else {
                int x = r - q * q;
                for (int y = 0; y < q; ++y) {
                    int pt = x * q + y;
                    if (keep[pt])
                        out.emplace_back(r + 1, pt + 1);
                }
            }
        }
    }
}

/* projective plane */
static long long countProjective(const Field &F, int R, int C) {
    int q = F.q;
    long long cnt = 0;
    for (int r = 0; r < R; ++r) {
        if (r < q * q) {                         // ordinary point (x,y)
            int x = r / q, y = r % q;
            for (int a = 0; a < q; ++a) {
                int ax = F.mul[a][x];
                int b  = F.add[y][F.neg[ax]];
                int col = a * q + b;                     // non‑vertical line
                if (col < C) ++cnt;
            }
            int col = q * q + x;                         // vertical line
            if (col < C) ++cnt;
        } else {                                  // point at infinity
            int d = r - q * q;
            if (d < q) {                     // direction ∞
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
static void emitProjective(const Field &F, int R, int C,
                           vector<pair<int,int>> &out) {
    int q = F.q;
    for (int r = 0; r < R; ++r) {
        if (r < q * q) {                     // ordinary point (x,y)
            int x = r / q, y = r % q;
            for (int a = 0; a < q; ++a) {
                int ax = F.mul[a][x];
                int b  = F.add[y][F.neg[ax]];
                int col = a * q + b;
                if (col < C) out.emplace_back(r + 1, col + 1);
            }
            int col = q * q + x;
            if (col < C) out.emplace_back(r + 1, col + 1);
        } else {                              // point at infinity
            int d = r - q * q;
            if (d < q) {
                for (int b = 0; b < q; ++b) {
                    int col = d * q + b;
                    if (col < C) out.emplace_back(r + 1, col + 1);
                }
            } else {
                for (int c = 0; c < q; ++c) {
                    int col = q * q + c;
                    if (col < C) out.emplace_back(r + 1, col + 1);
                }
            }
            int col = q * q + q;
            if (col < C) out.emplace_back(r + 1, col + 1);
        }
    }
}

/* Reed–Solomon polynomial design */
static long long countPolynomial(const Field &F, int p,
                                 int R, int C) {
    long long cnt = 0;
    int q = F.q;
    for (int i = 0; i < p; ++i) {
        for (int x = 0; x < q; ++x) {
            int r = i * q + x;
            if (r >= R) continue;
            for (int j = 0; j < p; ++j) {
                int y = F.add[F.mul[j][x]][i];
                int c = j * q + y;
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
            int r = i * q + x;
            if (r >= R) continue;
            for (int j = 0; j < p; ++j) {
                int y = F.add[F.mul[j][x]][i];
                int c = j * q + y;
                if (c < C) out.emplace_back(r + 1, c + 1);
            }
        }
    }
}

/*-------------------------------------------------------------
   3.   Candidate description
   -------------------------------------------------------------*/
struct Candidate {
    long long base = -1;                // edges before augmentation
    enum Kind { NONE, AFFINE, PROJECTIVE, POLYNOMIAL, PAIRWISE } kind = NONE;
    bool rowsArePoints = true;          // only for affine (orientation)
    bool transpose = false;             // apply a final (row ↔ col) swap
    int q = 0;                          // field order, 0 for pairwise
    int p = 0;                          // extra param (degree for polynomial,
                                        //     pair count for pairwise)
    int R = 0, C = 0;                   // rows / columns actually used
    const Field *F = nullptr;           // pointer to field (null for pairwise)
};

/*-------------------------------------------------------------
   4.   Greedy augmentation – maximal C4‑free extension
   -------------------------------------------------------------*/
static void greedyAugment(const Candidate &cand,
                          int n, int m,
                          vector<pair<int,int>> &answer,
                          mt19937 &rng) {
    const int WORDS = (m + 63) >> 6;
    vector<vector<uint64_t>> bits(n, vector<uint64_t>(WORDS, 0));
    vector<vector<int>> colRows(m);
    answer.clear();

    /* ---- insert the core design --------------------------------- */
    vector<pair<int,int>> core;
    if (cand.kind != Candidate::NONE) {
        if (cand.kind == Candidate::AFFINE) {
            emitAffineOptim(*cand.F, cand.rowsArePoints,
                            cand.R, cand.C, core);
        } else if (cand.kind == Candidate::PROJECTIVE) {
            emitProjective(*cand.F, cand.R, cand.C, core);
        } else if (cand.kind == Candidate::POLYNOMIAL) {
            emitPolynomial(*cand.F, cand.p, cand.R, cand.C, core);
        } else { // PAIRWISE
            long long pairCnt = cand.p;
            for (int i = 0; i < n && pairCnt > 0; ++i) {
                for (int j = i + 1; j < n && pairCnt > 0; ++j) {
                    int col = int(cand.p - pairCnt);
                    core.emplace_back(i + 1, col + 1);
                    core.emplace_back(j + 1, col + 1);
                    --pairCnt;
                }
            }
            for (int col = cand.p; col < m; ++col)
                core.emplace_back(1, col + 1);
        }
        if (cand.transpose) for (auto &pr : core) swap(pr.first, pr.second);
        for (auto &pr : core) {
            int r = pr.first - 1, c = pr.second - 1;
            if (r >= n || c >= m) continue;
            bits[r][c >> 6] |= 1ULL << (c & 63);
            colRows[c].push_back(r);
            answer.emplace_back(r + 1, c + 1);
        }
    }

    /* ---- diagonal of the untouched suffix (always safe) -------- */
    int rowPtr = cand.R, colPtr = cand.C;
    while (rowPtr < n && colPtr < m) {
        bits[rowPtr][colPtr >> 6] |= 1ULL << (colPtr & 63);
        colRows[colPtr].push_back(rowPtr);
        answer.emplace_back(rowPtr + 1, colPtr + 1);
        ++rowPtr; ++colPtr;
    }

    /* ---- random greedy fill -------------------------------------- */
    vector<pair<int,int>> empty;
    empty.reserve((size_t)n * m - answer.size());
    for (int r = 0; r < n; ++r)
        for (int c = 0; c < m; ++c)
            if ((bits[r][c >> 6] & (1ULL << (c & 63))) == 0)
                empty.emplace_back(r, c);
    shuffle(empty.begin(), empty.end(), rng);

    auto canAdd = [&](int r, int c) -> bool {
        const auto &rowsHere = colRows[c];
        for (int r2 : rowsHere) {
            const auto &A = bits[r];
            const auto &B = bits[r2];
            for (int w = 0; w < WORDS; ++w)
                if (A[w] & B[w]) return false;   // second common neighbour
        }
        return true;
    };

    for (auto &rc : empty) {
        int r = rc.first, c = rc.second;
        if (canAdd(r, c)) {
            bits[r][c >> 6] |= 1ULL << (c & 63);
            colRows[c].push_back(r);
            answer.emplace_back(r + 1, c + 1);
        }
    }
}

/*-------------------------------------------------------------
   5.   Main
   -------------------------------------------------------------*/
int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n, m;
    if (!(cin >> n >> m)) return 0;

    /* Degenerate grids have closed-form optima: with <=2 rows (or columns) no K2,2 can use
       both lines twice, so a full line plus one extra cell is optimal. Emit directly with
       buffered output — the generic pipeline below spends ~0.8s here for the same k. */
    if (n <= 2 || m <= 2) {
        string out;
        out.reserve(16u * (size_t)max(n, m) + 32);
        auto emit = [&out](int r, int c) {
            out += to_string(r); out += ' '; out += to_string(c); out += '\n';
        };
        if (n == 1 || m == 1) {
            out += to_string((long long)n * m); out += '\n';
            if (n == 1) for (int c = 1; c <= m; ++c) emit(1, c);
            else       for (int r = 1; r <= n; ++r) emit(r, 1);
        } else if (n == 2) {                 // full first row + one cell in the second
            out += to_string(m + 1); out += '\n';
            for (int c = 1; c <= m; ++c) emit(1, c);
            emit(2, 1);
        } else {                              // m == 2: full first column + one extra
            out += to_string(n + 1); out += '\n';
            for (int r = 1; r <= n; ++r) emit(r, 1);
            emit(1, 2);
        }
        cout << out;
        return 0;
    }

    /* ---- list of prime powers ≤ 317 (√100 000) ----------------- */
    const int LIM = 317;
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

    /* ---- lazy field cache ---------------------------------------- */
    unordered_map<int, Field> fieldCache;
    auto getField = [&](int q) -> const Field& {
        auto it = fieldCache.find(q);
        if (it == fieldCache.end()) {
            Field F = buildField(q);
            auto ins = fieldCache.emplace(q, std::move(F));
            return ins.first->second;
        }
        return it->second;
    };

    /* ---- generate candidates -------------------------------------- */
    vector<Candidate> candList;
    // trivial (empty) candidate – keeps the list non‑empty
    candList.push_back({0, Candidate::NONE, true, false, 0, 0, 0, 0, nullptr});

    for (int q : primePowers) {
        const Field &F = getField(q);
        /* ----- AFFINE (rows = points) ----- */
        {
            int maxRows = q * q;
            int maxCols = q * q + q;
            int R = min(n, maxRows);
            int C = min(m, maxCols);
            long long base = countAffineOptim(F, true, R, C);
            candList.push_back({base, Candidate::AFFINE, true, false,
                                q, 0, R, C, &F});
        }
        /* ----- AFFINE (rows = lines) ----- */
        {
            int maxRows = q * q + q;
            int maxCols = q * q;
            int R = min(n, maxRows);
            int C = min(m, maxCols);
            long long base = countAffineOptim(F, false, R, C);
            candList.push_back({base, Candidate::AFFINE, false, false,
                                q, 0, R, C, &F});
        }
        /* ----- PROJECTIVE (both orientations) ----- */
        {
            int tot = q * q + q + 1;
            int R = min(n, tot);
            int C = min(m, tot);
            long long base = countProjective(F, R, C);
            candList.push_back({base, Candidate::PROJECTIVE, false, false,
                                q, 0, R, C, &F});
            Candidate t = candList.back();
            t.transpose = true;
            candList.push_back(t);
        }
        /* ----- REED–SOLOMON (ordinary) ----- */
        {
            int p_needed = (int)ceil(sqrt((double)m));
            if (p_needed <= q) {
                int maxRows = p_needed * q;
                int maxCols = p_needed * q;
                int R = min(n, maxRows);
                int C = min(m, maxCols);
                long long base = countPolynomial(F, p_needed, R, C);
                candList.push_back({base, Candidate::POLYNOMIAL, false, false,
                                    q, p_needed, R, C, &F});
            }
        }
        /* ----- REED–SOLOMON (transposed) ----- */
        {
            int p_needed = (int)ceil(sqrt((double)n));
            if (p_needed <= q) {
                int maxRows = p_needed * q;
                int maxCols = p_needed * q;
                // rows ↔ columns
                int R = min(m, maxRows);
                int C = min(n, maxCols);
                long long base = countPolynomial(F, p_needed, R, C);
                candList.push_back({base, Candidate::POLYNOMIAL, false, true,
                                    q, p_needed, R, C, &F});
            }
        }
        /* ----- PAIRWISE (normal) ----- */
        {
            long long pairCnt = min<long long>(m, 1LL * n * (n - 1) / 2);
            long long edges = m + pairCnt;                 // 2*pairCnt + (m-pairCnt)
            candList.push_back({edges, Candidate::PAIRWISE, false, false,
                                0, (int)pairCnt, n, m, nullptr});
        }
        /* ----- PAIRWISE (transposed) ----- */
        {
            long long pairCnt = min<long long>(n, 1LL * m * (m - 1) / 2);
            long long edges = n + pairCnt;
            candList.push_back({edges, Candidate::PAIRWISE, false, true,
                                0, (int)pairCnt, m, n, nullptr});
        }
    }

    /* ---- keep only the most promising candidates ------------------- */
    const int KEEP = 30;
    sort(candList.begin(), candList.end(),
         [](const Candidate &a, const Candidate &b){ return a.base > b.base; });
    if ((int)candList.size() > KEEP) candList.resize(KEEP);

    /* ---- greedy augmentation -------------------------------------- */
    vector<pair<int,int>> bestAns;
    size_t bestSize = 0;
    mt19937 rng((uint64_t)chrono::steady_clock::now().time_since_epoch().count());

    const int ATTEMPTS = 4;
    vector<pair<int,int>> curAns;
    for (const Candidate &c : candList) {
        for (int it = 0; it < ATTEMPTS; ++it) {
            greedyAugment(c, n, m, curAns, rng);
            if (curAns.size() > bestSize) {
                bestSize = curAns.size();
                bestAns.swap(curAns);
            }
        }
    }

    /* ---- output --------------------------------------------------- */
    cout << bestAns.size() << '\n';
    for (auto &pr : bestAns) cout << pr.first << ' ' << pr.second << '\n';
    return 0;
}