// Traveling Santa with carrot penalty -- minimal, maximally-portable variant.
// Deliberately avoids: k-d tree recursion, <random>, <chrono>, heavy nested
// lambdas, malloc/realloc, and <string>. Uses only C++11-era constructs and a
// clock()-based CPU-time budget (the judge caps CPU time), so it compiles on any
// g++ >= 4.9 and cannot TLE (budget << limit). Always prints a valid tour.
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <vector>
#include <algorithm>
using namespace std;

typedef long long ll;
typedef unsigned int u32;
typedef unsigned long long u64;

static inline double cpuTime() { return (double)clock() / (double)CLOCKS_PER_SEC; }
static const double TLIM = 0.85; // CPU seconds; judge cpuLimit is >=1s, so huge margin

static int N;
static vector<ll> X, Y;
static vector<char> isP;

static inline double D(int a, int b) {
    double dx = (double)X[a] - (double)X[b];
    double dy = (double)Y[a] - (double)Y[b];
    return sqrt(dx * dx + dy * dy);
}

// ---- Hilbert order (n = 2^31) ----
static const int HB = 31;
static inline u64 hilb(u32 x, u32 y) {
    u64 d = 0;
    u32 n = 1u << HB;
    for (u32 s = n >> 1; s > 0; s >>= 1) {
        u32 rx = (x & s) ? 1u : 0u;
        u32 ry = (y & s) ? 1u : 0u;
        d += (u64)s * (u64)s * (u64)((3u * rx) ^ ry);
        if (ry == 0) {
            if (rx == 1) { x = (n - 1) - x; y = (n - 1) - y; }
            u32 t = x; x = y; y = t;
        }
    }
    return d;
}

static double penalized(const vector<int>& t) {
    double s = 0.0;
    for (int k = 1; k <= N; ++k) {
        int a = t[k - 1], b = t[k];
        double m = (k % 10 == 0 && !isP[a]) ? 1.1 : 1.0;
        s += m * D(a, b);
    }
    return s;
}
static double euclid(const vector<int>& t) {
    double s = 0.0;
    for (int k = 1; k <= N; ++k) s += D(t[k - 1], t[k]);
    return s;
}

// ---- simple buffered stdin reader (no malloc/realloc) ----
static vector<char> g_in;
static size_t g_pos;
static inline ll readLL() {
    size_t n = g_in.size();
    while (g_pos < n) { char c = g_in[g_pos]; if ((c >= '0' && c <= '9') || c == '-') break; g_pos++; }
    int sgn = 1;
    if (g_pos < n && g_in[g_pos] == '-') { sgn = -1; g_pos++; }
    ll v = 0;
    while (g_pos < n && g_in[g_pos] >= '0' && g_in[g_pos] <= '9') { v = v * 10 + (g_in[g_pos] - '0'); g_pos++; }
    return v * sgn;
}

int main() {
    // read all of stdin
    {
        char buf[1 << 16];
        size_t r;
        while ((r = fread(buf, 1, sizeof(buf), stdin)) > 0) g_in.insert(g_in.end(), buf, buf + r);
        g_pos = 0;
    }
    N = (int)readLL();
    if (N < 2) { printf("2\n0\n0\n"); return 0; }
    X.resize(N); Y.resize(N);
    for (int i = 0; i < N; i++) { X[i] = readLL(); Y[i] = readLL(); }

    // primes over city IDs (0 and 1 are not prime)
    {
        int M = N < 2 ? 2 : N;
        isP.assign(M, 1);
        isP[0] = 0; if (M > 1) isP[1] = 0;
        for (ll i = 2; i * i < (ll)M; ++i) if (isP[i])
            for (ll j = i * i; j < (ll)M; j += i) isP[j] = 0;
    }

    // ---- Hilbert-sorted order ----
    vector<pair<u64, int> > ord(N);
    {
        const ll OFF = 1000000000LL; // |x|,|y| <= 1e9 -> [0,2e9] < 2^31
        for (int i = 0; i < N; i++)
            ord[i] = make_pair(hilb((u32)(X[i] + OFF), (u32)(Y[i] + OFF)), i);
        sort(ord.begin(), ord.end());
    }

    // ---- seed tour: break the Hilbert cycle at city 0 ----
    vector<int> tour(N + 1);
    {
        int p = 0;
        for (int i = 0; i < N; i++) if (ord[i].second == 0) { p = i; break; }
        tour[0] = 0;
        for (int m = 1; m < N; ++m) tour[m] = ord[(p + m) % N].second;
        tour[N] = 0;
    }

    // ---- candidate lists from Hilbert-order neighbours (no k-d tree) ----
    int half = 8;
    if (N <= 2000) half = 20; else if (N <= 50000) half = 12;
    vector<int> candFlat;           // flattened neighbour lists
    vector<int> candBeg(N + 1, 0);   // candBeg[c]..candBeg[c+1] range in candFlat
    {
        vector<int> cnt(N, 0);
        for (int s = 0; s < N; ++s) {
            int lo = s - half < 0 ? 0 : s - half;
            int hi = s + half >= N ? N - 1 : s + half;
            cnt[ord[s].second] += (hi - lo); // excludes self
        }
        for (int c = 0; c < N; c++) candBeg[c + 1] = candBeg[c] + cnt[c];
        candFlat.assign(candBeg[N], 0);
        vector<int> cur(candBeg.begin(), candBeg.end() - 1);
        for (int s = 0; s < N; ++s) {
            int city = ord[s].second;
            int lo = s - half < 0 ? 0 : s - half;
            int hi = s + half >= N ? N - 1 : s + half;
            for (int t = lo; t <= hi; ++t) if (t != s) candFlat[cur[city]++] = ord[t].second;
        }
    }

    // position map (city -> index in tour); city 0 maps to index 0
    vector<int> pos(N + 1, 0);
    for (int i = N; i >= 0; --i) pos[tour[i]] = i;

    int MAXREV;
    if (N <= 2000) MAXREV = N + 5;
    else if (N <= 100000) MAXREV = 1000;
    else MAXREV = 400;
    const double eps = 1e-3;

    // active queue with in-queue flags (don't-look bits)
    vector<int> Q; Q.reserve(N + 2);
    vector<char> inq(N + 2, 0);
    for (int i = 1; i < N; i++) { inq[i] = 1; Q.push_back(i); }

    double cost = euclid(tour);
    int cnt = 0;
    while (!Q.empty()) {
        if (++cnt >= 512) { cnt = 0; if (cpuTime() > TLIM) break; }
        int i = Q.back(); Q.pop_back(); inq[i] = 0;
        int a = tour[i - 1], b = tour[i];
        double dab = D(a, b);
        bool imp = false;

        // 2-opt: reconnect a to a near city v (candidates of a)
        for (int ci = candBeg[a]; ci < candBeg[a + 1]; ++ci) {
            int v = candFlat[ci];
            int pv = pos[v];
            if (pv == i || pv == i - 1) continue;
            if (pv > i) {
                int len = pv - i + 1;
                if (len > MAXREV) continue;
                int d = tour[pv + 1];
                double delta = D(a, v) + D(b, d) - dab - D(v, d);
                if (delta < -eps) {
                    reverse(tour.begin() + i, tour.begin() + pv + 1);
                    for (int t = i; t <= pv; t++) pos[tour[t]] = t;
                    cost += delta;
                    int z;
                    z = i - 1; if (z >= 1 && z < N && !inq[z]) { inq[z] = 1; Q.push_back(z); }
                    z = i;     if (z >= 1 && z < N && !inq[z]) { inq[z] = 1; Q.push_back(z); }
                    z = pv;    if (z >= 1 && z < N && !inq[z]) { inq[z] = 1; Q.push_back(z); }
                    z = pv + 1;if (z >= 1 && z < N && !inq[z]) { inq[z] = 1; Q.push_back(z); }
                    imp = true; break;
                }
            } else { // pv < i-1
                int len = i - pv - 1;
                if (len > MAXREV) continue;
                int d = tour[pv + 1];
                double delta = D(a, v) + D(d, b) - D(v, d) - dab;
                if (delta < -eps) {
                    reverse(tour.begin() + pv + 1, tour.begin() + i);
                    for (int t = pv + 1; t < i; t++) pos[tour[t]] = t;
                    cost += delta;
                    int z;
                    z = pv;    if (z >= 1 && z < N && !inq[z]) { inq[z] = 1; Q.push_back(z); }
                    z = pv + 1;if (z >= 1 && z < N && !inq[z]) { inq[z] = 1; Q.push_back(z); }
                    z = i - 1; if (z >= 1 && z < N && !inq[z]) { inq[z] = 1; Q.push_back(z); }
                    z = i;     if (z >= 1 && z < N && !inq[z]) { inq[z] = 1; Q.push_back(z); }
                    imp = true; break;
                }
            }
        }

        // Or-opt: relocate segment [i..i+L-1] (L=1,2,3) next to a near city
        if (!imp) {
            for (int L = 1; L <= 3 && !imp; ++L) {
                if (i + L - 1 > N - 1) break;
                int segA = tour[i], segB = tour[i + L - 1];
                int pre = tour[i - 1], post = tour[i + L];
                double save = D(pre, segA) + D(segB, post) - D(pre, post);
                for (int ci = candBeg[segA]; ci < candBeg[segA + 1]; ++ci) {
                    int v = candFlat[ci];
                    int j = pos[v];
                    if (j >= i + L) {
                        if (j - i + 1 > MAXREV) continue;
                        int ja = tour[j], jb = tour[j + 1];
                        double delta = -save + D(ja, segA) + D(segB, jb) - D(ja, jb);
                        if (delta < -eps) {
                            reverse(tour.begin() + i, tour.begin() + i + L);
                            reverse(tour.begin() + i + L, tour.begin() + j + 1);
                            reverse(tour.begin() + i, tour.begin() + j + 1);
                            for (int t = i; t <= j; t++) pos[tour[t]] = t;
                            cost += delta;
                            int z;
                            z = i - 1; if (z >= 1 && z < N && !inq[z]) { inq[z] = 1; Q.push_back(z); }
                            z = i;     if (z >= 1 && z < N && !inq[z]) { inq[z] = 1; Q.push_back(z); }
                            z = j;     if (z >= 1 && z < N && !inq[z]) { inq[z] = 1; Q.push_back(z); }
                            z = j + 1; if (z >= 1 && z < N && !inq[z]) { inq[z] = 1; Q.push_back(z); }
                            imp = true; break;
                        }
                    } else if (j <= i - 2) {
                        if (j + L > N - 1) continue;
                        if (i + L - 1 - j > MAXREV) continue;
                        int ja = tour[j], jb = tour[j + 1];
                        double delta = -save + D(ja, segA) + D(segB, jb) - D(ja, jb);
                        if (delta < -eps) {
                            reverse(tour.begin() + j + 1, tour.begin() + i + L);
                            reverse(tour.begin() + j + 1, tour.begin() + j + 1 + L);
                            reverse(tour.begin() + j + 1 + L, tour.begin() + i + L);
                            for (int t = j + 1; t <= i + L - 1; t++) pos[tour[t]] = t;
                            cost += delta;
                            int z;
                            z = j;       if (z >= 1 && z < N && !inq[z]) { inq[z] = 1; Q.push_back(z); }
                            z = j + 1;   if (z >= 1 && z < N && !inq[z]) { inq[z] = 1; Q.push_back(z); }
                            z = i + L - 1; if (z >= 1 && z < N && !inq[z]) { inq[z] = 1; Q.push_back(z); }
                            z = i + L;   if (z >= 1 && z < N && !inq[z]) { inq[z] = 1; Q.push_back(z); }
                            imp = true; break;
                        }
                    }
                }
            }
        }

        if (imp) { if (i >= 1 && i < N && !inq[i]) { inq[i] = 1; Q.push_back(i); } }
    }

    // ---- validity check + monotonic fallback (guarantees a valid, non-worse tour) ----
    bool ok = ((int)tour.size() == N + 1) && tour[0] == 0 && tour[N] == 0;
    if (ok) {
        vector<char> seen(N, 0);
        for (int i = 1; i < N && ok; i++) {
            int c = tour[i];
            if (c < 0 || c >= N || seen[c]) ok = false; else seen[c] = 1;
        }
        for (int c = 1; c < N && ok; c++) if (!seen[c]) ok = false;
    }
    vector<int> mono(N + 1);
    for (int i = 0; i < N; i++) mono[i] = i; mono[N] = 0;

    bool useTour = ok;
    if (useTour) {
        double Lt = penalized(tour), Lm = penalized(mono);
        if (Lt > Lm) useTour = false;
    }
    const vector<int>& out = useTour ? tour : mono;

    // ---- fast output (manual int formatting, no <string>) ----
    {
        vector<char> ob;
        ob.reserve((size_t)(N + 2) * 8);
        char tmp[16];
        // first line: N+1
        {
            long val = (long)N + 1; int p = 0;
            if (val == 0) tmp[p++] = '0';
            while (val > 0) { tmp[p++] = char('0' + val % 10); val /= 10; }
            while (p > 0) ob.push_back(tmp[--p]);
            ob.push_back('\n');
        }
        for (int i = 0; i <= N; i++) {
            int val = out[i]; int p = 0;
            if (val == 0) tmp[p++] = '0';
            while (val > 0) { tmp[p++] = char('0' + val % 10); val /= 10; }
            while (p > 0) ob.push_back(tmp[--p]);
            ob.push_back('\n');
        }
        fwrite(ob.data(), 1, ob.size(), stdout);
    }
    return 0;
}
