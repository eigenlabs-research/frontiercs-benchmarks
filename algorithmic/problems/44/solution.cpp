#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include <ctime>
#include <cstring>
#include <string>
#include <cstdlib>
#include <random>
using namespace std;
using ll = long long;
using u32 = uint32_t;
using u64 = uint64_t;

static clock_t T0;
static inline double elapsed(){
    return (double)(clock() - T0) / (double)CLOCKS_PER_SEC;
}
// 2-opt / ILS budget; leave room for penalty pass + I/O under the 2.5s limit.
// (Wall = TL_PEN + output + process overhead. Kept conservative for safety.)
static double TL_MAIN = 1.55;
static double TL_PEN  = 1.78;

// ---------------- fast input ----------------
static char* inbuf; static size_t inpos, inlen;
static void initIn(){
    size_t cap = 1u<<20; inbuf = (char*)malloc(cap); inlen = 0; size_t r;
    while((r = fread(inbuf+inlen, 1, cap-inlen, stdin)) > 0){
        inlen += r;
        if(inlen == cap){ cap <<= 1; inbuf = (char*)realloc(inbuf, cap); }
    }
    inpos = 0;
}
static inline ll readLL(){
    while(inpos < inlen){ char c = inbuf[inpos]; if((c>='0'&&c<='9')||c=='-') break; inpos++; }
    int s = 1; if(inpos < inlen && inbuf[inpos]=='-'){ s = -1; inpos++; }
    ll v = 0;
    while(inpos < inlen && inbuf[inpos]>='0' && inbuf[inpos]<='9'){ v = v*10 + (inbuf[inpos]-'0'); inpos++; }
    return v * s;
}

int N;
vector<ll> X, Y;
vector<char> isP;

static inline double D(int a, int b){
    double dx = (double)X[a] - (double)X[b];
    double dy = (double)Y[a] - (double)Y[b];
    return sqrt(dx*dx + dy*dy);
}

// ---------------- Hilbert order (n = 2^31) ----------------
static const int HB = 31;
static inline u64 hilb(u32 x, u32 y){
    u64 d = 0;
    u32 n = 1u << HB;            // 2^31 = 2147483648
    for(u32 s = n>>1; s>0; s>>=1){
        u32 rx = (x & s) ? 1u : 0u;
        u32 ry = (y & s) ? 1u : 0u;
        d += (u64)s * (u64)s * ((3u*rx) ^ ry);
        if(ry == 0){
            if(rx == 1){ x = (n-1) - x; y = (n-1) - y; }
            u32 t = x; x = y; y = t;
        }
    }
    return d;
}

static double penalized(const vector<int>& t){
    double s = 0.0;
    for(int k=1; k<=N; ++k){
        int a = t[k-1], b = t[k];
        double m = (k%10==0 && !isP[a]) ? 1.1 : 1.0;
        s += m * D(a,b);
    }
    return s;
}
static double euclidCost(const vector<int>& t){
    double s = 0.0;
    for(int k=1; k<=N; ++k) s += D(t[k-1], t[k]);
    return s;
}

// ---------------- k-d tree (2D) for k-NN candidate lists ----------------
static vector<int> g_kd;     // point IDs reordered into an implicit k-d tree
static int g_KC;

static void kdBuild(int l, int r, int dep){
    if(r - l <= 1) return;
    int axis = dep & 1, mid = (l + r) >> 1;
    if(axis)
        nth_element(g_kd.begin()+l, g_kd.begin()+mid, g_kd.begin()+r,
                    [](int a, int b){ return Y[a] < Y[b]; });
    else
        nth_element(g_kd.begin()+l, g_kd.begin()+mid, g_kd.begin()+r,
                    [](int a, int b){ return X[a] < X[b]; });
    kdBuild(l, mid, dep+1);
    kdBuild(mid+1, r, dep+1);
}

struct NNHeap {
    int K, n; vector<double> hd; vector<int> hi;
    NNHeap(int k): K(k), n(0), hd(k), hi(k) {}
    void reset(){ n = 0; }
    double maxD() const { return n < K ? -1.0 : hd[0]; }
    void push(double d2, int id){
        if(n < K){
            hd[n] = d2; hi[n] = id; ++n;
            int i = n - 1;
            while(i > 0){ int p = (i-1)>>1; if(hd[p] < hd[i]){ swap(hd[p],hd[i]); swap(hi[p],hi[i]); i = p; } else break; }
        } else if(d2 < hd[0]){
            hd[0] = d2; hi[0] = id;
            int i = 0;
            while(2*i+1 < n){ int l = 2*i+1, r = l+1, b = l; if(r<n && hd[r]>hd[l]) b = r; if(hd[b]>hd[i]){ swap(hd[b],hd[i]); swap(hi[b],hi[i]); i = b; } else break; }
        }
    }
};

static void kdQuery(int l, int r, int dep, int p, NNHeap& h){
    if(l >= r) return;
    int axis = dep & 1, mid = (l + r) >> 1;
    int id = g_kd[mid];
    ll dx = X[p] - X[id], dy = Y[p] - Y[id];
    if(id != p) h.push((double)dx*dx + (double)dy*dy, id);
    ll diff = axis ? (Y[p] - Y[id]) : (X[p] - X[id]);
    int nearL, nearR, farL, farR;
    if(diff < 0){ nearL = l; nearR = mid; farL = mid+1; farR = r; }
    else        { nearL = mid+1; nearR = r; farL = l; farR = mid; }
    kdQuery(nearL, nearR, dep+1, p, h);
    if((double)diff*diff < h.maxD() || h.n < g_KC)
        kdQuery(farL, farR, dep+1, p, h);
}

// ---------------- exact Held-Karp for tiny N (penalized objective) ----------------
// dp[S][i] = min penalized cost of a path 0 -> ... -> city(i+1) visiting exactly
// the cities of S (bit i <-> city i+1). Steps used = popcount(S); the multiplier of
// the next transition is known because its step index is popcount(S)+1.
static bool heldKarp(vector<int>& tour){
    int M = N - 1;
    if(M > 18) return false;                    // 2^18 * 18^2 ~ 21M transitions
    int FULL = 1 << M;
    // distance matrix
    double dm[20][20];
    for(int i=0;i<N;i++) for(int j=0;j<N;j++) dm[i][j] = D(i,j);
    vector<double> dp((size_t)FULL * M, 1e30);
    vector<int8_t> par((size_t)FULL * M, -1);
    for(int i=0;i<M;i++) dp[((size_t)(1u<<i))*M + i] = dm[0][i+1];   // step 1, never penalized
    for(int S=1; S<FULL; ++S){
        int pc = __builtin_popcount((unsigned)S);
        if(pc == M) break;
        int t = pc + 1;                          // step index of the next transition
        bool penStep = (t % 10 == 0);
        const size_t rowS = (size_t)S * M;
        int rem = (~S) & (FULL - 1);
        for(int i=0;i<M;i++){
            if(!((S>>i)&1)) continue;
            double base = dp[rowS + i];
            if(base >= 1e29) continue;
            double m = (penStep && !isP[i+1]) ? 1.1 : 1.0;
            const double* di = dm[i+1];
            for(int j = rem; j; j &= j-1){
                int jb = __builtin_ctz((unsigned)j);
                double nd = base + m * di[jb+1];
                size_t idx = (size_t)(S | (1u<<jb)) * M + jb;
                if(nd < dp[idx]){ dp[idx] = nd; par[idx] = (int8_t)i; }
            }
        }
    }
    double bestC = 1e30; int bi = 0;
    for(int i=0;i<M;i++){
        double m = (N % 10 == 0 && !isP[i+1]) ? 1.1 : 1.0;
        double c = dp[(size_t)(FULL-1)*M + i] + m * dm[i+1][0];
        if(c < bestC){ bestC = c; bi = i; }
    }
    tour.assign(N+1, 0);
    int S = FULL - 1, cur = bi;
    for(int k = N-1; k >= 1; --k){
        tour[k] = cur + 1;
        int p = par[(size_t)S*M + cur];
        S ^= (1 << cur);
        cur = p;
    }
    return true;
}

int main(){
    T0 = clock();
    initIn();
    N = (int)readLL();
    if(N < 2){ printf("2\n0\n0\n"); return 0; }
    X.resize(N); Y.resize(N);
    for(int i=0;i<N;i++){ X[i] = readLL(); Y[i] = readLL(); }

    // primes over city IDs (0 and 1 are not prime)
    {
        int M = max(2, N);
        isP.assign(M, 1);
        isP[0] = 0; if(M > 1) isP[1] = 0;
        for(ll i=2; i*i < (ll)M; ++i) if(isP[i])
            for(ll j=i*i; j < (ll)M; j += i) isP[j] = 0;
    }

    vector<int> tour(N+1);
    bool solvedExact = (N <= 19) && heldKarp(tour);

    if(!solvedExact){

    // ---- Hilbert-sorted order of all cities ----
    vector<pair<u64,int>> ord(N);
    {
        const ll OFF = 1000000000LL;          // |x|,|y| <= 1e9 -> [0, 2e9] < 2^31
        for(int i=0;i<N;i++)
            ord[i] = { hilb((u32)(X[i]+OFF), (u32)(Y[i]+OFF)), i };
        sort(ord.begin(), ord.end());
    }

    // ---- seed tour: break the Hilbert cycle at city 0 ----
    {
        int p = 0; for(int i=0;i<N;i++) if(ord[i].second == 0){ p = i; break; }
        tour[0] = 0;
        for(int m=1; m<N; ++m) tour[m] = ord[(p + m) % N].second;
        tour[N] = 0;
    }

    // ---- candidate lists: true k-NN via k-d tree (fallback: Hilbert neighbours) ----
    int KC = 24;
    if(N > 150000) KC = 20; else if(N <= 2000) KC = 32;
    vector<vector<int>> cand(N);
    for(int i=0;i<N;i++) cand[i].reserve(KC + 4);

    bool useKD = (N >= 4);   // k-d tree k-NN; for tiny N Hilbert neighbours already cover all
    if(useKD){
        g_KC = KC;
        g_kd.resize(N); iota(g_kd.begin(), g_kd.end(), 0);
        kdBuild(0, N, 0);
        NNHeap heap(KC);
        for(int p = 0; p < N; ++p){
            if((p & 0x3FF) == 0 && elapsed() > 1.1){ useKD = false; break; }
            heap.reset();
            kdQuery(0, N, 0, p, heap);
            for(int j = 0; j < heap.n; ++j) cand[p].push_back(heap.hi[j]);
        }
        // discard partial k-d results if we aborted midway
        if(!useKD){ for(int i=0;i<N;i++) cand[i].clear(); }
    }
    if(!useKD){
        // Hilbert-order neighbour fallback
        int half = 12;
        if(N > 100000) half = 10; else if(N <= 2000) half = 24;
        for(int s = 0; s < N; ++s){
            int city = ord[s].second;
            for(int d = 1; d <= half; ++d){
                if(s-d >= 0) cand[city].push_back(ord[s-d].second);
                if(s+d <  N) cand[city].push_back(ord[s+d].second);
            }
        }
    }

    // ---- flatten candidates to CSR, sorted by distance (enables gain pruning) ----
    vector<int> cs(N+1, 0);
    vector<int> ca;                // candidate ids
    vector<double> cd;             // exact D(i, cand) cached
    {
        size_t tot = 0;
        for(int i=0;i<N;i++) tot += cand[i].size();
        ca.resize(tot); cd.resize(tot);
        vector<pair<double,int>> tmp;
        size_t w = 0;
        for(int i=0;i<N;i++){
            cs[i] = (int)w;
            tmp.clear();
            for(int v : cand[i]) tmp.push_back({D(i,v), v});
            sort(tmp.begin(), tmp.end());
            for(auto& pr : tmp){ cd[w] = pr.first; ca[w] = pr.second; ++w; }
        }
        cs[N] = (int)w;
        vector<vector<int>>().swap(cand);   // free the nested lists
    }

    // position map (city -> index in tour); city 0 maps to index 0
    vector<int> pos(N+1);
    auto rebuildPos = [&](vector<int>& t){
        for(int i=(int)t.size()-1; i>=0; --i) pos[t[i]] = i;   // i=0 last => pos[0]=0
    };
    rebuildPos(tour);

    // 2-opt / Or-opt local search (candidate lists + don't-look bits).
    int MAXREV;
    if(N <= 2000)       MAXREV = N + 5;      // full 2-opt for tiny N
    else if(N <= 100000) MAXREV = 1000;
    else                 MAXREV = 400;
    // Coordinates are ~1e9, so double distance error is ~1e-6 (dx*dx ~ 1e18 overflows
    // the 52-bit mantissa). Use an epsilon well above that noise floor but negligible
    // vs. any real edge length (>= ~1e4 here), to avoid infinite floating-point shuffling.
    const double eps = 1e-3;
    // candidate-scan pruning threshold (cd sorted ascending per node); tuned empirically.
    const double PRUNE = 8.0;

    vector<int> Q; Q.reserve(N+2);
    vector<char> inq(N+2, 0);
    auto enqueue = [&](int z){
        if(z >= 1 && z < N && !inq[z]){ inq[z] = 1; Q.push_back(z); }
    };

    // Runs until the active queue is empty or `deadline`. Mutates tour/pos, adds deltas to cost.
    auto localSearch = [&](double deadline, double& cost){
        int cnt = 0;
        while(!Q.empty()){
            if(++cnt >= 512){ cnt = 0; if(elapsed() > deadline) break; }
            int i = Q.back(); Q.pop_back(); inq[i] = 0;
            int a = tour[i-1], b = tour[i];
            double dab = D(a,b);
            bool imp = false;

            // candidates of a (source): reconnect a to a spatially-near city v
            for(int k = cs[a], ke = cs[a+1]; k < ke; ++k){
                double dav = cd[k];
                if(dav >= dab * PRUNE) break;  // sorted: gain criterion
                int v = ca[k];
                int pv = pos[v];
                if(pv == i || pv == i-1) continue;
                if(pv > i){
                    int len = pv - i + 1;
                    if(len > MAXREV) continue;
                    int d = tour[pv+1];
                    double delta = dav + D(b,d) - dab - D(v,d);
                    if(delta < -eps){
                        reverse(tour.begin()+i, tour.begin()+pv+1);
                        for(int t=i;t<=pv;t++) pos[tour[t]] = t;
                        cost += delta;
                        enqueue(i-1); enqueue(i); enqueue(pv); enqueue(pv+1);
                        imp = true; break;
                    }
                } else { // pv < i-1
                    int len = i - pv - 1;
                    if(len > MAXREV) continue;
                    int d = tour[pv+1];
                    double delta = dav + D(d,b) - D(v,d) - dab;
                    if(delta < -eps){
                        reverse(tour.begin()+pv+1, tour.begin()+i);
                        for(int t=pv+1;t<i;t++) pos[tour[t]] = t;
                        cost += delta;
                        enqueue(pv); enqueue(pv+1); enqueue(i-1); enqueue(i);
                        imp = true; break;
                    }
                }
            }

            if(!imp){
                // candidates of b: move b forward next to a near city v (Or-opt flavour)
                for(int k = cs[b], ke = cs[b+1]; k < ke; ++k){
                    double dbv = cd[k];
                    if(dbv >= dab * PRUNE) break;  // sorted: gain criterion
                    int v = ca[k];
                    int pv = pos[v];
                    if(pv <= i) continue;
                    int len = pv - i;
                    if(len > MAXREV) continue;
                    int c = tour[pv-1];
                    double delta = D(a,c) + dbv - dab - D(c,v);
                    if(delta < -eps){
                        reverse(tour.begin()+i, tour.begin()+pv);
                        for(int t=i;t<pv;t++) pos[tour[t]] = t;
                        cost += delta;
                        enqueue(i-1); enqueue(i); enqueue(pv-1); enqueue(pv);
                        imp = true; break;
                    }
                }
            }

            if(!imp){
                // Or-opt: relocate segment [i..i+L-1] (L=1,2,3), forward or reversed.
                // Insertion points come from candidates of segA and segB.
                for(int L = 1; L <= 3 && !imp; ++L){
                    if(i + L - 1 > N - 1) break;          // segment must stay in the interior
                    int segA = tour[i], segB = tour[i + L - 1];
                    int pre = tour[i - 1], post = tour[i + L];   // i+L <= N
                    double save = D(pre,segA) + D(segB,post) - D(pre,post);
                    auto tryJ = [&](int j) -> bool {
                        if(j < 0 || j > N - 1) return false;   // never insert before city 0 / past the end
                        if(j >= i + L){                   // move segment right, insert after j
                            if(j - i + 1 > MAXREV) return false;
                            int ja = tour[j], jb = tour[j + 1];
                            double base = -save - D(ja,jb);
                            double dF = base + D(ja,segA) + D(segB,jb);
                            double dR = (L > 1) ? (base + D(ja,segB) + D(segA,jb)) : 1e18;
                            bool rev = dR < dF;
                            double delta = rev ? dR : dF;
                            if(delta < -eps){
                                if(!rev) reverse(tour.begin()+i, tour.begin()+i+L);
                                reverse(tour.begin()+i+L, tour.begin()+j+1);
                                reverse(tour.begin()+i, tour.begin()+j+1);
                                for(int t=i;t<=j;t++) pos[tour[t]] = t;
                                cost += delta;
                                enqueue(i-1); enqueue(i); enqueue(j); enqueue(j+1);
                                return true;
                            }
                        } else if(j <= i - 2){            // move segment left, insert after j
                            if(j + L > N - 1) return false;
                            if(i + L - 1 - j > MAXREV) return false;
                            int ja = tour[j], jb = tour[j + 1];
                            double base = -save - D(ja,jb);
                            double dF = base + D(ja,segA) + D(segB,jb);
                            double dR = (L > 1) ? (base + D(ja,segB) + D(segA,jb)) : 1e18;
                            bool rev = dR < dF;
                            double delta = rev ? dR : dF;
                            if(delta < -eps){
                                reverse(tour.begin()+j+1, tour.begin()+i+L);
                                if(!rev) reverse(tour.begin()+j+1, tour.begin()+j+1+L);
                                reverse(tour.begin()+j+1+L, tour.begin()+i+L);
                                for(int t=j+1;t<=i+L-1;t++) pos[tour[t]] = t;
                                cost += delta;
                                enqueue(j); enqueue(j+1); enqueue(i+L-1); enqueue(i+L);
                                return true;
                            }
                        }
                        return false;
                    };
                    for(int k = cs[segA], ke = cs[segA+1]; k < ke; ++k){
                        if(cd[k] >= save * PRUNE) break;            // sorted: gain criterion
                        if(tryJ(pos[ca[k]])){ imp = true; break; }
                    }
                    if(!imp) for(int k = cs[segB], ke = cs[segB+1]; k < ke; ++k){
                        if(cd[k] >= save * PRUNE) break;            // sorted: gain criterion
                        if(tryJ(pos[ca[k]] - 1)){ imp = true; break; }
                    }
                }
            }

            if(imp) enqueue(i);
        }
    };

    auto resetActive = [&](){ Q.clear(); fill(inq.begin(), inq.end(), (char)0); };

    // penalized extra length: sum of D over penalized (non-prime-source) steps.
    // penalized(t) == euclid(t) + 0.1 * penExtra(t)
    auto penExtra = [&](const vector<int>& t) -> double {
        double s = 0.0;
        for(int k=10; k<=N; k+=10) if(!isP[t[k-1]]) s += D(t[k-1], t[k]);
        return s;
    };

    // time budgets: small N gets a dedicated penalized-objective phase
    const bool SMALL = (N <= 2000);
    const double TL_A       = SMALL ? 1.00 : TL_MAIN;   // euclidean ILS end
    const double PEN_ILS_END = 1.55;                    // penalized ILS end (small N)
    const double PEN_END     = SMALL ? 1.78 : TL_PEN;   // final passes end

    // ---- initial local search (all positions active) ----
    double cost = euclidCost(tour);
    for(int i=1;i<N;i++) enqueue(i);
    localSearch(TL_A, cost);

    vector<int> best = tour; double bestCost = cost;
    vector<int> current = tour; double curCost = cost;   // walking point for ILS

    const bool TRACK_PEN = (N <= 20000);
    vector<int> bestPenTour; double bestPenCost = 1e308;
    if(TRACK_PEN){ bestPenTour = tour; bestPenCost = cost + 0.1 * penExtra(tour); }

    // ---- Iterated Local Search with double-bridge kicks (threshold acceptance) ----
    // (needs >= 6 cities so a valid 3-cut double-bridge with gaps exists)
    if(N >= 6){
        mt19937 rng(33333u);
        uniform_int_distribution<int> u(1, N-1);
        // windowed kicks for large N: uniform random cuts on big tours produce huge,
        // almost-always-rejected perturbations; local double-bridges repair fast.
        const int GW = (N > 3000) ? min(MAXREV - 8, 512) : 0;
        vector<int> nt(N+1);                 // reusable kick buffer (avoid per-iter alloc)
        while(elapsed() < TL_A){
            tour = current; cost = curCost;     // kick from the walking point, not the best
            rebuildPos(tour);

            int a, b, c;
            if(GW){
                int g1 = 2 + (int)(rng() % (u32)GW), g2 = 2 + (int)(rng() % (u32)GW);
                int hi = N - 1 - g1 - g2;
                if(hi < 1) continue;
                a = 1 + (int)(rng() % (u32)hi);
                b = a + g1; c = b + g2;
            } else {
                int tries = 0;
                do{
                    a = u(rng); b = u(rng); c = u(rng);
                    if(a>b) swap(a,b); if(b>c) swap(b,c); if(a>b) swap(a,b);
                    if(++tries > 24) break;
                } while(b - a < 2 || c - b < 2);
                if(b - a < 2 || c - b < 2 || a < 1 || c > N-1) continue;
            }

            // double-bridge on the path: S0 S1 S2 S3 -> S0 S2 S1 S3
            double kd = D(tour[a],tour[b+1]) + D(tour[c],tour[a+1]) + D(tour[b],tour[c+1])
                      - D(tour[a],tour[a+1]) - D(tour[b],tour[b+1]) - D(tour[c],tour[c+1]);
            int p = 0;
            for(int i=0;i<=a;i++)        nt[p++] = tour[i];
            for(int i=b+1;i<=c;i++)      nt[p++] = tour[i];
            for(int i=a+1;i<=b;i++)      nt[p++] = tour[i];
            for(int i=c+1;i<=N;i++)      nt[p++] = tour[i];
            tour.swap(nt);
            cost += kd;
            rebuildPos(tour);

            int b3 = a + (c - b);          // end of S2 in new positions
            resetActive();
            enqueue(a);   enqueue(a+1);
            enqueue(b3);  enqueue(b3+1);
            enqueue(c);   enqueue(c+1);
            localSearch(TL_A, cost);
            resetActive();

            if(cost < bestCost - eps){ best = tour; bestCost = cost; }
            if(TRACK_PEN){
                double pc = cost + 0.1 * penExtra(tour);
                if(pc < bestPenCost){ bestPenCost = pc; bestPenTour = tour; }
            }
            // Threshold acceptance: keep walking if the result is within ~1% of the best
            if(cost < curCost - eps || cost < bestCost * 1.01){
                current = tour; curCost = cost;
            }
        }
        if(TRACK_PEN) tour = bestPenTour;   // continue with the best penalized tour
        else          tour = best;
        rebuildPos(tour);
    }

    // ================= penalized-objective refinement =================
    const double peps = 1e-3;
    auto mAt = [&](int k, int src) -> double {
        return (k%10==0 && !isP[src]) ? 1.1 : 1.0;
    };
    // penalised contribution of steps [a..b] of tour t (a>=1, b<=N)
    auto contrib = [&](const vector<int>& t, int a, int b) -> double {
        double s = 0.0;
        for(int k=a; k<=b; ++k){
            int src = t[k-1], dst = t[k];
            s += mAt(k, src) * D(src, dst);
        }
        return s;
    };

    // exact penalized delta of reversing tour[l..r], 1<=l<r<=N-1, in O(2 + (r-l)/10):
    // internal euclidean edge sums are reversal-invariant, so only the boundary edges
    // and the 10% surcharges of internal penalized steps change.
    auto penRevDelta = [&](int l, int r) -> double {
        int a = tour[l-1], b = tour[l], c = tour[r], d = tour[r+1];
        double delta = mAt(l, a) * (D(a,c) - D(a,b));
        delta += mAt(r+1, b) * D(b,d) - mAt(r+1, c) * D(c,d);
        int t0 = ((l + 10) / 10) * 10;               // first multiple of 10 >= l+1
        for(int t = t0; t <= r; t += 10){
            int oldSrc = tour[t-1], oldDst = tour[t];
            int p = l + r - t;                        // new step t covers old edge (p, p+1), reversed
            int newSrc = tour[p+1], newDst = tour[p];
            double add = 0.0;
            if(!isP[newSrc]) add += D(newSrc, newDst);
            if(!isP[oldSrc]) add -= D(oldSrc, oldDst);
            delta += 0.1 * add;
        }
        return delta;
    };

    // queue-driven candidate 2-opt on the TRUE penalized objective
    auto penLocalSearch = [&](double deadline, double& cost, vector<int>& pos_){
        int cnt = 0;
        while(!Q.empty()){
            if(++cnt >= 256){ cnt = 0; if(elapsed() > deadline) break; }
            int i = Q.back(); Q.pop_back(); inq[i] = 0;
            int a = tour[i-1], b = tour[i];
            double dab = D(a,b);
            bool imp = false;
            for(int k = cs[a], ke = cs[a+1]; k < ke; ++k){
                if(cd[k] >= dab * 1.25) break;   // penalty slack is at most ~10%
                int v = ca[k];
                int pv = pos_[v];
                int l, r;
                if(pv > i){ l = i; r = pv; }
                else if(pv < i-1){ l = pv+1; r = i-1; }
                else continue;
                if(l >= r || r - l + 1 > MAXREV) continue;
                double delta = penRevDelta(l, r);
                if(delta < -peps){
                    reverse(tour.begin()+l, tour.begin()+r+1);
                    for(int t=l;t<=r;t++) pos_[tour[t]] = t;
                    cost += delta;
                    enqueue(l-1); enqueue(l); enqueue(r); enqueue(r+1);
                    imp = true; break;
                }
            }
            if(!imp){
                for(int k = cs[b], ke = cs[b+1]; k < ke; ++k){
                    if(cd[k] >= dab * 1.25) break;
                    int v = ca[k];
                    int pv = pos_[v];
                    if(pv <= i + 1) continue;
                    int l = i, r = pv - 1;
                    if(l >= r || r - l + 1 > MAXREV) continue;
                    double delta = penRevDelta(l, r);
                    if(delta < -peps){
                        reverse(tour.begin()+l, tour.begin()+r+1);
                        for(int t=l;t<=r;t++) pos_[tour[t]] = t;
                        cost += delta;
                        enqueue(l-1); enqueue(l); enqueue(r); enqueue(r+1);
                        imp = true; break;
                    }
                }
            }
            if(imp) enqueue(i);
        }
    };

    // small-window reversal pass around penalized source positions
    auto windowRevPass = [&](double deadline) -> bool {
        bool any = false;
        int W = (N <= 2000) ? 10 : 6;
        for(int t = 10; t <= N; t += 10){
            if((t & 0x3F) == 0 && elapsed() > deadline) break;
            int s = t - 1;                       // source position of penalised step t
            if(isP[tour[s]]) continue;           // already prime -> no penalty

            int lo0 = max(1, s - W), hi0 = min(N - 1, s + W);
            for(int lo = lo0; lo <= s; ++lo){
                for(int hi = s; hi <= hi0; ++hi){
                    if(lo == hi) continue;        // no-op reversal
                    int newSrc = tour[lo + hi - s];
                    if(!isP[newSrc]) continue;    // only useful if a prime lands at s
                    int aa = lo, bb = hi + 1;     // affected steps: lo .. hi+1
                    double oldC = contrib(tour, aa, bb);
                    reverse(tour.begin() + lo, tour.begin() + hi + 1);
                    double newC = contrib(tour, aa, bb);
                    if(newC < oldC - peps){
                        for(int k=lo; k<=hi; ++k) pos[tour[k]] = k;
                        any = true;
                    } else {
                        reverse(tour.begin() + lo, tour.begin() + hi + 1);  // undo
                    }
                }
            }
        }
        return any;
    };

    // swap a non-prime penalized source with a spatially-near prime city
    auto primeSwapPass = [&](double deadline) -> bool {
        bool any = false;
        auto stepC = [&](int k) -> double {
            int src = tour[k-1], dst = tour[k];
            return mAt(k, src) * D(src, dst);
        };
        auto swapCost = [&](int u, int v) -> double {   // u < v; affected steps only
            double s = stepC(u) + stepC(u+1) + stepC(v+1);
            if(v > u + 1) s += stepC(v);
            return s;
        };
        for(int t = 10; t <= N; t += 10){
            if(elapsed() > deadline) break;
            int s = t - 1;
            int c = tour[s];
            if(isP[c]) continue;
            double bestD = -peps; int bestV = -1;
            for(int k = cs[c], ke = cs[c+1]; k < ke; ++k){
                int p = ca[k];
                if(!isP[p]) continue;
                int v = pos[p];
                if(v < 1 || v > N-1 || v == s) continue;
                int u2 = min(s, v), w2 = max(s, v);
                double before = swapCost(u2, w2);
                swap(tour[u2], tour[w2]);
                double after = swapCost(u2, w2);
                swap(tour[u2], tour[w2]);
                double d = after - before;
                if(d < bestD){ bestD = d; bestV = v; }
            }
            if(bestV >= 0){
                int u2 = min(s, bestV), w2 = max(s, bestV);
                swap(tour[u2], tour[w2]);
                pos[tour[u2]] = u2; pos[tour[w2]] = w2;
                any = true;
            }
        }
        return any;
    };

    if(SMALL){
        // ---- dedicated penalized phase: pen-2opt to convergence, then pen-ILS ----
        double pc = penalized(tour);
        resetActive();
        for(int i=1;i<N;i++) enqueue(i);
        penLocalSearch(PEN_ILS_END, pc, pos);
        resetActive();
        {   // quick structural penalty passes
            int r = 0;
            while(r++ < 4 && elapsed() < PEN_ILS_END){
                bool ch = windowRevPass(PEN_ILS_END);
                ch |= primeSwapPass(PEN_ILS_END);
                if(!ch) break;
            }
        }
        pc = penalized(tour);
        vector<int> bestP = tour; double bestPC = pc;
        vector<int> curP = tour;  double curPC = pc;

        if(N >= 6){
            mt19937 rng(33340u);
            uniform_int_distribution<int> u(1, N-1);
            vector<int> nt(N+1);
            while(elapsed() < PEN_ILS_END){
                tour = curP; rebuildPos(tour);
                int a, b, c, tries = 0;
                do{
                    a = u(rng); b = u(rng); c = u(rng);
                    if(a>b) swap(a,b); if(b>c) swap(b,c); if(a>b) swap(a,b);
                    if(++tries > 24) break;
                } while(b - a < 2 || c - b < 2);
                if(b - a < 2 || c - b < 2 || a < 1 || c > N-1) continue;
                int p = 0;
                for(int i=0;i<=a;i++)   nt[p++] = tour[i];
                for(int i=b+1;i<=c;i++) nt[p++] = tour[i];
                for(int i=a+1;i<=b;i++) nt[p++] = tour[i];
                for(int i=c+1;i<=N;i++) nt[p++] = tour[i];
                tour.swap(nt);
                rebuildPos(tour);
                pc = penalized(tour);          // positions shifted -> full recompute

                int b3 = a + (c - b);
                resetActive();
                enqueue(a);  enqueue(a+1);
                enqueue(b3); enqueue(b3+1);
                enqueue(c);  enqueue(c+1);
                penLocalSearch(PEN_ILS_END, pc, pos);
                resetActive();

                if(pc < bestPC - peps){ bestP = tour; bestPC = pc; }
                if(pc < curPC - peps || pc < bestPC * 1.003){ curP = tour; curPC = pc; }
            }
        }
        tour = bestP; rebuildPos(tour);
        // final structural passes
        int r = 0;
        while(r++ < 6 && elapsed() < PEN_END){
            bool ch = windowRevPass(PEN_END);
            ch |= primeSwapPass(PEN_END);
            if(!ch) break;
        }
    } else {
        // ---- large N: penalized 2-opt sweep (all large N), then structural passes ----
        if(elapsed() < TL_PEN - 0.02){
            double pc = penalized(tour);
            resetActive();
            for(int i=1;i<N;i++) enqueue(i);
            penLocalSearch(TL_PEN - 0.02, pc, pos);
            resetActive();
        }
        int r = 0;
        while(r++ < 8 && elapsed() < TL_PEN){
            bool ch = windowRevPass(TL_PEN);
            ch |= primeSwapPass(TL_PEN);
            if(!ch) break;
        }
    }

    } // !solvedExact

    // ---- validity check + monotonic fallback ----
    auto valid = [&](const vector<int>& t) -> bool {
        if((int)t.size() != N+1) return false;
        if(t[0] != 0 || t[N] != 0) return false;
        vector<char> seen(N, 0);
        for(int i=1;i<N;i++){
            int c = t[i];
            if(c < 0 || c >= N || seen[c]) return false;
            seen[c] = 1;
        }
        for(int c=1;c<N;c++) if(!seen[c]) return false;
        return true;
    };

    vector<int> mono(N+1);
    for(int i=0;i<N;i++) mono[i] = i; mono[N] = 0;

    bool useTour = valid(tour);
    double Lt = useTour ? penalized(tour) : 1e308;
    double Lm = penalized(mono);
    if(Lt > Lm) useTour = false;
    const vector<int>& bestOut = useTour ? tour : mono;

    // ---- fast output (manual int formatting; _Exit skips big-vector destructor teardown) ----
    {
        vector<char> ob; ob.reserve((size_t)(N+2) * 8);
        char tmp[16];
        long v0 = (long)N + 1; int p0 = 0;
        if(v0 == 0) tmp[p0++] = '0';
        while(v0 > 0){ tmp[p0++] = char('0' + v0 % 10); v0 /= 10; }
        while(p0 > 0) ob.push_back(tmp[--p0]);
        ob.push_back('\n');
        for(int i=0;i<=N;i++){
            int v = bestOut[i]; int p = 0;
            if(v == 0) tmp[p++] = '0';
            while(v > 0){ tmp[p++] = char('0' + v % 10); v /= 10; }
            while(p > 0) ob.push_back(tmp[--p]);
            ob.push_back('\n');
        }
        fwrite(ob.data(), 1, ob.size(), stdout);
        fflush(stdout);
    }
    _Exit(0);
}
