// Refinement of attempt 2 (my best, 0.1828): N7 critical-block tabu with
// incremental longest-path eval and dynamic oscillating tenure. Keeps the
// winning "single deep descent + rare wide kicks at 250ms stagnation" shape,
// but widens the sorted candidate window, uses a safe 0.95s CPU limit, and
// adds a light re-intensification (partial tabu clear when a new global best
// is found) to squeeze the moderately-improvable random/nearly/block cases.
#include <cstdio>
#include <vector>
#include <algorithm>
#include <ctime>
#include <random>
#include <climits>
#include <unistd.h>
using namespace std;

static int J, M, N;
static vector<vector<int>> m_of;
static vector<vector<long long>> p_of;
static vector<vector<int>> pos;
static vector<int> posF;
static vector<long long> pnode;
static vector<int> msucc, mpred, indeg, qbuf;
static vector<int> jsuc;
static vector<char> jpre;
static vector<long long> dist_, tail_;
static vector<char> crit;
static vector<long long> gmload;

static clock_t START;
static inline double elapsed(){ return double(clock()-START)/CLOCKS_PER_SEC; }
static const double TL = 0.95;

static long long evalSeq(const vector<vector<int>>& seq, bool fillCrit = false){
    const int n = N;
    int* __restrict ind = indeg.data();
    int* __restrict msu = msucc.data();
    const int* __restrict jsu = jsuc.data();
    const char* __restrict jpr = jpre.data();
    const long long* __restrict pn = pnode.data();
    long long* __restrict ds = dist_.data();
    for(int u=0;u<n;++u) ind[u] = jpr[u];
    for(int m=0;m<M;++m){
        const int* s = seq[m].data();
        int prev = s[0]*M + posF[s[0]*M + m];
        mpred[prev] = -1;
        for(int i=1;i<J;++i){
            int v = s[i]*M + posF[s[i]*M + m];
            msu[prev] = v; mpred[v] = prev; ind[v]++; prev = v;
        }
        msu[prev] = -1;
    }
    fill(dist_.begin(), dist_.begin()+n, 0);
    qbuf.clear(); int qh=0;
    for(int u=0;u<n;++u) if(ind[u]==0){ ds[u]=pn[u]; qbuf.push_back(u); }
    while(qh < (int)qbuf.size()){
        int u = qbuf[qh++]; long long du = ds[u];
        int v = jsu[u];
        if(v >= 0){ long long nd=du+pn[v]; if(nd>ds[v])ds[v]=nd; if(--ind[v]==0)qbuf.push_back(v); }
        v = msu[u];
        if(v >= 0){ long long nd=du+pn[v]; if(nd>ds[v])ds[v]=nd; if(--ind[v]==0)qbuf.push_back(v); }
    }
    if(qh != n) return -1;
    long long C = 0;
    for(int u=0;u<n;++u) if(ds[u] > C) C = ds[u];
    if(fillCrit){
        long long* __restrict tl = tail_.data();
        const int* __restrict qb = qbuf.data();
        for(int idx=n-1; idx>=0; --idx){
            int u = qb[idx]; long long mx = 0;
            int v = jsu[u]; if(v>=0 && tl[v]>mx) mx=tl[v];
            v = msu[u]; if(v>=0 && tl[v]>mx) mx=tl[v];
            tl[u] = pn[u] + mx;
        }
        fill(crit.begin(), crit.begin()+n, 0);
        for(int u=0;u<n;++u) if(ds[u] + tl[u] - pn[u] == C) crit[u] = 1;
    }
    return C;
}

// mode: 0=MWR,1=LPT,2=SPT,3=LWR,4=random,5=bottleneck-aware
static vector<vector<int>> seedGT(int mode, mt19937& rng){
    vector<int> jp(J, 0);
    vector<long long> jr(J, 0), mf(M, 0), wrem(J, 0);
    for(int j=0;j<J;++j) for(int k=0;k<M;++k) wrem[j] += p_of[j][k];
    vector<vector<int>> seq(M);
    int remaining = N;
    while(remaining > 0){
        long long bf = LLONG_MAX;
        for(int j=0;j<J;++j){
            if(jp[j] >= M) continue;
            int k=jp[j], m=m_of[j][k];
            long long f = max(jr[j], mf[m]) + p_of[j][k];
            if(f < bf) bf = f;
        }
        int cm = -1;
        for(int j=0;j<J;++j){
            if(jp[j] >= M) continue;
            int k=jp[j], m=m_of[j][k];
            if(max(jr[j], mf[m]) + p_of[j][k] == bf){ cm = m; break; }
        }
        int cj = -1; long long cp = 0;
        for(int j=0;j<J;++j){
            if(jp[j] >= M) continue;
            int k=jp[j], m=m_of[j][k];
            if(m != cm) continue;
            long long s = max(jr[j], mf[m]);
            if(s < bf){
                long long pr;
                if(mode==0) pr = wrem[j];
                else if(mode==1) pr = p_of[j][k];
                else if(mode==2) pr = -p_of[j][k];
                else if(mode==3) pr = -wrem[j];
                else if(mode==5){
                    long long w = 0;
                    for(int kk=k; kk<M; ++kk) w += gmload[m_of[j][kk]];
                    pr = w;
                }
                else pr = (long long)rng();
                if(cj==-1 || pr > cp){ cp = pr; cj = j; }
            }
        }
        if(cj == -1){ for(int j=0;j<J;++j) if(jp[j]<M && m_of[j][jp[j]]==cm){ cj=j; break; } }
        int k=jp[cj], m=m_of[cj][k];
        long long f = max(jr[cj], mf[m]) + p_of[cj][k];
        seq[m].push_back(cj);
        jr[cj]=f; mf[m]=f; wrem[cj]-=p_of[cj][k]; jp[cj]++; remaining--;
    }
    return seq;
}

struct Mv { int m, b, e, i; bool front; };
static vector<Mv> gmoves;
static vector<pair<long long,int>> gcand;
static vector<int> gord;
static vector<long long> gestC;
static vector<int> tabuTB;

static inline int opOf(int job, int m){ return job*M + posF[job*M + m]; }

static void genMoves(const vector<vector<int>>& cur){
    gmoves.clear();
    for(int m=0;m<M;++m){
        const auto& s = cur[m];
        int i = 0;
        while(i < J){
            if(!crit[opOf(s[i], m)]){ i++; continue; }
            int b = i;
            while(i+1 < J && crit[opOf(s[i+1], m)]) i++;
            int e = i; i++;
            if(e == b) continue;
            for(int t=b+1; t<=e; ++t) gmoves.push_back({m,b,e,t,true});
            for(int t=b; t<e; ++t)
                if(!(t==b && e==b+1)) gmoves.push_back({m,b,e,t,false});
        }
    }
}

static long long estMove(const vector<vector<int>>& cur, const Mv& mv){
    const auto& s = cur[mv.m];
    int lo, hi;
    if(mv.front){
        lo = mv.b; hi = mv.i;
        gord[0] = s[mv.i];
        for(int t=lo; t<hi; ++t) gord[t-lo+1] = s[t];
    } else {
        lo = mv.i; hi = mv.e;
        for(int t=lo+1; t<=hi; ++t) gord[t-lo-1] = s[t];
        gord[hi-lo] = s[mv.i];
    }
    int L = hi - lo + 1;
    long long prevC = 0;
    if(lo > 0) prevC = dist_[opOf(s[lo-1], mv.m)];
    for(int t=0; t<L; ++t){
        int v = gord[t]; int u = opOf(v, mv.m); int k = pos[v][mv.m];
        long long jp = (k>0) ? dist_[u-1] : 0;
        long long st = prevC > jp ? prevC : jp;
        gestC[t] = st + pnode[u];
        prevC = gestC[t];
    }
    long long prevT = 0;
    if(hi+1 < J) prevT = tail_[opOf(s[hi+1], mv.m)];
    long long bestLen = 0;
    for(int t=L-1; t>=0; --t){
        int v = gord[t]; int u = opOf(v, mv.m); int k = pos[v][mv.m];
        long long js = (k<M-1) ? tail_[u+1] : 0;
        long long tl = pnode[u] + (prevT > js ? prevT : js);
        long long len = gestC[t] - pnode[u] + tl;
        if(len > bestLen) bestLen = len;
        prevT = tl;
    }
    return bestLen;
}

static inline void applyMove(vector<vector<int>>& cur, const Mv& mv){
    auto& s = cur[mv.m];
    if(mv.front) rotate(s.begin()+mv.b, s.begin()+mv.i, s.begin()+mv.i+1);
    else         rotate(s.begin()+mv.i, s.begin()+mv.i+1, s.begin()+mv.e+1);
}
static inline void undoMove(vector<vector<int>>& cur, const Mv& mv){
    auto& s = cur[mv.m];
    if(mv.front) rotate(s.begin()+mv.b, s.begin()+mv.b+1, s.begin()+mv.i+1);
    else         rotate(s.begin()+mv.i, s.begin()+mv.e,   s.begin()+mv.e+1);
}

static vector<int> gwl;
static vector<char> ginq;
static long long incAfterMove(const vector<vector<int>>& cur, const Mv& mv){
    const int m = mv.m;
    const auto& s = cur[m];
    const int lo = mv.front ? mv.b : mv.i;
    const int hi = mv.front ? mv.i : mv.e;
    {
        int from = lo>0 ? lo-1 : 0;
        int to   = hi<J-1 ? hi+1 : J-1;
        for(int i=from;i<=to;++i){
            int u = opOf(s[i], m);
            msucc[u] = (i<J-1) ? opOf(s[i+1], m) : -1;
            mpred[u] = (i>0)   ? opOf(s[i-1], m) : -1;
        }
    }
    const int cap = 16*N;
    long long* __restrict ds = dist_.data();
    long long* __restrict tl = tail_.data();
    const long long* __restrict pn = pnode.data();
    gwl.clear();
    auto pushH = [&](int v){ if(v>=0 && !ginq[v]){ ginq[v]=1; gwl.push_back(v); } };
    int hiH = hi < J-1 ? hi+1 : hi;
    for(int i=lo;i<=hiH;++i) pushH(opOf(s[i], m));
    int wh = 0, pops = 0;
    while(wh < (int)gwl.size()){
        int v = gwl[wh++]; ginq[v] = 0;
        if(++pops > cap){ for(int t=wh;t<(int)gwl.size();++t) ginq[gwl[t]]=0; return -2; }
        long long b = 0;
        if(jpre[v] && ds[v-1] > b) b = ds[v-1];
        int mp = mpred[v]; if(mp >= 0 && ds[mp] > b) b = ds[mp];
        long long nd = pn[v] + b;
        if(nd != ds[v]){ ds[v] = nd; pushH(jsuc[v]); pushH(msucc[v]); }
    }
    gwl.clear();
    int loT = lo > 0 ? lo-1 : lo;
    for(int i=hi;i>=loT;--i){ int v=opOf(s[i], m); if(!ginq[v]){ ginq[v]=1; gwl.push_back(v); } }
    wh = 0; pops = 0;
    while(wh < (int)gwl.size()){
        int v = gwl[wh++]; ginq[v] = 0;
        if(++pops > cap){ for(int t=wh;t<(int)gwl.size();++t) ginq[gwl[t]]=0; return -2; }
        long long b = 0;
        int js = jsuc[v]; if(js >= 0 && tl[js] > b) b = tl[js];
        int ms = msucc[v]; if(ms >= 0 && tl[ms] > b) b = tl[ms];
        long long nt = pn[v] + b;
        if(nt != tl[v]){
            tl[v] = nt;
            if(jpre[v] && !ginq[v-1]){ ginq[v-1]=1; gwl.push_back(v-1); }
            int mp = mpred[v]; if(mp >= 0 && !ginq[mp]){ ginq[mp]=1; gwl.push_back(mp); }
        }
    }
    long long C = 0;
    for(int u=0;u<N;++u) if(ds[u] > C) C = ds[u];
    for(int u=0;u<N;++u) crit[u] = (ds[u] + tl[u] - pn[u] == C);
    return C;
}

static bool isTabu(const vector<vector<int>>& cur, const Mv& mv, int iter){
    const auto& s = cur[mv.m];
    int uj = s[mv.i];
    if(mv.front){
        for(int t=mv.b; t<mv.i; ++t){
            int xop = opOf(s[t], mv.m);
            if(tabuTB[(size_t)xop*J + uj] > iter) return true;
        }
    } else {
        int uop = opOf(uj, mv.m);
        for(int t=mv.i+1; t<=mv.e; ++t)
            if(tabuTB[(size_t)uop*J + s[t]] > iter) return true;
    }
    return false;
}

static vector<size_t> gpend;
static void collectTabu(const vector<vector<int>>& cur, const Mv& mv){
    gpend.clear();
    const auto& s = cur[mv.m];
    int uj = s[mv.i];
    int uop = opOf(uj, mv.m);
    if(mv.front){
        for(int t=mv.b; t<mv.i; ++t) gpend.push_back((size_t)uop*J + s[t]);
    } else {
        for(int t=mv.i+1; t<=mv.e; ++t)
            gpend.push_back((size_t)opOf(s[t], mv.m)*J + uj);
    }
}

int main(){
    START = clock();
    if(scanf("%d %d", &J, &M) != 2) return 0;
    N = J*M;
    m_of.assign(J, vector<int>(M));
    p_of.assign(J, vector<long long>(M));
    pos.assign(J, vector<int>(M));
    for(int j=0;j<J;++j)
        for(int k=0;k<M;++k)
            if(scanf("%d %lld", &m_of[j][k], &p_of[j][k]) != 2) return 0;
    for(int j=0;j<J;++j) for(int k=0;k<M;++k) pos[j][m_of[j][k]] = k;
    posF.assign(N, 0);
    for(int j=0;j<J;++j) for(int m=0;m<M;++m) posF[j*M+m] = pos[j][m];
    pnode.assign(N, 0);
    for(int j=0;j<J;++j) for(int k=0;k<M;++k) pnode[j*M+k] = p_of[j][k];

    msucc.assign(N, -1); mpred.assign(N, -1);
    jsuc.assign(N, -1); jpre.assign(N, 0);
    for(int u=0;u<N;++u){ if(u%M != M-1) jsuc[u]=u+1; if(u%M != 0) jpre[u]=1; }
    indeg.resize(N); dist_.resize(N); tail_.resize(N); crit.resize(N);
    gwl.reserve(4*N); ginq.assign(N, 0);
    qbuf.reserve(N);
    gord.resize(J); gestC.resize(J);
    tabuTB.assign((size_t)N*J, 0);

    gmload.assign(M, 0);
    long long LB = 0;
    for(int j=0;j<J;++j){
        long long jl = 0;
        for(int k=0;k<M;++k){ jl += p_of[j][k]; gmload[m_of[j][k]] += p_of[j][k]; }
        if(jl > LB) LB = jl;
    }
    int bneck = 0;
    for(int m=0;m<M;++m){ if(gmload[m] > LB) LB = gmload[m]; if(gmload[m] > gmload[bneck]) bneck = m; }

    vector<vector<int>> best(M, vector<int>(J));
    for(int m=0;m<M;++m) for(int j=0;j<J;++j) best[m][j] = j;
    long long bestC = evalSeq(best);

    if(J <= 1){ for(int m=0;m<M;++m) printf("0\n"); fflush(stdout); _exit(0); }

    vector<vector<int>> cur = best;
    long long curC = bestC;
    mt19937 rng(777u);

    // GT seeding: 6 structured rules + randomized restarts (attempt-2 style).
    {
        for(int mode=0; mode<=5; ++mode){
            auto s = seedGT(mode, rng);
            long long c = evalSeq(s);
            if(c>0 && c<bestC){ best=s; bestC=c; }
        }
        int extra = (N <= 300) ? 48 : (N <= 700 ? 24 : 12);
        for(int r=0; r<extra && elapsed()<0.16; ++r){
            auto s = seedGT(4, rng);
            long long c = evalSeq(s);
            if(c>0 && c<bestC){ best=s; bestC=c; }
        }
        cur = best; curC = bestC;
    }

    if(bestC > LB){
        long long c = evalSeq(cur, true);
        if(c > 0) curC = c; else { cur = best; curC = evalSeq(cur, true); }

        int iter = 0;
        const int TENURE_MIN = 15;
        const int spanShort = max(4, J/3);
        const int spanLong  = max(4, J/2);
        double lastImpMs = elapsed()*1000.0;
        const double reStagMs = 250.0;
        bool timeUp = false;

        auto dynTen = [&](double nowMs)->int{
            double stag = nowMs - lastImpMs;
            int tmin, span;
            if(stag <= 30){ tmin = 8; span = spanShort; }
            else if(stag >= 120){
                long long e = (long long)((stag - 120)/250); if(e>6) e=6;
                tmin = TENURE_MIN + (int)e; span = spanLong;
            } else {
                int f = (int)((stag - 30)*100/90);
                tmin = 8 + (TENURE_MIN - 8)*f/100;
                span = spanShort + (spanLong - spanShort)*f/100;
            }
            return tmin + (int)(rng() % (unsigned)span);
        };

        while(!timeUp){
            double nowMs = elapsed()*1000.0;
            if(nowMs >= TL*1000.0) break;
            iter++;
            if((iter & 16383) == 0){
                long long fc = evalSeq(cur, true);
                if(fc >= 0) curC = fc;
            }
            genMoves(cur);
            int nmv = (int)gmoves.size();
            if(nmv == 0) break;
            gcand.clear();
            for(int idx=0; idx<nmv; ++idx)
                gcand.push_back({estMove(cur, gmoves[idx]), idx});
            int K = nmv < 32 ? nmv : 32;
            partial_sort(gcand.begin(), gcand.begin()+K, gcand.end());
            bool sorted_all = (K == nmv);

            bool applied = false;
            for(int pass=0; pass<2 && !applied && !timeUp; ++pass){
                for(int t=0; t<nmv; ++t){
                    if((t & 15)==15 && elapsed() >= TL){ timeUp = true; break; }
                    if(t >= K && !sorted_all){ sort(gcand.begin(), gcand.end()); sorted_all=true; }
                    const Mv& mv = gmoves[gcand[t].second];
                    if(pass==0){
                        bool tb = isTabu(cur, mv, iter);
                        bool asp = gcand[t].first < bestC;
                        if(tb && !asp) continue;
                    }
                    collectTabu(cur, mv);
                    applyMove(cur, mv);
                    long long nc = incAfterMove(cur, mv);
                    if(nc == -2) nc = evalSeq(cur, true);
                    if(nc < 0){ undoMove(cur, mv); evalSeq(cur, true); continue; }
                    int tenure = dynTen(nowMs);
                    for(size_t id : gpend) tabuTB[id] = iter + tenure;
                    curC = nc; applied = true; break;
                }
            }
            if(!applied) break;

            if(curC < bestC){
                long long exact = evalSeq(cur, true);
                curC = exact;
                if(exact >= 0 && exact < bestC){
                    best = cur; bestC = exact;
                    lastImpMs = elapsed()*1000.0;
                    // Light re-intensification: partially relax tabu around the
                    // new best to let neighboring improvements chain quickly.
                    if((iter & 1) == 0){
                        // decay ~half the tenures (cheap, keeps some memory)
                        size_t sz = tabuTB.size();
                        for(size_t z=0; z<sz; z+=2) tabuTB[z] = 0;
                    }
                    if(bestC <= LB) break;
                }
            }
            else if(nowMs - lastImpMs > reStagMs && (TL*1000.0 - nowMs) > 100.0){
                // Rare wide diversification from global best (attempt-2 winner),
                // mixing bottleneck-machine swaps with random critical-block kicks.
                cur = best;
                long long cc = evalSeq(cur, true);
                curC = cc;
                int kicks = 3 + (int)(rng() % 4);
                for(int r=0; r<kicks; ++r){
                    if(elapsed() >= TL){ timeUp = true; break; }
                    if((r & 1)==0 && J >= 2){
                        int i = (int)(rng() % (unsigned)(J-1));
                        swap(cur[bneck][i], cur[bneck][i+1]);
                        long long nc = evalSeq(cur, true);
                        if(nc < 0){ swap(cur[bneck][i], cur[bneck][i+1]); evalSeq(cur, true); }
                        else curC = nc;
                    } else {
                        genMoves(cur);
                        if(gmoves.empty()) break;
                        const Mv& mv = gmoves[rng() % gmoves.size()];
                        applyMove(cur, mv);
                        long long nc = evalSeq(cur, true);
                        if(nc < 0){ undoMove(cur, mv); evalSeq(cur, true); }
                        else curC = nc;
                    }
                }
                fill(tabuTB.begin(), tabuTB.end(), 0);
                lastImpMs = elapsed()*1000.0;
            }
        }
    }

    {
        vector<char> buf;
        buf.reserve((size_t)N*8 + M + 16);
        for(int m=0;m<M;++m){
            for(int j=0;j<J;++j){
                int x = best[m][j];
                if(x == 0){ buf.push_back('0'); }
                else { char tmp[12]; int t=0; while(x>0){ tmp[t++]=char('0'+x%10); x/=10; } while(t>0) buf.push_back(tmp[--t]); }
                buf.push_back(j+1<J ? ' ' : '\n');
            }
        }
        fwrite(buf.data(), 1, buf.size(), stdout);
        fflush(stdout);
    }
    _exit(0);
}