#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC optimize("O3,unroll-loops")
#endif
#include <cstdio>
#include <vector>
#include <array>
#include <algorithm>
#include <chrono>
#include <random>
#include <climits>
#include <ctime>
#include <cmath>
#include <unistd.h>
using namespace std;
#ifdef NO_REOPT
#define NO_SWEEP
#define NO_TRIG
#endif
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
static long long gCmax = 0;
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
msu[prev] = v;
mpred[v] = prev;
ind[v]++;
prev = v;
}
msu[prev] = -1;
}
fill(dist_.begin(), dist_.begin()+n, 0);
qbuf.clear();
int qh=0;
for(int u=0;u<n;++u) if(ind[u]==0){ ds[u]=pn[u]; qbuf.push_back(u); }
while(qh < (int)qbuf.size()){
int u = qbuf[qh++];
long long du = ds[u];
int v = jsu[u];
if(v >= 0){
long long nd = du + pn[v];
if(nd > ds[v]) ds[v] = nd;
if(--ind[v]==0) qbuf.push_back(v);
}
v = msu[u];
if(v >= 0){
long long nd = du + pn[v];
if(nd > ds[v]) ds[v] = nd;
if(--ind[v]==0) qbuf.push_back(v);
}
}
if(qh != n) return -1;
long long C = 0;
for(int u=0;u<n;++u) if(ds[u] > C) C = ds[u];
gCmax = C;
if(fillCrit){
long long* __restrict tl = tail_.data();
const int* __restrict qb = qbuf.data();
for(int idx=n-1; idx>=0; --idx){
int u = qb[idx];
long long mx = 0;
int v = jsu[u];
if(v >= 0 && tl[v] > mx) mx = tl[v];
v = msu[u];
if(v >= 0 && tl[v] > mx) mx = tl[v];
tl[u] = pn[u] + mx;
}
fill(crit.begin(), crit.begin()+n, 0);
for(int u=0;u<n;++u) if(ds[u] + tl[u] - pn[u] == C) crit[u] = 1;
}
return C;
}
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
int k = jp[j], m = m_of[j][k];
long long s = max(jr[j], mf[m]);
long long f = s + p_of[j][k];
if(f < bf) bf = f;
}
int cm = -1;
for(int j=0;j<J;++j){
if(jp[j] >= M) continue;
int k = jp[j], m = m_of[j][k];
long long s = max(jr[j], mf[m]);
if(s + p_of[j][k] == bf){ cm = m; break; }
}
int cj = -1; long long cp = 0;
for(int j=0;j<J;++j){
if(jp[j] >= M) continue;
int k = jp[j], m = m_of[j][k];
if(m != cm) continue;
long long s = max(jr[j], mf[m]);
if(s < bf){
long long pr;
if(mode==0) pr = wrem[j];
else if(mode==1) pr = p_of[j][k];
else if(mode==2) pr = -p_of[j][k];
else pr = (long long)rng();
if(cj==-1 || pr > cp){ cp = pr; cj = j; }
}
}
if(cj == -1){
for(int j=0;j<J;++j) if(jp[j]<M && m_of[j][jp[j]]==cm){ cj=j; break; }
}
int k = jp[cj], m = m_of[cj][k];
long long s = max(jr[cj], mf[m]);
long long f = s + p_of[cj][k];
seq[m].push_back(cj);
jr[cj] = f; mf[m] = f; wrem[cj] -= p_of[cj][k];
jp[cj]++; remaining--;
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
// N6 critical-block neighborhood (Nowicki–Smutnicki / Zhang):
// any operation in a critical block may move to either end of that block.
static void genMoves(const vector<vector<int>>& cur){
gmoves.clear();
const long long* __restrict ds = dist_.data();
const long long* __restrict tl = tail_.data();
const long long* __restrict pn = pnode.data();
const long long C = gCmax;
#define ISCRIT_(u) (ds[u] + tl[u] - pn[u] == C)
for(int m=0;m<M;++m){
const auto& s = cur[m];
int i = 0;
while(i < J){
if(!ISCRIT_(opOf(s[i], m))){ i++; continue; }
int b = i;
while(i+1 < J && ISCRIT_(opOf(s[i+1], m))) i++;
int e = i; i++;
if(e == b) continue;
for(int t=b+1; t<=e; ++t) gmoves.push_back({m,b,e,t,true});
for(int t=b; t<e; ++t)
if(!(t==b && e==b+1))
gmoves.push_back({m,b,e,t,false});
}
}
#undef ISCRIT_
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
if(++pops > cap){
for(int t=wh;t<(int)gwl.size();++t) ginq[gwl[t]] = 0;
return -2;
}
long long b = 0;
if(jpre[v] && ds[v-1] > b) b = ds[v-1];
int mp = mpred[v];
if(mp >= 0 && ds[mp] > b) b = ds[mp];
long long nd = pn[v] + b;
if(nd != ds[v]){ ds[v] = nd; pushH(jsuc[v]); pushH(msucc[v]); }
}
gwl.clear();
int loT = lo > 0 ? lo-1 : lo;
for(int i=hi;i>=loT;--i){ int v = opOf(s[i], m); if(!ginq[v]){ ginq[v]=1; gwl.push_back(v); } }
int headPops = pops;
wh = 0; pops = 0;
while(wh < (int)gwl.size()){
int v = gwl[wh++]; ginq[v] = 0;
if(++pops > cap){
for(int t=wh;t<(int)gwl.size();++t) ginq[gwl[t]] = 0;
return -2;
}
long long b = 0;
int js = jsuc[v];
if(js >= 0 && tl[js] > b) b = tl[js];
int ms = msucc[v];
if(ms >= 0 && tl[ms] > b) b = tl[ms];
long long nt = pn[v] + b;
if(nt != tl[v]){
tl[v] = nt;
if(jpre[v] && !ginq[v-1]){ ginq[v-1]=1; gwl.push_back(v-1); }
int mp = mpred[v];
if(mp >= 0 && !ginq[mp]){ ginq[mp]=1; gwl.push_back(mp); }
}
}
long long C = 0;
for(int mm=0;mm<M;++mm){
int u = opOf(cur[mm][J-1], mm);
if(ds[u] > C) C = ds[u];
}
gCmax = C;
#ifdef DIAG
extern long long g_pops, g_calls;
g_pops += pops + headPops; g_calls++;
#endif
return C;
}
#ifdef DIAG
long long g_pops = 0, g_calls = 0;
#endif
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
static long long evalReduced(int mExcl, const vector<vector<int>>& seq){
const int n = N;
int* __restrict ind = indeg.data();
int* __restrict msu = msucc.data();
const int* __restrict jsu = jsuc.data();
const char* __restrict jpr = jpre.data();
const long long* __restrict pn = pnode.data();
long long* __restrict ds = dist_.data();
for(int u=0;u<n;++u) ind[u] = jpr[u];
for(int m=0;m<M;++m){
if(m == mExcl){
for(int j=0;j<J;++j){ int u = j*M + posF[j*M + m]; msu[u] = -1; mpred[u] = -1; }
continue;
}
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
qbuf.clear();
int qh=0;
for(int u=0;u<n;++u) if(ind[u]==0){ ds[u]=pn[u]; qbuf.push_back(u); }
while(qh < (int)qbuf.size()){
int u = qbuf[qh++];
long long du = ds[u];
int v = jsu[u];
if(v >= 0){
long long nd = du + pn[v];
if(nd > ds[v]) ds[v] = nd;
if(--ind[v]==0) qbuf.push_back(v);
}
v = msu[u];
if(v >= 0){
long long nd = du + pn[v];
if(nd > ds[v]) ds[v] = nd;
if(--ind[v]==0) qbuf.push_back(v);
}
}
if(qh != n) return -1;
long long C = 0;
for(int u=0;u<n;++u) if(ds[u] > C) C = ds[u];
{
long long* __restrict tl = tail_.data();
const int* __restrict qb = qbuf.data();
for(int idx=n-1; idx>=0; --idx){
int u = qb[idx];
long long mx = 0;
int v = jsu[u];
if(v >= 0 && tl[v] > mx) mx = tl[v];
v = msu[u];
if(v >= 0 && tl[v] > mx) mx = tl[v];
tl[u] = pn[u] + mx;
}
}
return C;
}
static long long carBestC;
static vector<int> carBestSeq;
static int carNodes, carNodeCap;
static chrono::steady_clock::time_point carDeadline;
static vector<long long> carR, carQ, carP;
#ifdef DIAG
long long g_roAtt = 0, g_roAcc = 0, g_roCyc = 0, g_roUs = 0;
#endif
static long long schrage1(int n, const long long* r, const long long* q, const long long* p, int* seqOut){
static vector<char> done; done.assign(n, 0);
long long t = 0, C = 0;
for(int i=0;i<n;++i){
long long rmin = LLONG_MAX;
for(int j=0;j<n;++j) if(!done[j] && r[j] < rmin) rmin = r[j];
if(t < rmin) t = rmin;
int b = -1;
for(int j=0;j<n;++j){
if(done[j] || r[j] > t) continue;
if(b < 0 || q[j] > q[b] || (q[j] == q[b] && p[j] > p[b])) b = j;
}
done[b] = 1; seqOut[i] = b;
t += p[b];
if(t + q[b] > C) C = t + q[b];
}
return C;
}
static void carlier(int n, long long* r, long long* q, const long long* p, int depth){
if(carNodes >= carNodeCap || depth > 400) return;
if((carNodes & 63) == 0 && chrono::steady_clock::now() >= carDeadline){ carNodes = carNodeCap; return; }
carNodes++;
int seq[64];
long long C = schrage1(n, r, q, p, seq);
if(C < carBestC){
carBestC = C;
for(int i=0;i<n;++i) carBestSeq[i] = seq[i];
}
int bpos = -1;
{
long long t = 0;
for(int i=0;i<n;++i){
int j = seq[i];
if(t < r[j]) t = r[j];
t += p[j];
if(t + q[j] == C) bpos = i;
}
}
if(bpos < 0) return;
int apos = -1;
{
long long sum = 0;
long long qb = q[seq[bpos]];
for(int i=bpos;i>=0;--i){
sum += p[seq[i]];
if(r[seq[i]] + sum + qb == C) apos = i;
}
}
if(apos < 0) return;
int cpos = -1;
for(int i=apos;i<bpos;++i) if(q[seq[i]] < q[seq[bpos]]) cpos = i;
if(cpos < 0) return;
int c = seq[cpos];
long long sump = 0, rmin = LLONG_MAX, qmin = LLONG_MAX;
for(int i=cpos+1;i<=bpos;++i){
int j = seq[i];
sump += p[j];
if(r[j] < rmin) rmin = r[j];
if(q[j] < qmin) qmin = q[j];
}
long long hJ = rmin + sump + qmin;
long long hJc = (rmin < r[c] ? rmin : r[c]) + sump + p[c] + (qmin < q[c] ? qmin : q[c]);
long long LB = hJ > hJc ? hJ : hJc;
if(LB >= carBestC) return;
{
long long old = q[c];
long long nq = sump + qmin;
if(nq > old){
q[c] = nq;
carlier(n, r, q, p, depth+1);
q[c] = old;
} else carlier(n, r, q, p, depth+1);
}
if(carNodes >= carNodeCap) return;
{
long long old = r[c];
long long nr = rmin + sump;
if(nr > old){
r[c] = nr;
carlier(n, r, q, p, depth+1);
r[c] = old;
}
}
}
static vector<int> roOld;
static bool reoptMachine(int m, vector<vector<int>>& cur, long long& curC, long long& bestC,
vector<vector<int>>& best, chrono::steady_clock::time_point T_end, chrono::steady_clock::time_point& lastImpT){
#ifdef DIAG
g_roAtt++;
auto diagT0 = chrono::steady_clock::now();
struct RoTimer { chrono::steady_clock::time_point t0; ~RoTimer(){ g_roUs += chrono::duration_cast<chrono::microseconds>(chrono::steady_clock::now() - t0).count(); } } roTimer{diagT0};
#endif
long long Cr = evalReduced(m, cur);
if(Cr < 0){ evalSeq(cur, true); return false; }
carR.resize(J); carQ.resize(J); carP.resize(J);
for(int j=0;j<J;++j){
int u = j*M + posF[j*M + m];
carP[j] = pnode[u];
carR[j] = dist_[u] - pnode[u];
carQ[j] = tail_[u] - pnode[u];
}
carBestC = LLONG_MAX;
carNodes = 0; carNodeCap = 3000;
carBestSeq.assign(J, 0);
auto nowT = chrono::steady_clock::now();
carDeadline = nowT + chrono::milliseconds(15);
if(carDeadline > T_end) carDeadline = T_end;
carlier(J, carR.data(), carQ.data(), carP.data(), 0);
if(carBestC >= curC || carBestSeq == cur[m]){
long long cc = evalSeq(cur, true);
if(cc >= 0) curC = cc;
return false;
}
roOld = cur[m];
cur[m] = carBestSeq;
long long nc = evalSeq(cur, true);
if(nc >= 0 && nc < curC){
curC = nc;
if(nc < bestC){ best = cur; bestC = nc; }
lastImpT = chrono::steady_clock::now();
#ifdef DIAG
g_roAcc++;
#endif
return true;
}
#ifdef DIAG
if(nc < 0) g_roCyc++;
#endif
cur[m] = roOld;
long long cc = evalSeq(cur, true);
if(cc >= 0) curC = cc;
return false;
}

static unsigned long long solHash(const vector<vector<int>>& s){
unsigned long long h = 1469598103934665603ULL;
for(int m=0;m<M;++m) for(int j=0;j<J;++j)
h = h*1099511628211ULL ^ (unsigned long long)(s[m][j]+1);
return h;
}
// Large-neighborhood: extract critical jobs and reinsert each job into the
// best position on every machine (greedy best-insertion). Escapes N6 plateaus
// by changing many machine orders at once while preserving feasibility search.
static bool reinsertCriticalJobs(vector<vector<int>>& cur, long long& curC,
long long& bestC, vector<vector<int>>& best,
mt19937& rng, chrono::steady_clock::time_point T_end, int maxJobs){
if(curC < 0) return false;
evalSeq(cur, true);
// Score jobs by total critical processing time
vector<long long> jcrit(J, 0);
for(int u=0;u<N;++u){
if(dist_[u] + tail_[u] - pnode[u] == gCmax)
jcrit[u/M] += pnode[u];
}
vector<pair<long long,int>> ord;
for(int j=0;j<J;++j) if(jcrit[j] > 0) ord.push_back({-jcrit[j], j});
if(ord.empty()) return false;
sort(ord.begin(), ord.end());
int nJ = min(maxJobs, (int)ord.size());
bool improved = false;
vector<vector<int>> work = cur;
long long workC = curC;
// For each critical job: on each machine, remove it and reinsert at the
// position that minimizes true makespan (other machines unchanged). This is a
// large neighborhood move that N6 cannot express in one step.
for(int t=0;t<nJ && chrono::steady_clock::now() < T_end; ++t){
int job = ord[t].second;
work = cur;
for(int m=0;m<M;++m){
if(chrono::steady_clock::now() >= T_end) break;
auto& s = work[m];
auto it = find(s.begin(), s.end(), job);
if(it == s.end()) continue;
s.erase(it);
int L = (int)s.size();
int bp = 0; long long bc = LLONG_MAX;
int step = L > 24 ? max(1, L/14) : 1;
for(int p=0;p<=L;p+=step){
s.insert(s.begin()+p, job);
long long nc = evalSeq(work, false);
if(nc >= 0 && nc < bc){ bc = nc; bp = p; }
s.erase(s.begin()+p);
}
for(int p=max(0,bp-3); p<=min(L,bp+3); ++p){
s.insert(s.begin()+p, job);
long long nc = evalSeq(work, false);
if(nc >= 0 && nc < bc){ bc = nc; bp = p; }
s.erase(s.begin()+p);
}
s.insert(s.begin()+bp, job);
}
long long nc = evalSeq(work, true);
if(nc < 0){ work = cur; continue; }
if(nc < bestC){ best = work; bestC = nc; improved = true; }
if(nc <= curC){ cur = work; curC = nc; improved = true; }
else if(nc < curC + curC/40){ cur = work; curC = nc; } // mild diversify
else { work = cur; }
}
return improved;
}
// Path relinking: walk from src toward dst by order-correcting adjacent swaps.
// Intermediate schedules on the geodesic often beat both endpoints (i-TSAB idea).
static void pathRelink(const vector<vector<int>>& src, const vector<vector<int>>& dst,
long long& bestC, vector<vector<int>>& best,
mt19937& rng, chrono::steady_clock::time_point T_end, int maxSteps){
vector<vector<int>> cur = src;
long long curC = evalSeq(cur, true);
if(curC < 0) return;
static vector<int> posT; posT.resize(J);
for(int step=0; step<maxSteps; ++step){
if(chrono::steady_clock::now() >= T_end) break;
// Collect correcting adjacent swaps with estimates
struct Cand { long long est; int m,i; };
vector<Cand> cands;
cands.reserve(M*J);
for(int m=0;m<M;++m){
for(int i=0;i<J;++i) posT[dst[m][i]] = i;
const auto& s = cur[m];
for(int i=0;i+1<J;++i){
if(posT[s[i]] > posT[s[i+1]]){
// estimate swap of i,i+1 via estMove on 2-block
Mv mv{m,i,i+1,i+1,true};
long long est = estMove(cur, mv);
cands.push_back({est, m, i});
}
}
}
if(cands.empty()) break;
// Pick best few estimates; take the true-best among top-K
int K = (int)min((size_t)12, cands.size());
partial_sort(cands.begin(), cands.begin()+K, cands.end(),
[](const Cand& a, const Cand& b){ return a.est < b.est; });
int bm=-1, bi=-1; long long bnc = LLONG_MAX;
int tryN = min(K, 4);
for(int t=0;t<tryN;++t){
int m=cands[t].m, i=cands[t].i;
swap(cur[m][i], cur[m][i+1]);
long long nc = evalSeq(cur, false);
swap(cur[m][i], cur[m][i+1]);
if(nc >= 0 && nc < bnc){ bnc = nc; bm = m; bi = i; }
}
if(bm < 0){
// force any correcting swap
int m=cands[0].m, i=cands[0].i;
swap(cur[m][i], cur[m][i+1]);
long long nc = evalSeq(cur, true);
if(nc < 0){ swap(cur[m][i], cur[m][i+1]); evalSeq(cur, true); break; }
curC = nc;
} else {
swap(cur[bm][bi], cur[bm][bi+1]);
curC = evalSeq(cur, true);
}
if(curC >= 0 && curC < bestC){ best = cur; bestC = curC; }
// Occasional random jump among remaining correcting swaps for diversity
if((step & 7) == 7 && cands.size() > 1){
int t = 1 + (int)(rng() % (unsigned)min((size_t)6, cands.size()-1));
int m=cands[t].m, i=cands[t].i;
// only if still inverted
for(int z=0;z<J;++z) posT[dst[m][z]] = z;
if(posT[cur[m][i]] > posT[cur[m][i+1]]){
swap(cur[m][i], cur[m][i+1]);
long long nc = evalSeq(cur, true);
if(nc < 0){ swap(cur[m][i], cur[m][i+1]); evalSeq(cur, true); }
else {
curC = nc;
if(nc < bestC){ best = cur; bestC = nc; }
}
}
}
}
// Short N6 polish on the relinked solution
if(curC > 0 && chrono::steady_clock::now() < T_end){
evalSeq(cur, true);
for(int polish=0; polish<200 && chrono::steady_clock::now() < T_end; ++polish){
genMoves(cur);
if(gmoves.empty()) break;
gcand.clear();
for(int i=0;i<(int)gmoves.size();++i) gcand.push_back({estMove(cur, gmoves[i]), i});
int K = min((int)gcand.size(), 8);
partial_sort(gcand.begin(), gcand.begin()+K, gcand.end());
bool ok=false;
for(int t=0;t<K;++t){
const Mv& mv = gmoves[gcand[t].second];
applyMove(cur, mv);
long long nc = evalSeq(cur, true);
if(nc >= 0 && nc <= curC){
curC = nc; ok=true;
if(nc < bestC){ best = cur; bestC = nc; }
break;
}
undoMove(cur, mv);
}
if(!ok) break;
}
evalSeq(cur, true);
}
}
int main(){
auto T0 = chrono::steady_clock::now();
const auto budget = chrono::milliseconds(985);
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
msucc.assign(N, -1);
mpred.assign(N, -1);
jsuc.assign(N, -1); jpre.assign(N, 0);
for(int u=0;u<N;++u){
if(u % M != M-1) jsuc[u] = u+1;
if(u % M != 0) jpre[u] = 1;
}
indeg.resize(N); dist_.resize(N); tail_.resize(N); crit.resize(N);
gwl.reserve(4*N); ginq.assign(N, 0);
qbuf.reserve(N);
gord.resize(J); gestC.resize(J);
tabuTB.assign((size_t)N*J, 0);
vector<vector<int>> best(M, vector<int>(J));
for(int m=0;m<M;++m) for(int j=0;j<J;++j) best[m][j] = j;
long long bestC = evalSeq(best);
if(J <= 1){
for(int m=0;m<M;++m){ printf("0\n"); }
fflush(stdout);
_exit(0);
}
vector<vector<int>> cur = best;
long long curC = bestC;
auto trySeed = [&](const vector<vector<int>>& s){
long long c = evalSeq(s);
if(c > 0 && c < curC){
cur = s; curC = c;
if(c < bestC){ best = s; bestC = c; }
}
};
mt19937 rng(777u);
// Multi-rule + randomized Giffler-Thompson seeds (correct model for JSSP).
// Note: permutation-flow-shop NEH/IG is intentionally omitted — different routes
// make a single job permutation a wrong model and wasted ~100ms for zero gain.
trySeed(seedGT(0, rng));
trySeed(seedGT(1, rng));
trySeed(seedGT(2, rng));
trySeed(seedGT(3, rng));
for(int r=0;r<50 && chrono::steady_clock::now()-T0 < budget;++r)
trySeed(seedGT(3, rng));

long long LB = 0;
{
vector<long long> mload(M, 0);
for(int j=0;j<J;++j){
long long jl = 0;
for(int k=0;k<M;++k){ jl += p_of[j][k]; mload[m_of[j][k]] += p_of[j][k]; }
if(jl > LB) LB = jl;
}
for(int m=0;m<M;++m) if(mload[m] > LB) LB = mload[m];
}
auto T_end = T0 + budget;
if(bestC > LB){
long long c = evalSeq(cur, true);
if(c > 0) curC = c; else { cur = best; curC = evalSeq(cur, true); }
int iter = 0, sinceImp = 0;
const int TENURE_MIN = 12;
const int TENURE_SPAN = max(4, J/2);
bool timeUp = false;
const int spanShort = max(4, J/3);
const int spanLong  = TENURE_SPAN;
auto lastImpT = chrono::steady_clock::now();
auto dynTen = [&](chrono::steady_clock::time_point nowT)->int{
long long stag = chrono::duration_cast<chrono::milliseconds>(nowT - lastImpT).count();
int tmin, span;
if(stag <= 30){ tmin = 8; span = spanShort; }
else if(stag >= 120){
tmin = TENURE_MIN + (int)min(8LL, (stag - 120)/200);
span = spanLong;
} else {
int f = (int)((stag - 30)*100/90);
tmin = 8 + (TENURE_MIN - 8)*f/100;
span = spanShort + (spanLong - spanShort)*f/100;
}
return tmin + (int)(rng() % (unsigned)span);
};
struct Elite { vector<vector<int>> seq; long long C; unsigned long long h; };
vector<Elite> pool;
const int E = 8; int poolHead = 0;
auto poolAdd = [&](const vector<vector<int>>& s, long long C){
unsigned long long h = solHash(s);
for(const Elite& e : pool) if(e.h == h) return;
if((int)pool.size() < E) pool.push_back({s, C, h});
else {
// replace worst
int wi=0; for(int i=1;i<(int)pool.size();++i) if(pool[i].C>pool[wi].C) wi=i;
if(C <= pool[wi].C + pool[wi].C/100) { pool[wi] = {s, C, h}; }
else { pool[poolHead] = {s, C, h}; poolHead = (poolHead+1)%E; }
}
};
vector<vector<int>> rb; long long rbC = LLONG_MAX;
int failCnt = 0;
auto lastReoptT = T0;
int lastReoptM = -1;
poolAdd(best, bestC);
// Initial shifting-bottleneck sweep
if(J > 2 && J <= 64){
vector<pair<long long,int>> mord(M);
vector<long long> mload(M, 0);
for(int j=0;j<J;++j) for(int k=0;k<M;++k) mload[m_of[j][k]] += p_of[j][k];
for(int m=0;m<M;++m) mord[m] = {-mload[m], m};
sort(mord.begin(), mord.end());
for(int t=0;t<M;++t){
if(chrono::steady_clock::now() >= T_end) break;
reoptMachine(mord[t].second, cur, curC, bestC, best, T_end, lastImpT);
}
poolAdd(best, bestC);
if(bestC <= LB) timeUp = true;
}
// Main N6 tabu loop with full re-eval (incremental was ~N pops = no win)
while(!timeUp){
auto nowT = chrono::steady_clock::now();
if(nowT >= T_end) break;
iter++;
// Periodic critical-machine Carlier reopt
if(J > 2 && J <= 64 && nowT - lastImpT > chrono::milliseconds(100)
&& nowT - lastReoptT > chrono::milliseconds(40)){
lastReoptT = nowT;
static vector<long long> critW, loadW;
critW.assign(M, 0); loadW.assign(M, 0);
for(int u=0;u<N;++u){
int mm = m_of[u/M][u%M];
loadW[mm] += pnode[u];
if(dist_[u] + tail_[u] - pnode[u] == gCmax) critW[mm] += pnode[u];
}
int pick = -1;
for(int mm=0;mm<M;++mm){
if(mm == lastReoptM) continue;
if(pick < 0 || critW[mm] > critW[pick] || (critW[mm]==critW[pick] && loadW[mm]>loadW[pick])) pick = mm;
}
if(pick >= 0){
lastReoptM = pick;
reoptMachine(pick, cur, curC, bestC, best, T_end, lastImpT);
poolAdd(best, bestC);
}
}
genMoves(cur);
int nmv = (int)gmoves.size();
if(nmv == 0){
// random kick
int m = (int)(rng()% (unsigned)M);
int i = (int)(rng()% (unsigned)max(1,J-1));
swap(cur[m][i], cur[m][i+1]);
long long nc = evalSeq(cur, true);
if(nc < 0){ swap(cur[m][i], cur[m][i+1]); evalSeq(cur, true); }
else curC = nc;
continue;
}
gcand.clear();
for(int idx=0; idx<nmv; ++idx)
gcand.push_back({estMove(cur, gmoves[idx]), idx});
int K = nmv < 28 ? nmv : 28;
partial_sort(gcand.begin(), gcand.begin()+K, gcand.end());
bool sorted_all = (K == nmv);
bool applied = false;
// Select by estimate (O(block) each); full re-eval only for the chosen move.
// Full verification of top-K was tried and lost iteration throughput without
// enough quality gain under the 1s limit.
for(int pass=0; pass<2 && !applied && !timeUp; ++pass){
for(int t=0; t<nmv; ++t){
if((t & 15)==15 && chrono::steady_clock::now() >= T_end){ timeUp = true; break; }
if(t >= K && !sorted_all){ sort(gcand.begin(), gcand.end()); sorted_all = true; }
const Mv& mv = gmoves[gcand[t].second];
if(pass==0){
bool tb = isTabu(cur, mv, iter);
bool asp = gcand[t].first < bestC;
if(tb && !asp) continue;
}
collectTabu(cur, mv);
applyMove(cur, mv);
long long nc = evalSeq(cur, true);
if(nc < 0){ undoMove(cur, mv); evalSeq(cur, true); continue; }
int tenure = dynTen(nowT);
for(size_t id : gpend) tabuTB[id] = iter + tenure;
curC = nc; applied = true; break;
}
}
if(!applied) break;
if(curC < rbC){ rbC = curC; rb = cur; }
if(curC < bestC){
long long exact = evalSeq(cur, true);
curC = exact;
if(exact >= 0 && exact < bestC){
best = cur; bestC = exact; sinceImp = 0;
lastImpT = chrono::steady_clock::now();
poolAdd(cur, exact);
failCnt = 0;
if(bestC <= LB) break;
}
} else {
++sinceImp;
long long stag = chrono::duration_cast<chrono::milliseconds>(nowT - lastImpT).count();
// Restart earlier (220ms) and use PATH RELINKING as diversification
if(stag > 180){
if(rbC <= bestC + bestC/40 && !rb.empty()) poolAdd(rb, rbC);
// Path-relink between elites / best (foundational i-TSAB step)
if((int)pool.size() >= 2 && chrono::steady_clock::now() < T_end){
// Pair best elite with a diverse (worse) elite for longer geodesic
int bi=0; for(int z=1;z<(int)pool.size();++z) if(pool[z].C<pool[bi].C) bi=z;
int wi=0; for(int z=1;z<(int)pool.size();++z) if(pool[z].C>pool[wi].C) wi=z;
pathRelink(pool[bi].seq, pool[wi].seq, bestC, best, rng, T_end, 70);
pathRelink(pool[wi].seq, pool[bi].seq, bestC, best, rng, T_end, 50);
if(pool.size()>=3){
int k = (int)(rng()%pool.size());
if(k!=bi) pathRelink(pool[bi].seq, pool[k].seq, bestC, best, rng, T_end, 40);
}
poolAdd(best, bestC);
} else if(!pool.empty() && chrono::steady_clock::now() < T_end){
pathRelink(cur, best, bestC, best, rng, T_end, 50);
poolAdd(best, bestC);
}
// Critical-job LNS (large neighborhood beyond N6)
if(chrono::steady_clock::now() < T_end){
cur = best; curC = bestC;
reinsertCriticalJobs(cur, curC, bestC, best, rng, T_end, 3);
poolAdd(best, bestC);
}
// Reseed trajectory
if(pool.empty() || (rng() & 1)) cur = best;
else cur = pool[rng()%pool.size()].seq;
curC = evalSeq(cur, true);
// Structural kicks on critical blocks
int kicks = 3 + failCnt + (int)(rng()%3); if(kicks > 10) kicks = 10;
for(int r=0;r<kicks;++r){
if(chrono::steady_clock::now() >= T_end){ timeUp=true; break; }
genMoves(cur);
if(gmoves.empty()) break;
const Mv& mv = gmoves[rng()%gmoves.size()];
applyMove(cur, mv);
long long nc = evalSeq(cur, true);
if(nc < 0){ undoMove(cur, mv); evalSeq(cur, true); }
else {
curC = nc;
if(nc < bestC){ best=cur; bestC=nc; poolAdd(best,bestC); }
}
}
fill(tabuTB.begin(), tabuTB.end(), 0);
sinceImp = 0; rbC = LLONG_MAX; rb.clear();
lastImpT = chrono::steady_clock::now();
++failCnt;
}
}
}
}
// Final SB polish
if(J<=64 && chrono::steady_clock::now() < T_end && bestC > LB){
cur = best; curC = bestC;
vector<pair<long long,int>> mord(M);
vector<long long> mload(M,0);
for(int j=0;j<J;++j) for(int k=0;k<M;++k) mload[m_of[j][k]]+=p_of[j][k];
for(int m=0;m<M;++m) mord[m]={-mload[m],m};
sort(mord.begin(), mord.end());
auto dummy = chrono::steady_clock::now();
for(int t=0;t<M && chrono::steady_clock::now()<T_end;++t)
reoptMachine(mord[t].second, cur, curC, bestC, best, T_end, dummy);
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
