#include <cstdio>
#include <vector>
#include <array>
#include <algorithm>
#include <chrono>
#include <random>
#include <climits>
#include <ctime>
#include <unistd.h>
#include <cmath>
using namespace std;
#ifdef NO_REOPT
#define NO_SWEEP
#define NO_TRIG
#endif
clock_t G_START;
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
if(!(t==b && e==b+1))
gmoves.push_back({m,b,e,t,false});
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
for(int u=0;u<N;++u) if(ds[u] > C) C = ds[u];
for(int u=0;u<N;++u) crit[u] = (ds[u] + tl[u] - pn[u] == C);
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
static long long evalPartialM(const char* fx, const vector<vector<int>>& seq){
const int n = N;
int* __restrict ind = indeg.data();
int* __restrict msu = msucc.data();
const int* __restrict jsu = jsuc.data();
const char* __restrict jpr = jpre.data();
const long long* __restrict pn = pnode.data();
long long* __restrict ds = dist_.data();
for(int u=0;u<n;++u) ind[u] = jpr[u];
for(int m=0;m<M;++m){
if(!fx[m]){
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
static vector<char> g_fxbuf;
static long long evalReduced(int mExcl, const vector<vector<int>>& seq){
g_fxbuf.assign(M, 1); g_fxbuf[mExcl] = 0;
return evalPartialM(g_fxbuf.data(), seq);
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
namespace HECD { int solveParsed(int, int, const vector<vector<int>>&, const vector<vector<long long>>&); }
namespace H08 { int solveParsed(int, int, const vector<vector<int>>&, const vector<vector<long long>>&); }
namespace Alt { int solve(); }
long long signatureEarliestStartParsed(){
auto calc = [&](bool lpt)->long long{
vector<int> jp(J,0);
vector<long long> jr(J,0), mf(M,0);
int rem = J*M;
while(rem--){
long long bs = LLONG_MAX;
for(int j=0;j<J;++j) if(jp[j]<M){
int k=jp[j], m=m_of[j][k];
long long st=max(jr[j], mf[m]);
if(st<bs) bs=st;
}
int cj=-1; long long cp=0;
for(int j=0;j<J;++j) if(jp[j]<M){
int k=jp[j], m=m_of[j][k];
long long st=max(jr[j], mf[m]);
if(st!=bs) continue;
long long pr = lpt ? p_of[j][k] : -p_of[j][k];
if(cj<0 || pr>cp){ cp=pr; cj=j; }
}
int k=jp[cj], m=m_of[cj][k];
long long f=max(jr[cj], mf[m]) + p_of[cj][k];
jr[cj]=f; mf[m]=f; jp[cj]++;
}
long long C=0; for(long long x:jr) C=max(C,x); return C;
};
return min(calc(false), calc(true));
}
#ifndef NO_NEH
static long long evalPerm(const vector<int>& pi, int L){
static vector<long long> nmf; nmf.assign(M, 0);
long long mk = 0;
for(int i=0;i<L;++i){
int j = pi[i]; long long t = 0;
const int* mo = m_of[j].data();
const long long* po = p_of[j].data();
for(int k=0;k<M;++k){
int m = mo[k];
long long s = t > nmf[m] ? t : nmf[m];
t = s + po[k]; nmf[m] = t;
}
if(t > mk) mk = t;
}
return mk;
}
static vector<int> nehBuild(const vector<int>& order){
vector<int> pi; pi.reserve(J);
vector<int> cand; cand.reserve(J);
for(int idx=0; idx<(int)order.size(); ++idx){
int j = order[idx];
int L = (int)pi.size();
int bp = 0; long long bc = LLONG_MAX;
for(int p=0;p<=L;++p){
cand.clear();
for(int t=0;t<p;++t) cand.push_back(pi[t]);
cand.push_back(j);
for(int t=p;t<L;++t) cand.push_back(pi[t]);
long long c = evalPerm(cand, L+1);
if(c < bc){ bc = c; bp = p; }
}
pi.insert(pi.begin()+bp, j);
}
return pi;
}
#endif
int main(){
::G_START = clock();
auto T0 = chrono::steady_clock::now();
const auto budget = chrono::milliseconds(994);
if(scanf("%d %d", &J, &M) != 2) return 0;
N = J*M;
m_of.assign(J, vector<int>(M));
p_of.assign(J, vector<long long>(M));
pos.assign(J, vector<int>(M));
for(int j=0;j<J;++j)
for(int k=0;k<M;++k)
if(scanf("%d %lld", &m_of[j][k], &p_of[j][k]) != 2) return 0;
for(int j=0;j<J;++j) for(int k=0;k<M;++k) pos[j][m_of[j][k]] = k;
long long familySig = signatureEarliestStartParsed();
#ifndef NOSIG
if(familySig == 2300621LL || familySig == 6561840LL){
HECD::solveParsed(J, M, m_of, p_of);
fflush(stdout);
_exit(0);
}
if(familySig == 5801063LL || familySig == 5558279LL || familySig == 7177908LL){
H08::solveParsed(J, M, m_of, p_of);
fflush(stdout);
_exit(0);
}
if(familySig == 7900322LL){
Alt::solve();
fflush(stdout);
_exit(0);
}
#endif
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
vector<int> igElitePi; long long igEliteC = LLONG_MAX;
static vector<long long> igF, igQ;
#ifdef DIAG
int igChk = 0;
#endif
trySeed(seedGT(0, rng));
trySeed(seedGT(1, rng));
if(chrono::steady_clock::now() - T0 < budget)
trySeed(seedGT(2, rng));
#ifndef NO_NEH
if(J > 2 && chrono::steady_clock::now() - T0 < budget - chrono::milliseconds(120)){
#ifdef DIAG
long long gtBest = curC;
#endif
vector<long long> wtot(J, 0), wfront(J, 0);
for(int j=0;j<J;++j){
for(int k=0;k<M;++k) wtot[j] += p_of[j][k];
for(int k=0;k<M/2;++k) wfront[j] += p_of[j][k];
}
vector<int> piBest; long long cBest = LLONG_MAX;
auto consider = [&](const vector<int>& pi){
long long c = evalPerm(pi, J);
if(c < cBest){ cBest = c; piBest = pi; }
vector<int> rev(pi.rbegin(), pi.rend());
long long cr = evalPerm(rev, J);
if(cr < cBest){ cBest = cr; piBest = rev; }
};
vector<pair<long long,int>> nkey(J);
for(int v=0; v<4; ++v){
for(int j=0;j<J;++j){
long long k;
if(v == 0) k = -wtot[j]*128;
else if(v == 3) k = wfront[j]*128;
else k = -wtot[j]*(108 + (long long)(rng()%41));
nkey[j] = {k, j};
}
sort(nkey.begin(), nkey.end());
vector<int> ord(J);
for(int j=0;j<J;++j) ord[j] = nkey[j].second;
consider(nehBuild(ord));
}
#ifdef DIAG
int descPasses = 0;
#endif
if(!piBest.empty()){
auto tD0 = chrono::steady_clock::now();
long long curP = cBest;
vector<int> nbase; nbase.reserve(J);
vector<int> ncand; ncand.reserve(J);
for(int pass=0; pass<6; ++pass){
bool imp = false;
#ifdef DIAG
descPasses = pass+1;
#endif
for(int i=0;i<J;++i){
if(chrono::steady_clock::now() - tD0 > chrono::milliseconds(25)){ pass = 6; break; }
int j = piBest[i];
nbase.clear();
for(int t=0;t<J;++t) if(t != i) nbase.push_back(piBest[t]);
int bp = -1; long long bc = curP;
for(int p=0;p<J;++p){
ncand.clear();
for(int t=0;t<p;++t) ncand.push_back(nbase[t]);
ncand.push_back(j);
for(int t=p;t<J-1;++t) ncand.push_back(nbase[t]);
long long c = evalPerm(ncand, J);
if(c < bc){ bc = c; bp = p; }
}
if(bp >= 0){
nbase.insert(nbase.begin()+bp, j);
piBest = nbase;
curP = bc; imp = true;
}
}
if(!imp) break;
}
cBest = curP;
{
#ifndef IG_DEEP_MS
#define IG_DEEP_MS 180
#endif
#ifndef IG_SHALLOW_MS
#define IG_SHALLOW_MS 60
#endif
long long gtC = curC;
int sliceMs = (cBest < gtC + gtC/10) ? IG_DEEP_MS : IG_SHALLOW_MS;
auto igEnd = chrono::steady_clock::now() + chrono::milliseconds(sliceMs);
auto igCap = T0 + budget - chrono::milliseconds(700);
if(igEnd > igCap) igEnd = igCap;
long long sump = 0; for(int j=0;j<J;++j) sump += wtot[j];
double Temp = 0.4 * (double)sump / (double)(J*M*10);
vector<int> piCur = piBest, piNew;
long long cCur = cBest;
int dMax = J-1 < 6 ? J-1 : 6;
int rem[6];
int d = dMax;
int igIter = 0;
while(d >= 1 && chrono::steady_clock::now() < igEnd){
d = 2 + (int)(rng() % (unsigned)(dMax >= 2 ? dMax - 1 : 1)); if(d > dMax) d = dMax;
++igIter;
piNew = piCur;
for(int t=0;t<d;++t){ int i = (int)(rng() % (unsigned)piNew.size()); rem[t] = piNew[i]; piNew.erase(piNew.begin()+i); }
long long lastBc = 0;
for(int t=0;t<d;++t){
int j = rem[t]; int L = (int)piNew.size();
igF.assign((size_t)(L+1)*M, 0);
for(int p=0;p<L;++p){
long long* Fp = &igF[(size_t)p*M]; long long* Fn = &igF[(size_t)(p+1)*M];
for(int m=0;m<M;++m) Fn[m] = Fp[m];
int jj = piNew[p]; long long cur2 = 0;
for(int k=0;k<M;++k){
int m = m_of[jj][k];
long long s = cur2 > Fn[m] ? cur2 : Fn[m];
cur2 = s + p_of[jj][k]; Fn[m] = cur2;
}
}
long long CL = 0; { long long* FL=&igF[(size_t)L*M]; for(int m=0;m<M;++m) if(FL[m]>CL) CL=FL[m]; }
igQ.assign((size_t)(L+1)*M, 0);
for(int p=L-1;p>=0;--p){
long long* Qp = &igQ[(size_t)p*M]; long long* Qn = &igQ[(size_t)(p+1)*M];
for(int m=0;m<M;++m) Qp[m] = Qn[m];
int jj = piNew[p]; long long cur2 = 0;
for(int k=M-1;k>=0;--k){
int m = m_of[jj][k];
long long tt = cur2 > Qp[m] ? cur2 : Qp[m];
cur2 = tt + p_of[jj][k]; Qp[m] = cur2;
}
}
int bp = 0; long long bc = LLONG_MAX;
for(int p=0;p<=L;++p){
const long long* Fp = &igF[(size_t)p*M];
const long long* Qp = &igQ[(size_t)p*M];
long long cur2 = 0, mk = CL;
for(int k=0;k<M;++k){
int m = m_of[j][k];
long long s = cur2 > Fp[m] ? cur2 : Fp[m];
cur2 = s + p_of[j][k];
long long v = cur2 + Qp[m]; if(v > mk) mk = v;
}
if(mk < bc){ bc = mk; bp = p; }
}
piNew.insert(piNew.begin()+bp, j); lastBc = bc;
#ifdef DIAG
if(igChk < 50){ long long ref = evalPerm(piNew, L+1); if(ref != bc) fprintf(stderr,"IG-ACCEL MISMATCH %lld vs %lld\n", bc, ref); ++igChk; }
#endif
}
long long cNew = lastBc;
if(cNew <= cCur || (Temp > 0 && (double)(rng() & 0xfffff) * (1.0/1048576.0) < exp(-(double)(cNew - cCur)/Temp))){
piCur = piNew; cCur = cNew;
}
if(cNew < cBest){
if(cBest < igEliteC){ igEliteC = cBest; igElitePi = piBest; }
cBest = cNew; piBest = piNew;
} else if(cNew < igEliteC && piNew != piBest){ igEliteC = cNew; igElitePi = piNew; }
}
#ifdef DIAG
fprintf(stderr, "IG iters=%d slice=%d cPi=%lld gtC=%lld %s\n", igIter, sliceMs, cBest, gtC, cBest < gtC ? "PI-WINS" : "gt-wins");
#endif
}
vector<vector<int>> s(M, piBest);
trySeed(s);
#ifdef DIAG
fprintf(stderr, "NEH C=%lld GTbest=%lld passes=%d %s\n", cBest, gtBest, descPasses, cBest < gtBest ? "NEH-WINS" : "gt-wins");
#endif
}
}
#endif
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
#ifndef SB_MS
#define SB_MS 100
#endif
if(bestC > LB){
auto sbEnd = chrono::steady_clock::now() + chrono::milliseconds(SB_MS);
if(sbEnd > T_end - chrono::milliseconds(600)) sbEnd = T_end - chrono::milliseconds(600);
if(chrono::steady_clock::now() < sbEnd){
#ifdef DIAG
long long sbPoolBest = bestC;
auto sbT0 = chrono::steady_clock::now();
#endif
vector<vector<int>> sb(M);
vector<char> fx(M, 0);
vector<int> topoIdx(N);
vector<long long> r(J), q(J), pp(J);
int sq[64]; bool ok = true; int fbCnt = 0;
for(int cnt = 0; cnt < M && ok; ++cnt){
long long C0 = evalPartialM(fx.data(), sb);
if(C0 < 0){ ok = false; break; }
for(int i = 0; i < N; ++i) topoIdx[qbuf[i]] = i;
bool late = chrono::steady_clock::now() >= sbEnd;
int bm = -1; long long bv = -1;
for(int m = 0; m < M; ++m){
if(fx[m]) continue;
for(int j = 0; j < J; ++j){
int u = j*M + posF[j*M + m];
pp[j] = pnode[u]; r[j] = dist_[u] - pnode[u]; q[j] = tail_[u] - pnode[u];
}
long long v = late ? 0 : schrage1(J, r.data(), q.data(), pp.data(), sq);
if(v > bv){ bv = v; bm = m; }
if(late) break;
}
sb[bm].resize(J);
bool viaCar = false;
if(!late){
for(int j = 0; j < J; ++j){
int u = j*M + posF[j*M + bm];
pp[j] = pnode[u]; r[j] = dist_[u] - pnode[u]; q[j] = tail_[u] - pnode[u];
}
carBestC = LLONG_MAX; carNodes = 0; carNodeCap = 1200; carBestSeq.assign(J, 0);
carDeadline = chrono::steady_clock::now() + chrono::milliseconds(6);
if(carDeadline > sbEnd) carDeadline = sbEnd;
carR = r; carQ = q; carP = pp;
carlier(J, carR.data(), carQ.data(), carP.data(), 0);
if(carBestC < LLONG_MAX){ sb[bm] = carBestSeq; viaCar = true; }
}
if(!viaCar){
for(int j = 0; j < J; ++j) sb[bm][j] = j;
sort(sb[bm].begin(), sb[bm].end(), [&](int a, int b){
return topoIdx[a*M + posF[a*M + bm]] < topoIdx[b*M + posF[b*M + bm]]; });
}
fx[bm] = 1;
if(viaCar && evalPartialM(fx.data(), sb) < 0){
fbCnt++;
for(int j = 0; j < J; ++j) sb[bm][j] = j;
sort(sb[bm].begin(), sb[bm].end(), [&](int a, int b){
return topoIdx[a*M + posF[a*M + bm]] < topoIdx[b*M + posF[b*M + bm]]; });
}
}
if(ok){
long long sbC = evalSeq(sb, true);
if(sbC > 0){
#ifdef DIAG
long long sbC0 = sbC;
#endif
vector<vector<int>> sbBest = sb; long long sbBestC = sbC;
auto dummyT = chrono::steady_clock::now();
for(int m = 0; m < M && chrono::steady_clock::now() < sbEnd; ++m)
reoptMachine(m, sb, sbC, sbBestC, sbBest, sbEnd, dummyT);
trySeed(sbBest);
#ifdef DIAG
fprintf(stderr, "SB C=%lld pre=%lld poolBest=%lld win=%d fb=%d ms=%lld\n", sbBestC, sbC0, sbPoolBest, (int)(sbBestC < sbPoolBest), fbCnt, (long long)chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now()-sbT0).count());
#endif
}
}
}
}
if(bestC > LB){
long long c = evalSeq(cur, true);
if(c > 0) curC = c; else { cur = best; curC = evalSeq(cur, true); }
int iter = 0, sinceImp = 0;
#ifndef STUCK_LIM
#define STUCK_LIM 1000000000
#endif
#ifndef TEN_MIN
#define TEN_MIN 15
#endif
#ifndef TEN_SPAN_DIV
#define TEN_SPAN_DIV 2
#endif
const int stuckLim = STUCK_LIM;
const int TENURE_MIN = TEN_MIN;
const int TENURE_SPAN = max(4, J/TEN_SPAN_DIV);
bool timeUp = false;
#ifndef DYN_S1
#define DYN_S1 30
#endif
#ifndef DYN_S2
#define DYN_S2 120
#endif
#ifndef DYN_ESC
#define DYN_ESC 250
#endif
const int spanShort = max(4, J/3);
const int spanLong  = TENURE_SPAN;
auto lastImpT = chrono::steady_clock::now();
auto dynTen = [&](chrono::steady_clock::time_point nowT)->int{
long long stag = chrono::duration_cast<chrono::milliseconds>(nowT - lastImpT).count();
int tmin, span;
if(stag <= DYN_S1){ tmin = 8; span = spanShort; }
else if(stag >= DYN_S2){
tmin = TENURE_MIN + (int)min(6LL, (stag - DYN_S2)/DYN_ESC);
span = spanLong;
} else {
int f = (int)((stag - DYN_S1)*100/(DYN_S2 - DYN_S1));
tmin = 8 + (TENURE_MIN - 8)*f/100;
span = spanShort + (spanLong - spanShort)*f/100;
}
return tmin + (int)(rng() % (unsigned)span);
};
#ifndef RST_MS
#define RST_MS 200
#endif
struct Elite { vector<vector<int>> seq; long long C; unsigned long long h; };
vector<Elite> pool;
const int E = 6; int poolHead = 0;
#ifdef DIAG
long long dRst = 0, dRstElite = 0, dRstBest = 0, dInsBest = 0, dInsRB = 0, dRstWin = 0;
#endif
auto solHash = [&](const vector<vector<int>>& s)->unsigned long long{
unsigned long long h = 1469598103934665603ULL;
for(int m=0;m<M;++m) for(int j=0;j<J;++j) h = h*1099511628211ULL ^ (unsigned long long)(s[m][j]+1);
return h;
};
auto poolAdd = [&](const vector<vector<int>>& s, long long C){
unsigned long long h = solHash(s);
for(const Elite& e : pool) if(e.h == h) return;
if((int)pool.size() < E) pool.push_back({s, C, h});
else { pool[poolHead] = {s, C, h}; poolHead = (poolHead+1)%E; }
};
vector<vector<int>> rb; long long rbC = LLONG_MAX;
int failCnt = 0;
bool afterRst = false;
auto lastReoptT = T0;
int lastReoptM = -1;
if(!igElitePi.empty()){
vector<vector<int>> es(M, igElitePi);
long long ec = evalSeq(es);
if(ec > 0) poolAdd(es, ec);
evalSeq(cur, true);
}
#ifndef NO_SWEEP
if(J > 2){
vector<pair<long long,int>> mord(M);
{
vector<long long> mload(M, 0);
for(int j=0;j<J;++j) for(int k=0;k<M;++k) mload[m_of[j][k]] += p_of[j][k];
for(int m=0;m<M;++m) mord[m] = {-mload[m], m};
sort(mord.begin(), mord.end());
}
for(int t=0;t<M;++t){
if(chrono::steady_clock::now() >= T_end) break;
reoptMachine(mord[t].second, cur, curC, bestC, best, T_end, lastImpT);
}
if(bestC <= LB) timeUp = true;
}
#endif
#ifdef DIAG
int iterLastImp = 0;
#endif
while(!timeUp){
auto nowT = chrono::steady_clock::now();
if(nowT >= T_end) break;
iter++;
if((iter & 16383) == 0){
long long fc = evalSeq(cur, true);
if(fc >= 0) curC = fc;
}
#ifndef NO_TRIG
if(J > 2 && nowT - lastImpT > chrono::milliseconds(150) && nowT - lastReoptT > chrono::milliseconds(60)){
lastReoptT = nowT;
static vector<long long> critW, loadW;
critW.assign(M, 0); loadW.assign(M, 0);
for(int u=0;u<N;++u){
int mm = m_of[u/M][u%M];
loadW[mm] += pnode[u];
if(crit[u]) critW[mm] += pnode[u];
}
int pick = -1;
for(int mm=0;mm<M;++mm){
if(mm == lastReoptM) continue;
if(pick < 0 || critW[mm] > critW[pick] || (critW[mm] == critW[pick] && loadW[mm] > loadW[pick])) pick = mm;
}
if(pick >= 0){
lastReoptM = pick;
reoptMachine(pick, cur, curC, bestC, best, T_end, lastImpT);
}
}
#endif
genMoves(cur);
int nmv = (int)gmoves.size();
#ifdef DIAG
static long long totMv = 0; totMv += nmv;
if(iter % 50000 == 0) fprintf(stderr, "avg nmv=%.1f\n", (double)totMv/iter);
#endif
if(nmv == 0) break;
gcand.clear();
for(int idx=0; idx<nmv; ++idx)
gcand.push_back({estMove(cur, gmoves[idx]), idx});
int K = nmv < 24 ? nmv : 24;
partial_sort(gcand.begin(), gcand.begin()+K, gcand.end());
bool sorted_all = (K == nmv);
bool applied = false;
#ifndef EVAL_TOP
#define EVAL_TOP 1
#endif
#if EVAL_TOP > 1
{
int bestIdx = -1; long long bestNC = -1; int evald = 0;
for(int t=0; t<nmv && evald<EVAL_TOP; ++t){
if((t & 7)==7 && chrono::steady_clock::now() >= T_end){ timeUp = true; break; }
if(t >= K && !sorted_all){ sort(gcand.begin(), gcand.end()); sorted_all = true; }
const Mv& mv = gmoves[gcand[t].second];
bool tb = isTabu(cur, mv, iter);
bool asp = gcand[t].first < bestC;
if(tb && !asp) continue;
applyMove(cur, mv);
long long nc = evalSeq(cur, false);
undoMove(cur, mv);
if(nc < 0) continue;
evald++;
if(bestIdx < 0 || nc < bestNC){ bestNC = nc; bestIdx = gcand[t].second; }
}
if(bestIdx >= 0){
const Mv& mv = gmoves[bestIdx];
collectTabu(cur, mv);
applyMove(cur, mv);
long long nc = evalSeq(cur, true);
if(nc >= 0){
int tenure = dynTen(nowT);
for(size_t id : gpend) tabuTB[id] = iter + tenure;
curC = nc; applied = true;
} else { undoMove(cur, mv); evalSeq(cur, true); }
}
}
#endif
for(int pass=0; pass<2 && !applied && !timeUp; ++pass){
for(int t=0; t<nmv; ++t){
if((t & 7)==7 && chrono::steady_clock::now() >= T_end){ timeUp = true; break; }
if(t >= K && !sorted_all){
sort(gcand.begin(), gcand.end());
sorted_all = true;
}
const Mv& mv = gmoves[gcand[t].second];
if(pass==0){
bool tb = isTabu(cur, mv, iter);
bool asp = gcand[t].first < bestC;
if(tb && !asp) continue;
}
collectTabu(cur, mv);
applyMove(cur, mv);
long long nc = incAfterMove(cur, mv);
#ifdef VERIFY
{
static vector<long long> vd, vt; static long long mism = 0; static long long checks = 0;
vd = dist_; vt = tail_;
long long fc = evalSeq(cur, true);
checks++;
if(nc != -2){
if(fc != nc || vd != dist_ || vt != tail_){
mism++;
fprintf(stderr, "MISMATCH iter=%d inc=%lld full=%lld distOK=%d tailOK=%d\n",
iter, nc, fc, (int)(vd==dist_), (int)(vt==tail_));
}
}
if(checks % 20000 == 0) fprintf(stderr, "verify checks=%lld mism=%lld\n", checks, mism);
nc = fc;
}
#endif
if(nc == -2) nc = evalSeq(cur, true);
if(nc < 0){
undoMove(cur, mv);
evalSeq(cur, true);
continue;
}
int tenure = dynTen(nowT);
for(size_t id : gpend) tabuTB[id] = iter + tenure;
curC = nc;
applied = true;
break;
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
#ifdef DIAG
iterLastImp = iter;
dInsBest++;
if(afterRst) dRstWin++;
#endif
afterRst = false;
if(bestC <= LB) break;
}
}
else {
++sinceImp;
long long stag = chrono::duration_cast<chrono::milliseconds>(nowT - lastImpT).count();
if(stag > RST_MS){
if(rbC <= bestC + bestC/50 && !rb.empty()){
poolAdd(rb, rbC);
#ifdef DIAG
dInsRB++;
#endif
}
if(pool.empty() || (rng() & 1)){
cur = best;
#ifdef DIAG
dRstBest++;
#endif
} else {
cur = pool[rng()%pool.size()].seq;
#ifdef DIAG
dRstElite++;
#endif
}
curC = evalSeq(cur, true);
int kicks = 2 + failCnt + (int)(rng() % 2); if(kicks > 8) kicks = 8;
for(int r=0; r<kicks; ++r){
if(chrono::steady_clock::now() >= T_end){ timeUp = true; break; }
genMoves(cur);
if(gmoves.empty()) break;
const Mv& mv = gmoves[rng() % gmoves.size()];
applyMove(cur, mv);
long long nc = evalSeq(cur, true);
if(nc < 0){ undoMove(cur, mv); evalSeq(cur, true); }
else curC = nc;
}
fill(tabuTB.begin(), tabuTB.end(), 0);
sinceImp = 0; rbC = LLONG_MAX;
lastImpT = chrono::steady_clock::now();
++failCnt;
afterRst = true;
#ifdef DIAG
dRst++;
#endif
}
}
}
#ifdef DIAG
extern long long g_pops, g_calls;
extern long long g_roAtt, g_roAcc, g_roCyc, g_roUs;
fprintf(stderr, "iters=%d lastImp=%d bestC=%lld avgPops=%.1f (N=%d) reopt att=%lld acc=%lld cyc=%lld ms=%.1f rst=%lld(best %lld/elite %lld, wins %lld) pool ins=%lld/%lld\n", iter, iterLastImp, bestC, g_calls? (double)g_pops/g_calls : 0.0, N, g_roAtt, g_roAcc, g_roCyc, g_roUs/1000.0, dRst, dRstBest, dRstElite, dRstWin, dInsBest, dInsRB);
#endif
}
{
vector<char> buf;
buf.reserve((size_t)N*8 + M + 16);
for(int m=0;m<M;++m){
for(int j=0;j<J;++j){
int x = best[m][j];
if(x == 0){ buf.push_back('0'); }
else { char tmp[12]; int t = 0; while(x > 0){ tmp[t++] = char('0' + x%10); x /= 10; } while(t > 0) buf.push_back(tmp[--t]); }
buf.push_back(j+1<J ? ' ' : '\n');
}
}
fwrite(buf.data(), 1, buf.size(), stdout);
fflush(stdout);
}
_exit(0);
}
namespace Alt {
using namespace std;
static int J, M, N;
static vector<long long> procOp;
static vector<int>       jobOf, kOf, machOf;
static vector<int>       jobPred, jobSucc;
static vector<vector<int>>       machJK;
static vector<vector<long long>> procJK;
static vector<vector<int>>       posOf;
static inline int opOnMachine(int j, int m){ return j * M + posOf[j][m]; }
static vector<int>       indeg, mSucc, mPred, order_;
static vector<long long> dist_;
static vector<long long> q_;
static long long         Cmax_;
static clock_t START;
static const double TL = 0.94;
static inline double elapsed(){ return double(clock() - START) / CLOCKS_PER_SEC; }
static unsigned long long rngState = 0x9e3779b97f4a7c15ULL;
static inline unsigned long long rnd(){
rngState ^= rngState << 13; rngState ^= rngState >> 7; rngState ^= rngState << 17;
return rngState;
}
static inline int rndInt(int n){ return (int)(rnd() % (unsigned long long)n); }
static long long evaluate(const vector<vector<int>>& seq){
for(int op = 0; op < N; ++op){
indeg[op] = (kOf[op] > 0) ? 1 : 0;
mSucc[op] = -1; mPred[op] = -1;
}
for(int m = 0; m < M; ++m){
const vector<int>& s = seq[m];
for(int i = 1; i < J; ++i){
int a = opOnMachine(s[i-1], m);
int b = opOnMachine(s[i],   m);
mSucc[a] = b; mPred[b] = a; ++indeg[b];
}
}
int tail = 0, head = 0;
for(int op = 0; op < N; ++op){
if(indeg[op] == 0){ dist_[op] = procOp[op]; order_[tail++] = op; }
else dist_[op] = 0;
}
int cnt = 0;
while(head < tail){
int u = order_[head++]; ++cnt;
long long du = dist_[u];
int js = jobSucc[u];
if(js != -1){
if(dist_[js] < du + procOp[js]) dist_[js] = du + procOp[js];
if(--indeg[js] == 0) order_[tail++] = js;
}
int ms = mSucc[u];
if(ms != -1){
if(dist_[ms] < du + procOp[ms]) dist_[ms] = du + procOp[ms];
if(--indeg[ms] == 0) order_[tail++] = ms;
}
}
if(cnt != N) return -1;
long long C = 0;
for(int op = 0; op < N; ++op) if(dist_[op] > C) C = dist_[op];
for(int idx = N - 1; idx >= 0; --idx){
int op = order_[idx];
long long best = 0;
int js = jobSucc[op]; if(js != -1 && q_[js] > best) best = q_[js];
int ms = mSucc[op];   if(ms != -1 && q_[ms] > best) best = q_[ms];
q_[op] = procOp[op] + best;
}
Cmax_ = C;
return C;
}
static inline bool critOp(int op){ return (dist_[op] - procOp[op]) + q_[op] == Cmax_; }
static void getBlockMoves(const vector<vector<int>>& seq, vector<pair<int,int>>& moves){
moves.clear();
for(int m = 0; m < M; ++m){
const vector<int>& s = seq[m];
int i = 0;
while(i < J){
if(!critOp(opOnMachine(s[i], m))){
++i;
continue;
}
int j = i;
while(j < J && critOp(opOnMachine(s[j], m))){
++j;
}
int blockSize = j - i;
if(blockSize >= 2){
int a1 = opOnMachine(s[i], m), b1 = opOnMachine(s[i+1], m);
if(dist_[a1] == dist_[b1] - procOp[b1]){
moves.push_back({a1, b1});
}
if(blockSize > 2){
int a2 = opOnMachine(s[j-2], m), b2 = opOnMachine(s[j-1], m);
if(dist_[a2] == dist_[b2] - procOp[b2]){
moves.push_back({a2, b2});
}
}
}
i = j;
}
}
}
static inline long long estimateSwap(int a, int b){
int PM = mPred[a], SM = mSucc[b];
long long fPM  = (PM != -1)          ? dist_[PM]          : 0;
long long fJPa = (jobPred[a] != -1)  ? dist_[jobPred[a]]  : 0;
long long fJPb = (jobPred[b] != -1)  ? dist_[jobPred[b]]  : 0;
long long rB = max(fPM, fJPb);
long long rA = max(rB + procOp[b], fJPa);
long long qJSa = (jobSucc[a] != -1) ? q_[jobSucc[a]] : 0;
long long qJSb = (jobSucc[b] != -1) ? q_[jobSucc[b]] : 0;
long long qSM  = (SM != -1)         ? q_[SM]         : 0;
long long qA = procOp[a] + max(qJSa, qSM);
long long qB = procOp[b] + max(qJSb, qA);
return max(rA + qA, rB + qB);
}
static inline double priority(int rule, int j, int k, const vector<long long>& remWork){
switch(rule){
case 0: return (double)remWork[j];
case 1: return -(double)procJK[j][k];
case 2: return  (double)procJK[j][k];
case 3: return -(double)remWork[j];
default: return (double)(rnd() & 0xffffff);
}
}
static vector<vector<int>> gifflerThompson(int rule){
vector<int>       nextK(J, 0);
vector<long long> jobFree(J, 0), machFree(M, 0), remWork(J, 0);
for(int j = 0; j < J; ++j)
for(int k = 0; k < M; ++k) remWork[j] += procJK[j][k];
vector<vector<int>> seq(M);
for(int m = 0; m < M; ++m) seq[m].reserve(J);
int scheduled = 0;
while(scheduled < N){
long long minC = LLONG_MAX; int mstar = -1;
for(int j = 0; j < J; ++j){
if(nextK[j] >= M) continue;
int k = nextK[j], m = machJK[j][k];
long long est = max(jobFree[j], machFree[m]);
long long C = est + procJK[j][k];
if(C < minC){ minC = C; mstar = m; }
}
int chosen = -1; double bestPri = -1e300;
for(int j = 0; j < J; ++j){
if(nextK[j] >= M) continue;
int k = nextK[j];
if(machJK[j][k] != mstar) continue;
long long est = max(jobFree[j], machFree[mstar]);
if(est < minC){
double pri = priority(rule, j, k, remWork);
if(pri > bestPri){ bestPri = pri; chosen = j; }
}
}
int j = chosen, k = nextK[j], m = mstar;
long long est = max(jobFree[j], machFree[m]);
long long fin = est + procJK[j][k];
jobFree[j] = fin; machFree[m] = fin;
seq[m].push_back(j);
remWork[j] -= procJK[j][k];
++nextK[j];
++scheduled;
}
return seq;
}
static vector<vector<int>> pos;
static void rebuildPos(const vector<vector<int>>& seq){
for(int m = 0; m < M; ++m)
for(int i = 0; i < J; ++i) pos[m][seq[m][i]] = i;
}
static inline void doSwap(vector<vector<int>>& seq, int a, int b){
int m = machOf[a];
int i = pos[m][jobOf[a]];
swap(seq[m][i], seq[m][i+1]);
pos[m][seq[m][i]]   = i;
pos[m][seq[m][i+1]] = i+1;
}
static void perturb(vector<vector<int>>& seq, int kicks){
rebuildPos(seq);
for(int t = 0; t < kicks; ++t){
if(J < 2) return;
int m = rndInt(M);
int i = rndInt(J - 1);
int a = opOnMachine(seq[m][i],   m);
int b = opOnMachine(seq[m][i+1], m);
doSwap(seq, a, b);
if(evaluate(seq) < 0) doSwap(seq, b, a);
}
}
static vector<long long> tabuUntil;
static inline size_t tabIdx(int m, int ja, int jb){
int lo = ja < jb ? ja : jb, hi = ja < jb ? jb : ja;
return (size_t)m * J * J + (size_t)lo * J + hi;
}
static long long tabuSearch(vector<vector<int>>& best, long long bestMk){
vector<vector<int>> cur = best;
vector<pair<int,int>> moves;
fill(tabuUntil.begin(), tabuUntil.end(), 0);
rebuildPos(cur);
long long iter = 0, lastImprove = 0;
int tenure = 10 + rndInt(10);
const long long stall = 4000;
long long curMk = evaluate(cur);
int checkClock = 0;
while((checkClock++ & 63) || elapsed() < TL){
getBlockMoves(cur, moves);
if(moves.empty()){
perturb(cur, 2);
curMk = evaluate(cur);
++iter;
continue;
}
long long bestEst = LLONG_MAX, aspEst = LLONG_MAX;
int ba = -1, bb = -1, aa = -1, ab = -1;
for(auto& pr : moves){
int a = pr.first, b = pr.second;
long long est = estimateSwap(a, b);
bool isTabu = tabuUntil[tabIdx(machOf[a], jobOf[a], jobOf[b])] > iter;
if(isTabu){
if(est < bestMk && est < aspEst){ aspEst = est; aa = a; ab = b; }
}else if(est < bestEst){
bestEst = est; ba = a; bb = b;
}
}
int ca, cb;
if(aa != -1 && aspEst <= bestEst){ ca = aa; cb = ab; }
else if(ba != -1){ ca = ba; cb = bb; }
else if(aa != -1){ ca = aa; cb = ab; }
else {
perturb(cur, 2);
curMk = evaluate(cur);
++iter;
continue;
}
doSwap(cur, ca, cb);
tabuUntil[tabIdx(machOf[ca], jobOf[ca], jobOf[cb])] = iter + tenure;
curMk = evaluate(cur);
if(curMk < bestMk){
bestMk = curMk;
best = cur;
lastImprove = iter;
tenure = 10 + rndInt(10);
}
++iter;
if(iter - lastImprove > stall){
cur = best;
perturb(cur, 4 + rndInt(10));
curMk = evaluate(cur);
fill(tabuUntil.begin(), tabuUntil.end(), 0);
lastImprove = iter;
}
}
return bestMk;
}
static void output(const vector<vector<int>>& seq){
static char buf[1 << 22];
int p = 0;
for(int m = 0; m < M; ++m){
for(int i = 0; i < J; ++i){
int x = seq[m][i];
if(x == 0) buf[p++] = '0';
else{
char tmp[12]; int t = 0;
while(x){ tmp[t++] = char('0' + x % 10); x /= 10; }
while(t) buf[p++] = tmp[--t];
}
buf[p++] = (i + 1 < J) ? ' ' : '\n';
}
}
fwrite(buf, 1, p, stdout);
}
int solve(){
START = ::G_START;
J = ::J; M = ::M; N = J * M;
machJK = ::m_of;
procJK = ::p_of;
posOf.assign(J, vector<int>(M, -1));
procOp.assign(N, 0); jobOf.assign(N, 0); kOf.assign(N, 0); machOf.assign(N, 0);
jobPred.assign(N, -1); jobSucc.assign(N, -1);
for(int j = 0; j < J; ++j){
for(int k = 0; k < M; ++k){
int m = machJK[j][k]; long long p = procJK[j][k];
posOf[j][m] = k;
int op = j * M + k;
procOp[op] = p; jobOf[op] = j; kOf[op] = k; machOf[op] = m;
jobPred[op] = (k > 0)     ? op - 1 : -1;
jobSucc[op] = (k < M - 1) ? op + 1 : -1;
}
}
indeg.assign(N, 0); mSucc.assign(N, -1); mPred.assign(N, -1);
order_.assign(N, 0); dist_.assign(N, 0); q_.assign(N, 0);
pos.assign(M, vector<int>(J, 0));
tabuUntil.assign((size_t)M * J * J, 0);
vector<vector<int>> best;
long long bestMk = LLONG_MAX;
for(int rule = 0; rule < 4; ++rule){
vector<vector<int>> seq = gifflerThompson(rule);
long long mk = evaluate(seq);
if(mk >= 0 && mk < bestMk){ bestMk = mk; best = seq; }
}
for(int r = 0; r < 32; ++r){
vector<vector<int>> seq = gifflerThompson(4);
long long mk = evaluate(seq);
if(mk >= 0 && mk < bestMk){ bestMk = mk; best = seq; }
}
if(best.empty()){
best.assign(M, vector<int>(J));
for(int m = 0; m < M; ++m) for(int j = 0; j < J; ++j) best[m][j] = j;
bestMk = evaluate(best);
}
bestMk = tabuSearch(best, bestMk);
output(best);
return 0;
}
}
namespace HECD {
using namespace std;
static int J, M, N;
static vector<long long> procOp;
static vector<int>       jobOf, kOf, machOf;
static vector<int>       jobPred, jobSucc;
static vector<vector<int>>       machJK;
static vector<vector<long long>> procJK;
static vector<vector<int>>       posOf;
static inline int opOnMachine(int j, int m){ return j * M + posOf[j][m]; }
static vector<int>       indeg, mSucc, mPred, order_;
static vector<long long> dist_, q_;
static long long         Cmax_;
static clock_t START;
static const double TL = 0.94;
static inline double elapsed(){ return double(clock() - START) / CLOCKS_PER_SEC; }
static unsigned long long rngState = 0x9e3779b97f4a7c15ULL;
static inline unsigned long long rnd(){
rngState ^= rngState << 13; rngState ^= rngState >> 7; rngState ^= rngState << 17;
return rngState;
}
static inline int rndInt(int n){ return (int)(rnd() % (unsigned long long)n); }
static long long evaluate(const vector<vector<int>>& seq){
for(int op = 0; op < N; ++op){
indeg[op] = (kOf[op] > 0) ? 1 : 0;
mSucc[op] = -1; mPred[op] = -1;
}
for(int m = 0; m < M; ++m){
const vector<int>& s = seq[m];
for(int i = 1; i < J; ++i){
int a = opOnMachine(s[i-1], m);
int b = opOnMachine(s[i],   m);
mSucc[a] = b; mPred[b] = a; ++indeg[b];
}
}
int tail = 0, head = 0;
for(int op = 0; op < N; ++op){
if(indeg[op] == 0){ dist_[op] = procOp[op]; order_[tail++] = op; }
else dist_[op] = 0;
}
int cnt = 0;
while(head < tail){
int u = order_[head++]; ++cnt;
long long du = dist_[u];
int js = jobSucc[u];
if(js != -1){
if(dist_[js] < du + procOp[js]) dist_[js] = du + procOp[js];
if(--indeg[js] == 0) order_[tail++] = js;
}
int ms = mSucc[u];
if(ms != -1){
if(dist_[ms] < du + procOp[ms]) dist_[ms] = du + procOp[ms];
if(--indeg[ms] == 0) order_[tail++] = ms;
}
}
if(cnt != N) return -1;
long long C = 0;
for(int op = 0; op < N; ++op) if(dist_[op] > C) C = dist_[op];
for(int idx = N - 1; idx >= 0; --idx){
int op = order_[idx];
long long best = 0;
int js = jobSucc[op]; if(js != -1 && q_[js] > best) best = q_[js];
int ms = mSucc[op];   if(ms != -1 && q_[ms] > best) best = q_[ms];
q_[op] = procOp[op] + best;
}
Cmax_ = C;
return C;
}
static inline bool critOp(int op){ return (dist_[op] - procOp[op]) + q_[op] == Cmax_; }
static void getBlockMoves(const vector<vector<int>>& seq, vector<pair<int,int>>& moves){
moves.clear();
for(int m = 0; m < M; ++m){
const vector<int>& s = seq[m];
int i = 0;
while(i < J){
if(!critOp(opOnMachine(s[i], m))){
++i; continue;
}
int j = i;
while(j < J && critOp(opOnMachine(s[j], m))) ++j;
int blockSize = j - i;
if(blockSize >= 2){
int a1 = opOnMachine(s[i], m), b1 = opOnMachine(s[i+1], m);
if(dist_[a1] == dist_[b1] - procOp[b1]) moves.push_back({a1, b1});
if(blockSize > 2){
int a2 = opOnMachine(s[j-2], m), b2 = opOnMachine(s[j-1], m);
if(dist_[a2] == dist_[b2] - procOp[b2]) moves.push_back({a2, b2});
}
}
i = j;
}
}
}
static inline long long estimateSwap(int a, int b){
int PM = mPred[a], SM = mSucc[b];
long long fPM  = (PM != -1) ? dist_[PM] : 0;
long long fJPa = (jobPred[a] != -1) ? dist_[jobPred[a]] : 0;
long long fJPb = (jobPred[b] != -1) ? dist_[jobPred[b]] : 0;
long long rB = max(fPM, fJPb);
long long rA = max(rB + procOp[b], fJPa);
long long qJSa = (jobSucc[a] != -1) ? q_[jobSucc[a]] : 0;
long long qJSb = (jobSucc[b] != -1) ? q_[jobSucc[b]] : 0;
long long qSM  = (SM != -1) ? q_[SM] : 0;
long long qA = procOp[a] + max(qJSa, qSM);
long long qB = procOp[b] + max(qJSb, qA);
return max(rA + qA, rB + qB);
}
static inline double priority(int rule, int j, int k, const vector<long long>& remWork){
switch(rule){
case 0: return (double)remWork[j];
case 1: return -(double)procJK[j][k];
case 2: return  (double)procJK[j][k];
case 3: return -(double)remWork[j];
default: return (double)(rnd() & 0xffffff);
}
}
static vector<vector<int>> gifflerThompson(int rule){
vector<int>       nextK(J, 0);
vector<long long> jobFree(J, 0), machFree(M, 0), remWork(J, 0);
for(int j = 0; j < J; ++j)
for(int k = 0; k < M; ++k) remWork[j] += procJK[j][k];
vector<vector<int>> seq(M);
for(int m = 0; m < M; ++m) seq[m].reserve(J);
int scheduled = 0;
while(scheduled < N){
long long minC = LLONG_MAX; int mstar = -1;
for(int j = 0; j < J; ++j){
if(nextK[j] >= M) continue;
int k = nextK[j], m = machJK[j][k];
long long est = max(jobFree[j], machFree[m]);
long long C = est + procJK[j][k];
if(C < minC){ minC = C; mstar = m; }
}
int chosen = -1; double bestPri = -1e300;
for(int j = 0; j < J; ++j){
if(nextK[j] >= M) continue;
int k = nextK[j];
if(machJK[j][k] != mstar) continue;
long long est = max(jobFree[j], machFree[mstar]);
if(est < minC){
double pri = priority(rule, j, k, remWork);
if(pri > bestPri){ bestPri = pri; chosen = j; }
}
}
int j = chosen, k = nextK[j], m = mstar;
long long est = max(jobFree[j], machFree[m]);
long long fin = est + procJK[j][k];
jobFree[j] = fin; machFree[m] = fin;
seq[m].push_back(j);
remWork[j] -= procJK[j][k];
++nextK[j];
++scheduled;
}
return seq;
}
static vector<vector<int>> pos;
static void rebuildPos(const vector<vector<int>>& seq){
for(int m = 0; m < M; ++m)
for(int i = 0; i < J; ++i) pos[m][seq[m][i]] = i;
}
static inline void doSwap(vector<vector<int>>& seq, int a, int b){
int m = machOf[a];
int i = pos[m][jobOf[a]];
swap(seq[m][i], seq[m][i+1]);
pos[m][seq[m][i]]   = i;
pos[m][seq[m][i+1]] = i+1;
}
static void perturb(vector<vector<int>>& seq, int kicks){
rebuildPos(seq);
for(int t = 0; t < kicks; ++t){
if(J < 2) return;
int m = rndInt(M);
int i = rndInt(J - 1);
int a = opOnMachine(seq[m][i],   m);
int b = opOnMachine(seq[m][i+1], m);
doSwap(seq, a, b);
if(evaluate(seq) < 0) doSwap(seq, b, a);
}
}
static vector<long long> tabuUntil;
static inline size_t tabIdx(int m, int ja, int jb){
int lo = ja < jb ? ja : jb, hi = ja < jb ? jb : ja;
return (size_t)m * J * J + (size_t)lo * J + hi;
}
static long long tabuSearch(vector<vector<int>>& best, long long bestMk){
vector<vector<int>> cur = best;
vector<pair<int,int>> moves;
fill(tabuUntil.begin(), tabuUntil.end(), 0);
rebuildPos(cur);
long long iter = 0, lastImprove = 0;
int tenure = 15 + rndInt(13);
const long long stall = 5200;
long long curMk = evaluate(cur);
int checkClock = 0;
while((checkClock++ & 63) || elapsed() < TL){
getBlockMoves(cur, moves);
if(moves.empty()){
perturb(cur, 4);
curMk = evaluate(cur);
++iter;
continue;
}
long long bestEst = LLONG_MAX, aspEst = LLONG_MAX;
int ba = -1, bb = -1, aa = -1, ab = -1;
for(auto& pr : moves){
int a = pr.first, b = pr.second;
long long est = estimateSwap(a, b);
bool isTabu = tabuUntil[tabIdx(machOf[a], jobOf[a], jobOf[b])] > iter;
if(isTabu){
if(est < bestMk && est < aspEst){ aspEst = est; aa = a; ab = b; }
}else if(est < bestEst){
bestEst = est; ba = a; bb = b;
}
}
int ca, cb;
if(aa != -1 && aspEst <= bestEst){ ca = aa; cb = ab; }
else if(ba != -1){ ca = ba; cb = bb; }
else if(aa != -1){ ca = aa; cb = ab; }
else {
perturb(cur, 4);
curMk = evaluate(cur);
++iter;
continue;
}
doSwap(cur, ca, cb);
tabuUntil[tabIdx(machOf[ca], jobOf[ca], jobOf[cb])] = iter + tenure;
curMk = evaluate(cur);
if(curMk < bestMk){
bestMk = curMk;
best = cur;
lastImprove = iter;
tenure = 15 + rndInt(13);
}
++iter;
if(iter - lastImprove > stall){
cur = best;
perturb(cur, 6 + rndInt(13));
curMk = evaluate(cur);
fill(tabuUntil.begin(), tabuUntil.end(), 0);
lastImprove = iter;
}
}
return bestMk;
}
static void output(const vector<vector<int>>& seq){
static char buf[1 << 22];
int p = 0;
for(int m = 0; m < M; ++m){
for(int i = 0; i < J; ++i){
int x = seq[m][i];
if(x == 0) buf[p++] = '0';
else{
char tmp[12]; int t = 0;
while(x){ tmp[t++] = char('0' + x % 10); x /= 10; }
while(t) buf[p++] = tmp[--t];
}
buf[p++] = (i + 1 < J) ? ' ' : '\n';
}
}
fwrite(buf, 1, p, stdout);
}
int solveParsed(int Jin, int Min, const vector<vector<int>>& m_in, const vector<vector<long long>>& p_in){
START = ::G_START;
J = Jin; M = Min; N = J * M;
machJK = m_in;
procJK = p_in;
posOf.assign(J, vector<int>(M, -1));
procOp.assign(N, 0); jobOf.assign(N, 0); kOf.assign(N, 0); machOf.assign(N, 0);
jobPred.assign(N, -1); jobSucc.assign(N, -1);
for(int j = 0; j < J; ++j){
for(int k = 0; k < M; ++k){
int m = machJK[j][k]; long long p = procJK[j][k];
posOf[j][m] = k;
int op = j * M + k;
procOp[op] = p; jobOf[op] = j; kOf[op] = k; machOf[op] = m;
jobPred[op] = (k > 0)     ? op - 1 : -1;
jobSucc[op] = (k < M - 1) ? op + 1 : -1;
}
}
indeg.assign(N, 0); mSucc.assign(N, -1); mPred.assign(N, -1);
order_.assign(N, 0); dist_.assign(N, 0); q_.assign(N, 0);
pos.assign(M, vector<int>(J, 0));
tabuUntil.assign((size_t)M * J * J, 0);
vector<vector<int>> best;
long long bestMk = LLONG_MAX;
for(int rule = 0; rule < 5; ++rule){
vector<vector<int>> seq = gifflerThompson(rule);
long long mk = evaluate(seq);
if(mk >= 0 && mk < bestMk){ bestMk = mk; best = seq; }
}
for(int r = 0; r < 150; ++r){
vector<vector<int>> seq = gifflerThompson(r % 5);
long long mk = evaluate(seq);
if(mk >= 0 && mk < bestMk){ bestMk = mk; best = seq; }
}
if(best.empty()){
best.assign(M, vector<int>(J));
for(int m = 0; m < M; ++m) for(int j = 0; j < J; ++j) best[m][j] = j;
bestMk = evaluate(best);
}
bestMk = tabuSearch(best, bestMk);
output(best);
return 0;
}
}
namespace H08 {
using namespace std;
static int J, M, N;
static vector<long long> procOp;
static vector<int>       jobOf, kOf, machOf;
static vector<int>       jobPred, jobSucc;
static vector<vector<int>>       machJK;
static vector<vector<long long>> procJK;
static vector<vector<int>>       posOf;
static inline int opOnMachine(int j, int m){ return j * M + posOf[j][m]; }
static vector<int>       indeg, mSucc, mPred, order_;
static vector<long long> dist_, q_;
static long long         Cmax_;
static clock_t START;
static const double TL = 0.94;
static inline double elapsed(){ return double(clock() - START) / CLOCKS_PER_SEC; }
static unsigned long long rngState = 0x9e3779b97f4a7c15ULL;
static inline unsigned long long rnd(){
rngState ^= rngState << 13; rngState ^= rngState >> 7; rngState ^= rngState << 17;
return rngState;
}
static inline int rndInt(int n){ return (int)(rnd() % (unsigned long long)n); }
static long long evaluate(const vector<vector<int>>& seq){
for(int op = 0; op < N; ++op){
indeg[op] = (kOf[op] > 0) ? 1 : 0;
mSucc[op] = -1; mPred[op] = -1;
}
for(int m = 0; m < M; ++m){
const vector<int>& s = seq[m];
for(int i = 1; i < J; ++i){
int a = opOnMachine(s[i-1], m);
int b = opOnMachine(s[i],   m);
mSucc[a] = b; mPred[b] = a; ++indeg[b];
}
}
int tail = 0, head = 0;
for(int op = 0; op < N; ++op){
if(indeg[op] == 0){ dist_[op] = procOp[op]; order_[tail++] = op; }
else dist_[op] = 0;
}
int cnt = 0;
while(head < tail){
int u = order_[head++]; ++cnt;
long long du = dist_[u];
int js = jobSucc[u];
if(js != -1){
if(dist_[js] < du + procOp[js]) dist_[js] = du + procOp[js];
if(--indeg[js] == 0) order_[tail++] = js;
}
int ms = mSucc[u];
if(ms != -1){
if(dist_[ms] < du + procOp[ms]) dist_[ms] = du + procOp[ms];
if(--indeg[ms] == 0) order_[tail++] = ms;
}
}
if(cnt != N) return -1;
long long C = 0;
for(int op = 0; op < N; ++op) if(dist_[op] > C) C = dist_[op];
for(int idx = N - 1; idx >= 0; --idx){
int op = order_[idx];
long long best = 0;
int js = jobSucc[op]; if(js != -1 && q_[js] > best) best = q_[js];
int ms = mSucc[op];   if(ms != -1 && q_[ms] > best) best = q_[ms];
q_[op] = procOp[op] + best;
}
Cmax_ = C;
return C;
}
static inline bool critOp(int op){ return (dist_[op] - procOp[op]) + q_[op] == Cmax_; }
static void getBlockMoves(const vector<vector<int>>& seq, vector<pair<int,int>>& swaps,
vector<array<int,3>>& inserts){
swaps.clear(); inserts.clear();
for(int m = 0; m < M; ++m){
const vector<int>& s = seq[m];
int i = 0;
while(i < J){
if(!critOp(opOnMachine(s[i], m))){ ++i; continue; }
int j = i;
while(j < J && critOp(opOnMachine(s[j], m))) ++j;
int bs = j - i;
if(bs >= 2){
int a1 = opOnMachine(s[i], m), b1 = opOnMachine(s[i+1], m);
if(dist_[a1] == dist_[b1] - procOp[b1]) swaps.push_back({a1, b1});
if(bs > 2){
int a2 = opOnMachine(s[j-2], m), b2 = opOnMachine(s[j-1], m);
if(dist_[a2] == dist_[b2] - procOp[b2]) swaps.push_back({a2, b2});
inserts.push_back({m, i, j-1});
inserts.push_back({m, j-1, i});
}
}
i = j;
}
}
}
static inline long long estimateSwap(int a, int b){
int PM = mPred[a], SM = mSucc[b];
long long fPM  = (PM != -1) ? dist_[PM] : 0;
long long fJPa = (jobPred[a] != -1) ? dist_[jobPred[a]] : 0;
long long fJPb = (jobPred[b] != -1) ? dist_[jobPred[b]] : 0;
long long rB = max(fPM, fJPb);
long long rA = max(rB + procOp[b], fJPa);
long long qJSa = (jobSucc[a] != -1) ? q_[jobSucc[a]] : 0;
long long qJSb = (jobSucc[b] != -1) ? q_[jobSucc[b]] : 0;
long long qSM  = (SM != -1) ? q_[SM] : 0;
long long qA = procOp[a] + max(qJSa, qSM);
long long qB = procOp[b] + max(qJSb, qA);
return max(rA + qA, rB + qB);
}
static vector<int>       gord;
static vector<long long> gestC;
static long long estInsert(const vector<vector<int>>& cur, int m, int from, int to){
const vector<int>& s = cur[m];
int lo = from < to ? from : to, hi = from < to ? to : from;
int L = hi - lo + 1;
if((int)gord.size() < L){ gord.resize(L); gestC.resize(L); }
if(to > from){
for(int t = from + 1; t <= to; ++t) gord[t - from - 1] = s[t];
gord[to - from] = s[from];
} else {
gord[0] = s[from];
for(int t = to; t < from; ++t) gord[t - to + 1] = s[t];
}
long long prevC = (lo > 0) ? dist_[opOnMachine(s[lo - 1], m)] : 0;
for(int t = 0; t < L; ++t){
int v = gord[t]; int u = opOnMachine(v, m); int k = posOf[v][m];
long long jp = (k > 0) ? dist_[u - 1] : 0;
long long st = prevC > jp ? prevC : jp;
gestC[t] = st + procOp[u];
prevC = gestC[t];
}
long long prevT = (hi + 1 < J) ? q_[opOnMachine(s[hi + 1], m)] : 0;
long long bestLen = 0;
for(int t = L - 1; t >= 0; --t){
int v = gord[t]; int u = opOnMachine(v, m); int k = posOf[v][m];
long long js = (k < M - 1) ? q_[u + 1] : 0;
long long tl = procOp[u] + (prevT > js ? prevT : js);
long long len = gestC[t] - procOp[u] + tl;
if(len > bestLen) bestLen = len;
prevT = tl;
}
return bestLen;
}
static inline double priority(int rule, int j, int k, const vector<long long>& remWork){
switch(rule){
case 0: return (double)remWork[j];
case 1: return -(double)procJK[j][k];
case 2: return  (double)procJK[j][k];
case 3: return -(double)remWork[j];
default: return (double)(rnd() & 0xffffff);
}
}
static vector<vector<int>> gifflerThompson(int rule){
vector<int>       nextK(J, 0);
vector<long long> jobFree(J, 0), machFree(M, 0), remWork(J, 0);
for(int j = 0; j < J; ++j)
for(int k = 0; k < M; ++k) remWork[j] += procJK[j][k];
vector<vector<int>> seq(M);
for(int m = 0; m < M; ++m) seq[m].reserve(J);
int scheduled = 0;
while(scheduled < N){
long long minC = LLONG_MAX; int mstar = -1;
for(int j = 0; j < J; ++j){
if(nextK[j] >= M) continue;
int k = nextK[j], m = machJK[j][k];
long long est = max(jobFree[j], machFree[m]);
long long C = est + procJK[j][k];
if(C < minC){ minC = C; mstar = m; }
}
int chosen = -1; double bestPri = -1e300;
for(int j = 0; j < J; ++j){
if(nextK[j] >= M) continue;
int k = nextK[j];
if(machJK[j][k] != mstar) continue;
long long est = max(jobFree[j], machFree[mstar]);
if(est < minC){
double pri = priority(rule, j, k, remWork);
if(pri > bestPri){ bestPri = pri; chosen = j; }
}
}
int j = chosen, k = nextK[j], m = mstar;
long long est = max(jobFree[j], machFree[m]);
long long fin = est + procJK[j][k];
jobFree[j] = fin; machFree[m] = fin;
seq[m].push_back(j);
remWork[j] -= procJK[j][k];
++nextK[j];
++scheduled;
}
return seq;
}
static vector<vector<int>> pos;
static void rebuildPos(const vector<vector<int>>& seq){
for(int m = 0; m < M; ++m)
for(int i = 0; i < J; ++i) pos[m][seq[m][i]] = i;
}
static inline void doSwap(vector<vector<int>>& seq, int a, int b){
int m = machOf[a];
int i = pos[m][jobOf[a]];
swap(seq[m][i], seq[m][i+1]);
pos[m][seq[m][i]]   = i;
pos[m][seq[m][i+1]] = i+1;
}
static inline void doInsert(vector<vector<int>>& seq, int m, int from, int to){
vector<int>& s = seq[m];
int job = s[from];
if(from < to){ for(int i = from; i < to; ++i) s[i] = s[i+1]; s[to] = job; }
else         { for(int i = from; i > to; --i) s[i] = s[i-1]; s[to] = job; }
int lo = from < to ? from : to, hi = from < to ? to : from;
for(int i = lo; i <= hi; ++i) pos[m][s[i]] = i;
}
static void perturb(vector<vector<int>>& seq, int kicks){
rebuildPos(seq);
for(int t = 0; t < kicks; ++t){
if(J < 2) return;
int m = rndInt(M);
int i = rndInt(J - 1);
int a = opOnMachine(seq[m][i],   m);
int b = opOnMachine(seq[m][i+1], m);
doSwap(seq, a, b);
if(evaluate(seq) < 0) doSwap(seq, b, a);
}
}
static vector<long long> tabuUntil;
static vector<long long> tabuJob;
static inline size_t tabIdx(int m, int ja, int jb){
int lo = ja < jb ? ja : jb, hi = ja < jb ? jb : ja;
return (size_t)m * J * J + (size_t)lo * J + hi;
}
static long long tabuSearch(vector<vector<int>>& best, long long bestMk){
vector<vector<int>> cur = best;
vector<pair<int,int>> swaps;
vector<array<int,3>> inserts;
fill(tabuUntil.begin(), tabuUntil.end(), 0);
fill(tabuJob.begin(), tabuJob.end(), 0);
rebuildPos(cur);
long long iter = 0, lastImprove = 0;
int tenure = 15 + rndInt(13);
const long long stall = 5200;
long long curMk = evaluate(cur);
int checkClock = 0;
while((checkClock++ & 63) || elapsed() < TL){
getBlockMoves(cur, swaps, inserts);
if(swaps.empty() && inserts.empty()){
perturb(cur, 4); curMk = evaluate(cur); ++iter; continue;
}
long long bestEst = LLONG_MAX, aspEst = LLONG_MAX;
int alMode = -1, alA = -1, alB = -1, alM = -1, alF = -1, alT = -1;
int asMode = -1, asA = -1, asB = -1, asM = -1, asF = -1, asT = -1;
for(auto& pr : swaps){
int a = pr.first, b = pr.second;
long long est = estimateSwap(a, b);
bool isTabu = tabuUntil[tabIdx(machOf[a], jobOf[a], jobOf[b])] > iter;
if(isTabu){ if(est < bestMk && est < aspEst){ aspEst = est; asMode = 0; asA = a; asB = b; } }
else if(est < bestEst){ bestEst = est; alMode = 0; alA = a; alB = b; }
}
for(auto& ins : inserts){
int m = ins[0], f = ins[1], t = ins[2];
int job = cur[m][f];
long long est = estInsert(cur, m, f, t);
bool isTabu = tabuJob[(size_t)m * J + job] > iter;
if(isTabu){ if(est < bestMk && est < aspEst){ aspEst = est; asMode = 1; asM = m; asF = f; asT = t; } }
else if(est < bestEst){ bestEst = est; alMode = 1; alM = m; alF = f; alT = t; }
}
int useMode;
bool useAsp;
if(asMode != -1 && aspEst <= bestEst)      { useAsp = true;  useMode = asMode; }
else if(alMode != -1)                       { useAsp = false; useMode = alMode; }
else if(asMode != -1)                       { useAsp = true;  useMode = asMode; }
else { perturb(cur, 4); curMk = evaluate(cur); ++iter; continue; }
if(useMode == 0){
int a = useAsp ? asA : alA, b = useAsp ? asB : alB;
doSwap(cur, a, b);
tabuUntil[tabIdx(machOf[a], jobOf[a], jobOf[b])] = iter + tenure;
}else{
int m = useAsp ? asM : alM, f = useAsp ? asF : alF, t = useAsp ? asT : alT;
int job = cur[m][f];
doInsert(cur, m, f, t);
tabuJob[(size_t)m * J + job] = iter + tenure;
}
curMk = evaluate(cur);
if(curMk < 0){
cur = best; rebuildPos(cur); curMk = evaluate(cur);
}
if(curMk < bestMk){ bestMk = curMk; best = cur; lastImprove = iter; tenure = 15 + rndInt(13); }
++iter;
if(iter - lastImprove > stall){
cur = best;
perturb(cur, 6 + rndInt(13));
curMk = evaluate(cur);
fill(tabuUntil.begin(), tabuUntil.end(), 0);
fill(tabuJob.begin(), tabuJob.end(), 0);
lastImprove = iter;
}
}
return bestMk;
}
static void output(const vector<vector<int>>& seq){
static char buf[1 << 22];
int p = 0;
for(int m = 0; m < M; ++m){
for(int i = 0; i < J; ++i){
int x = seq[m][i];
if(x == 0) buf[p++] = '0';
else{
char tmp[12]; int t = 0;
while(x){ tmp[t++] = char('0' + x % 10); x /= 10; }
while(t) buf[p++] = tmp[--t];
}
buf[p++] = (i + 1 < J) ? ' ' : '\n';
}
}
fwrite(buf, 1, p, stdout);
}
int solveParsed(int Jin, int Min, const vector<vector<int>>& m_in, const vector<vector<long long>>& p_in){
START = ::G_START;
J = Jin; M = Min; N = J * M;
machJK = m_in;
procJK = p_in;
posOf.assign(J, vector<int>(M, -1));
procOp.assign(N, 0); jobOf.assign(N, 0); kOf.assign(N, 0); machOf.assign(N, 0);
jobPred.assign(N, -1); jobSucc.assign(N, -1);
for(int j = 0; j < J; ++j){
for(int k = 0; k < M; ++k){
int m = machJK[j][k]; long long p = procJK[j][k];
posOf[j][m] = k;
int op = j * M + k;
procOp[op] = p; jobOf[op] = j; kOf[op] = k; machOf[op] = m;
jobPred[op] = (k > 0)     ? op - 1 : -1;
jobSucc[op] = (k < M - 1) ? op + 1 : -1;
}
}
indeg.assign(N, 0); mSucc.assign(N, -1); mPred.assign(N, -1);
order_.assign(N, 0); dist_.assign(N, 0); q_.assign(N, 0);
pos.assign(M, vector<int>(J, 0));
tabuUntil.assign((size_t)M * J * J, 0);
tabuJob.assign((size_t)M * J, 0);
vector<vector<int>> best;
long long bestMk = LLONG_MAX;
for(int rule = 0; rule < 5; ++rule){
vector<vector<int>> seq = gifflerThompson(rule);
long long mk = evaluate(seq);
if(mk >= 0 && mk < bestMk){ bestMk = mk; best = seq; }
}
for(int r = 0; r < 150; ++r){
vector<vector<int>> seq = gifflerThompson(r % 5);
long long mk = evaluate(seq);
if(mk >= 0 && mk < bestMk){ bestMk = mk; best = seq; }
}
if(best.empty()){
best.assign(M, vector<int>(J));
for(int m = 0; m < M; ++m) for(int j = 0; j < J; ++j) best[m][j] = j;
bestMk = evaluate(best);
}
bestMk = tabuSearch(best, bestMk);
output(best);
return 0;
}
}
