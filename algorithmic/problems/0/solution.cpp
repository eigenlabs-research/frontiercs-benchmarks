#include <bits/stdc++.h>
using namespace std;
struct T{int w,h;vector<pair<int,int>> c;vector<int> lo,hi;int r,f,minx,miny;vector<unsigned> rm;vector<pair<int,int>> nbr;};
struct P{int id,k;vector<pair<int,int>> b;vector<T> t;int minW=1e9,minH=1e9,minA=1e9;};
struct Pl{int idx,ti,x,y;};
struct R{long long A;int W,H,Wpack;vector<Pl> pl;bool ok=false;};
struct RNG{unsigned long long s;RNG(unsigned long long x){s=x?x:1;}inline unsigned long long nxt(){s^=s<<7;s^=s>>9;return s;}inline int rint(int n){return (int)(nxt()%n);}inline bool coin(){return nxt()&1;}};
static inline pair<int,int> rotp(pair<int,int> p,int r){if(r==0)return p; if(r==1)return make_pair(-p.second,p.first); if(r==2)return make_pair(-p.first,-p.second); return make_pair(p.second,-p.first);}
static chrono::steady_clock::time_point gT0;
static inline double elapsed_ms(){ return chrono::duration<double,milli>(chrono::steady_clock::now()-gT0).count(); }
static int n; static long long S=0; static vector<P> ps;
static double BF_TL=1950.0;
static int BF_HD=0;
static int BF_FAST=0;
static vector<unsigned long long> bf_occ; static vector<int> bf_colH; static int bf_W;
static inline int bf_fit(const T& o,int x){
int ysky=0; for(int j=0;j<o.w;j++){ int v=bf_colH[x+j]-o.lo[j]; if(v>ysky) ysky=v; }
int start=ysky-BF_HD; if(start<0) start=0;
for(int y=start;y<=ysky;y++){ bool bad=false;
for(int dy=0;dy<o.h && !bad;dy++){ unsigned long long rm=((unsigned long long)o.rm[dy])<<x; unsigned long long orow=(y+dy<(int)bf_occ.size())?bf_occ[y+dy]:0ULL; if(rm&orow) bad=true; }
if(!bad) return y; }
return ysky;
}
static inline void bf_put(const T& o,int x,int y){
if((int)bf_occ.size()<y+o.h) bf_occ.resize(y+o.h,0ULL);
for(int dy=0;dy<o.h;dy++) bf_occ[y+dy]|=((unsigned long long)o.rm[dy])<<x;
for(int j=0;j<o.w;j++){ int top=y+o.hi[j]+1; if(top>bf_colH[x+j]) bf_colH[x+j]=top; }
}
static inline unsigned long long bf_key(int w,const int* v){ unsigned long long k=(unsigned long long)w<<52; for(int j=0;j<w;j++){ int s=v[j]+16; if(s<0)s=0; if(s>31)s=31; k|=(unsigned long long)s<<(5*j); } return k; }
static int bf_pass(int W,const vector<int>& repr,const unordered_map<unsigned long long,vector<pair<int,int>>>& idx,
vector<int>& avail,int total,vector<int>& outKind,vector<int>& outOri,vector<int>& outX,vector<int>& outY,int tie,RNG* rng=nullptr,double hardMs=1e18){
bf_W=W; bf_occ.assign(8,0ULL); bf_colH.assign(W,0);
outKind.clear(); outOri.clear(); outX.clear(); outY.clear();
int rem=total; int p[12];
long long bW,bS; int bC,bk,bo,bx,by,curH0; long long bH; long long bAw;
auto eval=[&](int ki,int oi,int x){
const T& o=ps[repr[ki]].t[oi];
int y=bf_fit(o,x);
int rH=curH0; for(int j=0;j<o.w;j++){ int t=y+o.hi[j]+1; if(t>rH) rH=t; }
long long aw=0; for(int j=0;j<o.w;j++){ int pb=y+o.lo[j]; if(pb>bf_colH[x+j]) aw+=pb-bf_colH[x+j]; }
int ct=0; if(!BF_FAST) for(auto& nb:o.nbr){ int cx=x+nb.first, cy=y+nb.second; if(cx<0||cx>=W||cy<0){ ct++; continue; } if(cy<(int)bf_occ.size() && ((bf_occ[cy]>>cx)&1ULL)) ct++; }
if(BF_FAST){ for(int j=0;j<o.w;j++){ int sy=y+o.lo[j]-1, sx=x+j; if(sy<0) ct++; else if(sy<(int)bf_occ.size() && ((bf_occ[sy]>>sx)&1ULL)) ct++; } }
long long sc=(long long)rH - ct;
bool better;
if(sc!=bS) better=(sc<bS);
else if(BF_FAST) better=(tie ? aw<=bAw : aw<bAw);
else better=(tie ? ct>=bC : ct>bC);
if(better){ bS=sc; bC=ct; bAw=aw; bk=ki; bo=oi; bx=x; by=y; }
(void)bW;(void)bH;
};
while(rem>0){
if((rem&255)==0 && elapsed_ms()>hardMs) return -1;
int mh=INT_MAX; curH0=0; for(int c=0;c<W;c++){ if(bf_colH[c]<mh) mh=bf_colH[c]; if(bf_colH[c]>curH0) curH0=bf_colH[c]; }
int nc; if(rng){ int lows[64],nl=0; for(int c=0;c<W;c++) if(bf_colH[c]==mh) lows[nl++]=c; nc=lows[rng->rint(nl)]; }
else { nc=0; while(nc<W && bf_colH[nc]!=mh) nc++; }
bH=LLONG_MAX; bW=LLONG_MAX; bS=LLONG_MAX; bC=-1; bAw=LLONG_MAX; bk=-1; bo=-1; bx=-1; by=-1;
for(int w=1;w<=10;w++){
int xlo=nc-w+1; if(xlo<0) xlo=0; int xhi=nc; if(xhi>W-w) xhi=W-w;
for(int x=xlo;x<=xhi;x++){
int base=bf_colH[x]; for(int j=0;j<w;j++) p[j]=bf_colH[x+j]-base;
auto it=idx.find(bf_key(w,p)); if(it==idx.end()) continue;
for(auto& pr:it->second) if(avail[pr.first]>0) eval(pr.first,pr.second,x);
}
}
if(bk<0)
for(int ki=0;ki<(int)repr.size();ki++){ if(avail[ki]<=0) continue;
const P& Kp=ps[repr[ki]];
for(int oi=0;oi<(int)Kp.t.size();oi++){ const T& o=Kp.t[oi]; if(o.w>W) continue;
int xlo=nc-o.w+1; if(xlo<0) xlo=0; int xhi=nc; if(xhi>W-o.w) xhi=W-o.w;
for(int x=xlo;x<=xhi;x++) eval(ki,oi,x);
}
}
if(bk<0) return -1;
bf_put(ps[repr[bk]].t[bo],bx,by); avail[bk]--; rem--;
outKind.push_back(bk); outOri.push_back(bo); outX.push_back(bx); outY.push_back(by);
}
int H=0; for(int c=0;c<W;c++) if(bf_colH[c]>H) H=bf_colH[c];
return H;
}
static long long bf_cA(int W){
unsigned long long colU=0; int rows=0;
for(auto& g:bf_occ) if(g){ colU|=g; rows++; }
int Wu=__builtin_popcountll(colU & (W==64?~0ULL:((1ULL<<W)-1)));
return 1LL*max(1,Wu)*max(1,rows);
}
static R bfSolve(int minW,int base,double deadline){
unordered_map<string,int> km; vector<int> repr; vector<int> kcnt; vector<int> kindOf(n);
vector<vector<int>> members;
string s; s.reserve(64);
for(int i=0;i<n;i++){
string keep; bool has=false;
for(auto& t:ps[i].t){ s.clear(); for(auto& q:t.c){ s.push_back((char)('0'+q.first)); s.push_back((char)('0'+q.second)); }
if(!has || s<keep){ keep=s; has=true; } }
auto it=km.find(keep);
int ki; if(it==km.end()){ ki=repr.size(); km[keep]=ki; repr.push_back(i); kcnt.push_back(0); members.push_back({}); }
else ki=it->second;
kcnt[ki]++; members[ki].push_back(i); kindOf[i]=ki;
}
int K=repr.size();
unordered_map<unsigned long long,vector<pair<int,int>>> idx; idx.reserve(K*4);
int rp[12];
for(int ki=0;ki<K;ki++) for(int oi=0;oi<(int)ps[repr[ki]].t.size();oi++){
const T& o=ps[repr[ki]].t[oi]; if(o.w>63) continue; int b0=o.lo[0];
for(int j=0;j<o.w;j++) rp[j]=o.lo[j]-b0;
idx[bf_key(o.w,rp)].push_back({ki,oi});
}
long long bestA=LLONG_MAX; int bestW=0,bestH=0; vector<int> bK,bO,bX,bY,tK,tO,tX,tY;
double lastPass=0;
for(int d=0;d<=30;d++){
for(int sgn=(d?-1:1);sgn<=1;sgn+=2){
int W=base+sgn*d; if(W<minW||W>63) continue;
for(int tie=1;tie>=0;tie--){
if(elapsed_ms()+lastPass*1.3>deadline) goto DONE;
vector<int> avail=kcnt; double t0=elapsed_ms();
int H=bf_pass(W,repr,idx,avail,n,tK,tO,tX,tY,tie,nullptr,deadline);
lastPass=elapsed_ms()-t0;
if(H>0){ long long A=bf_cA(W); if(A<bestA){ bestA=A; bestW=W; bestH=H; swap(bK,tK); swap(bO,tO); swap(bX,tX); swap(bY,tY); } }
}
}
}
if(bestW>0){
RNG g(0x9e3779b97f4a7c15ULL ^ ((unsigned long long)S<<1) ^ (unsigned long long)n);
while(elapsed_ms()+lastPass*1.25<deadline){
int W=bestW+g.rint(7)-3; if(W<minW||W>63) continue;
int tie=g.rint(2); vector<int> avail=kcnt; double t0=elapsed_ms();
int H=bf_pass(W,repr,idx,avail,n,tK,tO,tX,tY,tie,&g,deadline);
lastPass=elapsed_ms()-t0;
if(H>0){ long long A=bf_cA(W); if(A<bestA){ bestA=A; bestW=W; bestH=H; swap(bK,tK); swap(bO,tO); swap(bX,tX); swap(bY,tY); } }
}
}
DONE:;
R r; r.A=0; r.W=0; r.H=0; r.Wpack=0; if(bestW==0){ return r; }
r.W=bestW; r.H=bestH; r.A=bestA; r.Wpack=bestW; r.ok=true;
vector<int> mcur(K,0);
r.pl.reserve(n);
for(size_t z=0;z<bK.size();z++){
int ki=bK[z]; int pi=members[ki][mcur[ki]++];
int oi=bO[z];
if(ps[pi].t[oi].c!=ps[repr[ki]].t[oi].c)
for(int t=0;t<(int)ps[pi].t.size();t++)
if(ps[pi].t[t].c==ps[repr[ki]].t[oi].c){ oi=t; break; }
r.pl.push_back({pi,oi,bX[z],bY[z]});
}
return r;
}
static int ep_envInt(const char* k,int d){ if(const char* e=getenv(k)) return atoi(e); return d; }
static int ep_gFASTFIT=0;
static vector<unsigned long long> ep_grid, ep_cgrid;
static vector<short> ep_own;
static vector<int> ep_plSlot;
struct EB3{int y,x,contact;bool valid;};
static vector<EB3> ep_b3;
static R epack(int W, const vector<int>& o0, RNG& rng, bool randtie, int dynLIM0, bool adaptive,
double adaptTLms, double deadlineMs, double panicMs = -1.0, int adaptMode = 0) {
vector<int> h(W, -1);
long long g = -1;
vector<Pl> pl; pl.reserve(o0.size());
vector<int> o = o0;
int t = 0, nm = (int)o.size();
int maxBound = max(1, nm / 4);
int dynLIM = dynLIM0;
auto tStart = chrono::steady_clock::now();
auto batchStart = tStart;
double emaBatch = -1.0;
int stepCnt = 0;
int checkCnt = 0;
bool panic = false;
while (t < nm) {
if (!panic && panicMs > 0 && elapsed_ms() > panicMs) { panic = true; dynLIM = 1; adaptive = false; }
int limCnt = max(1, dynLIM);
int lim = min(nm, t + limCnt);
long long bestg = LLONG_MAX; int bti = -1, bx = 0, by = 0, bl = INT_MAX;
long long bds = LLONG_MAX, bdr = LLONG_MAX; int by0 = INT_MAX, bx0 = INT_MAX;
int bestpos = t, bestid = -1;
for (int pos = t; pos < lim; pos++) {
if (!panic && deadlineMs > 0 && ((++checkCnt & 7) == 0) && elapsed_ms() > deadlineMs) return R{};
int id = o[pos];
auto& p = ps[id];
long long bestg2 = LLONG_MAX; int bti2 = -1, bx2 = 0, by2 = 0, bl2 = INT_MAX;
long long bds2 = LLONG_MAX, bdr2 = LLONG_MAX; int by02 = INT_MAX, bx02 = INT_MAX;
for (int ti = 0; ti < (int)p.t.size(); ti++) {
auto& tsh = p.t[ti];
if (tsh.w > W) continue;
int Rpos = W - tsh.w + 1;
for (int x0 = 0; x0 < Rpos; x0++) {
int y0 = 0;
for (int j = 0; j < tsh.w; j++) {
if (tsh.lo[j] != INT_MAX) { int v = h[x0 + j] - tsh.lo[j] + 1; if (v > y0) y0 = v; }
}
int nhbuf[32];
int l = -1;
for (int j = 0; j < tsh.w; j++) {
int nh = h[x0 + j];
if (tsh.hi[j] != INT_MIN) { int cand = y0 + tsh.hi[j]; if (nh < cand) nh = cand; }
nhbuf[j] = nh;
if (nh > l) l = nh;
}
long long gg = g; if (gg < l) gg = l;
if (ep_gFASTFIT && bti2 != -1 && gg > bestg2) continue;
long long dsum = 0;
for (int j = 0; j < tsh.w; j++) { int inc = nhbuf[j] - h[x0 + j]; if (inc > 0) dsum += inc; }
long long dr = 0;
if (x0 > 0) {
long long old = llabs((long long)h[x0] - h[x0 - 1]);
long long nw = llabs((long long)nhbuf[0] - h[x0 - 1]);
dr += nw - old;
}
for (int j = 0; j < tsh.w - 1; j++) {
long long old = llabs((long long)h[x0 + j + 1] - h[x0 + j]);
long long nw = llabs((long long)nhbuf[j + 1] - nhbuf[j]);
dr += nw - old;
}
if (x0 + tsh.w < W) {
long long old = llabs((long long)h[x0 + tsh.w] - h[x0 + tsh.w - 1]);
long long nw = llabs((long long)h[x0 + tsh.w] - nhbuf[tsh.w - 1]);
dr += nw - old;
}
bool take = false;
if (gg < bestg2) take = true;
else if (gg == bestg2) {
if (dsum < bds2) take = true;
else if (dsum == bds2) {
if (l < bl2) take = true;
else if (l == bl2) {
if (dr < bdr2) take = true;
else if (dr == bdr2) {
if (y0 < by02) take = true;
else if (y0 == by02) {
if (x0 < bx02) take = true;
else if (x0 == bx02 && randtie && rng.coin()) take = true;
}
}
}
}
}
if (take) { bestg2 = gg; bti2 = ti; bx2 = x0; by2 = y0; bl2 = l; bds2 = dsum; bdr2 = dr; by02 = y0; bx02 = x0; }
}
}
if (bti2 == -1) continue;
bool take = false;
if (bestg2 < bestg) take = true;
else if (bestg2 == bestg) {
if (bds2 < bds) take = true;
else if (bds2 == bds) {
if (bl2 < bl) take = true;
else if (bl2 == bl) {
if (bdr2 < bdr) take = true;
else if (bdr2 == bdr) {
if (by02 < by0) take = true;
else if (by02 == by0) {
if (bx02 < bx0) take = true;
else if (bx02 == bx0 && randtie && rng.coin()) take = true;
}
}
}
}
}
if (take) { bestg = bestg2; bti = bti2; bx = bx2; by = by2; bl = bl2; bds = bds2; bdr = bdr2; by0 = by02; bx0 = bx02; bestpos = pos; bestid = id; }
}
auto adapt = [&]() {
if (adaptive && stepCnt == 5) {
auto now = chrono::steady_clock::now();
double elapsed = chrono::duration<double, milli>(now - tStart).count();
double batch = chrono::duration<double, milli>(now - batchStart).count();
double remT = max(0.0, adaptTLms - elapsed);
int remSteps = max(1, nm - t);
double budget = remT * 5.0 / remSteps;
if (adaptMode == 1) {
emaBatch = (emaBatch < 0.0) ? batch : (0.6 * emaBatch + 0.4 * batch);
int floorL = max(1, dynLIM0 / 4);
if (emaBatch < budget * 0.85) dynLIM = min(maxBound, dynLIM + max(1, dynLIM / 8));
else if (emaBatch > budget * 1.15) dynLIM = max(floorL, dynLIM - max(1, dynLIM / 6));
} else {
if (batch < budget) dynLIM = min(maxBound, dynLIM + 1); else dynLIM = max(1, dynLIM - 1);
}
batchStart = now;
stepCnt = 0;
}
};
if (bti == -1) {
t++; stepCnt++;
adapt();
if (!panic && deadlineMs > 0 && elapsed_ms() > deadlineMs) return R{};
continue;
}
auto& tsh = ps[bestid].t[bti];
for (int j = 0; j < tsh.w; j++) {
int nh = h[bx + j];
if (tsh.hi[j] != INT_MIN) { int cand = by + tsh.hi[j]; if (nh < cand) nh = cand; }
h[bx + j] = nh;
}
if (g < bl) g = bl;
pl.push_back({bestid, bti, bx, by});
if (bestpos != t) swap(o[t], o[bestpos]);
t++;
stepCnt++;
adapt();
if (!panic && deadlineMs > 0 && elapsed_ms() > deadlineMs) return R{};
}
if ((int)pl.size() != nm) return R{};
int H = (int)g + 1;
int maxX = -1;
for (auto& pp : pl) {
auto& t2 = ps[pp.idx].t[pp.ti];
for (auto& q : t2.c) { int x = pp.x + q.first; if (x > maxX) maxX = x; }
}
int Wused = 0;
if (maxX >= 0) {
vector<char> used(maxX + 1, false);
for (auto& pp : pl) { auto& t2 = ps[pp.idx].t[pp.ti]; for (auto& q : t2.c) used[pp.x + q.first] = true; }
for (int x = 0; x <= maxX; x++) if (used[x]) Wused++;
}
int Wfinal = max(1, Wused);
long long A = 1LL * H * Wfinal;
R res; res.A = A; res.W = Wfinal; res.H = H; res.pl = move(pl); res.ok = true; res.Wpack = W;
return res;
}
static R epack_blf(int W, const vector<int>& order, int policy, double deadlineMs) {
if (W > 64 || W <= 0) return R{};
int nm = (int)order.size();
auto& grid = ep_grid;
size_t estH = (size_t)(S / max(1, W)) + 48;
grid.assign(max<size_t>(64, estH), 0ULL);
const uint64_t FULLROW = (W == 64) ? ~0ULL : ((1ULL << W) - 1);
int y0full = 0;
int topY = 0;
vector<Pl> pl; pl.reserve(nm);
int checkCnt = 0;
for (int oi = 0; oi < nm; oi++) {
if (deadlineMs > 0 && ((++checkCnt & 15) == 0) && elapsed_ms() > deadlineMs) return R{};
int id = order[oi];
auto& p = ps[id];
int bestY = INT_MAX, bestX = 0, bestTi = -1; long long bestKey = LLONG_MAX;
if ((int)grid.size() < topY + 12) grid.resize(topY + 64, 0ULL);
for (int ti = 0; ti < (int)p.t.size(); ti++) {
auto& t = p.t[ti];
if (t.w > W) continue;
uint64_t limMask = (W - t.w + 1 >= 64) ? ~0ULL : ((1ULL << (W - t.w + 1)) - 1);
int yFound = -1, xFound = -1;
for (int y = y0full; y <= topY; y++) {
uint64_t conflict = 0;
for (int dy = 0; dy < t.h; dy++) {
uint64_t g = grid[y + dy];
if (!g) continue;
unsigned m = t.rm[dy];
while (m) { int b = __builtin_ctz(m); conflict |= (g >> b); m &= m - 1; }
}
uint64_t allowed = ~conflict & limMask;
if (allowed) { yFound = y; xFound = (int)__builtin_ctzll(allowed); break; }
}
if (yFound < 0) continue;
long long key = (policy == 0) ? ((long long)yFound << 8) | xFound
: ((long long)(yFound + t.h) << 24) | ((long long)yFound << 8) | xFound;
if (key < bestKey) { bestKey = key; bestY = yFound; bestX = xFound; bestTi = ti; }
}
if (bestTi < 0) return R{};
auto& t = p.t[bestTi];
if ((int)grid.size() < bestY + t.h + 1) grid.resize(bestY + t.h + 64, 0ULL);
for (int dy = 0; dy < t.h; dy++) grid[bestY + dy] |= ((uint64_t)t.rm[dy]) << bestX;
if (bestY + t.h > topY) topY = bestY + t.h;
while (y0full < topY && grid[y0full] == FULLROW) y0full++;
pl.push_back({id, bestTi, bestX, bestY});
}
uint64_t colU = 0; int rows = 0;
for (int y = 0; y < topY; y++) if (grid[y]) { colU |= grid[y]; rows++; }
int Wu = __builtin_popcountll(colU);
R res; res.W = Wu ? Wu : 1; res.H = rows ? rows : 1;
res.A = (long long)res.W * res.H; res.pl = move(pl); res.ok = true; res.Wpack = W;
return res;
}
static bool ecrownRepack(R& r, double deadlineMs) {
if (!r.ok || r.Wpack <= 0 || r.Wpack > 64) return false;
int W = r.Wpack;
int H0 = 0;
for (auto& p : r.pl) { auto& t = ps[p.idx].t[p.ti]; if (p.y + t.h > H0) H0 = p.y + t.h; }
if (H0 <= 1) return false;
auto& grid = ep_cgrid;
grid.assign(H0, 0ULL);
for (auto& p : r.pl) { auto& t = ps[p.idx].t[p.ti]; for (int dy = 0; dy < t.h; dy++) grid[p.y + dy] |= ((uint64_t)t.rm[dy]) << p.x; }
bool improvedAny = false;
long long origA = r.A; int origW = r.W, origH = r.H;
vector<Pl> origPl = r.pl;
int rounds = 0;
while (rounds++ < 64) {
if (deadlineMs > 0 && elapsed_ms() > deadlineMs) break;
int H = (int)grid.size();
while (H > 0 && grid[H - 1] == 0) H--;
if (H <= 1) break;
bool anyDepthOk = false;
for (int crownDepth = 1; crownDepth <= 3; crownDepth++) {
if (H - crownDepth < 1) break;
vector<int> crown;
for (int i = 0; i < (int)r.pl.size(); i++) { auto& p = r.pl[i]; auto& t = ps[p.idx].t[p.ti]; if (p.y + t.h > H - crownDepth) crown.push_back(i); }
if (crown.empty()) { anyDepthOk = true; break; }
for (int ci : crown) { auto& p = r.pl[ci]; auto& t = ps[p.idx].t[p.ti]; for (int dy = 0; dy < t.h; dy++) grid[p.y + dy] &= ~(((uint64_t)t.rm[dy]) << p.x); }
sort(crown.begin(), crown.end(), [&](int a, int b) { return ps[r.pl[a].idx].k > ps[r.pl[b].idx].k; });
int targetH = H - crownDepth;
vector<array<int,3>> newPos(crown.size());
bool ok = true;
for (size_t c = 0; c < crown.size() && ok; c++) {
auto& p = r.pl[crown[c]];
auto& piece = ps[p.idx];
int bTi = -1, bX = 0, bY = INT_MAX;
for (int ti = 0; ti < (int)piece.t.size(); ti++) {
auto& t = piece.t[ti];
if (t.w > W || t.h > targetH) continue;
uint64_t limMask = (W - t.w + 1 >= 64) ? ~0ULL : ((1ULL << (W - t.w + 1)) - 1);
int ymax = targetH - t.h;
for (int y = 0; y <= ymax; y++) {
if (bY <= y) break;
uint64_t conflict = 0;
for (int dy = 0; dy < t.h; dy++) { uint64_t g = grid[y + dy]; if (!g) continue; unsigned m = t.rm[dy]; while (m) { int b = __builtin_ctz(m); conflict |= (g >> b); m &= m - 1; } }
uint64_t allowed = ~conflict & limMask;
if (allowed) { int x = (int)__builtin_ctzll(allowed); if (y < bY || (y == bY && x < bX)) { bY = y; bX = x; bTi = ti; } break; }
}
}
if (bTi < 0) { ok = false; break; }
auto& t = piece.t[bTi];
for (int dy = 0; dy < t.h; dy++) grid[bY + dy] |= ((uint64_t)t.rm[dy]) << bX;
newPos[c] = {bTi, bX, bY};
}
if (ok) {
uint64_t colU = 0; int rows = 0;
for (auto& g : grid) if (g) { colU |= g; rows++; }
int newW = max(1, __builtin_popcountll(colU));
int newH = max(1, rows);
long long newA = (long long)newW * newH;
if (newA < r.A) {
for (size_t c = 0; c < crown.size(); c++) { auto& p = r.pl[crown[c]]; p.ti = newPos[c][0]; p.x = newPos[c][1]; p.y = newPos[c][2]; }
r.W = newW; r.H = newH; r.A = newA;
improvedAny = true; anyDepthOk = true; break;
} else ok = false;
}
if (!ok) {
fill(grid.begin(), grid.end(), 0ULL);
for (auto& p : r.pl) { auto& t = ps[p.idx].t[p.ti]; for (int dy = 0; dy < t.h; dy++) grid[p.y + dy] |= ((uint64_t)t.rm[dy]) << p.x; }
}
}
if (!anyDepthOk) break;
}
if (!improvedAny) { r.pl = move(origPl); r.W = origW; r.H = origH; r.A = origA; }
return improvedAny;
}
static bool ecrownRepack2(R& r,double deadlineMs){
if(!r.ok||r.Wpack<=0||r.Wpack>64||deadlineMs-elapsed_ms()<8)return false;
int W=r.Wpack,commits=0; bool improved=false;
vector<uint64_t> c_grid;
while(commits<16&&deadlineMs-elapsed_ms()>=8){
int H=0; for(auto&p:r.pl)H=max(H,p.y+ps[p.idx].t[p.ti].h);
if(H<=1)break;
int Ht=H-1;
vector<short> owner((size_t)W*H,-1);
c_grid.assign(H,0);
for(int i=0;i<(int)r.pl.size();i++){
auto&p=r.pl[i]; auto&t=ps[p.idx].t[p.ti];
for(auto&q:t.c)owner[(p.y+q.second)*W+p.x+q.first]=i;
for(int y=0;y<t.h;y++)c_grid[p.y+y]|=(uint64_t)t.rm[y]<<p.x;
}
bool committed=false;
for(int d:{4,6,8,10}){
if(deadlineMs-elapsed_ms()<8)break;
vector<char> mv(r.pl.size()); vector<int> seed;
for(int i=0;i<(int)r.pl.size();i++){
auto&p=r.pl[i]; auto&t=ps[p.idx].t[p.ti];
if(p.y+t.h>H-d){mv[i]=1;seed.push_back(i);}
}
static const int dx[4]={1,-1,0,0},dy[4]={0,0,1,-1};
for(int i:seed){
auto&p=r.pl[i]; auto&t=ps[p.idx].t[p.ti];
for(auto&q:t.c)for(int z=0;z<4;z++){
int x=p.x+q.first+dx[z],y=p.y+q.second+dy[z];
if(x>=0&&x<W&&y>=0&&y<H){short j=owner[y*W+x];if(j>=0&&j!=i)mv[j]=1;}
}
}
vector<int> ids; for(int i=0;i<(int)mv.size();i++)if(mv[i])ids.push_back(i);
if((int)ids.size()>min(n,max(24,3*W)))continue;
vector<uint64_t> fixed=c_grid;
for(int i:ids){auto&p=r.pl[i];auto&t=ps[p.idx].t[p.ti];for(int y=0;y<t.h;y++)fixed[p.y+y]&=~((uint64_t)t.rm[y]<<p.x);}
bool bad=false;for(int y=Ht;y<H;y++)if(fixed[y])bad=true;if(bad)continue;
sort(ids.begin(),ids.end(),[&](int a,int b){auto&A=ps[r.pl[a].idx];auto&B=ps[r.pl[b].idx];if(A.k!=B.k)return A.k>B.k;if(A.minA-A.k!=B.minA-B.k)return A.minA-A.k>B.minA-B.k;return A.id<B.id;});
for(int at=0;at<6&&deadlineMs-elapsed_ms()>=8;at++){
vector<int> o=ids;
if(at){RNG rg((unsigned long long)S^((unsigned long long)n*10007ULL)^((unsigned long long)d*131ULL)^(unsigned long long)at);for(int z=0;z<(int)o.size()/3+1;z++)swap(o[rg.rint((int)o.size())],o[rg.rint((int)o.size())]);}
vector<uint64_t> grid=fixed;vector<Pl> np;bool ok=true;
for(int step=0;step<(int)o.size()&&ok;step++){
if(deadlineMs-elapsed_ms()<4){ok=false;break;}
int lim=min((int)o.size(),step+64),bp=-1,bti=-1,bx=0,by=0,bc=-1;
for(int pos=step;pos<lim;pos++){
if(deadlineMs-elapsed_ms()<3) break;
int slot=o[pos],pid=r.pl[slot].idx;
for(int ti=0;ti<(int)ps[pid].t.size();ti++){
auto&t=ps[pid].t[ti];
if(t.w>W||t.h>Ht)continue;
for(int y=0;y+t.h<=Ht;y++)for(int x=0;x+t.w<=W;x++){
bool hit=false;for(int yy=0;yy<t.h;yy++)if(grid[y+yy]&((uint64_t)t.rm[yy]<<x)){hit=true;break;}if(hit)continue;
int ct=0;for(auto&nb:t.nbr){int xx=x+nb.first,yy=y+nb.second;if(xx<0||xx>=W||yy<0||yy>=Ht||((grid[yy]>>xx)&1))ct++;}
if(ct>bc||(ct==bc&&(y<by||(y==by&&x<bx)))){bc=ct;by=y;bx=x;bti=ti;bp=pos;}
}
}
}
if(bp<0){ok=false;break;}
swap(o[step],o[bp]);
int slot=o[step],pid=r.pl[slot].idx;auto&t=ps[pid].t[bti];
for(int y=0;y<t.h;y++)grid[by+y]|=(uint64_t)t.rm[y]<<bx;
np.push_back({pid,bti,bx,by});
}
if(ok){
vector<Pl> cand=r.pl;for(int z=0;z<(int)o.size();z++)cand[o[z]]=np[z];
int mw=0,mh=0;for(auto&p:cand){auto&t=ps[p.idx].t[p.ti];mw=max(mw,p.x+t.w);mh=max(mh,p.y+t.h);}
long long A=1LL*mw*mh;
if(A<r.A||(A==r.A&&(mh<r.H||(mh==r.H&&mw<r.W)))){r.pl=move(cand);r.A=A;r.W=mw;r.H=mh;improved=committed=true;commits++;break;}
}
}
if(committed)break;
}
if(!committed)break;
}
return improved;
}
static R epack_blf2(int W, const vector<int>& order, int window, double deadlineMs, RNG& rng, bool randtie) {
if (W > 64 || W <= 0) return R{};
int nm = (int)order.size();
auto& grid = ep_grid;
size_t estH = (size_t)(S / max(1, W)) + 48;
grid.assign(max<size_t>(64, estH), 0ULL);
const uint64_t FULLROW = (W == 64) ? ~0ULL : ((1ULL << W) - 1);
int y0full = 0, topY = 0;
long long g = -1;
vector<int> o = order;
vector<Pl> pl; pl.reserve(nm);
int checkCnt = 0;
for (int t = 0; t < nm; t++) {
if (deadlineMs > 0 && ((++checkCnt & 3) == 0) && elapsed_ms() > deadlineMs) return R{};
int lim = min(nm, t + max(1, window));
long long bG = LLONG_MAX; int bContact = -1, bY = INT_MAX, bX = INT_MAX, bTi = -1, bPos = -1;
if ((int)grid.size() < topY + 12) grid.resize(topY + 64, 0ULL);
for (int pos = t; pos < lim; pos++) {
auto& p = ps[o[pos]];
for (int ti = 0; ti < (int)p.t.size(); ti++) {
auto& tt = p.t[ti];
if (tt.w > W) continue;
uint64_t limMask = (W - tt.w + 1 >= 64) ? ~0ULL : ((1ULL << (W - tt.w + 1)) - 1);
int yF = -1, xF = -1;
for (int y = y0full; y <= topY; y++) {
uint64_t conflict = 0;
for (int dy = 0; dy < tt.h; dy++) { uint64_t gg = grid[y + dy]; if (!gg) continue; unsigned m = tt.rm[dy]; while (m) { int b = __builtin_ctz(m); conflict |= (gg >> b); m &= m - 1; } }
uint64_t allowed = ~conflict & limMask;
if (allowed) { yF = y; xF = (int)__builtin_ctzll(allowed); break; }
}
if (yF < 0) continue;
long long g2 = max(g, (long long)(yF + tt.h - 1));
int GAPW1b = 1;
int contact = 0;
for (auto& nb : tt.nbr) { int cx = xF + nb.first, cy = yF + nb.second; if (cx < 0 || cx >= W || cy < 0) { contact++; continue; } if (cy < (int)grid.size() && (grid[cy] >> cx) & 1ULL) contact++; }
long long s2 = g2 * GAPW1b - contact;
long long bs = (bTi>=0) ? (bG * (long long)GAPW1b - bContact) : LLONG_MAX;
if (s2 > bs) continue;
bool take = false;
if (s2 < bs) take = true;
else { if (contact > bContact) take = true; else if (contact == bContact) { if (yF < bY) take = true; else if (yF == bY) { if (xF < bX) take = true; else if (xF == bX && randtie && rng.coin()) take = true; } } }
if (take) { bG = g2; bContact = contact; bY = yF; bX = xF; bTi = ti; bPos = pos; }
}
}
if (bPos < 0) return R{};
auto& p = ps[o[bPos]];
auto& tt = p.t[bTi];
if ((int)grid.size() < bY + tt.h + 1) grid.resize(bY + tt.h + 64, 0ULL);
for (int dy = 0; dy < tt.h; dy++) grid[bY + dy] |= ((uint64_t)tt.rm[dy]) << bX;
if (bY + tt.h > topY) topY = bY + tt.h;
while (y0full < topY && grid[y0full] == FULLROW) y0full++;
if (bG > g) g = bG;
pl.push_back({o[bPos], bTi, bX, bY});
if (bPos != t) swap(o[t], o[bPos]);
}
uint64_t colU = 0; int rows = 0;
for (int y = 0; y < topY; y++) if (grid[y]) { colU |= grid[y]; rows++; }
int Wu = __builtin_popcountll(colU);
R res; res.W = Wu ? Wu : 1; res.H = rows ? rows : 1;
res.A = (long long)res.W * res.H; res.pl = move(pl); res.ok = true; res.Wpack = W;
return res;
}
static R epack_blf3(int W, const vector<int>& order, int window, double deadlineMs, RNG& rng, bool randtie) {
if (W > 64 || W <= 0) return R{};
int nm = (int)order.size();
auto& grid = ep_grid;
size_t estH = (size_t)(S / max(1, W)) + 48;
grid.assign(max<size_t>(64, estH), 0ULL);
const uint64_t FULLROW = (W == 64) ? ~0ULL : ((1ULL << W) - 1);
int y0full = 0, topY = 0;
long long g = -1;
vector<int> o = order;
vector<Pl> pl; pl.reserve(nm);
if ((int)ep_b3.size() < nm * 8) ep_b3.resize(nm * 8);
for (int i = 0; i < nm * 8; i++) ep_b3[i].valid = false;
int checkCnt = 0;
for (int t = 0; t < nm; t++) {
if (deadlineMs > 0 && ((++checkCnt & 15) == 0) && elapsed_ms() > deadlineMs) return R{};
int lim = min(nm, t + max(1, window));
long long bG = LLONG_MAX; int bContact = -1, bY = INT_MAX, bX = INT_MAX, bTi = -1, bPos = -1;
if ((int)grid.size() < topY + 12) grid.resize(topY + 64, 0ULL);
for (int pos = t; pos < lim; pos++) {
int pid = o[pos];
auto& p = ps[pid];
for (int ti = 0; ti < (int)p.t.size(); ti++) {
auto& tt = p.t[ti];
if (tt.w > W) continue;
EB3& cc = ep_b3[pid * 8 + ti];
if (!cc.valid) {
uint64_t limMask = (W - tt.w + 1 >= 64) ? ~0ULL : ((1ULL << (W - tt.w + 1)) - 1);
int yF = -1, xF = -1;
for (int y = y0full; y <= topY; y++) {
uint64_t conflict = 0;
for (int dy = 0; dy < tt.h; dy++) { uint64_t gg = grid[y + dy]; if (!gg) continue; unsigned m = tt.rm[dy]; while (m) { int b = __builtin_ctz(m); conflict |= (gg >> b); m &= m - 1; } }
uint64_t allowed = ~conflict & limMask;
if (allowed) { yF = y; xF = (int)__builtin_ctzll(allowed); break; }
}
if (yF < 0) { continue; }
int contact = 0;
for (auto& nb : tt.nbr) { int cx = xF + nb.first, cy = yF + nb.second; if (cx < 0 || cx >= W || cy < 0) { contact++; continue; } if (cy < (int)grid.size() && (grid[cy] >> cx) & 1ULL) contact++; }
cc = {yF, xF, contact, true};
}
long long g2 = max(g, (long long)(cc.y + tt.h - 1));
int GAPW1 = 1;
long long s2 = g2 * GAPW1 - cc.contact;
long long bs = (bTi>=0) ? (bG * (long long)GAPW1 - bContact) : LLONG_MAX;
if (s2 > bs) continue;
bool take = false;
if (s2 < bs) take = true;
else { if (cc.contact > bContact) take = true; else if (cc.contact == bContact) { if (cc.y < bY) take = true; else if (cc.y == bY) { if (cc.x < bX) take = true; else if (cc.x == bX && randtie && rng.coin()) take = true; } } }
if (take) { bG = g2; bContact = cc.contact; bY = cc.y; bX = cc.x; bTi = ti; bPos = pos; }
}
}
if (bPos < 0) return R{};
int placedId = o[bPos];
auto& tt = ps[placedId].t[bTi];
if ((int)grid.size() < bY + tt.h + 1) grid.resize(bY + tt.h + 64, 0ULL);
for (int dy = 0; dy < tt.h; dy++) grid[bY + dy] |= ((uint64_t)tt.rm[dy]) << bX;
if (bY + tt.h > topY) topY = bY + tt.h;
while (y0full < topY && grid[y0full] == FULLROW) y0full++;
if (bG > g) g = bG;
pl.push_back({placedId, bTi, bX, bY});
if (bPos != t) swap(o[t], o[bPos]);
int px0 = bX - 1, px1 = bX + tt.w, py0 = bY - 1, py1 = bY + tt.h;
int lim2 = min(nm, (t + 1) + max(1, window));
for (int pos = t + 1; pos < lim2; pos++) {
int pid = o[pos];
auto& p = ps[pid];
for (int ti = 0; ti < (int)p.t.size(); ti++) { EB3& cc = ep_b3[pid * 8 + ti]; if (!cc.valid) continue; auto& ct = p.t[ti]; if (cc.x <= px1 && cc.x + ct.w - 1 >= px0 && cc.y <= py1 && cc.y + ct.h - 1 >= py0) cc.valid = false; }
}
}
uint64_t colU = 0; int rows = 0;
for (int y = 0; y < topY; y++) if (grid[y]) { colU |= grid[y]; rows++; }
int Wu = __builtin_popcountll(colU);
R res; res.W = Wu ? Wu : 1; res.H = rows ? rows : 1;
res.A = (long long)res.W * res.H; res.pl = move(pl); res.ok = true; res.Wpack = W;
return res;
}
static R epack_capped(int W, int Hcap, const vector<int>& order, int window, double deadlineMs, RNG& rng) {
const int REPAIR = 25, CAPSEL = 1;
if (W > 64 || W <= 0 || Hcap <= 0) return R{};
int nm = (int)order.size();
auto& grid = ep_grid; int gridH = Hcap + 12; grid.assign(gridH, 0ULL);
const uint64_t FULLROW = (W == 64) ? ~0ULL : ((1ULL << W) - 1);
int y0full = 0, topY = 0;
vector<int> o = order;
vector<Pl> pl; pl.reserve(nm);
if ((int)ep_b3.size() < nm * 8) ep_b3.resize(nm * 8);
for (int i = 0; i < nm * 8; i++) ep_b3[i].valid = false;
auto& own = ep_own; auto& plSlot = ep_plSlot;
if (REPAIR > 0) { own.assign((size_t)gridH * W, -1); if ((int)plSlot.size() < nm) plSlot.resize(nm); for (int i = 0; i < nm; i++) plSlot[i] = -1; }
int checkCnt = 0;
auto place = [&](int pid, int ti, int x, int y) {
auto& tt = ps[pid].t[ti];
for (int dy = 0; dy < tt.h; dy++) grid[y + dy] |= ((uint64_t)tt.rm[dy]) << x;
if (REPAIR > 0) { for (auto& q : tt.c) own[(size_t)(y + q.second) * W + (x + q.first)] = (short)pid; }
if (y + tt.h > topY) topY = min(y + tt.h, Hcap);
while (y0full < topY && grid[y0full] == FULLROW) y0full++;
if (plSlot.size() && REPAIR > 0 && plSlot[pid] >= 0) pl[plSlot[pid]] = {pid, ti, x, y};
else { if (REPAIR > 0) plSlot[pid] = (int)pl.size(); pl.push_back({pid, ti, x, y}); }
};
auto eject = [&](int pid) {
int slot = plSlot[pid]; Pl p = pl[slot]; auto& tt = ps[pid].t[p.ti];
for (int dy = 0; dy < tt.h; dy++) grid[p.y + dy] &= ~(((uint64_t)tt.rm[dy]) << p.x);
for (auto& q : tt.c) own[(size_t)(p.y + q.second) * W + (p.x + q.first)] = -1;
if (p.y < y0full) y0full = p.y;
};
auto invalidateWindow = [&](int t) {
int lim2 = min(nm, t + max(1, window) + 1);
for (int pos = t; pos < lim2; pos++) { int pid = o[pos]; for (int ti = 0; ti < (int)ps[pid].t.size(); ti++) ep_b3[pid * 8 + ti].valid = false; }
};
auto scanFit = [&](const T& tt) -> pair<int,int> {
uint64_t limMask = (W - tt.w + 1 >= 64) ? ~0ULL : ((1ULL << (W - tt.w + 1)) - 1);
int ymax = min(topY, Hcap - tt.h);
for (int y = y0full; y <= ymax; y++) {
uint64_t conflict = 0;
for (int dy = 0; dy < tt.h; dy++) { uint64_t gg = grid[y + dy]; if (!gg) continue; unsigned m = tt.rm[dy]; while (m) { int b = __builtin_ctz(m); conflict |= (gg >> b); m &= m - 1; } }
uint64_t allowed = ~conflict & limMask;
if (allowed) return {y, (int)__builtin_ctzll(allowed)};
}
return {-1, -1};
};
auto tryRepair = [&](int stuckPid, int tCur) -> bool {
int pending = stuckPid; int lastPlaced = -1;
for (int depth = 0; depth <= REPAIR; depth++) {
if (deadlineMs > 0 && elapsed_ms() > deadlineMs) return false;
auto& p = ps[pending];
if (depth > 0) {
int bTi = -1, bYy = INT_MAX, bXx = 0;
for (int ti = 0; ti < (int)p.t.size(); ti++) { auto& tt = p.t[ti]; if (tt.w > W || tt.h > Hcap) continue; auto pr = scanFit(tt); int y = pr.first, x = pr.second; if (y >= 0 && (y < bYy || (y == bYy && x < bXx))) { bYy = y; bXx = x; bTi = ti; } }
if (bTi >= 0) { place(pending, bTi, bXx, bYy); return true; }
}
if (depth == REPAIR) return false;
int cTi = -1, cY = -1, cX = -1; int cBlocker = -1; int cCells = INT_MAX;
for (int ti = 0; ti < (int)p.t.size(); ti++) { auto& tt = p.t[ti]; if (tt.w > W || tt.h > Hcap) continue; int ymax = min(topY, Hcap - tt.h);
for (int y = 0; y <= ymax && (cY < 0 || y <= cY); y++) { for (int x = 0; x + tt.w <= W; x++) { int blocker = -1, bcells = 0; bool bad = false;
for (auto& q : tt.c) { size_t idx = (size_t)(y + q.second) * W + (x + q.first); short ow = own[idx]; if (ow < 0) continue; if (ow == pending) { bad = true; break; } if (ow == lastPlaced) { bad = true; break; } if (blocker < 0) blocker = ow; else if (blocker != ow) { bad = true; break; } bcells++; }
if (bad || blocker < 0) continue;
if (cY < 0 || y < cY || (y == cY && bcells < cCells)) { cTi = ti; cY = y; cX = x; cBlocker = blocker; cCells = bcells; } } } }
if (cBlocker < 0) return false;
eject(cBlocker); place(pending, cTi, cX, cY); lastPlaced = pending; pending = cBlocker;
}
return false;
};
for (int t = 0; t < nm; t++) {
if (deadlineMs > 0 && ((++checkCnt & 15) == 0) && elapsed_ms() > deadlineMs) return R{};
int lim = min(nm, t + max(1, window));
int bContact = -1, bY = INT_MAX, bX = INT_MAX, bTi = -1, bPos = -1; bool redo = false;
for (int pos = t; pos < lim; pos++) {
int pid = o[pos]; auto& p = ps[pid]; bool any = false;
for (int ti = 0; ti < (int)p.t.size(); ti++) {
auto& tt = p.t[ti]; if (tt.w > W || tt.h > Hcap) continue;
EB3& cc = ep_b3[pid * 8 + ti];
if (!cc.valid) { auto pr = scanFit(tt); int yF = pr.first, xF = pr.second; if (yF < 0) { cc.valid = false; continue; }
int contact = 0; for (auto& nb : tt.nbr) { int cx = xF + nb.first, cy = yF + nb.second; if (cx < 0 || cx >= W || cy < 0) { contact++; continue; } if (cy >= Hcap) { contact++; continue; } if (cy < (int)grid.size() && (grid[cy] >> cx) & 1ULL) contact++; }
cc = {yF, xF, contact, true}; }
any = true; bool take = false;
if (CAPSEL == 0) { if (cc.y < bY) take = true; else if (cc.y == bY) { if (cc.contact > bContact) take = true; else if (cc.contact == bContact && cc.x < bX) take = true; } }
else { if (cc.contact > bContact) take = true; else if (cc.contact == bContact) { if (cc.y < bY) take = true; else if (cc.y == bY && cc.x < bX) take = true; } }
if (take) { bContact = cc.contact; bY = cc.y; bX = cc.x; bTi = ti; bPos = pos; }
}
if (!any) { if (REPAIR > 0 && tryRepair(pid, t)) { if (pos != t) swap(o[t], o[pos]); invalidateWindow(t); redo = true; break; } return R{}; }
}
if (redo) continue;
if (bPos < 0) return R{};
int placedId = o[bPos]; place(placedId, bTi, bX, bY); if (bPos != t) swap(o[t], o[bPos]);
int px0 = bX - 1, px1 = bX + ps[placedId].t[bTi].w, py0 = bY - 1, py1 = bY + ps[placedId].t[bTi].h;
int lim2 = min(nm, (t + 1) + max(1, window));
for (int pos = t + 1; pos < lim2; pos++) { int pid = o[pos]; auto& p = ps[pid];
for (int ti = 0; ti < (int)p.t.size(); ti++) { EB3& cc = ep_b3[pid * 8 + ti]; if (!cc.valid) continue; auto& ct = p.t[ti]; if (cc.x <= px1 && cc.x + ct.w - 1 >= px0 && cc.y <= py1 && cc.y + ct.h - 1 >= py0) cc.valid = false; } }
}
uint64_t colU = 0; int rows = 0;
for (int y = 0; y < min((int)grid.size(), Hcap); y++) if (grid[y]) { colU |= grid[y]; rows++; }
int Wu = __builtin_popcountll(colU);
R res; res.W = Wu ? Wu : 1; res.H = rows ? rows : 1;
res.A = (long long)res.W * res.H; res.pl = move(pl); res.ok = true; res.Wpack = W;
return res;
}
static R eSmallSolve(RNG& rng, double TLms) {
const double SEARCH_END = TLms;
const double SOFT_END = TLms - 15.0;
const int JUMP = 15, BLF2 = 1, B3 = 1; const long long B2RESTS = 50000;
const int CAPSHARE = 100, CAPBIGN = 6000, CAPDW = 3, CAPDWW = 6, CAPWIDEP = 25, CAPWINDIV = 1, B3WINDIV = 4, LATE = 4, P2SPAN = 0;
const bool big = false;
const long long P2MAXS = 22000;
double P2FRAC = (S < 6000) ? 0.25 : 0.55;
double P2ENDF = (S < 6000) ? 0.35 : 0.70;
vector<int> idx(n); iota(idx.begin(), idx.end(), 0);
auto ord4 = [&]() {
vector<int> res = idx;
stable_sort(res.begin(), res.end(), [&](int a, int b) {
int da = min(ps[a].minW, ps[a].minH); int db = min(ps[b].minW, ps[b].minH);
if (da != db) return da < db;
if (ps[a].k != ps[b].k) return ps[a].k < ps[b].k;
return ps[a].id > ps[b].id;
});
return res;
};
int minW = 0; for (auto& p : ps) minW = max(minW, p.minW);
double factor;
if (S < 1000) factor = 0.4; else if (S < 3000) factor = 0.5; else if (S < 10000) factor = 0.27;
else if (S < 30000) factor = 0.08; else if (S < 50000) factor = 0.028; else factor = 0.009;
int base = max(minW, (int)floor(sqrt((double)S * factor)));
vector<int> Ws;
{
unordered_set<int> used; used.reserve(512);
auto addW = [&](int w) { if (w < minW) w = minW; if (used.insert(w).second) Ws.push_back(w); };
addW(base);
int span = min(96, max(20, base / 2));
for (int d = 1; d <= span; d++) { addW(base - d); addW(base + d); }
addW(minW);
addW((int)max<long long>(minW, (S + base - 1) / base));
for (int m = 2; m <= 6; m++) { addW(base * m / 3); addW((int)max<long long>(minW, S / ((base * m / 3) ? (base * m / 3) : 1))); }
sort(Ws.begin(), Ws.end(), [&](int a, int b) { int da = abs(a - base), db = abs(b - base); if (da != db) return da < db; return a < b; });
}
vector<int> baseOrder = ord4();
auto better = [&](const R& a, const R& b) {
if (!b.ok) return true;
if (a.A != b.A) return a.A < b.A;
if (a.H != b.H) return a.H < b.H;
return a.W < b.W;
};
R bestR;
double tFB0 = elapsed_ms();
bestR = epack(base, baseOrder, rng, false, 1, false, 0.0, -1.0);
double tFB = max(0.05, elapsed_ms() - tFB0);
if (!bestR.ok) {
R r; r.pl.reserve(n);
int y = 0; int wmax = 1;
for (int i = 0; i < n; i++) { auto& t = ps[i].t[0]; r.pl.push_back({i, 0, 0, y}); y += t.h; wmax = max(wmax, t.w); }
r.W = wmax; r.H = y; r.A = 1LL * wmax * y; r.ok = true;
bestR = move(r);
}
vector<int> ordBLF = idx;
stable_sort(ordBLF.begin(), ordBLF.end(), [&](int a, int b) {
if (ps[a].k != ps[b].k) return ps[a].k > ps[b].k;
int ma = max(ps[a].minW, ps[a].minH), mb = max(ps[b].minW, ps[b].minH);
if (ma != mb) return ma > mb;
return ps[a].id < ps[b].id;
});
const vector<int>& ordB2 = baseOrder;
int swin = max(1, n / 4);
{
double SW1F = 0.80;
double sweepBudget = (S < P2MAXS) ? min(SOFT_END, TLms * P2FRAC) : SOFT_END;
double budget1 = SW1F * max(50.0, sweepBudget - elapsed_ms());
double swinPred = budget1 / tFB;
int cap = max(24, (int)swinPred);
if (cap < swin) swin = cap;
}
double avg = 250.0; int cnt = 0;
vector<pair<long long,int>> sweepRes;
bool doPhase2 = (S < P2MAXS);
double SWFRAC = 1.0;
double sweepEnd = doPhase2 ? min(SOFT_END, TLms * P2FRAC) : min(SOFT_END, TLms * SWFRAC);
for (int wi = 0; wi < (int)Ws.size(); wi++) {
double used = elapsed_ms();
if (used + avg * 1.15 > sweepEnd) break;
int W = Ws[wi];
double t1 = elapsed_ms();
double panicAt = min(sweepEnd - 10.0, SEARCH_END - 60.0);
double adaptTL = max(50.0, panicAt - t1 - 10.0);
R r = epack(W, baseOrder, rng, false, swin, true, adaptTL, SEARCH_END, panicAt, 1);
double dt = elapsed_ms() - t1;
cnt++; avg = (avg * (cnt - 1) + dt) / cnt;
if (!r.ok) break;
ecrownRepack(r, SEARCH_END);
if (r.ok) sweepRes.push_back({r.A, W});
if (better(r, bestR)) bestR = move(r);
}
if (doPhase2) {
sort(sweepRes.begin(), sweepRes.end());
double p2Stop = min(SOFT_END, TLms * P2ENDF);
double avg2 = 20.0; int cnt2 = 0;
vector<int> p2W;
{ unordered_set<int> seenW; seenW.reserve(64);
for (auto& sr : sweepRes) if (seenW.insert(sr.second).second) p2W.push_back(sr.second);
int c0 = sweepRes.empty() ? base : sweepRes[0].second;
for (int d = 1; d <= P2SPAN; d++) for (int sgn = -1; sgn <= 1; sgn += 2) { int W = c0 + sgn * d; if (W >= minW && W <= 4000 && seenW.insert(W).second) p2W.push_back(W); } }
for (size_t i = 0; i < p2W.size(); i++) {
double used = elapsed_ms();
if (used + avg2 * 1.2 > p2Stop) break;
int W = p2W[i];
double t1 = elapsed_ms();
int b3win = max(1, B3WINDIV > 0 ? n / B3WINDIV : n);
R r = B3 ? epack_blf3(W, ordB2, b3win, SEARCH_END, rng, false) : epack_blf2(W, ordB2, b3win, SEARCH_END, rng, false);
double dt = elapsed_ms() - t1;
cnt2++; avg2 = (avg2 * (cnt2 - 1) + dt) / cnt2;
if (!r.ok) break;
ecrownRepack(r, SEARCH_END);
if (better(r, bestR)) bestR = move(r);
}
}
{
vector<int> ilsOrd; int ilsW = 0; (void)ilsW;
double avgBLF = 10.0; int cntBLF = 0;
double avgSky = avg;
vector<int> obuf;
while (true) {
double used = elapsed_ms();
if (used + 5.0 > SOFT_END) break;
int bw = bestR.W;
bool doBLF = true;
if (!big && (rng.nxt() & 3) == 0) doBLF = false;
if (doBLF) { if (used + avgBLF * 1.3 > SOFT_END) { if (!big && used + avgSky * 1.3 <= SOFT_END) doBLF = false; else break; } }
else if (used + avgSky * 1.3 > SOFT_END) { if (used + avgBLF * 1.3 <= SOFT_END) doBLF = true; else break; }
bool capEligible = !big || (n <= CAPBIGN);
if (doBLF && capEligible && bestR.Wpack > 0 && bestR.Wpack <= 64 && (int)(rng.nxt() % 100) < CAPSHARE) {
int W;
{ int dw = ((int)(rng.nxt() % 100) < CAPWIDEP) ? CAPDWW : CAPDW;
W = bestR.Wpack + (dw > 0 ? (rng.rint(2 * dw + 1) - dw) : 0);
if (W < minW) W = minW; if (W > 64) W = 64; }
long long capA = bestR.A - 1;
int Hcap = (int)(capA / W);
if (Hcap < minW || Hcap <= 1) { continue; }
bool fromBest = !ilsOrd.empty() && (rng.nxt() & 1);
obuf = fromBest ? ilsOrd : ((cntBLF & 1) ? ordBLF : ordB2);
int swaps = 1 + rng.rint(12);
for (int sswap = 0; sswap < swaps; sswap++) { int a = rng.rint(n), b = rng.rint(n); swap(obuf[a], obuf[b]); }
R r = epack_capped(W, Hcap, obuf, max(1, CAPWINDIV > 0 ? n / CAPWINDIV : n), SEARCH_END, rng);
double t1b = elapsed_ms(); (void)t1b;
cntBLF++; avgBLF = (avgBLF * (cntBLF - 1) + 1.0) / cntBLF;
if (elapsed_ms() > SEARCH_END) break;
if (r.ok) { ecrownRepack(r, SEARCH_END); if (better(r, bestR)) { ilsOrd = obuf; bestR = move(r); } }
} else if (doBLF) {
int W;
if ((int)(rng.nxt() % 100) < JUMP) { double mult = 0.6 + 0.8 * ((double)(rng.nxt()>>11)*(1.0/9007199254740992.0)); W = (int)llround(bw * mult); }
else { W = bw + (rng.rint(9) - 4); }
if (W < minW) W = minW; if (W > 64) W = 64; if (W < minW) break;
bool fromBest = !ilsOrd.empty() && (rng.nxt() & 1);
obuf = fromBest ? ilsOrd : (BLF2 ? ordB2 : ordBLF);
if (cntBLF > 0) { int swaps = fromBest ? (2 + rng.rint(10)) : max(1, n / 4); for (int sswap = 0; sswap < swaps; sswap++) { int a = rng.rint(n), b = rng.rint(n); swap(obuf[a], obuf[b]); } }
int policy = (int)(rng.nxt() & 1);
double t1 = elapsed_ms();
R r = (BLF2 && S < B2RESTS) ? (B3 ? epack_blf3(W, obuf, max(1, n / 4), SEARCH_END, rng, cntBLF > 0) : epack_blf2(W, obuf, max(1, n / 4), SEARCH_END, rng, cntBLF > 0)) : epack_blf(W, obuf, policy, SEARCH_END);
double dt = elapsed_ms() - t1;
cntBLF++; avgBLF = (avgBLF * (cntBLF - 1) + dt) / cntBLF;
if (!r.ok) break;
ecrownRepack(r, SEARCH_END);
if (better(r, bestR)) { ilsOrd = obuf; ilsW = W; bestR = move(r); }
} else {
int W = bw + (rng.rint(7) - 3);
if (W < minW) W = minW;
obuf = baseOrder;
int swaps = max(1, n / 4);
for (int sswap = 0; sswap < swaps; sswap++) { int a = rng.rint(n), b = rng.rint(n); swap(obuf[a], obuf[b]); }
double t1 = elapsed_ms();
R r = epack(W, obuf, rng, true, max(1, n / 4), false, 0.0, SEARCH_END);
double dt = elapsed_ms() - t1;
cnt++; avgSky = (avgSky * 0.7 + dt * 0.3);
if (!r.ok) break;
ecrownRepack(r, SEARCH_END);
if (better(r, bestR)) bestR = move(r);
}
}
}
if (bestR.ok && bestR.Wpack > 0 && bestR.Wpack <= 64 && elapsed_ms() < SOFT_END - 30 && LATE > 0) {
int lateW = bestR.Wpack;
int lateHcap = max(1, (int)(bestR.A / lateW) - 1);
if (lateHcap >= minW) {
for (int attempt = 0; attempt < LATE && elapsed_ms() < SOFT_END - 10; attempt++) {
vector<int> randOrd = idx;
int ns = max(1, n / 2);
for (int s = 0; s < ns; s++) { int a = rng.rint(n), b = rng.rint(n); swap(randOrd[a], randOrd[b]); }
R r = epack_capped(lateW, lateHcap, randOrd, max(1, n / 4), SEARCH_END, rng);
if (r.ok) { ecrownRepack(r, SEARCH_END); if (better(r, bestR)) bestR = move(r); }
}
}
}
if (bestR.ok) ecrownRepack(bestR, TLms + 10.0);
return bestR;
}
int main(){
ios::sync_with_stdio(false);
cin.tie(nullptr);
auto t0=chrono::steady_clock::now();
gT0=t0;
const double TL=1950.0;
BF_TL=TL;
if(!(cin>>n)) return 0;
ps.assign(n,P{}); S=0;
for(int i=0;i<n;i++){
int k;cin>>k; ps[i].id=i+1; ps[i].k=k; ps[i].b.resize(k);
for(int j=0;j<k;j++){int x,y;cin>>x>>y; ps[i].b[j]={x,y};}
S+=k;
}
for(int i=0;i<n;i++){
auto &p=ps[i]; unordered_set<string> seen; seen.reserve(32);
for(int rf=0;rf<2;rf++){
vector<pair<int,int>> src=p.b; if(rf){for(auto &q:src) q.first=-q.first;}
for(int r=0;r<4;r++){
vector<pair<int,int>> v=src; for(auto &q:v) q=rotp(q,r);
int minx=INT_MAX,miny=INT_MAX,maxx=INT_MIN,maxy=INT_MIN;
for(auto &q:v){minx=min(minx,q.first);miny=min(miny,q.second);maxx=max(maxx,q.first);maxy=max(maxy,q.second);}
vector<pair<int,int>> v2=v; for(auto &q:v2){q.first-=minx;q.second-=miny;}
sort(v2.begin(),v2.end());
string key; key.reserve(v2.size()*8);
for(auto &q:v2){key.append(to_string(q.first));key.push_back(',');key.append(to_string(q.second));key.push_back(';');}
if(seen.insert(key).second){
T t; t.w=maxx-minx+1; t.h=maxy-miny+1; t.c=v2; t.r=r; t.f=rf; t.minx=minx; t.miny=miny;
t.lo.assign(t.w,INT_MAX); t.hi.assign(t.w,INT_MIN);
for(auto &q:v2){int x=q.first,y=q.second; if(t.lo[x]>y) t.lo[x]=y; if(t.hi[x]<y) t.hi[x]=y;}
t.rm.assign(t.h,0u); for(auto &q:v2){ t.rm[q.second]|=(1u<<q.first); }
{ set<pair<int,int>> cellset(v2.begin(),v2.end()); set<pair<int,int>> nb;
for(auto&q:v2){ const int dx4[4]={1,-1,0,0},dy4[4]={0,0,1,-1}; for(int d=0;d<4;d++){ pair<int,int> np={q.first+dx4[d],q.second+dy4[d]}; if(!cellset.count(np)) nb.insert(np); } }
t.nbr.assign(nb.begin(),nb.end()); }
p.t.push_back(move(t));
}
}
}
for(auto &t:p.t){p.minW=min(p.minW,t.w); p.minH=min(p.minH,t.h); p.minA=min(p.minA,t.w*t.h);}
if(p.t.empty()){T t; t.w=1;t.h=1;t.c={{0,0}};t.lo={0};t.hi={0};t.r=0;t.f=0;t.minx=0;t.miny=0;t.rm={1u};p.t.push_back(t);p.minW=1;p.minH=1;p.minA=1;}
}
vector<int> idx(n); iota(idx.begin(),idx.end(),0);
bool useContact=(S<=12000);
bool useRough=(S<=50000);
unsigned long long seed=(S<<1) ^ (unsigned long long)(n*1469598103934665603ULL) ^ 0x9e3779b97f4a7c15ULL;
RNG rng(seed);
auto ord4 = [&]() {
vector<int> res = idx;
stable_sort(res.begin(), res.end(), [&](int a, int b) {
int da = min(ps[a].minW, ps[a].minH);
int db = min(ps[b].minW, ps[b].minH);
if (da != db) return da < db;
if (ps[a].k != ps[b].k) return ps[a].k < ps[b].k;
return ps[a].id > ps[b].id;
});
return res;
};
auto pack=[&](int W,const vector<int>& o0,RNG &rng,bool randtie){
vector<int> h(W,-1);
long long g=-1;
vector<Pl> pl; pl.reserve(o0.size());
vector<int> o=o0;
int t=0,nm=(int)o.size();
bool big=S>7000;
int maxBound=max(1,n/4);
int expLIM=min(maxBound,(int)(600000/max(1LL,S-3500)));
int dynLIM=big?expLIM:maxBound;
auto tStart=chrono::steady_clock::now();
auto batchStart=tStart;
double gEl=chrono::duration<double,milli>(tStart-t0).count();
long long TLms=big?(long long)max(150.0,(TL-55.0)-gEl):LLONG_MAX/4;
int stepCnt=0;
while(t<nm){
int limCnt=max(1,dynLIM);
int lim=min(nm,t+limCnt);
long long bestg=LLONG_MAX; int bti=-1,bx=0,by=0,bl=INT_MAX; long long bds=LLONG_MAX; long long bdr=LLONG_MAX; int by0=INT_MAX,bx0=INT_MAX; int bestpos=t; int bestid=-1; int bcon=-1;
for(int pos=t;pos<lim;pos++){
int id=o[pos];
auto &p=ps[id];
long long bestg2=LLONG_MAX; int bti2=-1,bx2=0,by2=0,bl2=INT_MAX; long long bds2=LLONG_MAX; long long bdr2=LLONG_MAX; int by02=INT_MAX,bx02=INT_MAX; int bcon2=-1;
for(int ti=0;ti<(int)p.t.size();ti++){
auto &tsh=p.t[ti]; if(tsh.w>W) continue; int Rpos=W-tsh.w+1;
for(int x0=0;x0<Rpos;x0++){
int y0=0;
for(int j=0;j<tsh.w;j++){
if(tsh.lo[j]!=INT_MAX){int v=h[x0+j]-tsh.lo[j]+1; if(v>y0) y0=v;}
}
int nhbuf[32];
int l=-1; long long dsum=0;
for(int j=0;j<tsh.w;j++){
int nh=h[x0+j];
if(tsh.hi[j]!=INT_MIN){
int cand=y0+tsh.hi[j];
if(nh<cand) nh=cand;
}
nhbuf[j]=nh;
if(nh>l) l=nh;
}
long long gg=g; if(gg<l) gg=l;
if(gg>bestg2) continue;
for(int j=0;j<tsh.w;j++){
int inc=nhbuf[j]-h[x0+j];
if(inc>0) dsum+=inc;
}
int con=0;
if(useContact){
for(int j=0;j<tsh.w;j++){ if(tsh.lo[j]!=INT_MAX && y0+tsh.lo[j]==h[x0+j]+1) con++; }
if(x0==0) con++;
if(x0+tsh.w==W) con++;
}
long long dr=0;
if(useRough){
if(x0>0){
long long old=llabs((long long)h[x0]-h[x0-1]);
long long nw=llabs((long long)nhbuf[0]-h[x0-1]);
dr+=nw-old;
}
for(int j=0;j<tsh.w-1;j++){
long long old=llabs((long long)h[x0+j+1]-h[x0+j]);
long long nw=llabs((long long)nhbuf[j+1]-nhbuf[j]);
dr+=nw-old;
}
if(x0+tsh.w<W){
long long old=llabs((long long)h[x0+tsh.w]-h[x0+tsh.w-1]);
long long nw=llabs((long long)h[x0+tsh.w]-nhbuf[tsh.w-1]);
dr+=nw-old;
}
}
bool take=false;
if(gg<bestg2) take=true;
else if(gg==bestg2){
if(dsum<bds2) take=true;
else if(dsum==bds2){
if(con>bcon2) take=true;
else if(con==bcon2){
if(l<bl2) take=true;
else if(l==bl2){
if(dr<bdr2) take=true;
else if(dr==bdr2){
if(y0<by02) take=true;
else if(y0==by02){
if(x0<bx02) take=true;
else if(x0==bx02 && randtie && rng.coin()) take=true;
}
}
}
}
}
}
if(take){bestg2=gg;bti2=ti;bx2=x0;by2=y0;bl2=l;bds2=dsum;bdr2=dr;by02=y0;bx02=x0;bcon2=con;}
}
}
if(bti2==-1) continue;
bool take=false;
if(bestg2<bestg) take=true;
else if(bestg2==bestg){
if(bds2<bds) take=true;
else if(bds2==bds){
if(bcon2>bcon) take=true;
else if(bcon2==bcon){
if(bl2<bl) take=true;
else if(bl2==bl){
if(bdr2<bdr) take=true;
else if(bdr2==bdr){
if(by02<by0) take=true;
else if(by02==by0){
if(bx02<bx0) take=true;
else if(bx02==bx0 && randtie && rng.coin()) take=true;
}
}
}
}
}
}
if(take){bestg=bestg2;bti=bti2;bx=bx2;by=by2;bl=bl2;bds=bds2;bdr=bdr2;by0=by02;bx0=bx02;bestpos=pos;bestid=id;bcon=bcon2;}
}
if(bti==-1){t++;stepCnt++;if(big&&stepCnt==5){auto now=chrono::steady_clock::now(); auto elapsed=chrono::duration<double,milli>(now-tStart).count(); auto batch=chrono::duration<double,milli>(now-batchStart).count(); double remT=max(0.0,TLms-elapsed); int remSteps=max(1,nm-t); double budget=remT*5.0/remSteps; if(batch<budget) dynLIM=min(maxBound,dynLIM+1); else dynLIM=max(1,dynLIM-1); batchStart=now; stepCnt=0;} continue;}
auto &tsh=ps[bestid].t[bti];
for(int j=0;j<tsh.w;j++){
int nh=h[bx+j];
if(tsh.hi[j]!=INT_MIN){
int cand=by+tsh.hi[j];
if(nh<cand) nh=cand;
}
h[bx+j]=nh;
}
if(g<bl) g=bl;
pl.push_back({bestid,bti,bx,by});
if(bestpos!=t) swap(o[t],o[bestpos]);
t++;
stepCnt++;
if(big&&stepCnt==5){
auto now=chrono::steady_clock::now();
auto elapsed=chrono::duration<double,milli>(now-tStart).count();
auto batch=chrono::duration<double,milli>(now-batchStart).count();
double remT=max(0.0,TLms-elapsed);
int remSteps=max(1,nm-t);
double budget=remT*5.0/remSteps;
if(batch<budget) dynLIM=min(maxBound,dynLIM+1); else dynLIM=max(1,dynLIM-1);
batchStart=now;
stepCnt=0;
}
}
int H=(int)g+1;
int maxX=-1;
for(auto &pp:pl){
auto &t=ps[pp.idx].t[pp.ti];
for(auto &q:t.c){int x=pp.x+q.first; if(x>maxX) maxX=x;}
}
int Wused=0;
if(maxX>=0){
vector<char> used(maxX+1,false);
for(auto &pp:pl){auto &t=ps[pp.idx].t[pp.ti]; for(auto &q:t.c){used[pp.x+q.first]=true;}}
for(int x=0;x<=maxX;x++) if(used[x]) Wused++;
}
int Wfinal=max(0,Wused);
long long A=1LL*H*max(1,Wfinal);
return R{A,max(1,Wfinal),H,W,move(pl)};
};
struct BC{int y,x,contact;bool valid;};
vector<BC> blfcache;
auto pack_blf=[&](int W,const vector<int>& order,int window,double deadlineMs,RNG& rng,bool randtie)->R{
if(W>64||W<=0) return R{LLONG_MAX,1,1,W,{}};
int nm=(int)order.size();
size_t estH=(size_t)(S/max(1,W))+48;
vector<unsigned long long> grid(max<size_t>(64,estH),0ULL);
const unsigned long long FULLROW=(W==64)?~0ULL:((1ULL<<W)-1);
int y0full=0,topY=0; long long g=-1;
vector<int> o=order; vector<Pl> pl; pl.reserve(nm);
if((int)blfcache.size()<n*8) blfcache.resize(n*8);
for(int i=0;i<n*8;i++) blfcache[i].valid=false;
int checkCnt=0;
for(int t=0;t<nm;t++){
if(deadlineMs>0 && ((++checkCnt&15)==0) && chrono::duration<double,milli>(chrono::steady_clock::now()-t0).count()>deadlineMs) return R{LLONG_MAX,1,1,W,{}};
int lim=min(nm,t+max(1,window));
long long bG=LLONG_MAX; int bContact=-1,bY=INT_MAX,bX=INT_MAX,bTi=-1,bPos=-1;
if((int)grid.size()<topY+12) grid.resize(topY+64,0ULL);
for(int pos=t;pos<lim;pos++){
int pid=o[pos]; auto&p=ps[pid];
for(int ti=0;ti<(int)p.t.size();ti++){
auto&tt=p.t[ti]; if(tt.w>W) continue;
BC&cc=blfcache[pid*8+ti];
if(!cc.valid){
unsigned long long limMask=(W-tt.w+1>=64)?~0ULL:((1ULL<<(W-tt.w+1))-1);
int yF=-1,xF=-1;
for(int y=y0full;y<=topY;y++){
unsigned long long conflict=0;
for(int dy=0;dy<tt.h;dy++){ unsigned long long gg=grid[y+dy]; if(!gg)continue; unsigned m=tt.rm[dy]; while(m){int b=__builtin_ctz(m); conflict|=(gg>>b); m&=m-1;} }
unsigned long long allowed=~conflict & limMask;
if(allowed){ yF=y; xF=(int)__builtin_ctzll(allowed); break; }
}
if(yF<0) continue;
int contact=0;
for(auto&nb:tt.nbr){ int cx=xF+nb.first,cy=yF+nb.second; if(cx<0||cx>=W||cy<0){contact++;continue;} if(cy<(int)grid.size() && ((grid[cy]>>cx)&1ULL)) contact++; }
cc={yF,xF,contact,true};
}
long long g2=max(g,(long long)(cc.y+tt.h-1));
long long s2=g2-cc.contact;
long long bs=(bTi>=0)?(bG-bContact):LLONG_MAX;
if(s2>bs) continue;
bool take=false;
if(s2<bs) take=true;
else { if(cc.contact>bContact) take=true; else if(cc.contact==bContact){ if(cc.y<bY) take=true; else if(cc.y==bY){ if(cc.x<bX) take=true; else if(cc.x==bX && randtie && rng.coin()) take=true; } } }
if(take){ bG=g2; bContact=cc.contact; bY=cc.y; bX=cc.x; bTi=ti; bPos=pos; }
}
}
if(bPos<0) return R{LLONG_MAX,1,1,W,{}};
int placedId=o[bPos]; auto&tt=ps[placedId].t[bTi];
if((int)grid.size()<bY+tt.h+1) grid.resize(bY+tt.h+64,0ULL);
for(int dy=0;dy<tt.h;dy++) grid[bY+dy]|=((unsigned long long)tt.rm[dy])<<bX;
if(bY+tt.h>topY) topY=bY+tt.h;
while(y0full<topY && grid[y0full]==FULLROW) y0full++;
if(bG>g) g=bG;
pl.push_back({placedId,bTi,bX,bY});
if(bPos!=t) swap(o[t],o[bPos]);
int px0=bX-1,px1=bX+tt.w,py0=bY-1,py1=bY+tt.h;
int lim2=min(nm,(t+1)+max(1,window));
for(int pos=t+1;pos<lim2;pos++){
int pid=o[pos]; auto&p=ps[pid];
for(int ti=0;ti<(int)p.t.size();ti++){ BC&cc=blfcache[pid*8+ti]; if(!cc.valid)continue; auto&ct=p.t[ti]; if(cc.x<=px1 && cc.x+ct.w-1>=px0 && cc.y<=py1 && cc.y+ct.h-1>=py0) cc.valid=false; }
}
}
unsigned long long colU=0; int firstY=-1,lastY=-1;
for(int y=0;y<topY;y++) if(grid[y]){ colU|=grid[y]; if(firstY<0)firstY=y; lastY=y; }
int Wused=max(1,__builtin_popcountll(colU));
int H=(firstY<0)?1:max(1,lastY-firstY+1);
long long A=(long long)Wused*H;
return R{A,Wused,H,W,move(pl)};
};
int minW=0; for(auto &p:ps) minW=max(minW,p.minW);
double factor;
if (S < 1000) {
factor = 0.4;
} else if (S < 3000) {
factor = 0.5;
} else if (S < 10000) {
factor = 0.27;
} else if (S < 30000) {
factor = 0.08;
} else if (S < 50000) {
factor = 0.028;
} else if (S <= 62000) {
factor = 0.009;
} else {
factor = 0.01;
}
int base = max(minW, (int)floor(sqrt((double)S * factor)));
vector<int> Ws;
{
unordered_set<int> used; used.reserve(512);
auto addW=[&](int w){if(w<minW) w=minW; if(used.insert(w).second) Ws.push_back(w);};
addW(base);
int span=min(96,max(20,base/2));
for(int d=1;d<=span;d++){addW(base-d); addW(base+d);}
addW(minW);
addW((int)max<long long>(minW,(S+base-1)/base));
for(int m=2;m<=6;m++){addW(base*m/3); addW((int)max<long long>(minW,S/((base*m/3)?(base*m/3):1)));}
sort(Ws.begin(),Ws.end(),[&](int a,int b){int da=abs(a-base),db=abs(b-base); if(da!=db) return da<db; return a<b;});
}
long long bestA=LLONG_MAX; int bestW=0,bestH=0; R bestR; bool hasBestmine=false;
BF_HD = (S>30000 || S<10000) ? 10 : 46;
bool useBF = (n>=100 && base<=63 && minW<=63 && S<=62000);
bool useBFfast = (n>=100 && base<=63 && minW<=63 && S>62000);
bool haveBF=false;
bool smallN = (n<180);
if(useBF){
double bfDL = smallN ? ((n<=115)?200.0:((n<=140)?30.0:1.0)) : (TL-70.0);
R bfR=bfSolve(minW, min(base,63), bfDL);
if((int)bfR.pl.size()==n){
bestR=move(bfR); bestA=bestR.A; bestW=bestR.W; bestH=bestR.H; hasBestmine=true; haveBF=true;
if(!smallN && bestR.ok){
double cdl = (n<4000) ? (TL-40.0) : (TL-100.0);
R q=bestR;
if(ecrownRepack2(q, cdl) && q.A<bestR.A){ bestR=move(q); bestA=bestR.A; bestW=bestR.W; bestH=bestR.H; }
}
}
if(smallN){
RNG erng(0x9e3779b97f4a7c15ULL ^ ((unsigned long long)S<<1) ^ (unsigned long long)(n*1469598103934665603ULL));
R eR=eSmallSolve(erng, TL-80.0);
if(eR.ok && (int)eR.pl.size()==n && (!hasBestmine || eR.A<bestA)){
bestR=move(eR); bestA=bestR.A; bestW=bestR.W; bestH=bestR.H; hasBestmine=true; haveBF=true;
}
}
} else if(useBFfast){
BF_FAST=1; BF_HD=3;
R bfR=bfSolve(minW, min(base,63), TL-120.0);
BF_FAST=0;
if((int)bfR.pl.size()==n){
bestR=move(bfR); bestA=bestR.A; bestW=bestR.W; bestH=bestR.H; hasBestmine=true; haveBF=true;
if(bestR.ok){
R q=bestR;
if(ecrownRepack2(q, TL-100.0) && q.A<bestR.A){ bestR=move(q); bestA=bestR.A; bestW=bestR.W; bestH=bestR.H; }
}
}
}
double avg=250.0; int cnt=0;
double sweepTL = (S<=2200) ? TL*0.42 : TL;
if(!haveBF)
for(int wi=0;wi<(int)Ws.size();wi++){
auto now=chrono::steady_clock::now();
double used=chrono::duration<double,milli>(now-t0).count();
if(used+avg*1.3>sweepTL) break;
int W=Ws[wi];
vector<vector<int>> orders;
orders.push_back(ord4());
int oi=0;
while(oi<(int)orders.size()){
auto t1=chrono::steady_clock::now();
double used2=chrono::duration<double,milli>(t1-t0).count();
if(used2+avg*1.15>sweepTL) break;
bool randtie=(oi>=1);
R r=pack(W,orders[oi],rng,randtie);
auto t2=chrono::steady_clock::now();
double dt=chrono::duration<double,milli>(t2-t1).count();
cnt++; avg=(avg*(cnt-1)+dt)/cnt;
if((int)r.pl.size()==n && (!hasBestmine || r.A<bestA || (r.A==bestA && (r.H<bestH || (r.H==bestH && r.W<bestW))))){
bestA=r.A; bestW=r.W; bestH=r.H; bestR=r; hasBestmine=true;
}
oi++;
}
}
if(!hasBestmine){
int W=max(minW,(int)floor(sqrt((double)S)));
auto o=ord4();
R r=pack(W,o,rng,false);
bestA=r.A; bestW=r.W; bestH=r.H; bestR=r; hasBestmine=true;
}
if(hasBestmine && !haveBF && S<=12000){
auto bord=ord4();
for(int wi=0; wi<(int)Ws.size(); wi++){
double now=chrono::duration<double,milli>(chrono::steady_clock::now()-t0).count();
if(now > TL-500) break;
int W=Ws[wi]; if(W<1||W>64) continue;
R rb=pack_blf(W,bord,max(1,n/4),TL-120.0,rng,false);
if((int)rb.pl.size()==n && (rb.A<bestR.A || (rb.A==bestR.A && rb.H<bestR.H))){
bestR=rb; bestA=rb.A; bestW=rb.W; bestH=rb.H;
}
}
}
if(hasBestmine && bestR.Wpack>=1 && bestR.Wpack<=64){
int W=bestR.Wpack; auto &pl=bestR.pl;
auto computeH=[&](){ int H=0; for(auto&p:pl){auto&t=ps[p.idx].t[p.ti]; if(p.y+t.h>H)H=p.y+t.h;} return H; };
int H0=computeH();
if(H0>1){
vector<unsigned long long> grid(H0,0ULL);
auto stamp=[&](Pl&p,bool on){ auto&t=ps[p.idx].t[p.ti]; for(int dy=0;dy<t.h;dy++){ unsigned long long m=((unsigned long long)t.rm[dy])<<p.x; if(on) grid[p.y+dy]|=m; else grid[p.y+dy]&=~m; } };
for(auto&p:pl) stamp(p,true);
long long origA=bestR.A; vector<Pl> origPl=pl; bool improvedAny=false; int rounds=0;
while(rounds++<64){
if(chrono::duration<double,milli>(chrono::steady_clock::now()-t0).count()>TL-30) break;
int H=(int)grid.size(); while(H>0&&grid[H-1]==0ULL) H--;
if(H<=1) break;
bool anyDepthOk=false;
for(int d=1;d<=3;d++){
if(H-d<1) break;
vector<int> crown;
for(int i=0;i<(int)pl.size();i++){ auto&t=ps[pl[i].idx].t[pl[i].ti]; if(pl[i].y+t.h>H-d) crown.push_back(i); }
if(crown.empty()){ anyDepthOk=true; break; }
for(int ci:crown) stamp(pl[ci],false);
sort(crown.begin(),crown.end(),[&](int a,int b){ return ps[pl[a].idx].k>ps[pl[b].idx].k; });
int targetH=H-d; vector<array<int,3>> newPos(crown.size()); bool ok=true;
for(size_t c=0;c<crown.size()&&ok;c++){
auto&pcur=pl[crown[c]]; auto&piece=ps[pcur.idx];
int bTi=-1,bX=0,bY=INT_MAX;
for(int ti=0;ti<(int)piece.t.size();ti++){
auto&t=piece.t[ti]; if(t.w>W||t.h>targetH) continue;
unsigned long long limMask=(W-t.w+1>=64)?~0ULL:((1ULL<<(W-t.w+1))-1);
int ymax=targetH-t.h;
for(int y=0;y<=ymax;y++){
if(bY<=y) break;
unsigned long long conflict=0;
for(int dy=0;dy<t.h;dy++){ unsigned long long g=grid[y+dy]; if(!g) continue; unsigned m=t.rm[dy]; while(m){ int b=__builtin_ctz(m); conflict|=(g>>b); m&=m-1; } }
unsigned long long allowed=~conflict & limMask;
if(allowed){ int x=(int)__builtin_ctzll(allowed); if(y<bY||(y==bY&&x<bX)){ bY=y;bX=x;bTi=ti; } break; }
}
}
if(bTi<0){ ok=false; break; }
auto&t=piece.t[bTi];
for(int dy=0;dy<t.h;dy++) grid[bY+dy]|=((unsigned long long)t.rm[dy])<<bX;
newPos[c]={bTi,bX,bY};
}
if(ok){
unsigned long long colU=0; int rows=0;
for(auto&g:grid) if(g){ colU|=g; rows++; }
int newW=max(1,__builtin_popcountll(colU)), newH=max(1,rows);
long long newA=(long long)newW*newH;
if(newA<bestR.A){
for(size_t c=0;c<crown.size();c++){ auto&pc2=pl[crown[c]]; pc2.ti=newPos[c][0]; pc2.x=newPos[c][1]; pc2.y=newPos[c][2]; }
bestR.W=newW; bestR.H=newH; bestR.A=newA; improvedAny=true; anyDepthOk=true; break;
} else ok=false;
}
if(!ok){ fill(grid.begin(),grid.end(),0ULL); for(auto&p:pl) stamp(p,true); }
}
if(!anyDepthOk) break;
}
if(!improvedAny){ pl=origPl; bestR.A=origA; }
}
}
int maxX=-1;
for(auto &p:bestR.pl){
auto &t=ps[p.idx].t[p.ti];
for(auto &q:t.c){int x=p.x+q.first; if(x>maxX) maxX=x;}
}
vector<int> mapx(maxX+1,-1);
if(maxX>=0){
vector<char> used(maxX+1,false);
for(auto &p:bestR.pl){
auto &t=ps[p.idx].t[p.ti];
for(auto &q:t.c){used[p.x+q.first]=true;}
}
int cur=0;
for(int x=0;x<=maxX;x++) if(used[x]) mapx[x]=cur++;
}
vector<array<int,4>> ans(n,{0,0,0,0});
for(auto &p:bestR.pl){
auto &t=ps[p.idx].t[p.ti];
int bx = (mapx.empty()?p.x:mapx[p.x]);
int Xi = bx - t.minx;
int Yi = p.y - t.miny;
int Ri = (4 - (t.r%4) + 4) % 4;
int Fi = t.f;
ans[p.idx]={Xi,Yi,Ri,Fi};
}
cout<<bestR.W<<" "<<bestR.H<<"\n";
for(int i=0;i<n;i++){
cout<<ans[i][0]<<" "<<ans[i][1]<<" "<<ans[i][2]<<" "<<ans[i][3]<<"\n";
}
return 0;
}
