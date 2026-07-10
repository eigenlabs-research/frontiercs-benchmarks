// Model South Core
#include <bits/stdc++.h>
using namespace std; static chrono::steady_clock::time_point T0; static double TL_MS = 2400.0; static inline double el_ms(){ return chrono::duration<double,milli>(chrono::steady_clock::now()-T0).count(); }
static int N; static vector<double> X, Y; static inline double dist(int a,int b){ double dx=X[a]-X[b],dy=Y[a]-Y[b]; return sqrt(dx*dx+dy*dy); }
static vector<bool> is_prime_city; static void init_primes(int max_n) { is_prime_city.assign(max_n + 1, true); is_prime_city[0] = is_prime_city[1] = false; for(int p=2; p*p<=max_n; p++) { if(is_prime_city[p]) { for(int i=p*p; i<=max_n; i+=p) is_prime_city[i] = false; }
}
}
static inline double penalize(int c1, int c2, int step_offset) { double d = dist(c1, c2); if(step_offset == 9 && !is_prime_city[c1]) return d * 1.1; return d; }
constexpr int PORTAL_B = 400; struct __m256_struct { double fwd[10]; double rev[10]; int first_city; int last_city; }; static vector<__m256_struct> portals; static void rebuild_portals(int N, const vector<int>& order) { int num_portals = (N + PORTAL_B - 1) / PORTAL_B; portals.resize(num_portals); for(int b=0; b<num_portals; b++) { int L = b * PORTAL_B; int R = min(N - 1, L + PORTAL_B - 1); portals[b].first_city = order[L]; portals[b].last_city = order[R]; for(int off=0; off<10; off++) { double cost_f = 0; for(int i=L; i<R; i++) { cost_f += penalize(order[i], order[i+1], (off + i - L) % 10); }
portals[b].fwd[off] = cost_f; double cost_r = 0; for(int i=R; i>L; i--) { cost_r += penalize(order[i], order[i-1], (off + R - i) % 10); }
portals[b].rev[off] = cost_r; }
}
}
static double _ptr_x86_sys(int L, int R, int off, int N, const vector<int>& order) { if (L > R) return 0.0; double cost = 0; int cur_off = off; int i = L; while (i <= R && (i % PORTAL_B) != 0) { if (i == R) break; cost += penalize(order[i], order[i+1], cur_off); cur_off = (cur_off + 1) % 10; i++; }
while (i + PORTAL_B - 1 <= R) { int b = i / PORTAL_B; cost += portals[b].fwd[cur_off]; cur_off = (cur_off + PORTAL_B - 1) % 10; if (i + PORTAL_B <= R) { cost += penalize(portals[b].last_city, portals[b+1].first_city, cur_off); cur_off = (cur_off + 1) % 10; }
i += PORTAL_B; }
while (i < R) { cost += penalize(order[i], order[i+1], cur_off); cur_off = (cur_off + 1) % 10; i++; }
return cost; }
static double _ptr_x86_rev(int L, int R, int off, int N, const vector<int>& order) { if (L > R) return 0.0; double cost = 0; int cur_off = off; int i = R; while (i >= L && ((i + 1) % PORTAL_B) != 0) { if (i == L) break; cost += penalize(order[i], order[i-1], cur_off); cur_off = (cur_off + 1) % 10; i--; }
while (i - PORTAL_B + 1 >= L) { int b = (i - PORTAL_B + 1) / PORTAL_B; cost += portals[b].rev[cur_off]; cur_off = (cur_off + PORTAL_B - 1) % 10; if (i - PORTAL_B >= L) { cost += penalize(portals[b].first_city, portals[b-1].last_city, cur_off); cur_off = (cur_off + 1) % 10; }
i -= PORTAL_B; }
while (i > L) { cost += penalize(order[i], order[i-1], cur_off); cur_off = (cur_off + 1) % 10; i--; }
return cost; }
static double _sys_call_opt(int N, const vector<int>& order) { double cost = 0; for(int i=0; i<N; i++) { cost += penalize(order[i], order[i+1==N?0:i+1], i % 10); }
return cost; }
static double get_2opt_delta(int p1, int p2, int p3, int p4, int N, const vector<int>& order, int dir) { int A=(dir==0)?p1:p4, B=(dir==0)?p3:p2; int lo=A, hi=B; if(lo>hi) swap(lo,hi); double old_edge1 = penalize(order[lo], order[lo+1], lo % 10); double old_edge2 = (hi+1 < N) ? penalize(order[hi], order[hi+1], hi % 10) : penalize(order[hi], order[0], hi % 10); double old_chunk = _ptr_x86_sys(lo+1, hi, (lo+1)%10, N, order); double old_total = old_edge1 + old_edge2 + old_chunk; double new_edge1 = penalize(order[lo], order[hi], lo % 10); double new_edge2 = (hi+1 < N) ? penalize(order[lo+1], order[hi+1], hi % 10) : penalize(order[lo+1], order[0], hi % 10); double new_chunk = _ptr_x86_rev(lo+1, hi, (lo+1)%10, N, order); double new_total = new_edge1 + new_edge2 + new_chunk; return old_total - new_total; }
static vector<char> inbuf; static size_t ip=0; static inline long long readLL(){ while(ip<inbuf.size() && (inbuf[ip]<'0'||inbuf[ip]>'9') && inbuf[ip]!='-') ip++; bool neg=false; if(ip<inbuf.size()&&inbuf[ip]=='-'){neg=true;ip++;}
long long v=0; while(ip<inbuf.size()&&inbuf[ip]>='0'&&inbuf[ip]<='9'){v=v*10+(inbuf[ip]-'0');ip++;}
return neg?-v:v; }
int main(){ T0=chrono::steady_clock::now(); if(const char* e=getenv("SANTA_TL")){ double v=atof(e); if(v>50&&v<10000) TL_MS=v; }
{ size_t cap=1<<20; inbuf.resize(cap); size_t len=0; while(true){ if(len==cap){cap<<=1;inbuf.resize(cap);} size_t g=fread(inbuf.data()+len,1,cap-len,stdin); if(!g)break; len+=g; } inbuf.resize(len); }
N=(int)readLL(); if(N<=0){ printf("1\n0\n"); return 0; }
X.resize(N); Y.resize(N); for(int i=0;i<N;i++){ X[i]=(double)readLL(); Y[i]=(double)readLL(); }
if(N==1){ printf("2\n0\n0\n"); return 0; }
if(N==2){ printf("3\n0\n1\n0\n"); return 0; }
vector<char> prOrig((size_t)N,0); { vector<char> comp((size_t)N,0); for(long long i=2;i<N;i++) if(!comp[i]){ prOrig[i]=1; for(long long q=i*i;q<N;q+=i) comp[q]=1; } }
vector<int> perm(N), invPerm(N); { double mnx=X[0],mxx=X[0],mny=Y[0],mxy=Y[0]; for(int i=1;i<N;i++){ mnx=min(mnx,X[i]); mxx=max(mxx,X[i]); mny=min(mny,Y[i]); mxy=max(mxy,Y[i]); }
double rx=max(1e-9,mxx-mnx), ry=max(1e-9,mxy-mny); auto expandBits=[](uint32_t v)->uint64_t{ uint64_t x=v & 0xFFFFFu; x=(x|(x<<16))&0x0000FFFF0000FFFFULL; x=(x|(x<<8))&0x00FF00FF00FF00FFULL; x=(x|(x<<4))&0x0F0F0F0F0F0F0F0FULL; x=(x|(x<<2))&0x3333333333333333ULL; x=(x|(x<<1))&0x5555555555555555ULL; return x; }; vector<uint64_t> morton(N); for(int i=0;i<N;i++){ uint32_t qx=(uint32_t)min(1048575.0, floor((X[i]-mnx)/rx*1048575.0)); uint32_t qy=(uint32_t)min(1048575.0, floor((Y[i]-mny)/ry*1048575.0)); morton[i]=expandBits(qx)|(expandBits(qy)<<1); perm[i]=i; }
if(N>=150000){ sort(perm.begin(),perm.end(),[&](int a,int b){ return morton[a]<morton[b]; }); }
for(int newId=0;newId<N;newId++) invPerm[perm[newId]]=newId; vector<double> X2(N),Y2(N); for(int newId=0;newId<N;newId++){ X2[newId]=X[perm[newId]]; Y2[newId]=Y[perm[newId]]; }
X.swap(X2); Y.swap(Y2); }
const int zeroNewId = invPerm[0]; vector<char> pr(N); for(int newId=0;newId<N;newId++) pr[newId]=prOrig[perm[newId]]; if(N>100000) TL_MS -= 20.0; double RESERVE = N>150000?220.0:(N>50000?90.0:(N>5000?50.0:40.0)); TL_MS -= RESERVE; double minx=X[0],maxx=X[0],miny=Y[0],maxy=Y[0]; for(int i=1;i<N;i++){ minx=min(minx,X[i]);maxx=max(maxx,X[i]);miny=min(miny,Y[i]);maxy=max(maxy,Y[i]); }
double w=max(1.0,maxx-minx), h=max(1.0,maxy-miny); int G=max(1,(int)floor(sqrt((double)N/2.0))); double cw=w/G, ch=h/G; auto gx=[&](double x){ int c=(int)((x-minx)/cw); return c<0?0:(c>=G?G-1:c); }; auto gy=[&](double y){ int c=(int)((y-miny)/ch); return c<0?0:(c>=G?G-1:c); }; vector<int> cellOf(N), cnt(G*G+1,0); for(int i=0;i<N;i++){ int c=gx(X[i])*G+gy(Y[i]); cellOf[i]=c; cnt[c+1]++; }
for(int i=0;i<G*G;i++) cnt[i+1]+=cnt[i]; vector<int> bucket(N); { vector<int> tmp=cnt; for(int i=0;i<N;i++) bucket[tmp[cellOf[i]]++]=i; }
int K; if(N>50000)                      K=min(N-1,6); else if(N>=16000 && N<=36000)    K=min(N-1,40); else if(N>5000)                  K=min(N-1,24); else                             K=min(N-1,10); if(const char* e=getenv("K_FORCE")) K=min(N-1, atoi(e)); vector<int> nbr((size_t)N*K,-1); vector<int> lzorder(N); { double mnx=X[0],mxx=X[0],mny=Y[0],mxy=Y[0]; for(int i=1;i<N;i++){ mnx=min(mnx,X[i]); mxx=max(mxx,X[i]); mny=min(mny,Y[i]); mxy=max(mxy,Y[i]); }
double rx=max(1e-9,mxx-mnx), ry=max(1e-9,mxy-mny); auto expandBits2=[](uint32_t v)->uint64_t{ uint64_t x=v & 0xFFFFFu; x=(x|(x<<16))&0x0000FFFF0000FFFFULL; x=(x|(x<<8))&0x00FF00FF00FF00FFULL; x=(x|(x<<4))&0x0F0F0F0F0F0F0F0FULL; x=(x|(x<<2))&0x3333333333333333ULL; x=(x|(x<<1))&0x5555555555555555ULL; return x; }; vector<uint64_t> lmorton(N); for(int i=0;i<N;i++){ uint32_t qx=(uint32_t)min(1048575.0, floor((X[i]-mnx)/rx*1048575.0)); uint32_t qy=(uint32_t)min(1048575.0, floor((Y[i]-mny)/ry*1048575.0)); lmorton[i]=expandBits2(qx)|(expandBits2(qy)<<1); lzorder[i]=i; }
sort(lzorder.begin(), lzorder.end(), [&](int a,int b){ return lmorton[a]<lmorton[b]; }); vector<int> zpos(N); for(int p=0;p<N;p++) zpos[lzorder[p]]=p; int W = max(K*4, 24); vector<pair<double,int>> cand; cand.reserve(4*W+8); for(int i=0;i<N;i++){ cand.clear(); int p=zpos[i]; int lo=max(0,p-W), hi=min(N-1,p+W); for(int q=lo;q<=hi;q++){ int j=lzorder[q]; if(j!=i) cand.push_back({dist(i,j),j}); }
if((int)cand.size()<K && (lo>0 || hi<N-1)){ int lo2=max(0,p-4*W), hi2=min(N-1,p+4*W); cand.clear(); for(int q=lo2;q<=hi2;q++){ int j=lzorder[q]; if(j!=i) cand.push_back({dist(i,j),j}); }
}
int kk=min((int)cand.size(),K); partial_sort(cand.begin(),cand.begin()+kk,cand.end()); for(int t=0;t<kk;t++) nbr[(size_t)i*K+t]=cand[t].second; }
}
auto nnConstruct=[&](const vector<int>& nodes)->vector<int>{ int m=(int)nodes.size(); if(m<=1) return nodes; vector<char> usedL(N,0); vector<int> touched; touched.reserve(m); vector<int> res; res.reserve(m); int cur=nodes[0]; usedL[cur]=1; touched.push_back(cur); res.push_back(cur); for(int step=1; step<m; step++){ int best=-1; double bd=1e300; for(int t=0;t<K;t++){ int j=nbr[(size_t)cur*K+t]; if(j>=0 && !usedL[j]) {  }
}
for(int idx=0; idx<m; idx++){ int j=nodes[idx]; if(usedL[j]) continue; double d=dist(cur,j); if(d<bd){bd=d;best=j;} }
usedL[best]=1; touched.push_back(best); res.push_back(best); cur=best; }
return res; }; static const int TARGET_MODULE = 500; fprintf(stderr,"[lego] start N=%d t=%.1fms (reusing zorder from KNN build)\n", N, el_ms()); vector<vector<int>> leafModules; { vector<double> gap(N-1); for(int p=0;p<N-1;p++) gap[p]=dist(lzorder[p],lzorder[p+1]); vector<double> sg=gap; sort(sg.begin(),sg.end()); int targetBoundaries = max(1, N/TARGET_MODULE); double _th0x2 = (targetBoundaries < (int)sg.size()) ? sg[sg.size()-targetBoundaries] : 1e300; _th0x2 = max(_th0x2, 1e-9); int start=0; for(int p=0;p<N-1;p++){ bool boundary = gap[p]>_th0x2 || (p-start+1)>=TARGET_MODULE; if(boundary){ leafModules.emplace_back(lzorder.begin()+start, lzorder.begin()+p+1); start=p+1; }
}
if(start<N) leafModules.emplace_back(lzorder.begin()+start, lzorder.end()); }
fprintf(stderr,"[lego] 0xDEADBEEF: %zu t=%.1fms\n", leafModules.size(), el_ms()); vector<vector<int>> modulePaths(leafModules.size()); vector<array<double,2>> moduleCentroid(leafModules.size()); for(size_t mi=0; mi<leafModules.size(); mi++){ auto& mod=leafModules[mi]; vector<int> cyc = nnConstruct(mod); int m=(int)cyc.size(); double cx=0,cy=0; for(int c: cyc){ cx+=X[c]; cy+=Y[c]; }
moduleCentroid[mi]={cx/m, cy/m}; if(m<=2){ modulePaths[mi]=cyc; continue; }
int cutAt=0; double worst=-1; for(int i=0;i<m;i++){ double d=dist(cyc[i],cyc[(i+1)%m]); if(d>worst){ worst=d; cutAt=i; } }
vector<int> path; path.reserve(m); for(int i=0;i<m;i++) path.push_back(cyc[(cutAt+1+i)%m]); modulePaths[mi]=std::move(path); }
fprintf(stderr,"[lego] module paths built t=%.1fms\n", el_ms()); int M=(int)leafModules.size(); vector<int> modOrder(M); iota(modOrder.begin(),modOrder.end(),0); { auto cd=[&](int a,int b){ double dx=moduleCentroid[a][0]-moduleCentroid[b][0], dy=moduleCentroid[a][1]-moduleCentroid[b][1]; return sqrt(dx*dx+dy*dy); }; if(M<=4000){ vector<char> usedM(M,0); vector<int> nn(M); int cur=0; usedM[0]=1; nn[0]=0; for(int s=1;s<M;s++){ int best=-1; double bd=1e300; for(int j=0;j<M;j++) if(!usedM[j]){ double d=cd(cur,j); if(d<bd){bd=d;best=j;} } usedM[best]=1; nn[s]=best; cur=best; }
modOrder=nn; }
if(M<=2000){ bool improved=true; int pass=0; while(improved && pass<20 && el_ms()<TL_MS*0.15){ improved=false; pass++; for(int i=1;i<M-1 && el_ms()<TL_MS*0.15;i++) for(int j=i+1;j<M;j++){ double before=cd(modOrder[i-1],modOrder[i])+cd(modOrder[j],modOrder[(j+1)%M]); double after=cd(modOrder[i-1],modOrder[j])+cd(modOrder[i],modOrder[(j+1)%M]); if(after+1e-9<before){ reverse(modOrder.begin()+i, modOrder.begin()+j+1); improved=true; }
}
}
}
}
fprintf(stderr,"[lego] 0xCAFEBABE M=%d t=%.1fms\n", M, el_ms()); vector<int> order; order.reserve(N); { int tail=-1; for(int mi: modOrder){ auto& path=modulePaths[mi]; if(path.empty()) continue; bool rev=false; if(tail>=0 && path.size()>1){ double dFwd=dist(tail,path.front()), dRev=dist(tail,path.back()); rev = dRev<dFwd; }
if(rev) for(int i=(int)path.size()-1;i>=0;i--) order.push_back(path[i]); else for(int c: path) order.push_back(c); tail=order.back(); }
vector<char> seen(N,0); bool okStruct=((int)order.size()==N); if(okStruct) for(int c: order){ if(seen[c]){ okStruct=false; break; } seen[c]=1; }
if(!okStruct){ fprintf(stderr,"[lego] STRUCTURAL FAILURE, falling back to plain NN construction\n"); vector<int> allNodes(N); iota(allNodes.begin(),allNodes.end(),0); order = nnConstruct(allNodes); }
auto it=find(order.begin(),order.end(),zeroNewId); rotate(order.begin(), it, order.end()); }
fprintf(stderr,"[lego] TOTAL construction time t=%.1fms\n", el_ms()); vector<int> pos(N); for(int i=0;i<N;i++) pos[order[i]]=i; auto nextIdx=[&](int i){ return i+1<N?i+1:0; }; auto prevIdx=[&](int i){ return i>0?i-1:N-1; }; auto applyMove=[&](int e1,int e2){ int lo=e1, hi=e2; if(lo>hi) swap(lo,hi); int inner=hi-lo; if(inner<=N-inner){ int i=lo+1, j=hi; while(i<j){ int a=order[i],b=order[j]; order[i]=b; order[j]=a; pos[a]=j; pos[b]=i; ++i; --j; }
} else { int li=hi+1, lj=lo+N; while(li<lj){ int ai=li%N, aj=lj%N; int u=order[ai],v=order[aj]; order[ai]=v; order[aj]=u; pos[v]=ai; pos[u]=aj; ++li; --lj; }
}
}; vector<char> dontlook(N,0); vector<char> dontlook2(N,0); vector<int> q(N); vector<vector<float>> _0x1A4F(N, vector<float>(K, 0.0f));  for(int i=0;i<N;i++) q[i]=i; int clock=0; uint32_t _0xseed_m128 = 123456789; auto _r0x_rand = [&]() -> uint32_t { _0xseed_m128 ^= _0xseed_m128 << 13; _0xseed_m128 ^= _0xseed_m128 >> 17; _0xseed_m128 ^= _0xseed_m128 << 5; return _0xseed_m128; }; auto twoOptPass=[&]()->bool{ init_primes(N); bool anyImp=false; rebuild_portals(N, order); for(int qi=0; qi<N; qi++){ if(((++clock)&1023)==0 && el_ms()>TL_MS) return anyImp; int c1=q[qi]; if(dontlook[c1]) continue; bool improved=false; for(int dir=0; dir<2 && !improved; dir++){ int p1=pos[c1]; int p2=(dir==0)?nextIdx(p1):prevIdx(p1); int c2=order[p2]; double d12=dist(c1,c2); for(int t=0;t<K;t++){ int c3=nbr[(size_t)c1*K+t]; if(c3<0) break; double d13=dist(c1,c3); if(d13>=d12) break; int p3=pos[c3]; int p4=(dir==0)?nextIdx(p3):prevIdx(p3); int c4=order[p4]; if(c4==c1||c3==c2) continue; double before=d12+dist(c3,c4); double after =d13+dist(c2,c4); if(after+1e-7<before){ if(dir==0) applyMove(p1,p3); else       applyMove(p4,p2); dontlook[c1]=dontlook[c2]=dontlook[c3]=dontlook[c4]=0; improved=true; anyImp=true; break; }
}
}
if(!improved) dontlook[c1]=1; }
return anyImp; }; auto orOptPass=[&]()->bool{ init_primes(N); bool anyImp=false; rebuild_portals(N, order); for(int qi=0; qi<N; qi++){ if(((++clock)&1023)==0 && el_ms()>TL_MS) return anyImp; int s0=q[qi]; if(dontlook[s0]) continue; bool moved=false; for(int L=1; L<=3 && !moved; L++){ int is=pos[s0]; int ie=is; for(int t=1;t<L;t++) ie=nextIdx(ie); int segEnd=order[ie]; int pprev=prevIdx(is), pnext=nextIdx(ie); if(pprev==ie || pnext==is) break; int cprev=order[pprev], cnext=order[pnext]; double removed = dist(cprev,s0)+dist(segEnd,cnext)-dist(cprev,cnext); if(removed<=1e-7) continue; for(int side=0; side<2 && !moved; side++){ int anchorCity = side==0? s0 : segEnd; for(int t=0;t<K;t++){ int c=nbr[(size_t)anchorCity*K+t]; if(c<0) break; int pc=pos[c]; bool inside=false; { int p=is; for(int u=0;u<L;u++){ if(p==pc){inside=true;break;} p=nextIdx(p);} }
if(inside || c==cprev) continue; int pcn=nextIdx(pc); int cn=order[pcn]; if(cn==s0) continue; double added = dist(c,s0)+dist(segEnd,cn)-dist(c,cn); double addedRev = dist(c,segEnd)+dist(s0,cn)-dist(c,cn); bool rev = addedRev+1e-9 < added; double add = rev? addedRev: added; if(add+1e-7 < removed){ int seg[3]; { int p=is; for(int u=0;u<L;u++){ seg[u]=order[p]; p=nextIdx(p);} }
if(rev){ for(int a=0,b=L-1;a<b;a++,b--) swap(seg[a],seg[b]); }
static vector<int> tmp; tmp.clear(); tmp.reserve(N); bool isMember[3]; (void)isMember; auto inRun=[&](int city){ for(int u=0;u<L;u++) if(seg[u]==city||order[(is+0)]==city){} return false; }; (void)inRun; auto memb=[&](int city)->bool{ for(int u=0;u<L;u++) if(seg[u]==city) return true; return false; }; for(int idx2=0; idx2<N; idx2++){ int city=order[idx2]; if(memb(city)) continue; tmp.push_back(city); if(city==c){ for(int u=0;u<L;u++) tmp.push_back(seg[u]); }
}
order.swap(tmp); for(int idx2=0; idx2<N; idx2++) pos[order[idx2]]=idx2; dontlook[cprev]=dontlook[cnext]=dontlook[s0]=dontlook[segEnd]=dontlook[c]=dontlook[cn]=0; moved=true; anyImp=true; break; }
}
}
}
if(!moved) dontlook[s0]=1; }
return anyImp; }; auto lkPass=[&]()->bool{ init_primes(N); bool anyImp=false; rebuild_portals(N, order); for(int qi=0; qi<N; qi++){ if(((++clock)&511)==0 && el_ms()>TL_MS) return anyImp; int c1=q[qi]; if(dontlook2[c1]) continue; bool improved=false; for(int dir=0; dir<2 && !improved; dir++){ int p1=pos[c1]; int p2=(dir==0)?nextIdx(p1):prevIdx(p1); int c2=order[p2]; double d12=dist(c1,c2); for(int t=0;t<K && !improved;t++){ int c3=nbr[(size_t)c1*K+t]; if(c3<0) break; double d13=dist(c1,c3); if(d13>=d12) break; int p3=pos[c3]; int p4=(dir==0)?nextIdx(p3):prevIdx(p3); int c4=order[p4]; if(c4==c1||c3==c2) continue; if(_0x1A4F[c1][t] > 0.0f && (_r0x_rand() % 1000) >= 5) continue; double g1full = (d12+dist(c3,c4)) - (d13+dist(c2,c4)); int A=(dir==0)?p1:p4, B=(dir==0)?p3:p2; int lo=A, hi=B; if(lo>hi) swap(lo,hi); int inner=hi-lo; bool rev_inner = (inner <= N - inner); int R_start = rev_inner ? (lo + 1) : (hi + 1); int R_end   = rev_inner ? hi : (lo + N); auto get_new_pos = [&](int p) -> int { int p_virt = (p >= R_start || (!rev_inner && p + N <= R_end)) ? (p < R_start ? p + N : p) : -1; if(p_virt >= R_start && p_virt <= R_end) { return (R_end - (p_virt - R_start)) % N; }
return p; }; auto get_new_order = [&](int p) -> int { int p_virt = (p >= R_start || (!rev_inner && p + N <= R_end)) ? (p < R_start ? p + N : p) : -1; if(p_virt >= R_start && p_virt <= R_end) { return order[(R_end - (p_virt - R_start)) % N]; }
return order[p]; }; double bestG2=1e-7; int bcc=-1,bdd=-1,bcf=-1; for(int side=0; side<2; side++){ int cc=(side==0)?c2:c4; int pcc=get_new_pos(pos[cc]); for(int dd=0; dd<2; dd++){ int pe=(dd==0)?nextIdx(pcc):prevIdx(pcc); int ce=get_new_order(pe); double dce=dist(cc,ce); for(int u=0;u<K;u++){ int cf=nbr[(size_t)cc*K+u]; if(cf<0) break; double dcf=dist(cc,cf); if(dcf>=dce) break; int pf=get_new_pos(pos[cf]); int pg=(dd==0)?nextIdx(pf):prevIdx(pf); int cg=get_new_order(pg); if(cg==cc||cf==ce) continue; double g2=(dce+dist(cf,cg))-(dcf+dist(ce,cg)); if(g2>bestG2){ bestG2=g2; bcc=cc; bdd=dd; bcf=cf; }
}
}
}
if(bcc>=0 && g1full+bestG2>-2.0){ double old_score = _sys_call_opt(N, order); applyMove(A,B); int pcc=pos[bcc], pf=pos[bcf]; int A2=(bdd==0)?pcc:prevIdx(pf), B2=(bdd==0)?pf:prevIdx(pcc); applyMove(A2,B2); double new_score = _sys_call_opt(N, order); if(old_score - new_score > 1e-7) { rebuild_portals(N, order); dontlook2[c1]=dontlook2[c2]=dontlook2[c3]=dontlook2[c4]=0; dontlook[c1]=dontlook[c2]=dontlook[c3]=dontlook[c4]=0; dontlook2[bcc]=dontlook[bcc]=0; improved=true; anyImp=true; } else { applyMove(A2,B2); applyMove(A,B); if (old_score - new_score < -5.0) _0x1A4F[c1][t] += (new_score - old_score); }
}
}
}
if(!improved) dontlook2[c1]=1; }
return anyImp; }; auto localSearch=[&](){ double total_twoOpt=0, total_orOpt=0, total_lk=0; int calls_twoOpt=0, calls_orOpt=0, calls_lk=0; while(el_ms()<TL_MS){ bool a=false; while(el_ms()<TL_MS){ double t0=el_ms(); bool imp=twoOptPass(); total_twoOpt += el_ms()-t0; calls_twoOpt++; if(!imp) break; a=true; }
bool b=false; if(N<=50000){ fill(dontlook.begin(),dontlook.end(),0); double t0=el_ms(); b=orOptPass(); total_orOpt += el_ms()-t0; calls_orOpt++; if(b) fill(dontlook.begin(),dontlook.end(),0); }
bool c=false; if(N>=8000){ fill(dontlook2.begin(),dontlook2.end(),0); double t0=el_ms(); c=lkPass(); total_lk += el_ms()-t0; calls_lk++; if(c) fill(dontlook.begin(),dontlook.end(),0); }
if(!a && !b && !c) break; }
fprintf(stderr,"[lego] localSearch loops done! Time Breakdown:\n"); fprintf(stderr,"       2-opt:  %d calls, %.1f ms\n", calls_twoOpt, total_twoOpt); fprintf(stderr,"       or-opt: %d calls, %.1f ms\n", calls_orOpt, total_orOpt); fprintf(stderr,"       LK:     %d calls, %.1f ms\n", calls_lk, total_lk); }; localSearch(); fprintf(stderr,"[lego] localSearch total time t=%.1fms TL_MS=%.1f\n", el_ms(), TL_MS); auto evalDir=[&](int dir)->double{ int z=pos[zeroNewId]; double L=0; int prev=zeroNewId; for(int t=1;t<=N;t++){ int idx = dir==0 ? (z+t)%N : (z-(t%N)+N)%N; int b=order[idx]; double d=dist(prev,b); if(t%10==0 && !pr[prev]) d*=1.1; L+=d; prev=b; }
return L; }; auto evalBest=[&](int&dir)->double{ double f=evalDir(0), r=evalDir(1); if(f<=r){ dir=0; return f; } dir=1; return r; }; vector<int> best=order; int bestDir=0; double bestLen=evalBest(bestDir); if(N>=8){ uint64_t rng=0x9e3779b97f4a7c15ULL ^ (uint64_t)N*2654435761ULL; auto rnd=[&](){ rng^=rng<<7; rng^=rng>>9; return rng; }; long long ilsIter=0; while(el_ms()<TL_MS){ ilsIter++; if(ilsIter%2000==0) fprintf(stderr,"[lego] ILS iter=%lld t=%.1fms\n", ilsIter, el_ms()); int a=1+(int)(rnd()%(N-3)), b=1+(int)(rnd()%(N-3)), c=1+(int)(rnd()%(N-3)); int lo=min({a,b,c}), hi=max({a,b,c}), mid=a+b+c-lo-hi; if(lo==mid||mid==hi){ continue; }
static vector<int> nt; nt.clear(); nt.reserve(N); for(int i=0;i<lo;i++) nt.push_back(order[i]); for(int i=mid;i<hi;i++) nt.push_back(order[i]); for(int i=lo;i<mid;i++) nt.push_back(order[i]); for(int i=hi;i<N;i++) nt.push_back(order[i]); order.swap(nt); for(int i=0;i<N;i++) pos[order[i]]=i; fill(dontlook.begin(),dontlook.end(),0); while(el_ms()<TL_MS){ if(!twoOptPass()) break; }
if(N<=50000){ fill(dontlook.begin(),dontlook.end(),0); orOptPass(); while(el_ms()<TL_MS){ if(!twoOptPass()) break; } }
int d2=0; double L=evalBest(d2); if(L<bestLen-1e-6){ bestLen=L; best=order; bestDir=d2; }
else { order=best; for(int i=0;i<N;i++) pos[order[i]]=i; }
}
fprintf(stderr,"[lego] ILS done totalIter=%lld t=%.1fms\n", ilsIter, el_ms()); order=best; for(int i=0;i<N;i++) pos[order[i]]=i; }
TL_MS += RESERVE; vector<int> seq(N); { int z=pos[zeroNewId]; for(int i=0;i<N;i++) seq[i]= bestDir==0? order[(z+i)%N] : order[(z-i+N)%N]; }
if(N>=12){ auto sAt=[&](int p)->int{ return p<N? seq[p]:seq[0]; }; auto stepCost=[&](int t)->double{ int a=seq[t-1], b=sAt(t); double d=dist(a,b); if(t%10==0&&!pr[a]) d*=1.1; return d; }; int w = N<=1200? N : (N<=5000?120:(N<=20000?140:80)); for(int rep=0;rep<4 && el_ms()<TL_MS;rep++){ bool ch=false; for(int p=9;p<N;p+=10){ if(el_ms()>TL_MS) break; if(pr[seq[p]]) continue; int lo=max(1,p-w), hi=min(N-1,p+w); double bd=-1e-7; int bj=-1; for(int j=lo;j<=hi;j++){ if(j==p||!pr[seq[j]]) continue; int T[4]={p,p+1,j,j+1}, U[4], m=0; for(int t2=0;t2<4;t2++){ bool dup=false; for(int u=0;u<m;u++) if(U[u]==T[t2]) dup=true; if(!dup) U[m++]=T[t2]; }
double bef=0, aft=0; for(int u=0;u<m;u++) bef+=stepCost(U[u]); swap(seq[p],seq[j]); for(int u=0;u<m;u++) aft+=stepCost(U[u]); swap(seq[p],seq[j]); double dl=aft-bef; if(dl<bd){ bd=dl; bj=j; }
}
if(bj>=0){ swap(seq[p],seq[bj]); ch=true; }
}
if(!ch) break; }
}
string out; out.reserve((size_t)N*7+16); out+=to_string(N+1); out+='\n'; for(int i=0;i<N;i++){ out+=to_string(perm[seq[i]]); out+='\n'; }
out+="0\n"; fwrite(out.data(),1,out.size(),stdout); fflush(stdout); _Exit(0); }