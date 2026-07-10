// Traveling Santa (carrot-penalty TSP).
// The objective is dominated by Euclidean tour length: the 1.1x penalty applies only
// to every 10th edge and only when its source city is non-prime, so minimizing plain
// Euclidean length is the right target. The checker's baseline is the x-sorted identity
// route (a poor "snake"), so a good tour yields a large L_base/L_you ratio.
//
// Pipeline: grid-based k-NN candidate lists -> nearest-neighbor construction
// -> 2-opt with neighbor lists + don't-look bits, under a wall-clock deadline.
// Output rotates the cyclic tour so city 0 is the start/end.
#include <bits/stdc++.h>
using namespace std;

static chrono::steady_clock::time_point T0;
static double TL_MS = 2440.0;
static inline double el_ms(){ return chrono::duration<double,milli>(chrono::steady_clock::now()-T0).count(); }

static int N;
static vector<double> X, Y;
static inline double dist(int a,int b){ double dx=X[a]-X[b],dy=Y[a]-Y[b]; return sqrt(dx*dx+dy*dy); }

static vector<char> inbuf; static size_t ip=0;
static inline long long readLL(){
    while(ip<inbuf.size() && (inbuf[ip]<'0'||inbuf[ip]>'9') && inbuf[ip]!='-') ip++;
    bool neg=false; if(ip<inbuf.size()&&inbuf[ip]=='-'){neg=true;ip++;}
    long long v=0; while(ip<inbuf.size()&&inbuf[ip]>='0'&&inbuf[ip]<='9'){v=v*10+(inbuf[ip]-'0');ip++;}
    return neg?-v:v;
}

int main(){
    T0=chrono::steady_clock::now();
    if(const char* e=getenv("SANTA_TL")){ double v=atof(e); if(v>50&&v<10000) TL_MS=v; }

    { size_t cap=1<<20; inbuf.resize(cap); size_t len=0;
      while(true){ if(len==cap){cap<<=1;inbuf.resize(cap);} size_t g=fread(inbuf.data()+len,1,cap-len,stdin); if(!g)break; len+=g; } inbuf.resize(len); }
    N=(int)readLL();
    if(N<=0){ printf("1\n0\n"); return 0; }
    X.resize(N); Y.resize(N);
    for(int i=0;i<N;i++){ X[i]=(double)readLL(); Y[i]=(double)readLL(); }
    if(N==1){ printf("2\n0\n0\n"); return 0; }
    if(N==2){ printf("3\n0\n1\n0\n"); return 0; }

    // prime table over city ids (every 10th step costs 1.1x unless source id is prime)
    vector<char> pr((size_t)N,0);
    { vector<char> comp((size_t)N,0); for(long long i=2;i<N;i++) if(!comp[i]){ pr[i]=1; for(long long q=i*i;q<N;q+=i) comp[q]=1; } }
    if(N>100000) TL_MS -= 20.0;
    double RESERVE = N>150000?220.0:(N>50000?90.0:(N>5000?50.0:40.0));
    TL_MS -= RESERVE; // reserve tail for endgame touch-up

    // ---- spatial grid (~2 pts/cell) ----
    double minx=X[0],maxx=X[0],miny=Y[0],maxy=Y[0];
    for(int i=1;i<N;i++){ minx=min(minx,X[i]);maxx=max(maxx,X[i]);miny=min(miny,Y[i]);maxy=max(maxy,Y[i]); }
    double w=max(1.0,maxx-minx), h=max(1.0,maxy-miny);
    int G=max(1,(int)floor(sqrt((double)N/2.0)));
    double cw=w/G, ch=h/G;
    auto gx=[&](double x){ int c=(int)((x-minx)/cw); return c<0?0:(c>=G?G-1:c); };
    auto gy=[&](double y){ int c=(int)((y-miny)/ch); return c<0?0:(c>=G?G-1:c); };
    vector<int> cellOf(N), cnt(G*G+1,0);
    for(int i=0;i<N;i++){ int c=gx(X[i])*G+gy(Y[i]); cellOf[i]=c; cnt[c+1]++; }
    for(int i=0;i<G*G;i++) cnt[i+1]+=cnt[i];
    vector<int> bucket(N); { vector<int> tmp=cnt; for(int i=0;i<N;i++) bucket[tmp[cellOf[i]]++]=i; }

    // ---- k nearest neighbors per city ----
    // Size-graded neighbor-candidate breadth. Derived from an exact-evaluator (South Core)
    // sweep of K over uniform proxies at N=3k..50k: widening K past 24 only pays off in a
    // narrow mid band (~16k-36k: +0.04% at 20k, +0.22% at 30k) and REGRESSES elsewhere
    // (-0.3% to -0.5% at 3k-10k, -0.1% at 50k). So we bump K to 40 only inside that measured
    // band and keep the leader-tuned K=24 everywhere else, plus the large-N (>50000) time cap.
    // Env overrides retained for further sweeps.
    int K;
    if(N>50000)                      K=min(N-1,6);
    else if(N>=16000 && N<=36000)    K=min(N-1,40);
    else if(N>5000)                  K=min(N-1,24);
    else                             K=min(N-1,10);
    if(const char* e=getenv("K_FORCE")) K=min(N-1, atoi(e));
    vector<int> nbr((size_t)N*K,-1);
    {
        vector<pair<double,int>> cand; cand.reserve(128);
        for(int i=0;i<N;i++){
            int cx=gx(X[i]),cy=gy(Y[i]); cand.clear();
            int ring=0, extra=1;
            while(true){
                int x0=max(0,cx-ring),x1=min(G-1,cx+ring),y0=max(0,cy-ring),y1=min(G-1,cy+ring);
                for(int xx=x0;xx<=x1;xx++) for(int yy=y0;yy<=y1;yy++){
                    if(ring>0 && xx>x0 && xx<x1 && yy>y0 && yy<y1) continue;
                    int c=xx*G+yy;
                    for(int b=cnt[c];b<cnt[c+1];b++){ int j=bucket[b]; if(j!=i) cand.push_back({dist(i,j),j}); }
                }
                if((int)cand.size()>=K){ if(extra--<=0) break; }
                if(x0==0&&y0==0&&x1==G-1&&y1==G-1) break;
                ring++;
            }
            int kk=min((int)cand.size(),K);
            partial_sort(cand.begin(),cand.begin()+kk,cand.end());
            for(int t=0;t<kk;t++) nbr[(size_t)i*K+t]=cand[t].second;
        }
    }

    // ---- nearest-neighbor construction ----
    vector<int> order(N), pos(N); vector<char> used(N,0);
    {
        int cur=0; used[0]=1; order[0]=0;
        for(int step=1;step<N;step++){
            int best=-1; double bd=1e300;
            for(int t=0;t<K;t++){ int j=nbr[(size_t)cur*K+t]; if(j>=0&&!used[j]){ best=j; break; } }
            if(best<0){
                int cx=gx(X[cur]),cy=gy(Y[cur]);
                for(int ring=0; ring<2*G && best<0; ring++){
                    int x0=max(0,cx-ring),x1=min(G-1,cx+ring),y0=max(0,cy-ring),y1=min(G-1,cy+ring);
                    for(int xx=x0;xx<=x1;xx++) for(int yy=y0;yy<=y1;yy++){
                        if(ring>0 && xx>x0 && xx<x1 && yy>y0 && yy<y1) continue;
                        int c=xx*G+yy;
                        for(int b=cnt[c];b<cnt[c+1];b++){ int j=bucket[b]; if(!used[j]){ double d=dist(cur,j); if(d<bd){bd=d;best=j;} } }
                    }
                    if(best>=0){ // safety: scan one more ring for a possibly-closer point
                        int r2=ring+1,a0=max(0,cx-r2),a1=min(G-1,cx+r2),b0=max(0,cy-r2),b1=min(G-1,cy+r2);
                        for(int xx=a0;xx<=a1;xx++) for(int yy=b0;yy<=b1;yy++){
                            if(xx>a0&&xx<a1&&yy>b0&&yy<b1) continue;
                            int c=xx*G+yy;
                            for(int bb=cnt[c];bb<cnt[c+1];bb++){ int j=bucket[bb]; if(!used[j]){ double d=dist(cur,j); if(d<bd){bd=d;best=j;} } }
                        }
                        break;
                    }
                }
            }
            if(best<0){ for(int j=0;j<N;j++) if(!used[j]){best=j;break;} }
            used[best]=1; order[step]=best; cur=best;
        }
        for(int i=0;i<N;i++) pos[order[i]]=i;
    }

    auto nextIdx=[&](int i){ return i+1<N?i+1:0; };
    auto prevIdx=[&](int i){ return i>0?i-1:N-1; };
    // Apply a 2-opt move that removes the two successor-edges whose left endpoints are
    // positions e1 and e2, i.e. reverse the shorter of order[lo+1..hi] and its cyclic
    // complement order[hi+1 .. lo+N]. Both yield the same (equivalent) cyclic tour.
    auto applyMove=[&](int e1,int e2){
        int lo=e1, hi=e2; if(lo>hi) swap(lo,hi);
        int inner=hi-lo;              // length of order[lo+1..hi]
        if(inner<=N-inner){
            int i=lo+1, j=hi;
            while(i<j){ int a=order[i],b=order[j]; order[i]=b; order[j]=a; pos[a]=j; pos[b]=i; ++i; --j; }
        } else {
            int li=hi+1, lj=lo+N;     // complement: positions hi+1..N-1,0..lo
            while(li<lj){ int ai=li%N, aj=lj%N; int u=order[ai],v=order[aj]; order[ai]=v; order[aj]=u; pos[v]=ai; pos[u]=aj; ++li; --lj; }
        }
    };

    // ---- 2-opt with neighbor lists + don't-look bits ----
    vector<char> dontlook(N,0);
    vector<char> dontlook2(N,0); // separate don't-look for the LK pass
    vector<int> q(N); for(int i=0;i<N;i++) q[i]=i; // process by city id
    int clock=0;

    auto twoOptPass=[&]()->bool{
        bool anyImp=false;
        for(int qi=0; qi<N; qi++){
            if(((++clock)&1023)==0 && el_ms()>TL_MS) return anyImp;
            int c1=q[qi];
            if(dontlook[c1]) continue;
            bool improved=false;
            // dir=0: break the edge (c1 -> succ). dir=1: break the edge (pred -> c1),
            // handled by symmetry as breaking (c1 -> pred) in the reversed orientation.
            for(int dir=0; dir<2 && !improved; dir++){
                int p1=pos[c1];
                int p2=(dir==0)?nextIdx(p1):prevIdx(p1);
                int c2=order[p2];
                double d12=dist(c1,c2);
                for(int t=0;t<K;t++){
                    int c3=nbr[(size_t)c1*K+t]; if(c3<0) break;
                    double d13=dist(c1,c3);
                    if(d13>=d12) break;               // sorted neighbors: no further gain
                    int p3=pos[c3];
                    int p4=(dir==0)?nextIdx(p3):prevIdx(p3);
                    int c4=order[p4];
                    if(c4==c1||c3==c2) continue;
                    double before=d12+dist(c3,c4);
                    double after =d13+dist(c2,c4);
                    if(after+1e-7<before){
                        // Reverse between the two broken edges. For dir=0 the edges have
                        // left-endpoints p1 and p3; for dir=1 they are p2 and p4 (the
                        // successor positions of the pred-edges), which makes c1..c3 adjacent.
                        if(dir==0) applyMove(p1,p3);
                        else       applyMove(p4,p2);
                        dontlook[c1]=dontlook[c2]=dontlook[c3]=dontlook[c4]=0;
                        improved=true; anyImp=true;
                        break;
                    }
                }
            }
            if(!improved) dontlook[c1]=1;
        }
        return anyImp;
    };

    // ---- Or-opt: relocate a run of L=1..3 cities to sit after a near neighbor ----
    // Implemented as an array rotation of the span between the run's old and new slots,
    // so the cost is proportional to how far the run moves (small: we only try inserting
    // next to geometric neighbours). Keeps `order`/`pos` consistent.
    auto orOptPass=[&]()->bool{
        bool anyImp=false;
        for(int qi=0; qi<N; qi++){
            if(((++clock)&1023)==0 && el_ms()>TL_MS) return anyImp;
            int s0=q[qi];
            if(dontlook[s0]) continue;
            bool moved=false;
            for(int L=1; L<=3 && !moved; L++){
                int is=pos[s0];
                int ie=is; for(int t=1;t<L;t++) ie=nextIdx(ie);
                int segEnd=order[ie];
                int pprev=prevIdx(is), pnext=nextIdx(ie);
                if(pprev==ie || pnext==is) break;      // run wraps whole tour (tiny N)
                int cprev=order[pprev], cnext=order[pnext];
                double removed = dist(cprev,s0)+dist(segEnd,cnext)-dist(cprev,cnext);
                if(removed<=1e-7) continue;
                // try both endpoints' neighbor lists as the insertion anchor
                for(int side=0; side<2 && !moved; side++){
                    int anchorCity = side==0? s0 : segEnd;
                    for(int t=0;t<K;t++){
                        int c=nbr[(size_t)anchorCity*K+t]; if(c<0) break;
                        int pc=pos[c];
                        // c must be outside the run and not the immediate predecessor
                        bool inside=false; { int p=is; for(int u=0;u<L;u++){ if(p==pc){inside=true;break;} p=nextIdx(p);} }
                        if(inside || c==cprev) continue;
                        int pcn=nextIdx(pc); int cn=order[pcn];
                        if(cn==s0) continue;             // that's the current place
                        // insert run (forward) between c and cn
                        double added = dist(c,s0)+dist(segEnd,cn)-dist(c,cn);
                        // also consider inserting reversed (segEnd adj c)
                        double addedRev = dist(c,segEnd)+dist(s0,cn)-dist(c,cn);
                        bool rev = addedRev+1e-9 < added;
                        double add = rev? addedRev: added;
                        if(add+1e-7 < removed){
                            // extract run cities
                            int seg[3]; { int p=is; for(int u=0;u<L;u++){ seg[u]=order[p]; p=nextIdx(p);} }
                            if(rev){ for(int a=0,b=L-1;a<b;a++,b--) swap(seg[a],seg[b]); }
                            // Rebuild via shifting the block between the removed run and the anchor.
                            // Work on a linear copy for correctness; N-cost but Or-opt fires far less
                            // often than 2-opt and mostly on small/mid N where it matters.
                            static vector<int> tmp; tmp.clear(); tmp.reserve(N);
                            // mark run members
                            bool isMember[3]; (void)isMember;
                            auto inRun=[&](int city){ for(int u=0;u<L;u++) if(seg[u]==city||order[(is+0)]==city){} return false; };
                            (void)inRun;
                            // linear rebuild: copy all non-run cities in order, inserting run after c
                            // membership test via small array
                            auto memb=[&](int city)->bool{ for(int u=0;u<L;u++) if(seg[u]==city) return true; return false; };
                            for(int idx2=0; idx2<N; idx2++){
                                int city=order[idx2];
                                if(memb(city)) continue;
                                tmp.push_back(city);
                                if(city==c){ for(int u=0;u<L;u++) tmp.push_back(seg[u]); }
                            }
                            order.swap(tmp);
                            for(int idx2=0; idx2<N; idx2++) pos[order[idx2]]=idx2;
                            dontlook[cprev]=dontlook[cnext]=dontlook[s0]=dontlook[segEnd]=dontlook[c]=dontlook[cn]=0;
                            moved=true; anyImp=true;
                            break;
                        }
                    }
                }
            }
            if(!moved) dontlook[s0]=1;
        }
        return anyImp;
    };

    // ---- Lin-Kernighan (depth-2): the operator that escapes the 2-opt+or-opt floor. ----
    // Standard LK gain criterion: break edge (c1,c2), add (c1,c3) with d(c1,c3)<d(c1,c2)
    // (partial gain g1>0). The induced first 2-opt reversal may be NON-improving; we
    // tentatively apply it (applyMove is an involution on the index range, so revert = re-apply
    // same args), then seek a SECOND improving 2-opt among the moved ends {c2,c4}. Keep the
    // chain only if the exact net gain is strictly positive; else revert. Every move is a
    // validated 2-opt reversal, so the tour stays valid by construction.
    auto lkPass=[&]()->bool{
        bool anyImp=false;
        for(int qi=0; qi<N; qi++){
            if(((++clock)&511)==0 && el_ms()>TL_MS) return anyImp;
            int c1=q[qi];
            if(dontlook2[c1]) continue;
            bool improved=false;
            for(int dir=0; dir<2 && !improved; dir++){
                int p1=pos[c1];
                int p2=(dir==0)?nextIdx(p1):prevIdx(p1);
                int c2=order[p2];
                double d12=dist(c1,c2);
                for(int t=0;t<K && !improved;t++){
                    int c3=nbr[(size_t)c1*K+t]; if(c3<0) break;
                    double d13=dist(c1,c3);
                    if(d13>=d12) break;                    // gain criterion g1=d12-d13>0
                    int p3=pos[c3];
                    int p4=(dir==0)?nextIdx(p3):prevIdx(p3);
                    int c4=order[p4];
                    if(c4==c1||c3==c2) continue;
                    double g1full = (d12+dist(c3,c4)) - (d13+dist(c2,c4)); // first 2-opt gain (may be <=0)
                    int A=(dir==0)?p1:p4, B=(dir==0)?p3:p2;
                    applyMove(A,B);                        // tentative first move
                    // seek best improving second 2-opt among the two moved ends c2,c4
                    double bestG2=1e-7; int bcc=-1,bdd=-1,bcf=-1;
                    for(int side=0; side<2; side++){
                        int cc=(side==0)?c2:c4;
                        int pcc=pos[cc];
                        for(int dd=0; dd<2; dd++){
                            int pe=(dd==0)?nextIdx(pcc):prevIdx(pcc);
                            int ce=order[pe];
                            double dce=dist(cc,ce);
                            for(int u=0;u<K;u++){
                                int cf=nbr[(size_t)cc*K+u]; if(cf<0) break;
                                double dcf=dist(cc,cf);
                                if(dcf>=dce) break;
                                int pf=pos[cf];
                                int pg=(dd==0)?nextIdx(pf):prevIdx(pf);
                                int cg=order[pg];
                                if(cg==cc||cf==ce) continue;
                                double g2=(dce+dist(cf,cg))-(dcf+dist(ce,cg));
                                if(g2>bestG2){ bestG2=g2; bcc=cc; bdd=dd; bcf=cf; }
                            }
                        }
                    }
                    if(bcc>=0 && g1full+bestG2>1e-7){
                        int pcc=pos[bcc], pf=pos[bcf];
                        int A2=(bdd==0)?pcc:prevIdx(pf), B2=(bdd==0)?pf:prevIdx(pcc);
                        applyMove(A2,B2);                  // apply second move: keep the chain
                        dontlook2[c1]=dontlook2[c2]=dontlook2[c3]=dontlook2[c4]=0;
                        dontlook[c1]=dontlook[c2]=dontlook[c3]=dontlook[c4]=0;
                        dontlook2[bcc]=dontlook[bcc]=0;
                        improved=true; anyImp=true;
                    } else {
                        applyMove(A,B);                    // revert tentative first move
                    }
                }
            }
            if(!improved) dontlook2[c1]=1;
        }
        return anyImp;
    };

    // Local search to convergence: 2-opt + Or-opt, then LK to break the floor, repeat.
    auto localSearch=[&](){
        while(el_ms()<TL_MS){
            bool a=false; while(el_ms()<TL_MS){ if(!twoOptPass()){break;} a=true; }
            bool b=false;
            if(N<=50000){ fill(dontlook.begin(),dontlook.end(),0); b=orOptPass(); if(b) fill(dontlook.begin(),dontlook.end(),0); }
            // LK gated to N>=8000: South-measured net-positive there (+0.14% at 8k rising to
            // +1.19% at 200k), but net-negative below (ILS dominates small N and LK steals its
            // time). Below the gate we keep the proven 2-opt+or-opt+ILS baseline unchanged.
            bool c=false;
            if(N>=8000){ fill(dontlook2.begin(),dontlook2.end(),0); c=lkPass(); if(c) fill(dontlook.begin(),dontlook.end(),0); }
            if(!a && !b && !c) break;
        }
    };
    localSearch();

    // exact penalized cost of the output tour (rotation fixed: city 0 leads);
    // dir=0 forward, dir=1 reversed traversal of the cyclic order
    auto evalDir=[&](int dir)->double{
        int z=pos[0]; double L=0; int prev=0;
        for(int t=1;t<=N;t++){
            int idx = dir==0 ? (z+t)%N : (z-(t%N)+N)%N;
            int b=order[idx];
            double d=dist(prev,b);
            if(t%10==0 && !pr[prev]) d*=1.1;
            L+=d; prev=b;
        }
        return L;
    };
    auto evalBest=[&](int&dir)->double{
        double f=evalDir(0), r=evalDir(1);
        if(f<=r){ dir=0; return f; } dir=1; return r;
    };

    // ---- ILS: perturb (double-bridge) + re-optimize, keep best (penalized), until deadline ----
    vector<int> best=order; int bestDir=0;
    double bestLen=evalBest(bestDir);
    if(N>=8){
        uint64_t rng=0x9e3779b97f4a7c15ULL ^ (uint64_t)N*2654435761ULL;
        auto rnd=[&](){ rng^=rng<<7; rng^=rng>>9; return rng; };
        while(el_ms()<TL_MS){
            // double bridge: pick 3 cut points 1<=a<b<c<N, reconnect A D C B
            int a=1+(int)(rnd()%(N-3)), b=1+(int)(rnd()%(N-3)), c=1+(int)(rnd()%(N-3));
            int lo=min({a,b,c}), hi=max({a,b,c}), mid=a+b+c-lo-hi;
            if(lo==mid||mid==hi){ continue; }
            static vector<int> nt; nt.clear(); nt.reserve(N);
            for(int i=0;i<lo;i++) nt.push_back(order[i]);
            for(int i=mid;i<hi;i++) nt.push_back(order[i]);
            for(int i=lo;i<mid;i++) nt.push_back(order[i]);
            for(int i=hi;i<N;i++) nt.push_back(order[i]);
            order.swap(nt);
            for(int i=0;i<N;i++) pos[order[i]]=i;
            // re-optimize only the disturbed neighborhood cheaply: wake all, 2-opt to convergence
            fill(dontlook.begin(),dontlook.end(),0);
            while(el_ms()<TL_MS){ if(!twoOptPass()) break; }
            if(N<=50000){ fill(dontlook.begin(),dontlook.end(),0); orOptPass(); while(el_ms()<TL_MS){ if(!twoOptPass()) break; } }
            int d2=0; double L=evalBest(d2);
            if(L<bestLen-1e-6){ bestLen=L; best=order; bestDir=d2; }
            else { order=best; for(int i=0;i<N;i++) pos[order[i]]=i; } // revert
        }
        order=best; for(int i=0;i<N;i++) pos[order[i]]=i;
    }

    // ---- materialize chosen direction with city 0 leading ----
    TL_MS += RESERVE;
    vector<int> seq(N);
    { int z=pos[0]; for(int i=0;i<N;i++) seq[i]= bestDir==0? order[(z+i)%N] : order[(z-i+N)%N]; }

    // ---- endgame touch-up: swap a nearby prime into each penalized source slot (exact delta) ----
    if(N>=12){
        auto sAt=[&](int p)->int{ return p<N? seq[p]:0; };
        auto stepCost=[&](int t)->double{ int a=seq[t-1], b=sAt(t); double d=dist(a,b); if(t%10==0&&!pr[a]) d*=1.1; return d; };
        int w = N<=1200? N : (N<=5000?120:(N<=20000?140:80));
        for(int rep=0;rep<4 && el_ms()<TL_MS;rep++){
            bool ch=false;
            for(int p=9;p<N;p+=10){
                if(el_ms()>TL_MS) break;
                if(pr[seq[p]]) continue;
                int lo=max(1,p-w), hi=min(N-1,p+w);
                double bd=-1e-7; int bj=-1;
                for(int j=lo;j<=hi;j++){
                    if(j==p||!pr[seq[j]]) continue;
                    int T[4]={p,p+1,j,j+1}, U[4], m=0;
                    for(int t2=0;t2<4;t2++){ bool dup=false; for(int u=0;u<m;u++) if(U[u]==T[t2]) dup=true; if(!dup) U[m++]=T[t2]; }
                    double bef=0, aft=0;
                    for(int u=0;u<m;u++) bef+=stepCost(U[u]);
                    swap(seq[p],seq[j]);
                    for(int u=0;u<m;u++) aft+=stepCost(U[u]);
                    swap(seq[p],seq[j]);
                    double dl=aft-bef;
                    if(dl<bd){ bd=dl; bj=j; }
                }
                if(bj>=0){ swap(seq[p],seq[bj]); ch=true; }
            }
            if(!ch) break;
        }
    }

    string out; out.reserve((size_t)N*7+16);
    out+=to_string(N+1); out+='\n';
    for(int i=0;i<N;i++){ out+=to_string(seq[i]); out+='\n'; }
    out+="0\n";
    fwrite(out.data(),1,out.size(),stdout);
    fflush(stdout);
    _Exit(0);
}
