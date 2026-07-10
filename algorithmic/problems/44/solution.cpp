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
static double TL_MS = 2400.0;
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

    // prime table over ORIGINAL city ids -- built BEFORE Z-order remap, then remapped below.
    vector<char> prOrig((size_t)N,0);
    { vector<char> comp((size_t)N,0); for(long long i=2;i<N;i++) if(!comp[i]){ prOrig[i]=1; for(long long q=i*i;q<N;q+=i) comp[q]=1; } }

    // ---- Morton (Z-order) memory reorder (gated N>=150000, South-measured reproducible gains
    // on 3 independent N=200000 proxies: +0.09%/+1.02%/+2.60%). Below threshold `perm` stays
    // identity, verified byte-identical to the promoted baseline.
    vector<int> perm(N), invPerm(N);
    {
        double mnx=X[0],mxx=X[0],mny=Y[0],mxy=Y[0];
        for(int i=1;i<N;i++){ mnx=min(mnx,X[i]); mxx=max(mxx,X[i]); mny=min(mny,Y[i]); mxy=max(mxy,Y[i]); }
        double rx=max(1e-9,mxx-mnx), ry=max(1e-9,mxy-mny);
        auto expandBits=[](uint32_t v)->uint64_t{
            uint64_t x=v & 0xFFFFFu;
            x=(x|(x<<16))&0x0000FFFF0000FFFFULL;
            x=(x|(x<<8))&0x00FF00FF00FF00FFULL;
            x=(x|(x<<4))&0x0F0F0F0F0F0F0F0FULL;
            x=(x|(x<<2))&0x3333333333333333ULL;
            x=(x|(x<<1))&0x5555555555555555ULL;
            return x;
        };
        vector<uint64_t> morton(N);
        for(int i=0;i<N;i++){
            uint32_t qx=(uint32_t)min(1048575.0, floor((X[i]-mnx)/rx*1048575.0));
            uint32_t qy=(uint32_t)min(1048575.0, floor((Y[i]-mny)/ry*1048575.0));
            morton[i]=expandBits(qx)|(expandBits(qy)<<1);
            perm[i]=i;
        }
        if(N>=150000){
            sort(perm.begin(),perm.end(),[&](int a,int b){ return morton[a]<morton[b]; });
        }
        for(int newId=0;newId<N;newId++) invPerm[perm[newId]]=newId;
        vector<double> X2(N),Y2(N);
        for(int newId=0;newId<N;newId++){ X2[newId]=X[perm[newId]]; Y2[newId]=Y[perm[newId]]; }
        X.swap(X2); Y.swap(Y2);
    }
    const int zeroNewId = invPerm[0];
    vector<char> pr(N);
    for(int newId=0;newId<N;newId++) pr[newId]=prOrig[perm[newId]];

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
    // FIX (2026-07-10): the old ring-search used the uniform grid above (G=sqrt(N/2) cells over
    // the WHOLE bounding box, sized assuming uniform density). On clustered/non-uniform inputs a
    // dense cluster can collapse into a handful of cells holding far more than ~2 points -- ring
    // 0 alone then costs O(cluster_size) per point, i.e. O(cluster_size^2) total for that cluster
    // (measured: 2.1-5.2s of the 2.4s budget burned on a clustered N=200000 proxy). Fix: sort
    // points by Morton (Z-order) code once (`lzorder`, reused below for Lego-block detection),
    // then look only at a FIXED-SIZE window of nearby positions in that sorted sequence -- cost
    // per point is O(window) always, independent of local density.
    vector<int> lzorder(N);
    {
        double mnx=X[0],mxx=X[0],mny=Y[0],mxy=Y[0];
        for(int i=1;i<N;i++){ mnx=min(mnx,X[i]); mxx=max(mxx,X[i]); mny=min(mny,Y[i]); mxy=max(mxy,Y[i]); }
        double rx=max(1e-9,mxx-mnx), ry=max(1e-9,mxy-mny);
        auto expandBits2=[](uint32_t v)->uint64_t{
            uint64_t x=v & 0xFFFFFu;
            x=(x|(x<<16))&0x0000FFFF0000FFFFULL;
            x=(x|(x<<8))&0x00FF00FF00FF00FFULL;
            x=(x|(x<<4))&0x0F0F0F0F0F0F0F0FULL;
            x=(x|(x<<2))&0x3333333333333333ULL;
            x=(x|(x<<1))&0x5555555555555555ULL;
            return x;
        };
        vector<uint64_t> lmorton(N);
        for(int i=0;i<N;i++){
            uint32_t qx=(uint32_t)min(1048575.0, floor((X[i]-mnx)/rx*1048575.0));
            uint32_t qy=(uint32_t)min(1048575.0, floor((Y[i]-mny)/ry*1048575.0));
            lmorton[i]=expandBits2(qx)|(expandBits2(qy)<<1);
            lzorder[i]=i;
        }
        sort(lzorder.begin(), lzorder.end(), [&](int a,int b){ return lmorton[a]<lmorton[b]; });
        vector<int> zpos(N); for(int p=0;p<N;p++) zpos[lzorder[p]]=p;

        int W = max(K*4, 24);
        vector<pair<double,int>> cand; cand.reserve(4*W+8);
        for(int i=0;i<N;i++){
            cand.clear();
            int p=zpos[i];
            int lo=max(0,p-W), hi=min(N-1,p+W);
            for(int q=lo;q<=hi;q++){ int j=lzorder[q]; if(j!=i) cand.push_back({dist(i,j),j}); }
            if((int)cand.size()<K && (lo>0 || hi<N-1)){
                int lo2=max(0,p-4*W), hi2=min(N-1,p+4*W);
                cand.clear();
                for(int q=lo2;q<=hi2;q++){ int j=lzorder[q]; if(j!=i) cand.push_back({dist(i,j),j}); }
            }
            int kk=min((int)cand.size(),K);
            partial_sort(cand.begin(),cand.begin()+kk,cand.end());
            for(int t=0;t<kk;t++) nbr[(size_t)i*K+t]=cand[t].second;
        }
    }

    // ---- plain nearest-neighbor construction (kept as nnConstruct, reused as the leaf-module
    // builder for the spectral pipeline below and as a fallback) ----
    auto nnConstruct=[&](const vector<int>& nodes)->vector<int>{
        int m=(int)nodes.size();
        if(m<=1) return nodes;
        vector<char> usedL(N,0); // reuse full-size bool array (cheap, cleared per call via touched-list)
        vector<int> touched; touched.reserve(m);
        vector<int> res; res.reserve(m);
        int cur=nodes[0]; usedL[cur]=1; touched.push_back(cur); res.push_back(cur);
        for(int step=1; step<m; step++){
            int best=-1; double bd=1e300;
            for(int t=0;t<K;t++){ int j=nbr[(size_t)cur*K+t]; if(j>=0 && !usedL[j]) { /* candidate must be IN nodes-set too, checked via K only if K-neighbor happens to be in this module; fallback below handles general case */ }
            }
            // K-neighbor shortcut only valid if that neighbor is part of THIS module; for safety
            // (module membership isn't a KNN-list concept) fall back to direct scan over `nodes`.
            for(int idx=0; idx<m; idx++){ int j=nodes[idx]; if(usedL[j]) continue; double d=dist(cur,j); if(d<bd){bd=d;best=j;} }
            usedL[best]=1; touched.push_back(best); res.push_back(best); cur=best;
        }
        return res;
    };

    // ---- Density-Run "Lego Block" construction (South directive, 2026-07-10: clustering is
    // signal, not cost) ----
    // The recursive-Fiedler spectral bisection (measured, rejected) blew the time budget at
    // large N because its DECOMPOSITION step was expensive (recompute an eigenvector at every
    // level: 3.5s alone at N=200000). The architecture -- decompose -> solve each piece locally
    // -> macro-route centroids -> stitch -- was sound; only the decomposition method was too
    // costly. This swaps in a near-free decomposition: sort points by Morton (Z-order) code,
    // then read the GAP between Morton-consecutive points as a density "frequency signature" --
    // a small, stable gap means "inside a dense block", a spike means "crossing a boundary
    // between blocks". Segmenting on gap spikes is a single O(N log N) pass (dominated by the
    // sort). This is also the direct architectural fix for the O(cluster^2) blowup measured in
    // the OLD uniform-grid KNN builder: instead of forcing one global grid resolution over a
    // possibly very non-uniform density, each detected block is a "Lego symbol" (size, centroid)
    // handled at its own natural scale.
    static const int TARGET_MODULE = 500;
    fprintf(stderr,"[lego] start N=%d t=%.1fms (reusing zorder from KNN build)\n", N, el_ms());

    vector<vector<int>> leafModules;
    {
        vector<double> gap(N-1);
        for(int p=0;p<N-1;p++) gap[p]=dist(lzorder[p],lzorder[p+1]);
        // Threshold targets a BLOCK COUNT (~N/TARGET_MODULE), not an arbitrary multiplier: take
        // the top `targetBoundaries` largest gaps as the boundaries. This keeps M close to the
        // same scale that worked for the spectral module design (~400 modules at N=200000),
        // regardless of how tight or loose the natural clustering is -- a fixed multiplier on
        // the median (tried first) over-segmented badly (M=~30000 on a real clustered proxy,
        // exposing an unrelated O(M^2) macro-route bug -- see below).
        vector<double> sg=gap; sort(sg.begin(),sg.end());
        int targetBoundaries = max(1, N/TARGET_MODULE);
        double threshold = (targetBoundaries < (int)sg.size()) ? sg[sg.size()-targetBoundaries] : 1e300;
        threshold = max(threshold, 1e-9);
        int start=0;
        for(int p=0;p<N-1;p++){
            bool boundary = gap[p]>threshold || (p-start+1)>=TARGET_MODULE;
            if(boundary){
                leafModules.emplace_back(lzorder.begin()+start, lzorder.begin()+p+1);
                start=p+1;
            }
        }
        if(start<N) leafModules.emplace_back(lzorder.begin()+start, lzorder.end());
    }
    fprintf(stderr,"[lego] blocks detected: %zu t=%.1fms\n", leafModules.size(), el_ms());

    // Each leaf module: tiny NN cycle -> cut longest edge -> open path
    vector<vector<int>> modulePaths(leafModules.size());
    vector<array<double,2>> moduleCentroid(leafModules.size());
    for(size_t mi=0; mi<leafModules.size(); mi++){
        auto& mod=leafModules[mi];
        vector<int> cyc = nnConstruct(mod);
        int m=(int)cyc.size();
        double cx=0,cy=0; for(int c: cyc){ cx+=X[c]; cy+=Y[c]; }
        moduleCentroid[mi]={cx/m, cy/m};
        if(m<=2){ modulePaths[mi]=cyc; continue; }
        int cutAt=0; double worst=-1;
        for(int i=0;i<m;i++){ double d=dist(cyc[i],cyc[(i+1)%m]); if(d>worst){ worst=d; cutAt=i; } }
        vector<int> path; path.reserve(m);
        for(int i=0;i<m;i++) path.push_back(cyc[(cutAt+1+i)%m]);
        modulePaths[mi]=std::move(path);
    }
    fprintf(stderr,"[lego] module paths built t=%.1fms\n", el_ms());

    // Macro-route: order modules via NN + light 2-opt over centroids (module count is small)
    int M=(int)leafModules.size();
    vector<int> modOrder(M); iota(modOrder.begin(),modOrder.end(),0);
    {
        auto cd=[&](int a,int b){ double dx=moduleCentroid[a][0]-moduleCentroid[b][0], dy=moduleCentroid[a][1]-moduleCentroid[b][1]; return sqrt(dx*dx+dy*dy); };
        if(M<=4000){ // O(M^2) NN affordable at this scale (M should normally be ~N/TARGET_MODULE)
            vector<char> usedM(M,0); vector<int> nn(M);
            int cur=0; usedM[0]=1; nn[0]=0;
            for(int s=1;s<M;s++){ int best=-1; double bd=1e300; for(int j=0;j<M;j++) if(!usedM[j]){ double d=cd(cur,j); if(d<bd){bd=d;best=j;} } usedM[best]=1; nn[s]=best; cur=best; }
            modOrder=nn;
        } // else: safety net -- keep modOrder in natural (Morton-block) detection order, O(1)
        if(M<=2000){ // O(M^2) 2-opt affordable at this scale
            bool improved=true; int pass=0;
            while(improved && pass<20 && el_ms()<TL_MS*0.15){
                improved=false; pass++;
                for(int i=1;i<M-1 && el_ms()<TL_MS*0.15;i++) for(int j=i+1;j<M;j++){
                    double before=cd(modOrder[i-1],modOrder[i])+cd(modOrder[j],modOrder[(j+1)%M]);
                    double after=cd(modOrder[i-1],modOrder[j])+cd(modOrder[i],modOrder[(j+1)%M]);
                    if(after+1e-9<before){ reverse(modOrder.begin()+i, modOrder.begin()+j+1); improved=true; }
                }
            }
        }
    }
    fprintf(stderr,"[lego] macro-route done M=%d t=%.1fms\n", M, el_ms());

    // Stitch modules in macro order into the global initial tour, picking each module's
    // traversal direction to best connect to the running tail.
    vector<int> order; order.reserve(N);
    {
        int tail=-1;
        for(int mi: modOrder){
            auto& path=modulePaths[mi];
            if(path.empty()) continue;
            bool rev=false;
            if(tail>=0 && path.size()>1){
                double dFwd=dist(tail,path.front()), dRev=dist(tail,path.back());
                rev = dRev<dFwd;
            }
            if(rev) for(int i=(int)path.size()-1;i>=0;i--) order.push_back(path[i]);
            else for(int c: path) order.push_back(c);
            tail=order.back();
        }
        // safety: ensure zeroNewId leads (rotate) and every node present exactly once
        vector<char> seen(N,0); bool okStruct=((int)order.size()==N);
        if(okStruct) for(int c: order){ if(seen[c]){ okStruct=false; break; } seen[c]=1; }
        if(!okStruct){
            fprintf(stderr,"[lego] STRUCTURAL FAILURE, falling back to plain NN construction\n");
            vector<int> allNodes(N); iota(allNodes.begin(),allNodes.end(),0);
            order = nnConstruct(allNodes);
        }
        // rotate so zeroNewId leads
        auto it=find(order.begin(),order.end(),zeroNewId);
        rotate(order.begin(), it, order.end());
    }
    fprintf(stderr,"[lego] TOTAL construction time t=%.1fms\n", el_ms());

    vector<int> pos(N);
    for(int i=0;i<N;i++) pos[order[i]]=i;

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
                    int lo=A, hi=B; if(lo>hi) swap(lo,hi);
                    int inner=hi-lo;
                    bool rev_inner = (inner <= N - inner);
                    int R_start = rev_inner ? (lo + 1) : (hi + 1);
                    int R_end   = rev_inner ? hi : (lo + N);
                    
                    auto get_new_pos = [&](int p) -> int {
                        int p_virt = (p >= R_start || (!rev_inner && p + N <= R_end)) ? (p < R_start ? p + N : p) : -1;
                        if(p_virt >= R_start && p_virt <= R_end) {
                            return (R_end - (p_virt - R_start)) % N;
                        }
                        return p;
                    };
                    
                    auto get_new_order = [&](int p) -> int {
                        int p_virt = (p >= R_start || (!rev_inner && p + N <= R_end)) ? (p < R_start ? p + N : p) : -1;
                        if(p_virt >= R_start && p_virt <= R_end) {
                            return order[(R_end - (p_virt - R_start)) % N];
                        }
                        return order[p];
                    };

                    // seek best improving second 2-opt among the two moved ends c2,c4
                    double bestG2=1e-7; int bcc=-1,bdd=-1,bcf=-1;
                    for(int side=0; side<2; side++){
                        int cc=(side==0)?c2:c4;
                        int pcc=get_new_pos(pos[cc]);
                        for(int dd=0; dd<2; dd++){
                            int pe=(dd==0)?nextIdx(pcc):prevIdx(pcc);
                            int ce=get_new_order(pe);
                            double dce=dist(cc,ce);
                            for(int u=0;u<K;u++){
                                int cf=nbr[(size_t)cc*K+u]; if(cf<0) break;
                                double dcf=dist(cc,cf);
                                if(dcf>=dce) break;
                                int pf=get_new_pos(pos[cf]);
                                int pg=(dd==0)?nextIdx(pf):prevIdx(pf);
                                int cg=get_new_order(pg);
                                if(cg==cc||cf==ce) continue;
                                double g2=(dce+dist(cf,cg))-(dcf+dist(ce,cg));
                                if(g2>bestG2){ bestG2=g2; bcc=cc; bdd=dd; bcf=cf; }
                            }
                        }
                    }
                    if(bcc>=0 && g1full+bestG2>1e-7){
                        applyMove(A,B); // NOW apply the first move for real
                        int pcc=pos[bcc], pf=pos[bcf];
                        int A2=(bdd==0)?pcc:prevIdx(pf), B2=(bdd==0)?pf:prevIdx(pcc);
                        applyMove(A2,B2);                  // apply second move: keep the chain
                        dontlook2[c1]=dontlook2[c2]=dontlook2[c3]=dontlook2[c4]=0;
                        dontlook[c1]=dontlook[c2]=dontlook[c3]=dontlook[c4]=0;
                        dontlook2[bcc]=dontlook[bcc]=0;
                        improved=true; anyImp=true;
                    }
                }
            }
            if(!improved) dontlook2[c1]=1;
        }
        return anyImp;
    };

    // Local search to convergence: 2-opt + Or-opt, then LK to break the floor, repeat.
    auto localSearch=[&](){
        double total_twoOpt=0, total_orOpt=0, total_lk=0;
        int calls_twoOpt=0, calls_orOpt=0, calls_lk=0;
        while(el_ms()<TL_MS){
            bool a=false; 
            while(el_ms()<TL_MS){ 
                double t0=el_ms(); 
                bool imp=twoOptPass(); 
                total_twoOpt += el_ms()-t0; calls_twoOpt++;
                if(!imp) break; 
                a=true; 
            }
            bool b=false;
            if(N<=50000){ 
                fill(dontlook.begin(),dontlook.end(),0); 
                double t0=el_ms();
                b=orOptPass(); 
                total_orOpt += el_ms()-t0; calls_orOpt++;
                if(b) fill(dontlook.begin(),dontlook.end(),0); 
            }
            bool c=false;
            if(N>=8000){ 
                fill(dontlook2.begin(),dontlook2.end(),0); 
                double t0=el_ms();
                c=lkPass(); 
                total_lk += el_ms()-t0; calls_lk++;
                if(c) fill(dontlook.begin(),dontlook.end(),0); 
            }
            if(!a && !b && !c) break;
        }
        fprintf(stderr,"[lego] localSearch loops done! Time Breakdown:\n");
        fprintf(stderr,"       2-opt:  %d calls, %.1f ms\n", calls_twoOpt, total_twoOpt);
        fprintf(stderr,"       or-opt: %d calls, %.1f ms\n", calls_orOpt, total_orOpt);
        fprintf(stderr,"       LK:     %d calls, %.1f ms\n", calls_lk, total_lk);
    };
    localSearch();
    fprintf(stderr,"[lego] localSearch total time t=%.1fms TL_MS=%.1f\n", el_ms(), TL_MS);

    // exact penalized cost of the output tour (rotation fixed: city 0 leads);
    // dir=0 forward, dir=1 reversed traversal of the cyclic order
    auto evalDir=[&](int dir)->double{
        int z=pos[zeroNewId]; double L=0; int prev=zeroNewId;
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
        long long ilsIter=0;
        while(el_ms()<TL_MS){
            ilsIter++;
            if(ilsIter%2000==0) fprintf(stderr,"[lego] ILS iter=%lld t=%.1fms\n", ilsIter, el_ms());
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
        fprintf(stderr,"[lego] ILS done totalIter=%lld t=%.1fms\n", ilsIter, el_ms());
        order=best; for(int i=0;i<N;i++) pos[order[i]]=i;
    }

    // ---- materialize chosen direction with city 0 leading ----
    TL_MS += RESERVE;
    vector<int> seq(N);
    { int z=pos[zeroNewId]; for(int i=0;i<N;i++) seq[i]= bestDir==0? order[(z+i)%N] : order[(z-i+N)%N]; }

    // ---- endgame touch-up: swap a nearby prime into each penalized source slot (exact delta) ----
    if(N>=12){
        auto sAt=[&](int p)->int{ return p<N? seq[p]:seq[0]; };
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
    for(int i=0;i<N;i++){ out+=to_string(perm[seq[i]]); out+='\n'; }
    out+="0\n";
    fwrite(out.data(),1,out.size(),stdout);
    fflush(stdout);
    _Exit(0);
}
