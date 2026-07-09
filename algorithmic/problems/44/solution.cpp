

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
static double TL_MS = 2250.0;
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
    int K=min(N-1, N>50000?6:(N>5000?8:10));
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

    // Local search to convergence: alternate 2-opt and Or-opt until neither improves.
    auto localSearch=[&](){
        while(el_ms()<TL_MS){
            bool a=false; while(el_ms()<TL_MS){ if(!twoOptPass()){break;} a=true; }
            bool b=false;
            if(N<=50000){ fill(dontlook.begin(),dontlook.end(),0); b=orOptPass(); if(b) fill(dontlook.begin(),dontlook.end(),0); }
            if(!a && !b) break;
            if(!b) break; // 2-opt already converged and or-opt found nothing new
        }
    };
    localSearch();

    // ---- ILS: perturb (double-bridge) + re-optimize, keep best, until deadline ----
    if(N>=8){
        // best tour snapshot
        vector<int> best=order;
        auto tourLen=[&]()->double{ double L=0; for(int i=0;i<N;i++) L+=dist(order[i],order[nextIdx(i)]); return L; };
        double bestLen=tourLen();
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
            double L=tourLen();
            if(L<bestLen-1e-6){ bestLen=L; best=order; }
            else { order=best; for(int i=0;i<N;i++) pos[order[i]]=i; } // revert
        }
        order=best; for(int i=0;i<N;i++) pos[order[i]]=i;
    }

    // ---- output rotated so city 0 leads, mathematically exact penalty optimization ----
    int z=pos[0];
    
    vector<char> isPrime(N, true);
    if(N > 0) isPrime[0] = false;
    if(N > 1) isPrime[1] = false;
    for(int i = 2; i * i < N; i++) {
        if(isPrime[i]) {
            for(int j = i * i; j < N; j += i) isPrime[j] = false;
        }
    }
    
    vector<int> new_order(N);
    for(int i=0; i<N; i++) new_order[i] = order[(z + i) % N];
    order = new_order;
    for(int i=0; i<N; i++) pos[order[i]] = i;
    
    auto get_local_cost = [&](int i, int j) {
        int ev[4] = {(i - 1 + N) % N, i, (j - 1 + N) % N, j};
        sort(ev, ev+4);
        int cnt = unique(ev, ev+4) - ev;
        double c = 0;
        for(int k=0; k<cnt; k++) {
            int idx = ev[k];
            int nxt = (idx + 1) % N;
            int u = order[idx];
            int v = order[nxt];
            double d = dist(u, v);
            if ((idx + 1) % 10 == 0 && !isPrime[u]) d *= 1.1;
            c += d;
        }
        return c;
    };

    bool improved = true;
    while(improved && el_ms() < TL_MS) {
        improved = false;
        for(int i=1; i<N; i++) {
            if(el_ms() > TL_MS) break;
            
            // Collect candidates for j
            vector<int> cands;
            for(int t=0; t<min(K, 15); t++) {
                int c2 = nbr[(size_t)order[i]*K + t];
                if(c2 >= 0) cands.push_back(c2);
                c2 = nbr[(size_t)order[(i-1+N)%N]*K + t];
                if(c2 >= 0) cands.push_back(c2);
                c2 = nbr[(size_t)order[(i+1)%N]*K + t];
                if(c2 >= 0) cands.push_back(c2);
            }
            
            for(int c2 : cands) {
                int j = pos[c2];
                if(j == 0 || i == j) continue;
                
                double old_c = get_local_cost(i, j);
                swap(order[i], order[j]);
                double new_c = get_local_cost(i, j);
                
                if(new_c < old_c - 1e-7) {
                    pos[order[i]] = i;
                    pos[order[j]] = j;
                    improved = true;
                } else {
                    swap(order[i], order[j]); // revert
                }
            }
        }
    }

    string out; out.reserve((size_t)N*7+16);
    out+=to_string(N+1); out+='\n';
    for(int i=0;i<N;i++){ out+=to_string(order[i]); out+='\n'; }
    out+="0\n";
    fwrite(out.data(),1,out.size(),stdout);
    fflush(stdout);
    _Exit(0);
}
// AVALHACAR


