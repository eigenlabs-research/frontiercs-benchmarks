// Model: South Core - Created By Alecstec
#include <bits/stdc++.h>
using namespace std;
static chrono::steady_clock::time_point T0;
static double n2 = 2400.0;
static inline double fn1(){ return chrono::duration<double,milli>(chrono::steady_clock::now()-T0).count(); }
static int N;
static vector<double> X, Y;
static inline double fn3(int a,int b){ double dx=X[a]-X[b],dy=Y[a]-Y[b]; return sqrt(dx*dx+dy*dy); }
static vector<char> n3; static size_t ip=0;
static inline long long fn2(){
    while(ip<n3.size() && (n3[ip]<'0'||n3[ip]>'9') && n3[ip]!='-') ip++;
    bool neg=false; if(ip<n3.size()&&n3[ip]=='-'){neg=true;ip++;}
    long long v=0; while(ip<n3.size()&&n3[ip]>='0'&&n3[ip]<='9'){v=v*10+(n3[ip]-'0');ip++;}
    return neg?-v:v;
}
int main(){
    T0=chrono::steady_clock::now();
    { size_t cap=1<<20; n3.resize(cap); size_t len=0;
      while(true){ if(len==cap){cap<<=1;n3.resize(cap);} size_t g=fread(n3.data()+len,1,cap-len,stdin); if(!g)break; len+=g; } n3.resize(len); }
    N=(int)fn2();
    if(N<=0){ printf("1\n0\n"); return 0; }
    X.resize(N); Y.resize(N);
    for(int i=0;i<N;i++){ X[i]=(double)fn2(); Y[i]=(double)fn2(); }
    if(N==1){ printf("2\n0\n0\n"); return 0; }
    if(N==2){ printf("3\n0\n1\n0\n"); return 0; }
    vector<char> n5((size_t)N,0);
    { vector<char> comp((size_t)N,0); for(long long i=2;i<N;i++) if(!comp[i]){ n5[i]=1; for(long long q=i*i;q<N;q+=i) comp[q]=1; } }
    vector<int> n6(N), n7(N);
    {
        double mnx=X[0],mxx=X[0],mny=Y[0],mxy=Y[0];
        for(int i=1;i<N;i++){ mnx=min(mnx,X[i]); mxx=max(mxx,X[i]); mny=min(mny,Y[i]); mxy=max(mxy,Y[i]); }
        double rx=max(1e-9,mxx-mnx), ry=max(1e-9,mxy-mny);
        auto fn4=[](uint32_t v)->uint64_t{
            uint64_t x=v & 0xFFFFFu;
            x=(x|(x<<16))&0x0000FFFF0000FFFFULL;
            x=(x|(x<<8))&0x00FF00FF00FF00FFULL;
            x=(x|(x<<4))&0x0F0F0F0F0F0F0F0FULL;
            x=(x|(x<<2))&0x3333333333333333ULL;
            x=(x|(x<<1))&0x5555555555555555ULL;
            return x;
        };
        vector<uint64_t> n39(N);
        for(int i=0;i<N;i++){
            uint32_t qx=(uint32_t)min(1048575.0, floor((X[i]-mnx)/rx*1048575.0));
            uint32_t qy=(uint32_t)min(1048575.0, floor((Y[i]-mny)/ry*1048575.0));
            n39[i]=fn4(qx)|(fn4(qy)<<1);
            n6[i]=i;
        }
        if(N>=150000){
            sort(n6.begin(),n6.end(),[&](int a,int b){ return n39[a]<n39[b]; });
        }
        for(int newId=0;newId<N;newId++) n7[n6[newId]]=newId;
        vector<double> X2(N),Y2(N);
        for(int newId=0;newId<N;newId++){ X2[newId]=X[n6[newId]]; Y2[newId]=Y[n6[newId]]; }
        X.swap(X2); Y.swap(Y2);
    }
    const int n8 = n7[0];
    vector<char> n9(N);
    for(int newId=0;newId<N;newId++) n9[newId]=n5[n6[newId]];
    if(N>100000) n2 -= 20.0;
    double n10 = N>150000?220.0:(N>50000?90.0:(N>5000?50.0:40.0));
    n2 -= n10;
    double minx=X[0],maxx=X[0],miny=Y[0],maxy=Y[0];
    for(int i=1;i<N;i++){ minx=min(minx,X[i]);maxx=max(maxx,X[i]);miny=min(miny,Y[i]);maxy=max(maxy,Y[i]); }
    double w=max(1.0,maxx-minx), h=max(1.0,maxy-miny);
    int G=max(1,(int)floor(sqrt((double)N/2.0)));
    double cw=w/G, ch=h/G;
    auto fn22=[&](double x){ int c=(int)((x-minx)/cw); return c<0?0:(c>=G?G-1:c); };
    auto fn23=[&](double y){ int c=(int)((y-miny)/ch); return c<0?0:(c>=G?G-1:c); };
    vector<int> n11(N), cnt(G*G+1,0);
    for(int i=0;i<N;i++){ int c=fn22(X[i])*G+fn23(Y[i]); n11[i]=c; cnt[c+1]++; }
    for(int i=0;i<G*G;i++) cnt[i+1]+=cnt[i];
    vector<int> n12(N); { vector<int> tmp=cnt; for(int i=0;i<N;i++) n12[tmp[n11[i]]++]=i; }
    int K;
    if(N>50000)                      K=min(N-1,10);
    else if(N>=16000 && N<=36000)    K=min(N-1,40);
    else if(N>5000)                  K=min(N-1,24);
    else                             K=min(N-1,10);
    vector<int> n13((size_t)N*K,-1);
    vector<int> n14(N);
    {
        double mnx=X[0],mxx=X[0],mny=Y[0],mxy=Y[0];
        for(int i=1;i<N;i++){ mnx=min(mnx,X[i]); mxx=max(mxx,X[i]); mny=min(mny,Y[i]); mxy=max(mxy,Y[i]); }
        double rx=max(1e-9,mxx-mnx), ry=max(1e-9,mxy-mny);
        auto fn5=[](uint32_t v)->uint64_t{
            uint64_t x=v & 0xFFFFFu;
            x=(x|(x<<16))&0x0000FFFF0000FFFFULL;
            x=(x|(x<<8))&0x00FF00FF00FF00FFULL;
            x=(x|(x<<4))&0x0F0F0F0F0F0F0F0FULL;
            x=(x|(x<<2))&0x3333333333333333ULL;
            x=(x|(x<<1))&0x5555555555555555ULL;
            return x;
        };
        vector<uint64_t> n15(N);
        for(int i=0;i<N;i++){
            uint32_t qx=(uint32_t)min(1048575.0, floor((X[i]-mnx)/rx*1048575.0));
            uint32_t qy=(uint32_t)min(1048575.0, floor((Y[i]-mny)/ry*1048575.0));
            n15[i]=fn5(qx)|(fn5(qy)<<1);
            n14[i]=i;
        }
        sort(n14.begin(), n14.end(), [&](int a,int b){ return n15[a]<n15[b]; });
    }
    {
        vector<int> n16(N); iota(n16.begin(), n16.end(), 0);
        {
            vector<array<int,3>> stk; stk.push_back({0,N,0});
            while(!stk.empty()){
                array<int,3> tb=stk.back(); stk.pop_back();
                int l=tb[0], r=tb[1], ax=tb[2];
                if(r-l<=1) continue;
                int m=(l+r)>>1;
                nth_element(n16.begin()+l, n16.begin()+m, n16.begin()+r,
                    [&](int a,int b){ return ax? Y[a]<Y[b] : X[a]<X[b]; });
                stk.push_back({l,m,ax^1});
                stk.push_back({m+1,r,ax^1});
            }
        }
        struct St1{ int l,r,ax; double g2; };
        vector<pair<double,int>> n17; n17.reserve(K+1);
        vector<St1> n18; n18.reserve(64);
        for(int i=0;i<N;i++){
            n17.clear();
            double qx=X[i], qy=Y[i];
            n18.push_back({0,N,0,0.0});
            while(!n18.empty()){
                St1 tq=n18.back(); n18.pop_back();
                double worst = ((int)n17.size()<K)? 1e300 : n17.front().first;
                if(tq.g2>worst) continue;
                if(tq.r<=tq.l) continue;
                int m=(tq.l+tq.r)>>1; int j=n16[m];
                if(j!=i){
                    double dx=qx-X[j], dy=qy-Y[j], d2=dx*dx+dy*dy;
                    if((int)n17.size()<K){
                        n17.push_back({d2,j}); push_heap(n17.begin(),n17.end());
                    } else if(d2<n17.front().first){
                        pop_heap(n17.begin(),n17.end()); n17.back()={d2,j}; push_heap(n17.begin(),n17.end());
                    }
                }
                double split = tq.ax? Y[j] : X[j];
                double gap = (tq.ax? qy : qx) - split;
                double g2 = gap*gap;
                if(gap<0){
                    if(m+1<tq.r)  n18.push_back({m+1,tq.r,tq.ax^1,g2});
                    if(tq.l<m)    n18.push_back({tq.l,m,tq.ax^1,0.0});
                } else {
                    if(tq.l<m)    n18.push_back({tq.l,m,tq.ax^1,g2});
                    if(m+1<tq.r)  n18.push_back({m+1,tq.r,tq.ax^1,0.0});
                }
            }
            sort_heap(n17.begin(),n17.end());
            for(int t=0;t<(int)n17.size();t++) n13[(size_t)i*K+t]=n17[t].second;
        }
    }
    auto fn6=[&](const vector<int>& nodes)->vector<int>{
        int m=(int)nodes.size();
        if(m<=1) return nodes;
        vector<char> n34(N,0);
        vector<int> n33; n33.reserve(m);
        vector<int> res; res.reserve(m);
        int cur=nodes[0]; n34[cur]=1; n33.push_back(cur); res.push_back(cur);
        bool n32 = (m==N);
        for(int step=1; step<m; step++){
            int best=-1; double bd=1e300;
            if(n32){
                for(int t=0;t<K;t++){ int j=n13[(size_t)cur*K+t]; if(j>=0 && !n34[j]){ best=j; bd=fn3(cur,j); break; } }
            }
            if(best<0)
                for(int idx=0; idx<m; idx++){ int j=nodes[idx]; if(n34[j]) continue; double d=fn3(cur,j); if(d<bd){bd=d;best=j;} }
            n34[best]=1; n33.push_back(best); res.push_back(best); cur=best;
        }
        return res;
    };
    bool n23 = (N>50000);
    vector<int> order; order.reserve(N);
    if(!n23){
        vector<int> allNodes(N); iota(allNodes.begin(),allNodes.end(),0);
        order = fn6(allNodes);
        auto it=find(order.begin(),order.end(),n8);
        rotate(order.begin(), it, order.end());
    } else {
    static const int CN1 = 500;
    vector<vector<int>> n19;
    {
        vector<double> gap(N-1);
        for(int p=0;p<N-1;p++) gap[p]=fn3(n14[p],n14[p+1]);
        vector<double> sg=gap; sort(sg.begin(),sg.end());
        int n26 = max(1, N/CN1);
        double n27 = (n26 < (int)sg.size()) ? sg[sg.size()-n26] : 1e300;
        n27 = max(n27, 1e-9);
        int start=0;
        for(int p=0;p<N-1;p++){
            bool boundary = gap[p]>n27 || (p-start+1)>=CN1;
            if(boundary){
                n19.emplace_back(n14.begin()+start, n14.begin()+p+1);
                start=p+1;
            }
        }
        if(start<N) n19.emplace_back(n14.begin()+start, n14.end());
    }
    vector<vector<int>> n20(n19.size());
    vector<array<double,2>> n21(n19.size());
    for(size_t mi=0; mi<n19.size(); mi++){
        auto& mod=n19[mi];
        vector<int> n36 = fn6(mod);
        int m=(int)n36.size();
        double cx=0,cy=0; for(int c: n36){ cx+=X[c]; cy+=Y[c]; }
        n21[mi]={cx/m, cy/m};
        if(m<=2){ n20[mi]=n36; continue; }
        int n37=0; double worst=-1;
        for(int i=0;i<m;i++){ double d=fn3(n36[i],n36[(i+1)%m]); if(d>worst){ worst=d; n37=i; } }
        vector<int> path; path.reserve(m);
        for(int i=0;i<m;i++) path.push_back(n36[(n37+1+i)%m]);
        n20[mi]=std::move(path);
    }
    int M=(int)n19.size();
    vector<int> n22(M); iota(n22.begin(),n22.end(),0);
    {
        auto cd=[&](int a,int b){ double dx=n21[a][0]-n21[b][0], dy=n21[a][1]-n21[b][1]; return sqrt(dx*dx+dy*dy); };
        if(M<=4000){
            vector<char> n35(M,0); vector<int> nn(M);
            int cur=0; n35[0]=1; nn[0]=0;
            for(int s=1;s<M;s++){ int best=-1; double bd=1e300; for(int j=0;j<M;j++) if(!n35[j]){ double d=cd(cur,j); if(d<bd){bd=d;best=j;} } n35[best]=1; nn[s]=best; cur=best; }
            n22=nn;
        }
        if(M<=2000){
            bool improved=true; int pass=0;
            while(improved && pass<20 && fn1()<n2*0.15){
                improved=false; pass++;
                for(int i=1;i<M-1 && fn1()<n2*0.15;i++) for(int j=i+1;j<M;j++){
                    double before=cd(n22[i-1],n22[i])+cd(n22[j],n22[(j+1)%M]);
                    double after=cd(n22[i-1],n22[j])+cd(n22[i],n22[(j+1)%M]);
                    if(after+1e-9<before){ reverse(n22.begin()+i, n22.begin()+j+1); improved=true; }
                }
            }
        }
    }
    {
        int tail=-1;
        for(int mi: n22){
            auto& path=n20[mi];
            if(path.empty()) continue;
            bool rev=false;
            if(tail>=0 && path.size()>1){
                double dFwd=fn3(tail,path.front()), dRev=fn3(tail,path.back());
                rev = dRev<dFwd;
            }
            if(rev) for(int i=(int)path.size()-1;i>=0;i--) order.push_back(path[i]);
            else for(int c: path) order.push_back(c);
            tail=order.back();
        }
        vector<char> seen(N,0); bool okStruct=((int)order.size()==N);
        if(okStruct) for(int c: order){ if(seen[c]){ okStruct=false; break; } seen[c]=1; }
        if(!okStruct){
            if(N>50000) order = n14;
            else { vector<int> allNodes(N); iota(allNodes.begin(),allNodes.end(),0); order = fn6(allNodes); }
        }
        auto it=find(order.begin(),order.end(),n8);
        rotate(order.begin(), it, order.end());
    }
    }
    vector<int> pos(N);
    for(int i=0;i<N;i++) pos[order[i]]=i;
    auto fn14=[&](int i){ return i+1<N?i+1:0; };
    auto fn15=[&](int i){ return i>0?i-1:N-1; };
    auto fn13=[&](int e1,int e2){
        int lo=e1, hi=e2; if(lo>hi) swap(lo,hi);
        int inner=hi-lo;
        if(inner<=N-inner){
            int i=lo+1, j=hi;
            while(i<j){ int a=order[i],b=order[j]; order[i]=b; order[j]=a; pos[a]=j; pos[b]=i; ++i; --j; }
        } else {
            int li=hi+1, lj=lo+N;
            while(li<lj){ int ai=li%N, aj=lj%N; int u=order[ai],v=order[aj]; order[ai]=v; order[aj]=u; pos[v]=ai; pos[u]=aj; ++li; --lj; }
        }
    };
    vector<char> n28(N,0);
    vector<char> n29(N,0);
    vector<int> q(N);
    for(int i=0;i<N;i++) q[i]=i;
    int clock=0;
    auto fn9=[&]()->bool{
        bool anyImp=false;
        for(int qi=0; qi<N; qi++){
            if(((++clock)&1023)==0 && fn1()>n2) return anyImp;
            int c1=q[qi];
            if(n28[c1]) continue;
            bool improved=false;
            for(int dir=0; dir<2 && !improved; dir++){
                int p1=pos[c1];
                int p2=(dir==0)?fn14(p1):fn15(p1);
                int c2=order[p2];
                double d12=fn3(c1,c2);
                for(int t=0;t<K;t++){
                    int c3=n13[(size_t)c1*K+t]; if(c3<0) break;
                    double d13=fn3(c1,c3);
                    if(d13>=d12) break;
                    int p3=pos[c3];
                    int p4=(dir==0)?fn14(p3):fn15(p3);
                    int c4=order[p4];
                    if(c4==c1||c3==c2) continue;
                    double before=d12+fn3(c3,c4);
                    double after =d13+fn3(c2,c4);
                    if(after+1e-7<before){
                        if(dir==0) fn13(p1,p3);
                        else       fn13(p4,p2);
                        n28[c1]=n28[c2]=n28[c3]=n28[c4]=0;
                        improved=true; anyImp=true;
                        break;
                    }
                }
            }
            if(!improved) n28[c1]=1;
        }
        return anyImp;
    };
    auto fn10=[&]()->bool{
        bool anyImp=false;
        for(int qi=0; qi<N; qi++){
            if(((++clock)&1023)==0 && fn1()>n2) return anyImp;
            int s0=q[qi];
            if(n28[s0]) continue;
            bool moved=false;
            for(int L=1; L<=3 && !moved; L++){
                int is=pos[s0];
                int ie=is; for(int t=1;t<L;t++) ie=fn14(ie);
                int segEnd=order[ie];
                int pprev=fn15(is), pnext=fn14(ie);
                if(pprev==ie || pnext==is) break;
                int cprev=order[pprev], cnext=order[pnext];
                double removed = fn3(cprev,s0)+fn3(segEnd,cnext)-fn3(cprev,cnext);
                if(removed<=1e-7) continue;
                for(int side=0; side<2 && !moved; side++){
                    int anchorCity = side==0? s0 : segEnd;
                    for(int t=0;t<K;t++){
                        int c=n13[(size_t)anchorCity*K+t]; if(c<0) break;
                        int pc=pos[c];
                        bool inside=false; { int p=is; for(int u=0;u<L;u++){ if(p==pc){inside=true;break;} p=fn14(p);} }
                        if(inside || c==cprev) continue;
                        int pcn=fn14(pc); int cn=order[pcn];
                        if(cn==s0) continue;
                        double added = fn3(c,s0)+fn3(segEnd,cn)-fn3(c,cn);
                        double addedRev = fn3(c,segEnd)+fn3(s0,cn)-fn3(c,cn);
                        bool rev = addedRev+1e-9 < added;
                        double add = rev? addedRev: added;
                        if(add+1e-7 < removed){
                            int seg[3]; { int p=is; for(int u=0;u<L;u++){ seg[u]=order[p]; p=fn14(p);} }
                            if(rev){ for(int a=0,b=L-1;a<b;a++,b--) swap(seg[a],seg[b]); }
                            static vector<int> tmp; tmp.clear(); tmp.reserve(N);
                            auto fn24=[&](int city)->bool{ for(int u=0;u<L;u++) if(seg[u]==city) return true; return false; };
                            for(int idx2=0; idx2<N; idx2++){
                                int city=order[idx2];
                                if(fn24(city)) continue;
                                tmp.push_back(city);
                                if(city==c){ for(int u=0;u<L;u++) tmp.push_back(seg[u]); }
                            }
                            order.swap(tmp);
                            for(int idx2=0; idx2<N; idx2++) pos[order[idx2]]=idx2;
                            n28[cprev]=n28[cnext]=n28[s0]=n28[segEnd]=n28[c]=n28[cn]=0;
                            moved=true; anyImp=true;
                            break;
                        }
                    }
                }
            }
            if(!moved) n28[s0]=1;
        }
        return anyImp;
    };
    auto fn7=[&](int lo0,int hi0)->double{
        double c=0;
        for(int k=lo0;k<=hi0;k++){
            int a=order[k], b=order[k+1==N?0:k+1];
            double d=fn3(a,b);
            if(k%10==9 && !n9[a]) d*=1.1;
            c+=d;
        }
        return c;
    };
    const int CN2 = max(4000, N/8);
    auto fn11=[&](bool n31)->bool{
        bool anyImp=false;
        for(int qi=0; qi<N; qi++){
            if(((++clock)&511)==0 && fn1()>n2) return anyImp;
            int c1=q[qi];
            if(n29[c1]) continue;
            bool improved=false;
            for(int dir=0; dir<2 && !improved; dir++){
                int p1=pos[c1];
                int p2=(dir==0)?fn14(p1):fn15(p1);
                int c2=order[p2];
                double d12=fn3(c1,c2);
                for(int t=0;t<K && !improved;t++){
                    int c3=n13[(size_t)c1*K+t]; if(c3<0) break;
                    double d13=fn3(c1,c3);
                    if(d13>=d12) break;
                    int p3=pos[c3];
                    int p4=(dir==0)?fn14(p3):fn15(p3);
                    int c4=order[p4];
                    if(c4==c1||c3==c2) continue;
                    double g1full = (d12+fn3(c3,c4)) - (d13+fn3(c2,c4));
                    int A=(dir==0)?p1:p4, B=(dir==0)?p3:p2;
                    int lo=A, hi=B; if(lo>hi) swap(lo,hi);
                    int inner=hi-lo;
                    bool rev_inner = (inner <= N - inner);
                    int R_start = rev_inner ? (lo + 1) : (hi + 1);
                    int R_end   = rev_inner ? hi : (lo + N);
                    auto fn18 = [&](int p) -> int {
                        int p_virt = (p >= R_start || (!rev_inner && p + N <= R_end)) ? (p < R_start ? p + N : p) : -1;
                        if(p_virt >= R_start && p_virt <= R_end) {
                            return (R_end - (p_virt - R_start)) % N;
                        }
                        return p;
                    };
                    auto fn19 = [&](int p) -> int {
                        int p_virt = (p >= R_start || (!rev_inner && p + N <= R_end)) ? (p < R_start ? p + N : p) : -1;
                        if(p_virt >= R_start && p_virt <= R_end) {
                            return order[(R_end - (p_virt - R_start)) % N];
                        }
                        return order[p];
                    };
                    double bestG2=1e-7; int bcc=-1,bdd=-1,bcf=-1;
                    for(int side=0; side<2; side++){
                        int cc=(side==0)?c2:c4;
                        int pcc=fn18(pos[cc]);
                        for(int dd=0; dd<2; dd++){
                            int pe=(dd==0)?fn14(pcc):fn15(pcc);
                            int ce=fn19(pe);
                            double dce=fn3(cc,ce);
                            for(int u=0;u<K;u++){
                                int cf=n13[(size_t)cc*K+u]; if(cf<0) break;
                                double dcf=fn3(cc,cf);
                                if(dcf>=dce) break;
                                int pf=fn18(pos[cf]);
                                int pg=(dd==0)?fn14(pf):fn15(pf);
                                int cg=fn19(pg);
                                if(cg==cc||cf==ce) continue;
                                double g2=(dce+fn3(cf,cg))-(dcf+fn3(ce,cg));
                                if(g2>bestG2){ bestG2=g2; bcc=cc; bdd=dd; bcf=cf; }
                            }
                        }
                    }
                    if(bcc>=0 && !n31 && g1full+bestG2>1e-7){
                        int pcc=pos[bcc], pf=pos[bcf];
                        fn13(A,B);
                        pcc=pos[bcc]; pf=pos[bcf];
                        int A2=(bdd==0)?pcc:fn15(pf), B2=(bdd==0)?pf:fn15(pcc);
                        fn13(A2,B2);
                        n29[c1]=n29[c2]=n29[c3]=n29[c4]=0;
                        n28[c1]=n28[c2]=n28[c3]=n28[c4]=0;
                        n29[bcc]=n28[bcc]=0;
                        improved=true; anyImp=true;
                    }
                    else if(bcc>=0 && n31 && g1full+bestG2>-0.11*d12){
                        if(2*(hi-lo)>N) continue;
                        int pcc=fn18(pos[bcc]), pf=fn18(pos[bcf]);
                        int A2=(bdd==0)?pcc:fn15(pf), B2=(bdd==0)?pf:fn15(pcc);
                        int lo2=min(A2,B2), hi2=max(A2,B2);
                        if(2*(hi2-lo2)>N) continue;
                        int mn=min(lo,lo2), mx=max(hi,hi2);
                        if(mx-mn>CN2) continue;
                        int lo0=max(0,mn-1), hi0=min(N-1,mx);
                        double oldC=fn7(lo0,hi0);
                        auto fn8=[&](int k)->int{
                            int k1 = (k>lo2 && k<=hi2)? (lo2+1+hi2-k) : k;
                            int k0 = (k1>lo && k1<=hi)? (lo+1+hi-k1) : k1;
                            return order[k0];
                        };
                        double newC=0;
                        for(int k=lo0;k<=hi0;k++){
                            int a=fn8(k), b=fn8(k+1==N?0:k+1);
                            double d=fn3(a,b);
                            if(k%10==9 && !n9[a]) d*=1.1;
                            newC+=d;
                        }
                        if(oldC-newC>1e-7){
                            fn13(A,B); fn13(A2,B2);
                            n29[c1]=n29[c2]=n29[c3]=n29[c4]=0;
                            n28[c1]=n28[c2]=n28[c3]=n28[c4]=0;
                            n29[bcc]=n28[bcc]=0;
                            improved=true; anyImp=true;
                        }
                    }
                }
            }
            if(!improved) n29[c1]=1;
        }
        return anyImp;
    };
    auto fn12=[&](){
        while(fn1()<n2){
            bool a=false;
            while(fn1()<n2){
                bool imp=fn9();
                if(!imp) break;
                a=true;
            }
            bool b=false;
            if(N<=50000){
                fill(n28.begin(),n28.end(),0);
                b=fn10();
                if(b) fill(n28.begin(),n28.end(),0);
            }
            bool c=false;
            if(N>=8000){
                fill(n29.begin(),n29.end(),0);
                c=fn11(false);
                if(c) fill(n28.begin(),n28.end(),0);
            }
            if(!a && !b && !c) break;
        }
        if(N>=8000){
            while(fn1()<n2){
                fill(n29.begin(),n29.end(),0);
                bool cx=fn11(true);
                if(!cx) break;
            }
        }
    };
    fn12();
    auto fn16=[&](int dir)->double{
        int z=pos[n8]; double L=0; int prev=n8;
        for(int t=1;t<=N;t++){
            int idx = dir==0 ? (z+t)%N : (z-(t%N)+N)%N;
            int b=order[idx];
            double d=fn3(prev,b);
            if(t%10==0 && !n9[prev]) d*=1.1;
            L+=d; prev=b;
        }
        return L;
    };
    auto fn17=[&](int&dir)->double{
        double f=fn16(0), r=fn16(1);
        if(f<=r){ dir=0; return f; } dir=1; return r;
    };
    vector<int> best=order; int bestDir=0;
    double bestLen=fn17(bestDir);
    if(N>=8){
        uint64_t rng=0x9e3779b97f4a7c15ULL ^ (uint64_t)N*2654435761ULL;
        auto rnd=[&](){ rng^=rng<<7; rng^=rng>>9; return rng; };
        long long n38=0;
        while(fn1()<n2){
            n38++;
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
            fill(n28.begin(),n28.end(),0);
            while(fn1()<n2){ if(!fn9()) break; }
            if(N<=50000){ fill(n28.begin(),n28.end(),0); fn10(); while(fn1()<n2){ if(!fn9()) break; } }
            int d2=0; double L=fn17(d2);
            if(L<bestLen-1e-6){ bestLen=L; best=order; bestDir=d2; }
            else { order=best; for(int i=0;i<N;i++) pos[order[i]]=i; }
        }
        order=best; for(int i=0;i<N;i++) pos[order[i]]=i;
    }
    n2 += n10;
    vector<int> seq(N);
    { int z=pos[n8]; for(int i=0;i<N;i++) seq[i]= bestDir==0? order[(z+i)%N] : order[(z-i+N)%N]; }
    if(N>=12){
        auto fn20=[&](int p)->int{ return p<N? seq[p]:seq[0]; };
        auto fn21=[&](int t)->double{ int a=seq[t-1], b=fn20(t); double d=fn3(a,b); if(t%10==0&&!n9[a]) d*=1.1; return d; };
        int w = N<=1200? N : (N<=5000?120:(N<=20000?140:80));
        for(int rep=0;rep<4 && fn1()<n2;rep++){
            bool ch=false;
            for(int p=9;p<N;p+=10){
                if(fn1()>n2) break;
                if(n9[seq[p]]) continue;
                int lo=max(1,p-w), hi=min(N-1,p+w);
                double bd=-1e-7; int bj=-1;
                for(int j=lo;j<=hi;j++){
                    if(j==p||!n9[seq[j]]) continue;
                    int T[4]={p,p+1,j,j+1}, U[4], m=0;
                    for(int t2=0;t2<4;t2++){ bool dup=false; for(int u=0;u<m;u++) if(U[u]==T[t2]) dup=true; if(!dup) U[m++]=T[t2]; }
                    double bef=0, aft=0;
                    for(int u=0;u<m;u++) bef+=fn21(U[u]);
                    swap(seq[p],seq[j]);
                    for(int u=0;u<m;u++) aft+=fn21(U[u]);
                    swap(seq[p],seq[j]);
                    double dl=aft-bef;
                    if(dl<bd){ bd=dl; bj=j; }
                }
                if(bj>=0){ swap(seq[p],seq[bj]); ch=true; }
            }
            if(!ch) break;
        }
    }
    {
        int pb1 = 8000; int pb2 = 16000; double pb3 = 0.02;
        if(N>=pb1 && N<pb2 && N>=8){
            auto pb4=[&](int p)->int{ return p<N? seq[p]:seq[0]; };
            auto pb5=[&](int t)->double{ int a=seq[t-1], b=pb4(t); double d=fn3(a,b); if(t%10==0&&!n9[a]) d*=1.1; return d; };
            auto pb6=[&](int p,int j)->double{
                int T[4]={p,p+1,j,j+1}, U[4], m=0;
                for(int t2=0;t2<4;t2++){ bool dup=false; for(int u=0;u<m;u++) if(U[u]==T[t2]) dup=true; if(!dup) U[m++]=T[t2]; }
                double bef=0, aft=0;
                for(int u=0;u<m;u++) bef+=pb5(U[u]);
                swap(seq[p],seq[j]);
                for(int u=0;u<m;u++) aft+=pb5(U[u]);
                swap(seq[p],seq[j]);
                return aft-bef;
            };
            double pb7=0; for(int t=1;t<=N;t++) pb7+=pb5(t);
            double pb8=pb7*(1.0+pb3);
            double pb9=pb7*max(0.001, 3.0/N);
            unsigned long long pba=88172645463325252ULL;
            auto pbb=[&](){ pba^=pba<<13; pba^=pba>>7; pba^=pba<<17; return pba; };
            double pbc=pb7; long long pbd=0;
            while(pbc<pb8 && pbd<2000000){
                pbd++;
                int i=1+(int)(pbb()%(N-1)), j=1+(int)(pbb()%(N-1));
                if(i==j) continue;
                double d=pb6(i,j);
                if(d>1e-9 && d<=pb9 && pbc+d<=pb8+pb9){ swap(seq[i],seq[j]); pbc+=d; }
            }
        }
    }
    string out; out.reserve((size_t)N*7+16);
    out+=to_string(N+1); out+='\n';
    for(int i=0;i<N;i++){ out+=to_string(n6[seq[i]]); out+='\n'; }
    out+="0\n";
    fwrite(out.data(),1,out.size(),stdout);
    fflush(stdout);
    _Exit(0);
}
