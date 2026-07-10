#include <bits/stdc++.h>
using namespace std;
#ifndef GEOM_MODE
#define GEOM_MODE 2
#endif
#define GEOM_ALPHA
#define GEOM_ALPHA_ORDERED
#define GEOM_KD_SEED
#define GEOM_GROUP_DUPLICATES

struct GeomPoint { long long x,y; int city; };
struct GeomResult { vector<pair<int,int>> edges,mstEdges; vector<int> mstTour; };
struct GeomDSU { vector<int> p,s;GeomDSU(int n):p(n),s(n,1){iota(p.begin(),p.end(),0);}int f(int x){return p[x]==x?x:p[x]=f(p[x]);}bool join(int a,int b){a=f(a);b=f(b);if(a==b)return false;if(s[a]<s[b])swap(a,b);p[b]=a;s[a]+=s[b];return true;} };

struct DelaunayCore {
    struct Q { int org=-1,rot=-1,nxt=-1; };
    vector<GeomPoint> p;vector<Q> q;vector<char> alive;vector<int> freeBlocks;
    static __int128 orient128(const GeomPoint&a,const GeomPoint&b,const GeomPoint&c){return (__int128)(b.x-a.x)*(c.y-a.y)-(__int128)(b.y-a.y)*(c.x-a.x);}
    int orient(int a,int b,int c)const{__int128 v=orient128(p[a],p[b],p[c]);return(v>0)-(v<0);}
    bool inCircle(int a,int b,int c,int d)const{
        auto norm=[](const GeomPoint&v){return (__int128)v.x*v.x+(__int128)v.y*v.y;};
        __int128 na=norm(p[a])-norm(p[d]),nb=norm(p[b])-norm(p[d]),nc=norm(p[c])-norm(p[d]);
        return orient128(p[d],p[a],p[b])*nc+orient128(p[d],p[b],p[c])*na+orient128(p[d],p[c],p[a])*nb>0;
    }
    int rev(int e)const{return q[q[e].rot].rot;}int dest(int e)const{return q[rev(e)].org;}
    int lnext(int e)const{int r=q[e].rot;return q[q[rev(r)].nxt].rot;}int oprev(int e)const{return q[q[q[e].rot].nxt].rot;}
    void splice(int a,int b){int an=q[a].nxt,bn=q[b].nxt,aa=q[an].rot,bb=q[bn].rot;swap(q[aa].nxt,q[bb].nxt);swap(q[a].nxt,q[b].nxt);}
    int makeEdge(int u,int v){
        int b;if(freeBlocks.empty()){b=alive.size();alive.push_back(1);q.resize(q.size()+4);}else{b=freeBlocks.back();freeBlocks.pop_back();alive[b]=1;}
        int e=4*b;q[e]={u,e+2,e};q[e+1]={v,e+3,e+1};q[e+2]={-1,e+1,e+3};q[e+3]={-1,e,e+2};return e;
    }
    int connect(int a,int b){int e=makeEdge(dest(a),q[b].org);splice(e,lnext(a));splice(rev(e),b);return e;}
    void deleteEdge(int e){int r=rev(e);splice(e,oprev(e));splice(r,oprev(r));int b=e/4;alive[b]=0;freeBlocks.push_back(b);}
    bool leftOf(int x,int e)const{return orient(x,q[e].org,dest(e))>0;}bool rightOf(int x,int e)const{return orient(x,q[e].org,dest(e))<0;}
    pair<int,int> build(int l,int r){
        int n=r-l+1;if(n==2){int a=makeEdge(l,r);return{a,rev(a)};}
        if(n==3){int a=makeEdge(l,l+1),b=makeEdge(l+1,r);splice(rev(a),b);int s=orient(l,l+1,r);if(!s)return{a,rev(b)};int c=connect(b,a);return s>0?make_pair(a,rev(b)):make_pair(rev(c),c);}
        int m=(l+r)/2;auto [ldo,ldi]=build(l,m);auto [rdi,rdo]=build(m+1,r);
        while(true){if(leftOf(q[rdi].org,ldi)){ldi=lnext(ldi);continue;}if(rightOf(q[ldi].org,rdi)){rdi=q[rev(rdi)].nxt;continue;}break;}
        int base=connect(rev(rdi),ldi);if(q[ldi].org==q[ldo].org)ldo=rev(base);if(q[rdi].org==q[rdo].org)rdo=base;auto valid=[&](int e){return rightOf(dest(e),base);};
        while(true){
            int lc=q[rev(base)].nxt;if(valid(lc))while(inCircle(dest(base),q[base].org,dest(lc),dest(q[lc].nxt))){int t=q[lc].nxt;deleteEdge(lc);lc=t;}
            int rc=oprev(base);if(valid(rc))while(inCircle(dest(base),q[base].org,dest(rc),dest(oprev(rc)))){int t=oprev(rc);deleteEdge(rc);rc=t;}
            bool lv=valid(lc),rv=valid(rc);if(!lv&&!rv)break;if(!lv||(rv&&inCircle(dest(lc),q[lc].org,q[rc].org,dest(rc))))base=connect(rc,rev(base));else base=connect(rev(base),rev(lc));
        }
        return{ldo,rdo};
    }
    vector<pair<int,int>> run(vector<GeomPoint> pts){
        p=move(pts);int n=p.size();if(n<2)return{};bool col=true;for(int i=2;i<n;i++)if(orient(0,1,i)){col=false;break;}
        if(col){vector<pair<int,int>> e;for(int i=1;i<n;i++)e.push_back({i-1,i});return e;}
        q.reserve(24*n);alive.reserve(6*n);build(0,n-1);vector<unsigned long long> key;key.reserve(3*n);
        for(int b=0;b<(int)alive.size();b++)if(alive[b]){int u=q[4*b].org,v=q[4*b+1].org;if(u>=0&&v>=0&&u!=v){if(u>v)swap(u,v);key.push_back(((unsigned long long)(unsigned)u<<32)|(unsigned)v);}}
        vector<unsigned long long> tmp(key.size());vector<int> count(1<<16);for(int shift=0;shift<64;shift+=16){fill(count.begin(),count.end(),0);for(auto x:key)count[(x>>shift)&65535]++;int sum=0;for(int&i:count){int c=i;i=sum;sum+=c;}for(auto x:key)tmp[count[(x>>shift)&65535]++]=x;key.swap(tmp);}
        vector<pair<int,int>> e;e.reserve(key.size());for(int i=0;i<(int)key.size();i++)if(!i||key[i]!=key[i-1])e.push_back({(int)(key[i]>>32),(int)(key[i]&0xffffffffu)});return e;
    }
};

static GeomResult buildGeometry(const vector<double>&X,const vector<double>&Y){
    int n=X.size();vector<GeomPoint> raw(n);for(int i=0;i<n;i++)raw[i]={(long long)X[i],(long long)Y[i],i};vector<int> ids(n);iota(ids.begin(),ids.end(),0);
    bool sortedX=true;for(int i=1;i<n;i++)if(raw[i].x<raw[i-1].x)sortedX=false;
    if(sortedX){for(int l=0;l<n;){int r=l+1;while(r<n&&raw[r].x==raw[l].x)r++;sort(ids.begin()+l,ids.begin()+r,[&](int a,int b){return raw[a].y<raw[b].y||(raw[a].y==raw[b].y&&a<b);});l=r;}}
    else sort(ids.begin(),ids.end(),[&](int a,int b){return tie(raw[a].x,raw[a].y,a)<tie(raw[b].x,raw[b].y,b);});
    vector<GeomPoint> uniq;vector<vector<int>> groups;for(int id:ids){if(uniq.empty()||uniq.back().x!=raw[id].x||uniq.back().y!=raw[id].y){uniq.push_back(raw[id]);groups.push_back({id});}else groups.back().push_back(id);}
    DelaunayCore dt;auto ue=dt.run(uniq);GeomResult out;out.edges.reserve(ue.size()+n);for(auto [a,b]:ue)out.edges.push_back({groups[a][0],groups[b][0]});for(auto&g:groups)for(int i=1;i<(int)g.size();i++)out.edges.push_back({g[0],g[i]});
    struct E{int a,b;long long d;};vector<E> es;es.reserve(out.edges.size());for(auto [a,b]:out.edges){long long dx=raw[a].x-raw[b].x,dy=raw[a].y-raw[b].y;es.push_back({a,b,dx*dx+dy*dy});}
    vector<E> tmp(es.size());vector<int> count(1<<16);for(int shift=0;shift<64;shift+=16){fill(count.begin(),count.end(),0);for(auto&e:es)count[((unsigned long long)e.d>>shift)&65535]++;int sum=0;for(int&i:count){int c=i;i=sum;sum+=c;}for(auto&e:es)tmp[count[((unsigned long long)e.d>>shift)&65535]++]=e;es.swap(tmp);}
    GeomDSU ds(n);vector<vector<int>> tr(n);for(auto&e:es)if(ds.join(e.a,e.b)){tr[e.a].push_back(e.b);tr[e.b].push_back(e.a);out.mstEdges.push_back({e.a,e.b});}
    for(int u=0;u<n;u++)sort(tr[u].begin(),tr[u].end(),[&](int a,int b){long long ax=raw[a].x-raw[u].x,ay=raw[a].y-raw[u].y,bx=raw[b].x-raw[u].x,by=raw[b].y-raw[u].y;int ah=ay>0||(ay==0&&ax>=0),bh=by>0||(by==0&&bx>=0);if(ah!=bh)return ah>bh;__int128 cr=(__int128)ax*by-(__int128)ay*bx;if(cr)return cr>0;__int128 ad=(__int128)ax*ax+(__int128)ay*ay,bd=(__int128)bx*bx+(__int128)by*by;return ad<bd||(ad==bd&&a<b);});
    auto contour=[&](int dir){vector<int> tour;tour.reserve(n);vector<pair<int,int>> st={{0,-1}};while(!st.empty()){auto [u,pa]=st.back();st.pop_back();tour.push_back(u);int d=tr[u].size(),pivot=-1;for(int i=0;i<d;i++)if(tr[u][i]==pa){pivot=i;break;}vector<int> kids;kids.reserve(d-(pa>=0));for(int z=1;z<=d;z++){int at=pivot<0?(dir>0?z-1:d-z):(pivot+dir*z+d)%d,v=tr[u][at];if(v!=pa)kids.push_back(v);}for(int i=(int)kids.size()-1;i>=0;i--)st.push_back({kids[i],u});}return tour;};
    auto a=contour(1),b=contour(-1);auto length=[&](const vector<int>&v){long double z=0;for(int i=0;i<n;i++){int j=(i+1)%n;long long dx=raw[v[i]].x-raw[v[j]].x,dy=raw[v[i]].y-raw[v[j]].y;z+=sqrt((long double)dx*dx+(long double)dy*dy);}return z;};out.mstTour=length(a)<=length(b)?move(a):move(b);return out;
}

struct DynamicKD {
    struct Node{int id=-1,left=-1,right=-1,parent=-1,alive=0;long long minx=0,maxx=0,miny=0,maxy=0;};
    const vector<double>&X,&Y;int n,root=-1;vector<Node> node;vector<int> where;vector<char> live;
    DynamicKD(const vector<double>&x,const vector<double>&y):X(x),Y(y),n(x.size()),where(n,-1),live(n){vector<int> ids;for(int i=1;i<n;i++)ids.push_back(i);node.reserve(ids.size());root=build(ids,0,ids.size(),0,-1);for(int i=1;i<n;i++)live[i]=1;}
    int build(vector<int>&ids,int l,int r,int dep,int par){if(l>=r)return-1;int m=(l+r)/2,axis=dep&1;nth_element(ids.begin()+l,ids.begin()+m,ids.begin()+r,[&](int a,int b){if(!axis)return X[a]<X[b]||(X[a]==X[b]&&(Y[a]<Y[b]||(Y[a]==Y[b]&&a<b)));return Y[a]<Y[b]||(Y[a]==Y[b]&&(X[a]<X[b]||(X[a]==X[b]&&a<b)));});int at=node.size(),id=ids[m];node.push_back(Node());auto&z=node[at];z.id=id;z.parent=par;z.alive=1;z.minx=z.maxx=(long long)X[id];z.miny=z.maxy=(long long)Y[id];where[id]=at;int a=build(ids,l,m,dep+1,at),b=build(ids,m+1,r,dep+1,at);node[at].left=a;node[at].right=b;for(int c:{a,b})if(c>=0){node[at].alive+=node[c].alive;node[at].minx=min(node[at].minx,node[c].minx);node[at].maxx=max(node[at].maxx,node[c].maxx);node[at].miny=min(node[at].miny,node[c].miny);node[at].maxy=max(node[at].maxy,node[c].maxy);}return at;}
    long double box(int at,int cur)const{long double dx=0,dy=0;auto&z=node[at];if(X[cur]<z.minx)dx=z.minx-X[cur];else if(X[cur]>z.maxx)dx=X[cur]-z.maxx;if(Y[cur]<z.miny)dy=z.miny-Y[cur];else if(Y[cur]>z.maxy)dy=Y[cur]-z.maxy;return dx*dx+dy*dy;}
    void query(int at,int cur,int&best,long double&bd)const{if(at<0||!node[at].alive)return;long double lo=box(at,cur);if(best>=0&&lo>=bd)return;auto&z=node[at];if(live[z.id]){long double dx=X[cur]-X[z.id],dy=Y[cur]-Y[z.id],d=dx*dx+dy*dy;if(best<0||d<bd){best=z.id;bd=d;}}int a=z.left,b=z.right;long double da=a<0||!node[a].alive?numeric_limits<long double>::infinity():box(a,cur),db=b<0||!node[b].alive?numeric_limits<long double>::infinity():box(b,cur);if(db<da){swap(a,b);swap(da,db);}if(best<0||da<bd)query(a,cur,best,bd);if(best<0||db<bd)query(b,cur,best,bd);}
    void erase(int id){if(id<=0||!live[id])return;live[id]=0;for(int p=where[id];p>=0;p=node[p].parent)--node[p].alive;}
    vector<int> complete(vector<int> out){out.reserve(n);for(int v:out)erase(v);int cur=out.empty()?0:out.back();if(out.empty())out.push_back(0);while((int)out.size()<n){int best=-1;long double bd=numeric_limits<long double>::infinity();query(root,cur,best,bd);if(best<0)break;out.push_back(best);erase(best);cur=best;}return out;}
    vector<int> route(){return complete({0});}
};

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

    vector<char> pr((size_t)N,0);
    { vector<char> comp((size_t)N,0); for(long long i=2;i<N;i++) if(!comp[i]){ pr[i]=1; for(long long q=i*i;q<N;q+=i) comp[q]=1; } }
    bool legacyMid = false;
    GeomResult geom;
    if(GEOM_MODE && !legacyMid) geom=buildGeometry(X,Y);
#ifdef GEOM_RESET_TIMER
    if(GEOM_MODE) T0=chrono::steady_clock::now();
#endif
    if(N>100000) TL_MS -= 20.0;
    double RESERVE = N>150000?220.0:(N>50000?90.0:(N>5000?50.0:40.0));
    TL_MS -= RESERVE; // reserve tail for endgame touch-up

    double minx=X[0],maxx=X[0],miny=Y[0],maxy=Y[0];
    for(int i=1;i<N;i++){ minx=min(minx,X[i]);maxx=max(maxx,X[i]);miny=min(miny,Y[i]);maxy=max(maxy,Y[i]); }
    double w=max(1.0,maxx-minx), h=max(1.0,maxy-miny);
    int G=max(1,(int)floor(sqrt((double)N/2.0)));
    int GX=G, GY=G;
    if(N<=5000){ int cells=max(1,N/2); GX=max(1,min(cells,(int)llround(sqrt((double)cells*w/h)))); GY=max(1,cells/GX); }
    double cw=w/GX, ch=h/GY;
    auto gx=[&](double x){ int c=(int)((x-minx)/cw); return c<0?0:(c>=GX?GX-1:c); };
    auto gy=[&](double y){ int c=(int)((y-miny)/ch); return c<0?0:(c>=GY?GY-1:c); };
    vector<int> cellOf(N), cnt(GX*GY+1,0);
    for(int i=0;i<N;i++){ int c=gx(X[i])*GY+gy(Y[i]); cellOf[i]=c; cnt[c+1]++; }
    for(int i=0;i<GX*GY;i++) cnt[i+1]+=cnt[i];
    vector<int> bucket(N); { vector<int> tmp=cnt; for(int i=0;i<N;i++) bucket[tmp[cellOf[i]]++]=i; }

    int K;
    if(N>50000)                      K=min(N-1,6);
    else if(N>=16000 && N<=36000)    K=min(N-1,40);
    else if(N>5000)                  K=min(N-1,24);
    else                             K=min(N-1,10);
    if(const char* e=getenv("K_FORCE")) K=min(N-1, atoi(e));
#ifdef GEOM_NO_GRID
    K=0;
#endif
    vector<int> nbr((size_t)N*K,-1);
    bool gridValid=K>0;long long gridOps=0,gridCap=legacyMid?LLONG_MAX:1LL*N*max(512,8*K);
    if(K>0){
        vector<pair<double,int>> cand; cand.reserve(128);
        for(int i=0;i<N&&gridValid;i++){
            int cx=gx(X[i]),cy=gy(Y[i]); cand.clear();
            int ring=0, extra=1;
            while(true){
                int x0=max(0,cx-ring),x1=min(GX-1,cx+ring),y0=max(0,cy-ring),y1=min(GY-1,cy+ring);
                for(int xx=x0;xx<=x1;xx++) for(int yy=y0;yy<=y1;yy++){
                    if(legacyMid ? (ring>0 && xx>x0 && xx<x1 && yy>y0 && yy<y1) : (max(abs(xx-cx),abs(yy-cy))!=ring))continue;
                    int c=xx*GY+yy;
                    for(int b=cnt[c];b<cnt[c+1];b++){if(++gridOps>gridCap){gridValid=false;break;}int j=bucket[b];if(j!=i)cand.push_back({dist(i,j),j});}
                    if(!gridValid)break;
                }
                if(!gridValid)break;
                if((int)cand.size()>=K){ if(extra--<=0) break; }
                if(x0==0&&y0==0&&x1==GX-1&&y1==GY-1) break;
                ring++;
            }
            int kk=min((int)cand.size(),K);
            partial_sort(cand.begin(),cand.begin()+kk,cand.end());
            for(int t=0;t<kk;t++) nbr[(size_t)i*K+t]=cand[t].second;
        }
    }
#ifdef GRID_MST_SEED
    vector<int> gridMstTour;
    {
        struct E{int a,b;long long d;};vector<E> es;es.reserve((size_t)N*(K+1));
        auto d2=[&](int a,int b){long long dx=(long long)X[a]-(long long)X[b],dy=(long long)Y[a]-(long long)Y[b];return dx*dx+dy*dy;};
        for(int i=0;i<N;i++)for(int t=0;t<K;t++){int j=nbr[(size_t)i*K+t];if(j>=0&&i<j)es.push_back({i,j,d2(i,j)});}
        for(int i=1;i<N;i++)es.push_back({i-1,i,d2(i-1,i)});
        sort(es.begin(),es.end(),[](auto&a,auto&b){return a.d<b.d;});GeomDSU ds(N);vector<vector<int>> tr(N);int usedEdges=0;
        for(auto&e:es)if(ds.join(e.a,e.b)){tr[e.a].push_back(e.b);tr[e.b].push_back(e.a);if(++usedEdges==N-1)break;}
        if(usedEdges==N-1){
            for(int u=0;u<N;u++)sort(tr[u].begin(),tr[u].end(),[&](int a,int b){
                long long ax=(long long)X[a]-(long long)X[u],ay=(long long)Y[a]-(long long)Y[u],bx=(long long)X[b]-(long long)X[u],by=(long long)Y[b]-(long long)Y[u];
                int ah=ay>0||(ay==0&&ax>=0),bh=by>0||(by==0&&bx>=0);if(ah!=bh)return ah>bh;__int128 cr=(__int128)ax*by-(__int128)ay*bx;if(cr)return cr>0;return d2(u,a)<d2(u,b);
            });
            auto contour=[&](int dir){vector<int> v;v.reserve(N);vector<pair<int,int>> st={{0,-1}};while(!st.empty()){auto[u,pa]=st.back();st.pop_back();v.push_back(u);int d=tr[u].size(),pivot=-1;for(int i=0;i<d;i++)if(tr[u][i]==pa)pivot=i;vector<int> kids;for(int z=1;z<=d;z++){int at=pivot<0?(dir>0?z-1:d-z):(pivot+dir*z+d)%d,w=tr[u][at];if(w!=pa)kids.push_back(w);}for(int i=(int)kids.size()-1;i>=0;i--)st.push_back({kids[i],u});}return v;};
            auto a=contour(1),b=contour(-1);auto len=[&](const vector<int>&v){double z=0;for(int i=0;i<N;i++)z+=dist(v[i],v[(i+1)%N]);return z;};gridMstTour=len(a)<=len(b)?move(a):move(b);
        }
    }
#endif
    vector<vector<int>> candidates(N);
    for(int i=0;i<N;i++) for(int t=0;t<K;t++) if(nbr[(size_t)i*K+t]>=0) candidates[i].push_back(nbr[(size_t)i*K+t]);
    bool augmented=(!legacyMid && (GEOM_MODE&1));
#ifdef GEOM_ALPHA_ORDERED
    const bool preserveCandidatePrefix=true;
#else
    const bool preserveCandidatePrefix=false;
#endif
    if(!legacyMid && (GEOM_MODE&1)){
        for(auto [a,b]:geom.edges){candidates[a].push_back(b);candidates[b].push_back(a);}
    }
#ifdef GEOM_ALPHA
    if(!legacyMid && (int)geom.mstEdges.size()==N-1){
        augmented=true;int LOG=1;while((1<<LOG)<=N)LOG++;
        vector<vector<pair<int,long long>>> tr(N);auto d2=[&](int a,int b){long long dx=(long long)X[a]-(long long)X[b],dy=(long long)Y[a]-(long long)Y[b];return dx*dx+dy*dy;};
        for(auto [a,b]:geom.mstEdges){long long edgeSq=d2(a,b);tr[a].push_back({b,edgeSq});tr[b].push_back({a,edgeSq});}
        vector<int> depth(N),parent(N,-1),stack={0};vector<long long> parentW(N);
        while(!stack.empty()){int u=stack.back();stack.pop_back();for(auto [v,edgeSq]:tr[u])if(v!=parent[u]){parent[v]=u;parentW[v]=edgeSq;depth[v]=depth[u]+1;stack.push_back(v);}}
        vector<int> up((size_t)LOG*N);vector<long long> mx((size_t)LOG*N);
        for(int i=0;i<N;i++){up[i]=parent[i]<0?i:parent[i];mx[i]=parentW[i];}
        for(int l=1;l<LOG;l++)for(int i=0;i<N;i++){int m=up[(size_t)(l-1)*N+i];up[(size_t)l*N+i]=up[(size_t)(l-1)*N+m];mx[(size_t)l*N+i]=max(mx[(size_t)(l-1)*N+i],mx[(size_t)(l-1)*N+m]);}
        auto bottleneck=[&](int a,int b){long long z=0;if(depth[a]<depth[b])swap(a,b);int d=depth[a]-depth[b];for(int l=0;l<LOG;l++)if(d>>l&1){z=max(z,mx[(size_t)l*N+a]);a=up[(size_t)l*N+a];}if(a==b)return z;for(int l=LOG-1;l>=0;l--)if(up[(size_t)l*N+a]!=up[(size_t)l*N+b]){z=max(z,max(mx[(size_t)l*N+a],mx[(size_t)l*N+b]));a=up[(size_t)l*N+a];b=up[(size_t)l*N+b];}return max(z,max(mx[a],mx[b]));};
        vector<vector<pair<long double,int>>> alpha(N);
        for(auto [a,b]:geom.edges){long long edgeSq=d2(a,b),pathMax=bottleneck(a,b);long double score=sqrt((long double)edgeSq)-sqrt((long double)pathMax);alpha[a].push_back({score,b});alpha[b].push_back({score,a});}
        for(int i=0;i<N;i++){auto &v=alpha[i];sort(v.begin(),v.end(),[&](auto&a,auto&b){return a.first<b.first||(a.first==b.first&&d2(i,a.second)<d2(i,b.second));});int take=min(8,(int)v.size());for(int j=0;j<take;j++)candidates[i].push_back(v[j].second);}
    }
#endif
    if(augmented&&!preserveCandidatePrefix){
        if(!preserveCandidatePrefix)for(int i=0;i<N;i++){
            auto &v=candidates[i];sort(v.begin(),v.end());v.erase(unique(v.begin(),v.end()),v.end());
            sort(v.begin(),v.end(),[&](int a,int b){double da=dist(i,a),db=dist(i,b);return da<db||(da==db&&a<b);});
        }
    }

    vector<int> order(N), pos(N); vector<char> used(N,0);
    bool haveOrder=false;
#ifdef GEOM_MST_ALWAYS
    if(!legacyMid && (GEOM_MODE&2) && (int)geom.mstTour.size()==N){order=geom.mstTour;haveOrder=true;}
#else
    if(gridValid){
        vector<int> trial(N);vector<char> trialUsed(N);int cur=0,filled=1;trialUsed[0]=1;trial[0]=0;bool ok=true;long long routeOps=0,routeCap=legacyMid?LLONG_MAX:1LL*N*max(512,8*K);
        for(int step=1;step<N&&ok;step++){
            int best=-1;double bd=1e300;for(int t=0;t<K;t++){int j=nbr[(size_t)cur*K+t];if(j>=0&&!trialUsed[j]){best=j;break;}}
            if(best<0){
                int cx=gx(X[cur]),cy=gy(Y[cur]);
                for(int ring=0;ring<2*max(GX,GY)&&best<0&&ok;ring++){
                    int x0=max(0,cx-ring),x1=min(GX-1,cx+ring),y0=max(0,cy-ring),y1=min(GY-1,cy+ring);
                    for(int xx=x0;xx<=x1&&ok;xx++)for(int yy=y0;yy<=y1;yy++){
                        if(legacyMid ? (ring>0 && xx>x0 && xx<x1 && yy>y0 && yy<y1) : (max(abs(xx-cx),abs(yy-cy))!=ring))continue;
                        int c=xx*GY+yy;
                        for(int b=cnt[c];b<cnt[c+1];b++){if(++routeOps>routeCap){ok=false;break;}int j=bucket[b];if(!trialUsed[j]){double d=dist(cur,j);if(d<bd){bd=d;best=j;}}}
                    }
                    if(best>=0&&ok){int r2=ring+1,a0=max(0,cx-r2),a1=min(GX-1,cx+r2),b0=max(0,cy-r2),b1=min(GY-1,cy+r2);for(int xx=a0;xx<=a1&&ok;xx++)for(int yy=b0;yy<=b1;yy++){if(legacyMid ? (xx>a0&&xx<a1&&yy>b0&&yy<b1) : (max(abs(xx-cx),abs(yy-cy))!=r2))continue;int c=xx*GY+yy;for(int b=cnt[c];b<cnt[c+1];b++){if(++routeOps>routeCap){ok=false;break;}int j=bucket[b];if(!trialUsed[j]){double d=dist(cur,j);if(d<bd){bd=d;best=j;}}}}break;}
                }
            }
            if(best<0)ok=false;else{trialUsed[best]=1;trial[step]=best;cur=best;filled=step+1;}
        }
        if(!ok&&filled>1){trial.resize(filled);trial=DynamicKD(X,Y).complete(move(trial));ok=(int)trial.size()==N;}
        if(ok){order.swap(trial);haveOrder=true;}
    }
#endif
#ifdef GEOM_KD_SEED
    if(!legacyMid && el_ms() < TL_MS - 300.0){
        vector<int> kd=DynamicKD(X,Y).route();auto cycleLen=[&](const vector<int>&v){double z=0;for(int i=0;i<N;i++)z+=dist(v[i],v[(i+1)%N]);return z;};
        if(getenv("GEOM_DEBUG"))cerr<<"grid="<<(haveOrder?cycleLen(order):-1)<<" kd="<<((int)kd.size()==N?cycleLen(kd):-1)<<"\n";
        if((int)kd.size()==N&&(!haveOrder||cycleLen(kd)<cycleLen(order))){order.swap(kd);haveOrder=true;}
    }
#endif
    if(!haveOrder){order.resize(N);iota(order.begin(),order.end(),0);}
    for(int i=0;i<N;i++)pos[order[i]]=i;
#ifndef GEOM_MST_ALWAYS
    if(!legacyMid && (GEOM_MODE&2) && (int)geom.mstTour.size()==N){
        auto cycleLen=[&](const vector<int>&v){double z=0;for(int i=0;i<N;i++)z+=dist(v[i],v[(i+1)%N]);return z;};
        double nnLen=cycleLen(order),mstLen=cycleLen(geom.mstTour);
        if(getenv("GEOM_DEBUG"))cerr<<"nn="<<nnLen<<" mst="<<mstLen<<" ratio="<<mstLen/nnLen<<"\n";
        if(mstLen<nnLen){order=geom.mstTour;for(int i=0;i<N;i++)pos[order[i]]=i;}
    }
#endif
#ifdef GRID_MST_SEED
    if((int)gridMstTour.size()==N){
        auto cycleLen=[&](const vector<int>&v){double z=0;for(int i=0;i<N;i++)z+=dist(v[i],v[(i+1)%N]);return z;};
        if(cycleLen(gridMstTour)<cycleLen(order)){order=gridMstTour;for(int i=0;i<N;i++)pos[order[i]]=i;}
    }
#endif
    if(augmented){
        for(int i=0;i<N;i++){
            int a=order[i],b=order[(i+1)%N];candidates[a].push_back(b);candidates[b].push_back(a);
        }
        for(int i=0;i<N;i++){
            auto &v=candidates[i];sort(v.begin(),v.end());v.erase(unique(v.begin(),v.end()),v.end());
            sort(v.begin(),v.end(),[&](int a,int b){double da=dist(i,a),db=dist(i,b);return da<db||(da==db&&a<b);});
        }
    }

    auto nextIdx=[&](int i){ return i+1<N?i+1:0; };
    auto prevIdx=[&](int i){ return i>0?i-1:N-1; };
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
            for(int dir=0; dir<2 && !improved; dir++){
                int p1=pos[c1];
                int p2=(dir==0)?nextIdx(p1):prevIdx(p1);
                int c2=order[p2];
                double d12=dist(c1,c2);
                for(int c3:candidates[c1]){
                    double d13=dist(c1,c3);
                    if(d13>=d12){if(preserveCandidatePrefix && augmented)continue;break;} // sorted prefix, then alpha extras
                    int p3=pos[c3];
                    int p4=(dir==0)?nextIdx(p3):prevIdx(p3);
                    int c4=order[p4];
                    if(c4==c1||c3==c2) continue;
                    double before=d12+dist(c3,c4);
                    double after =d13+dist(c2,c4);
                    if(after+1e-7<before){
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
                for(int side=0; side<2 && !moved; side++){
                    int anchorCity = side==0? s0 : segEnd;
                    for(int c:candidates[anchorCity]){
                        int pc=pos[c];
                        bool inside=false; { int p=is; for(int u=0;u<L;u++){ if(p==pc){inside=true;break;} p=nextIdx(p);} }
                        if(inside || c==cprev) continue;
                        int pcn=nextIdx(pc); int cn=order[pcn];
                        if(cn==s0) continue;             // that's the current place
                        double added = dist(c,s0)+dist(segEnd,cn)-dist(c,cn);
                        double addedRev = dist(c,segEnd)+dist(s0,cn)-dist(c,cn);
                        bool rev = addedRev+1e-9 < added;
                        double add = rev? addedRev: added;
                        if(add+1e-7 < removed){
                            int seg[3]; { int p=is; for(int u=0;u<L;u++){ seg[u]=order[p]; p=nextIdx(p);} }
                            if(rev){ for(int a=0,b=L-1;a<b;a++,b--) swap(seg[a],seg[b]); }
                            static vector<int> tmp; tmp.clear(); tmp.reserve(N);
                            bool isMember[3]; (void)isMember;
                            auto inRun=[&](int city){ for(int u=0;u<L;u++) if(seg[u]==city||order[(is+0)]==city){} return false; };
                            (void)inRun;
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
                for(int c3:candidates[c1]){
                    double d13=dist(c1,c3);
                    if(d13>=d12){if(preserveCandidatePrefix && augmented)continue;break;} // gain criterion g1=d12-d13>0
                    int p3=pos[c3];
                    int p4=(dir==0)?nextIdx(p3):prevIdx(p3);
                    int c4=order[p4];
                    if(c4==c1||c3==c2) continue;
                    double g1full = (d12+dist(c3,c4)) - (d13+dist(c2,c4)); // first 2-opt gain (may be <=0)
                    int A=(dir==0)?p1:p4, B=(dir==0)?p3:p2;
                    applyMove(A,B);                        // tentative first move
                    double bestG2=1e-7; int bcc=-1,bdd=-1,bcf=-1;
                    for(int side=0; side<2; side++){
                        int cc=(side==0)?c2:c4;
                        int pcc=pos[cc];
                        for(int dd=0; dd<2; dd++){
                            int pe=(dd==0)?nextIdx(pcc):prevIdx(pcc);
                            int ce=order[pe];
                            double dce=dist(cc,ce);
                            for(int cf:candidates[cc]){
                                double dcf=dist(cc,cf);
                                if(dcf>=dce){if(preserveCandidatePrefix && augmented)continue;break;}
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
#ifdef GEOM_SPLAY_LK
    auto fastLkPass=[&](){
        SplayTour tour(order);fill(dontlook2.begin(),dontlook2.end(),0);bool any=false;
        auto cycPrev=[&](int p){return p>0?p-1:N-1;};
        for(int qi=0;qi<N;qi++){
            if((qi&255)==0&&el_ms()>TL_MS)break;
            int c1=q[qi];if(dontlook2[c1])continue;bool improved=false;
            for(int dir=0;dir<2&&!improved;dir++){
                int p1=tour.rank(c1),p2=dir==0?(p1+1)%N:cycPrev(p1),c2=tour.cityAt(p2);double d12=dist(c1,c2);
                for(int c3:candidates[c1]){
                    double d13=dist(c1,c3);if(d13>=d12)break;
                    int p3=tour.rank(c3),p4=dir==0?(p3+1)%N:cycPrev(p3),c4=tour.cityAt(p4);
                    if(c4==c1||c3==c2)continue;
                    double g1=(d12+dist(c3,c4))-(d13+dist(c2,c4));
                    int A=dir==0?p1:p4,B=dir==0?p3:p2;tour.twoOpt(A,B);
                    double bestG=1e-7;int bc=-1,bd=-1,bf=-1;
                    for(int side=0;side<2;side++){
                        int cc=side==0?c2:c4,pcc=tour.rank(cc);
                        for(int dd=0;dd<2;dd++){
                            int pe=dd==0?(pcc+1)%N:cycPrev(pcc),ce=tour.cityAt(pe);double dce=dist(cc,ce);
                            for(int cf:candidates[cc]){
                                double dcf=dist(cc,cf);if(dcf>=dce)break;
                                int pf=tour.rank(cf),pg=dd==0?(pf+1)%N:cycPrev(pf),cg=tour.cityAt(pg);
                                if(cg==cc||cf==ce)continue;
                                double g=(dce+dist(cf,cg))-(dcf+dist(ce,cg));
                                if(g>bestG){bestG=g;bc=cc;bd=dd;bf=cf;}
                            }
                        }
                    }
                    if(bc>=0&&g1+bestG>1e-7){
                        int pc=tour.rank(bc),pf=tour.rank(bf),A2=bd==0?pc:cycPrev(pf),B2=bd==0?pf:cycPrev(pc);tour.twoOpt(A2,B2);
                        dontlook2[c1]=dontlook2[c2]=dontlook2[c3]=dontlook2[c4]=dontlook2[bc]=0;improved=any=true;
                    }else tour.twoOpt(A,B);
                    if(improved)break;
                }
            }
            if(!improved)dontlook2[c1]=1;
        }
        order=tour.materialize();for(int i=0;i<N;i++)pos[order[i]]=i;return any;
    };
#endif

    auto localSearch=[&](){
        while(el_ms()<TL_MS){
            bool a=false; while(el_ms()<TL_MS){ if(!twoOptPass()){break;} a=true; }
            bool b=false;
            if(N<=50000){ fill(dontlook.begin(),dontlook.end(),0); b=orOptPass(); if(b) fill(dontlook.begin(),dontlook.end(),0); }
            bool c=false;
#ifndef GEOM_SPLAY_LK
            if(N>=64){ fill(dontlook2.begin(),dontlook2.end(),0); c=lkPass(); if(c) fill(dontlook.begin(),dontlook.end(),0); }
#else
            if(N>=64){c=fastLkPass();if(c)fill(dontlook.begin(),dontlook.end(),0);}
#endif
            if(!a && !b && !c) break;
        }
    };
    localSearch();

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

    vector<int> best=order; int bestDir=0;
    double bestLen=evalBest(bestDir);
    if(N>=8){
        uint64_t rng=0x9e3779b97f4a7c15ULL ^ (uint64_t)N*2654435761ULL;
        auto rnd=[&](){ rng^=rng<<7; rng^=rng>>9; return rng; };
        int win = N<=60?N:max(24,min(N,N<=2000?N/4:(N<=50000?400:N)));
        int stale=0;
        while(el_ms()<TL_MS){
            int span = min(win, N-1);
            if(span<4) break;
            int base = win>=N? 0 : (int)(rnd()%N);
            auto W=[&](int off){ return (base+off)%N; };
            int a=1+(int)(rnd()%(span-3)), b=1+(int)(rnd()%(span-3)), c=1+(int)(rnd()%(span-3));
            int lo=min({a,b,c}), hi=max({a,b,c}), mid=a+b+c-lo-hi;
            if(lo==mid||mid==hi) continue;
            static vector<int> wbuf; wbuf.clear(); wbuf.reserve(span);
            for(int i=0;i<span;i++) wbuf.push_back(order[W(i)]);
            static vector<int> nt; nt.clear(); nt.reserve(span);
            for(int i=0;i<lo;i++) nt.push_back(wbuf[i]);
            for(int i=mid;i<hi;i++) nt.push_back(wbuf[i]);
            for(int i=lo;i<mid;i++) nt.push_back(wbuf[i]);
            for(int i=hi;i<span;i++) nt.push_back(wbuf[i]);
            for(int i=0;i<span;i++){ int city=nt[i]; order[W(i)]=city; pos[city]=W(i); }
            for(int i=0;i<span;i++) dontlook[wbuf[i]]=0;
            { int before=(base==0?N-1:base-1), after=W(span-1)+1<N?W(span-1)+1:0;
              dontlook[order[before]]=0; dontlook[order[(after)%N]]=0; }
            while(el_ms()<TL_MS){ if(!twoOptPass()) break; }
            if(N<=50000){ orOptPass(); while(el_ms()<TL_MS){ if(!twoOptPass()) break; } }
            int d2=0; double L=evalBest(d2);
            if(L<bestLen-1e-6){ bestLen=L; best=order; bestDir=d2; stale=0; }
            else { order=best; for(int i=0;i<N;i++) pos[order[i]]=i; if(++stale>200 && N<64) break; } // revert
        }
        order=best; for(int i=0;i<N;i++) pos[order[i]]=i;
    }

    TL_MS += RESERVE;
    vector<int> seq(N);
    { int z=pos[0]; for(int i=0;i<N;i++) seq[i]= bestDir==0? order[(z+i)%N] : order[(z-i+N)%N]; }

#ifdef GEOM_GROUP_DUPLICATES
    if(!legacyMid){
        vector<int> ids(N);iota(ids.begin(),ids.end(),0);sort(ids.begin(),ids.end(),[&](int a,int b){return X[a]<X[b]||(X[a]==X[b]&&(Y[a]<Y[b]||(Y[a]==Y[b]&&a<b)));});
        vector<int> groupOf(N),first;vector<vector<int>> members;
        for(int id:ids){if(first.empty()||X[id]!=X[first.back()]||Y[id]!=Y[first.back()]){first.push_back(id);members.push_back({});}groupOf[id]=(int)members.size()-1;members.back().push_back(id);}
        if((int)members.size()<N){
            vector<char> seen(members.size());vector<int> grouped;grouped.reserve(N);
            for(int city:seq){int g=groupOf[city];if(seen[g])continue;seen[g]=1;auto block=members[g];
                if(g==groupOf[0]){auto it=find(block.begin(),block.end(),0);rotate(block.begin(),it,block.end());}
                int endPos=(int)grouped.size()+(int)block.size()-1;
                if((endPos+1)%10==0){auto it=find_if(block.begin(),block.end(),[&](int v){return pr[v];});if(it!=block.end())iter_swap(it,block.end()-1);}
                grouped.insert(grouped.end(),block.begin(),block.end());
            }
            auto cost=[&](const vector<int>&v){double z=0;for(int t=1;t<=N;t++){int a=v[t-1],b=t<N?v[t]:0;double d=dist(a,b);if(t%10==0&&!pr[a])d*=1.1;z+=d;}return z;};
            if((int)grouped.size()==N&&grouped[0]==0&&cost(grouped)<cost(seq))seq.swap(grouped);
        }
    }
#endif

    if(N>=10){
        auto sAt=[&](int p)->int{ return p<N? seq[p]:0; };
        auto stepCost=[&](int t)->double{ int a=seq[t-1], b=sAt(t); double d=dist(a,b); if(t%10==0&&!pr[a]) d*=1.1; return d; };
        int window = N<=1200? N : (N<=5000?120:(N<=20000?140:80));
        for(int rep=0;rep<4 && el_ms()<TL_MS;rep++){
            bool changed=false;
            for(int p=9;p<N;p+=10){
                if(el_ms()>TL_MS) break;
                if(pr[seq[p]]) continue;
                int lo=max(1,p-window), hi=min(N-1,p+window);
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
                if(bj>=0){ swap(seq[p],seq[bj]); changed=true; }
            }
            if(!changed) break;
        }
    }
    {
        auto exact=[&](const vector<int>&v){double z=0;for(int t=1;t<=N;t++){int a=v[t-1],b=t<N?v[t]:0;double d=dist(a,b);if(t%10==0&&!pr[a])d*=1.1;z+=d;}return z;};
        vector<int> rev(N);rev[0]=0;for(int i=1;i<N;i++)rev[i]=seq[N-i];
        if(exact(rev)<exact(seq))seq.swap(rev);
    }

    string out; out.reserve((size_t)N*7+16);
    out+=to_string(N+1); out+='\n';
    for(int i=0;i<N;i++){ out+=to_string(seq[i]); out+='\n'; }
    out+="0\n";
    fwrite(out.data(),1,out.size(),stdout);
    fflush(stdout);
    _Exit(0);
}
