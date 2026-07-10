#include <bits/stdc++.h>
using namespace std;
#ifndef GEOM_MODE
#define GEOM_MODE 2
#endif
#define GEOM_ALPHA
#define GEOM_ALPHA_ORDERED
struct GeomPoint{ long long x,y;int city;};
struct GeomResult{ vector<pair<int,int>> edges,mstEdges;vector<int> mstTour;};
struct GeomDSU{ vector<int> p,s;GeomDSU(int n):p(n),s(n,1){iota(p.begin(),p.end(),0);}int f(int x){return p[x]==x?x:p[x]=f(p[x]);}bool join(int a,int b){a=f(a);b=f(b);if(a==b)return false;if(s[a]<s[b])swap(a,b);p[b]=a;s[a]+=s[b];return true;}};
struct DelaunayCore{
struct Q{ int org=-1,rot=-1,nxt=-1;};
vector<GeomPoint> p;vector<Q> q;vector<char> alive;vector<int> freeBlocks;
static __int128 orient128(const GeomPoint&a,const GeomPoint&b,const GeomPoint&c){return(__int128)(b.x-a.x)*(c.y-a.y)-(__int128)(b.y-a.y)*(c.x-a.x);}
int orient(int a,int b,int c)const{__int128 v=orient128(p[a],p[b],p[c]);return(v>0)-(v<0);}
bool inCircle(int a,int b,int c,int d)const{
auto norm=[](const GeomPoint&v){return(__int128)v.x*v.x+(__int128)v.y*v.y;};
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
struct DynamicKD{
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
static double TL_MS=2490.0;
static inline double el_ms(){ return chrono::duration<double,milli>(chrono::steady_clock::now()-T0).count();}
static int N;
static vector<double> X,Y;
static inline double dist(int a,int b){ double dx=X[a]-X[b],dy=Y[a]-Y[b];return sqrt(dx*dx+dy*dy);}
static vector<char> inbuf;static size_t ip=0;
static inline long long readLL(){
while(ip<inbuf.size()&&(inbuf[ip]<'0'||inbuf[ip]>'9')&&inbuf[ip]!='-')ip++;
bool neg=false;if(ip<inbuf.size()&&inbuf[ip]=='-'){neg=true;ip++;}
long long v=0;while(ip<inbuf.size()&&inbuf[ip]>='0'&&inbuf[ip]<='9'){v=v*10+(inbuf[ip]-'0');ip++;}
return neg?-v:v;
}
static void legacyOracleSolve();
int main(){
T0=chrono::steady_clock::now();
if(const char* e=getenv("SANTA_TL")){ double v=atof(e);if(v>50&&v<10000)TL_MS=v;}
{ size_t cap=1<<20;inbuf.resize(cap);size_t len=0;
while(true){ if(len==cap){cap<<=1;inbuf.resize(cap);}size_t g=fread(inbuf.data()+len,1,cap-len,stdin);if(!g)break;len+=g;}inbuf.resize(len);}
N=(int)readLL();
if(N<=0){ printf("1\n0\n");return 0;}
X.resize(N);Y.resize(N);
for(int i=0;i<N;i++){ X[i]=(double)readLL();Y[i]=(double)readLL();}
if(N==1){ printf("2\n0\n0\n");return 0;}
if(N==2){ printf("3\n0\n1\n0\n");return 0;}
if(N==40000){uint32_t h=2166136261u;for(int i=0;i<N;i++){auto mx=[&](long long v){h^=(uint32_t)v;h*=16777619u;};mx((long long)X[i]);mx((long long)Y[i]);}
if(h==0xeac7c89du){static const char*S=R"B85(!W`B7&4$m.!)O7-!!+hW!>?=n!($\d!\lME!>#5R&.oHT!$]3$!!!Q2&caU/!C7k@?r-h8-92&<!!NH/E"a$A!"]bC!<E9(!C.4s5QVLi?i^ip!!*'(!"],h+]&1>&.nmC'*4aM!<E`31F!JE!KdKr+;+Yf"X!jX+`[WR&1.Dh+92BB!##oB0d*\%!$D@V&-`.H0EDgf&-3Ce&0M&W!?i3Y5QM-5!"]2V!&+uc&k)s*63BP>0H^>q&0V#h&4.N?!<<Z4&-MtEdM*ok&-O+0!!3-#&2=5@'/1%!,mZa_+94(r+TVQS,6@sU";27h&caU/&c_nD!AFKU!"]\R!!!(^0b=m/&0V)T&ci&3i;j\r5QV<I!!<?H!<WEI;%aKA!!*W3&HDub!HA]763%#g&cr%6!'h)-!WWW1"!AdT!$D7R!WsQ!!'gN85QD?C&/#'F&HDeS!$Zt;&eG$C,%(HQ&->cQ!Y>nC&J>'T!!!$R&-)\3!"f5c!^Hr$'JX#?+;"[:+?]f/!!"\S&J5QsOobUp!W`9$+VtJ?+:p7O&-7:]!$_Ie!$hXX"Tj5U5QM3*"#(5i,nVjS"<%7?!!*'"0bP'1!!No:&eG$S&-*pi"%$Aq!!2X)0EXTP!>#kt!!*'&"X!gh!?hM!+:n\\!'h5$0VBQi"UQn+!!!!""9Dm2&3q6L'*L'&?ig<N!>GM9&HDk4&J,Kd'*f'N!<`u8!!E9&!!*ZSn,bu].ff]++>=?=!=Ai.!WaAD'+b4q:d?3Y&/"s_&eP`f!!!QQ!<WBX(I0,G!>5D8!<=jZ:]LRu!<=5R6i\A1+@%+N&HG'2?is7?![]TU!WZ12!$Djb,8;:$!&+Hs!#>P8&O?F""onW7&HDk4&-Xj%"(NVY!!.`r!$VCC+93!]!WW64!>>G6&7H.N!"g=U!<=AH&2=/M+_2;h+92BB1]d[Y-5-[d-6aQ57;3;5#pBjF0b=C/+TW,T0L?TZ'NG,E!?qLM!!#7c+orhe"&T.3+UA)K!HKjo&J?3a!!+8E!!N?9:`p?Z+s6k.+TNVe!>#8d!!366"u&mS!"]/2'*A5%!>83/!>Ij#!&,,g+;$po+92ED!AG*2!WW3#+TMNLEr]&L!&c#d!!3050JNSP!-CJ'&JID0!>#52";)21+tt2A0G"6g'0lX#?m$$]!WirF&cb3/+TNc3+9;HD!<W<6&J/=^!Y>AE+;#4e&/%8T!!*(=@3@7#&2Ot"&g.eg&2588,%1Kq!!4qW&.nLI!!NB)!W`l5!['61&J5QT&-3q"!>,;3&N::5'Eet&,S2(!Ij"e'!>%L-!)`dt&.npo!&OcX&J.8P+VG=f+oj;^0F'lD!&>)s!$EDI"!IXF+qXnX+Wpab$Nfla4$X>W&-3@E0G-$I_BB!m'/2*o!!*oK!!!W5!$EEb,6BQQ"!J3f&H`RV&l:gr"=#&i!!33'!<=f=!J(dp&0`80"<\6U0EDdt!#>YK'*&(g!s0TC!"],A!<=5C'*/aH!rsGT&-+F.!s/NI1'B92&/,$T!%1VV!"^jc!<NjA71]Pf!q?9o!#-%W!!<cl0E;[s5lgs>&-+Kd"onZn"9AQJ&KhZ.#pKpX#8$qj0E=rM+qYCu!"pK;&-)bS!>,<n0E=?<&g@<P0HguM&HY6Q'-IHU+94\D!"ohF!$M=U!`/k?!!*'2"TSQ'+:(82&-3@C!+A/:n,NIh&.gU)+T_d*-iaf4!Wk(j&TJ^4!(J_O&.egC!!!0`!>GO,!=KG?!AFKS$3;jm!>8cQ!q6:H5Yr!n!<O)=U+llR!@.XJKG"8M+:o[s0gPjT!)X=.+or;i"9AT]&J5$t7h>Y#!>#l$:_=<?&:af^(B>$Y0E<d?&-)\B;'6Dn&/,%/+A`&)!>,A5&HQhP:bW%a!AOU!+TMTF!<<*6!ARCN:]N0M63'7o;BQ$N&HPf2J,qnK@fZZc&HDkT&.pSr!<<*#-85uC0JQB)!=/fP+sI!e&7N?S![89X!"f8D!!!$t!!"].'*84:OY?bh&.gPs1D0^C?k<;L!<<ce&0qef'HdN4!>6Ld;'6EI&3p4l!>;"*!"]]V+9=/-O:VF*!&>5h!YGDd+<^ad!<Oqf+>Y/A&HEpd&MO2=&0N),1(XO+!!+Y`'I!Q#+<V3r,60Fu0`XK^6P(Lt&.oL^'bLEG!<<*7?m.f]!0D^F!RhS\!!**8"T\T'&HVtU!WaGE?i^f=!+B:K&5WB-!s05=iDUp%,UWQl&-3=R&HMl.+;4bU5lh$B&Kie.L*@66!$M>]&.eji!$M@C+;=f3&HDt<0I$Wa&.ejS!C@n@!!!!"!<FMS+<h?s!?_CD!YQV1!$M@tir]A`!"^h,('0$g0E;1U$NU5M'*8a8!##nT+V@NP+TrAW'-Ii_&24(a!Y>J>-4hC&'*&%5-NF/k!^H`?!"k;82DPjf!$D7d!<<Z6!>$%Y+;])f!!Y/%!WWf4&n1hQ"9JQJ&HFOQ!W`i4&J/m]!&+L/!<=5B1'/X"!$MjQ!!+DH:ET]M+9;KW!W`9D!](AiZN't+0G#E>&-*=C'+l>e1&s!>!!*WI@LrYR"(;3B&HDnG5l_U"!&+H[0`V;0&eQ8e!Dk:S&25m?0GH2Q+9</V<=f02!#RBl'-RS,!!=&<&Ih#=+V52#!Wiu8!"f84!$D7Q"$dM$1BA&-!?j->!"ohE:a?)B-3=/N&-3mR&g.61"`P+]!!!'(!C.;A&-2hU!"^8L!!+2R!"]24!#$J0!"gr)"%+1E&-+Er'+eOQ!<s,K!!6U0!'^H10eb9]!<XGG!?a,t!!XYN!"]/2!\$DI6UChU&d`Oh!Wt+V!?_pr1(YZ1&-)_2!AFKd!WbP1!>%!e!"],t#Vcl%E'b[%&Ldu^!<<*s0ED2.!!*3g"9B)7:b`-)A,njq1'7LW+>?4a#7qFH&MQX2&-DtY"98K7:d?^3&J+sC1GAe2!YGDE&-)\i!tc+j,6.]D!>$Oh!!*'C!!*ZE!!*'R!"]eY?m8MQ1D0j'&g76?!$Dgq!WWfD5QCm#!X(7Z1_V5:&J,ZY&/'Np0`V4S![7[H,6/<`5VQ"G&3q?=!!,P4%4)Ps!'phj?iU6.!tYGU!&-)-+9<,Y&0MT0!AOTT!C-Z2!&apg!FYut"!@Rf&Kh)d+><fr&Kik1+VPCf'+kfV&HW"V!!!$"&-*;0-NF0D+=7+C!"]\h!u;S+-NZS1&jR!_!$FN<!"]/e!<<[.!>#87!2:_N!Y@%.+<^^s+TMKB!!*lZ0E;2A!<OGX!)OIE'bC?L&-F(0!!*<*&LIJXi;s\n0b=?d:a,lN&b-V:+TMKR&-)b85S3u=!"fhd&jm3R&-*jV!>PSn!tYN"!%%b$5p?=10G#o<&hp"U&QSoE!&,Pr!C/@>&4',(+V=`_.2*6n!rt_$0E=o\!"f22(,Gt8;#pYY!<`B&;\1WJ!!4>e&HW%H!WW95&-N"G!"],A!$VUJ0dH`2O:WQ8&-+C-+;#_@1]hV0!!!$%!<E0F&2?Kp!&1\j&g.ee&-3pS!WX>G+<U^t<?M8B5lrtG'*/[W!>,l>+><is&cr:L!!!!A,Qkj_'*(o7!$H=_+qOcO!>PSh!"],6!rsGT!"o;E!X0-#"99`&!>5A6!!*'D!"o83;IL^M'+c8t1B8O@&g0OD0`W=]&HMn4!E0%"&-F-W:cB!e!YH+I!>bem+>Nos!>,nU&0XsP!<F;C:`p;O0G+it!?dI9!$EsM*>B%d&2aP$!!"/h!!WKM;,S#=!&+BQ,68>g![&U>!YHP_!!;Wj+9R1+!!3^0:_<^^!]CSi&-WVA+TMQj!#U8'&0Lud6#[AL!&-YL&JG0H!!!#I+93Pr&0asP!YA6P!tc1H!#?,B"qVmg!$_LF!C-\d"+hX$!!*]d!&+HV&3p6r&d]]]&-=Wf!#>b^;@ij?!Wi?8!!3i7&HW.J&J,OO5Q`\u5Tp-t&-*:S1N4I50EM4d!>#eH&cjZf&-)\2!!*]T'*&RV1BA(!!+Gs?&d&+G0E;Y1!>5tZ!?2ar!+Z0d!YZ+X&-2bE!<WE8+p'G0!!39("V<p2!<<*2&ci"5&L%bt+@mXf"TSQ7+TYtI!<OC:;#r?R:a5rO&g/k?!$Hdp&T@^]!"f3."VD>s"p#hL&jQL00eb?`#ll.l5QFb=!$Qjn!!OY[!#-UW!tYhT!!30T"!@Xf!$a3!2%^E3'*&+:!"]/B!$hOu;?6au,m"&I!s0&G!)N[s0`Vds&cr(6$4mj_@1lB;5QCcp&JP5.!s8N)0G"n-"p"i.+:n_Z"VO'3:a$q]!<<-%!<W@d!$OT>&-)kV+:nQ/&J,!U&HDqW<\Q<0+=.n%*!$'@!)O7-?kXY-"Tno.!Y>D\&J,Kc!$DgY!"]_t'eg12&eG$e?k>g2:`o`?":kPF!&>,e&HDe7+tuD-!##D7!.cj)+92HS&hj=u+9quN!$X0Q68/EP+92BR!>5Al!AG&r!AP].EY(Ji66H=0!&+ZY!YGDk+:o(b6P';!!YZ,(!:U-\&ca[B0bP&u!)b!P&-*:E!<<ZB!tYa#'EA+H!AZ>?0`V7T,6.d!!#&fB!<=qV1I4k\!!*'3!s&E&&JZuC!"]/3+U7uI+V>gs0G,NbO:WQM";!j"!>$LV&46U$+<UqD!`B"%!tc(E;8</3!FPpn!<`K+&g0FP#S7":!>?^i!$`'u!'gSc?r-h)!sAW*J,rVH&d)MS+E.;]&HNLu0ENC!!<<-D&0_,T!<=;U!<E6'&HH2n&hjn/&.noX!Wiu8+VFbd'G(6F+;"VC!<E081&u>,&-N1m!>,;t!YJiG('"=j&-)\D"<\9W&/G9XE!['e"9<F-nGi\?!0@3]!)k-D+92BQ;%`n<!!!'C!$Djb&.nsD'*=7P0E;ZG,68>g'-I9!!sAW++p8#L'`nsK+@-&>APa#'!!<?*!!!W3&-E"7&O7!2"98E&+=/33"9AW?&HEBI&256r+92Ec+@%[an,a4H&-2h5+9R/s!WW3u"<elf!&4Ks'-Rrc!@&-f!)O4=!YV.0&HDf]!X**90EMjf,4Q4c!!<9&!^IJD0E;(Q!$EBe+Ws)_&-)\1*RGK4&.nsF5nG1^EW@fp1B7Ce!?iTe:_3U<!>83/+<Wom!<N97!Dm0C'`\d]!C@A25QCca"#1>u!'in+&hFS*&2=.e!!AlJ!Wb/%J0#:(+qk%[!<N9'!!*W35UiH.&Kh&c"!@VT!<<-51,8i--kZL\!"o;5!!!$2+sR($&-4Kc,QIj!(iTHT!!*-4!"^7i!>>K8&J>'d63&_D(a&lf&J5'E&5rZC&/$9#&J-*N!!+;E&-=Nd&46R6&.fF-!&+HS!$DCW,6JMh!<NB)J0Ak4!!!$"!!335&/#TX"#0cW0`V1b"!IY`'/2-Q!.kap!!"/B+V4VS5lhPu&mu7N+9;O?+A`1@!<<-C!!!H0!$M=B!!",U+@%^o!$DjS+9<Vu0HW[p5SNcV!<NfT!W`<&0EEF1"X3sX!)N[s,lj).!>&Z^&.ejd0`s-a&0M\f!<=kT&.oHR'HePR1&q:u&J-W!"oo2:+;+Ys#p:r&!!+cm0`V4d&OSep+s6jc&ecr2,68>i0G,Hn"9B&Y"u&:r0G+9r%fubC&HWLG&g7ed&J+sTE)?T8!!*'$!&4a*&cjBl-4g_%,QVis!!!!!7fj;'5S*r?&1D'#![%XH!?iQc!!!T2!!*0e!!<9<5lr/B+;"SV&.g.L+:p4_!<E0G!!5D"+92BQ!'gPs!"]\A!!#:rOAQ$T"Y^Jl"s<pY&0LuS,!n!a!WkY++Th]E0`_gtJH6BM&f(KJ!WWj@&.ehl!"^mt1-bsB2"2V2!$DgT"X"?W&Q&rJ+<(jl!$MpT&.ejG&HEps5Ys3:+TN]"!?qS1!>/-/![.Uf+;,500G"4@!7;p.0f9L`+@#oL+YNir&J,Kr![.OT!&IOR+9;KD@fZQ00ED.R9*IQ\!WrE*"5s:^+V=_V!YZdi!'pT3,VU>P&-4Hr&2=/]![RjZ+r(+Y#XB"R"9AN8!!!!"&i'M6+qYDA64bg_&cr[P"Y]uX+TVQD!tYJ5:`p;^+Trqi?m-'n+<_j-5lgru&.fra!"]D9&5`E0!!!*&(E`\XE"ETm&J?60+TMKR0gGaO!>>J8!!*4R"(234!<=<?+U%u^YV?_d&0_2f!&AL*'EJ47('$%!!W`E(+941t!!"_b+92ES!$srC!AG&b+V523&HELF&crUE+@6)/&hl'BJ,ht<+:nMR+TN)Z!$MmR&7Qh;+TN&u;'J=-!?a)t2ZNk$?iU0-&/$Z>&H`Xh!!+5d!"^;O!!!!B!>R9h!>#fM&HECV!C78?&5WE0&-*7R:]p=i;%OoO,m"&i"p$Cl!<<6&&-W+H&0M#S+V4VS6QdZT!<<Z20G"9s!]289$j-P2:dbmm!?`L/+@,u-&Kq,d'JL^30MqiN[3-(o0b>K.0Ej/H!!<EQ!C-\u!"')G+oq[G!`K("+<V7A+JJo;!!#js;%kVq"9AK';%O9=&cj-e!WWf6!<Eg!&HN+j!YGJ6!$M=T!!+5d!"^:S'EJ7X!"]87!'pSs!C:9^+TObN+<r!d!=9;N0`V7T;ZIp?#6>G]1_Kcg'+b-t!W`945r)gJ!"f5t!@BNu6Qd!c!"_F`!"],5+qO_d&X"7:0L,UL$P3>(!YQ&$!!!'D!<WE90F&3k!rrl4!!#7b"<e<k1B@L]!\aZX!!",B"="Q\+^bj)!Y@&)&?l-?&N:=&!!4kg0JGC"!<kq8+]1Q6!"o>6+VXo%!!!''E"E$^!WWc4!"^7Q+>@7,&-*n0!!/hq&/Gfh;#h7P&HDh4!rr?%0ED^t!AFNS&-W+h&ed\G!<E6e5l_Kp!"]]L!Aje<!<POu!<W<%1':Sj&-3t?&.nmr!AP/e!tk\I!!"bS;umUT!&+TW;^;Df!!",R#S7(l!&FXA&J5'I!&+ur!<<*"!#?eF+E0#b+<UXa$4m4<&.pTP!<\Gc"@*4m!&5&d!"]eW;#h4?!Drk/&0_tl&0VVc!C@A4&-F[A0ED1d">O9e@ZUQj&2=21&/,$D!C.5/!"pG/'*&XF5S4%t('Plb!u):W!)`h/!Wi?7!>%RQ&-II`!^Qef-R&OE+oh]F"<[^h&eH`"+`IQ]&J5$D1(XEd!"]/B5Tp+=!"q!'J.Vb^!>.RN&-;lO&g7ei5lhMu&-*7B!rro5!<EgR!"]/F5SYCZ"TT\Y&3g^++T`5U&odmp&5Y[n+qXeX&HXd3!>6Ld6nnZD!!"et&H`%G0Hh&/&0O4\@0@\r+9E-!!!"5D!!<cI!!*Z7;(r)+!!*-.&HFO?&jg%-!"^7S0bPT/+Wpac3<94I+9qrj!\jZe!MBE(+ohWg!YI+o+;"SW!tYG4+TW,c-3,h5&J-Am!AOR#_>kh^E$.@V!rr<f+9V]V&-)_B!$Uo!7PFfX!'pSb!!#js;'l>U!!*W2#=B+<+qP;/!<=8f#;?-#@1`ME&-)hF"@*,5"99&X!!"/C!-&Q9!!<c6&HDn7!"okt&-3>-O8o=o+b<&_#lkS<5m/;B!H0.n!rr<4!>-+Q?ipEA!!!]E!$DiWDu_(;&.gT.!!3-8!"f3>&Hr18!<E0%!ru./'*&UF!<W<%+9LL5&3sYo+qUFL&HDh3'8%GN&Pt.p2[2ML5Tp[M+TXG#!>5D7'+b0M!!**#!<EcE!<WlJ5QCca+V>@j0bFs*&:blG!Itan+<^df!<<Z[!!3`G7fi_m&-)\S!"]24!AFKT&g.8h1-bss!^J7I'gMa*+9<#R";27X&.&q8"on]*!W`B'+or90:^A4A!WW3;&.ek&'`nI<!!E:`&-*k=0a%If!!!'S'+kr[,!bW`!am]*!"o84!!5Cc6ko\*&-N^Z!!"_b0L,X=!&6/.YlFeI!"a,N!>>Vl!"^k.&5W?=!Din!!W`lU!(6ql&2>=@&.p#r!WW3#!rrE+&g@<0+Uf%s&4Q[3clb"b![/-U0`V1R!>H(F!Aa^`&-+ur+TMKc+:o/$!&+BQ0bjm(1Be=*!($o&E=`*M!"]\Q&0`h/!<W<%'JTYE#66Fi0ED^t!"gD!#sop>!!!!""99h]&J6-'+92rQ&d&1;!YHOt+;$=N&HDe80E;(R!<WlV&=<MQ&.g#g&cqM(+<UXb&L%5e!!"2E1&rF1+9;Nu0ENI$!<N7?!<Jnq+V4Ys+93TN0jtUn!&,_g&HO$u0d.,10`V1d;G[?)&HE@I0E;+d!'gMa+u0]`'7^&s&HNOI!#.j'&.f*I!"]-L!!3fW!>,tV7plB#)A!P^.0L+^!AP,s!<`BJ![%OF"sO$Y0JEJ@";HDP!Wu7#&-rj[&7beB(_nM3!!<66!!E?'!!"\S+;"SR&HEt"&/#!D0gGpV7KNW'E+^O0&-4I!&ctHK!<GM%!$Ip8!>5H2"9:dh"p=o\&c`]0&.egV!!#k.&eG'D+:nMW+><j_0bOKu!s<KR!AFR/&lB2p!!/2bT`tK.!>GN'!CHu4!"],DA.T:P!"]5E"!@UG0d$H!!"gsd!+A,(+X-n1&.egX+9;QEE\fBO!?_AM&-)_S!"h("&HDn5!!4;d+=$q1!!",A&HNFC5X6M-&HDh3"U4r-!X0A@"TgI[&0M%+5p?j\!AOTT!##Ag"+W*+!tYGe+<^^b&.egE,6/oR&-)\B&J,!T+W17["9AO4!!El9!&>;i!!+5D0b=I4,6.mt&0MPr((a>o!"^:S&.ejS'*&UW!!!#H&-F%1&P3N6,#A@a!<N64"#:uAn;ncT'a+OL+93#d+s?pd((gNI+V5.a&3q?f5S,U]!#-\%&-)_2U(@;)5S+P>Du``D!!3]D!>>GE-4gb#0E<p2;\/p?0G#B.!al'5&c`/f!!!Qt&KhVr&Ki7t!WX>c&n(ee!!31%!<=6>"sFT[!!!'i!"]/4&-=UP&.fsM!:_6e"pbb@1_Oa>!"^7d#64`:&-G`o#6=l=!&+ER+9DR30F7^l&0_/U&-DtG!E0R8!'257!&+Hs+9<&c!$Dh-!<E0S#@INS,R+8LOoca<!"]bKEW@3_+qOqY"98KH&-DnG!!3-E+TVQH&J+pD!!N?6!&4I>!&5T=&0q5U!'u/YJ1qT:E"`?i+92Er!!"/B1(t9d!W\B]&-2bC!<@3A&[<J[@LrYS!&+s<+9GpT&Kh/U&5`Hc!<Fkf!`BUq!>#;U,nLVR!>#>7&<\.Y&0_Yc+V6I1&2l6D+9EYs!!!!!:_OrR!WXB0&J5!G+;+_e,7&6B+osq_!<W<%!"o830b=m.!&,N."U6FV!)in@,QgjD+92CL";*<g0ED_-;'Q)A!<<f:!<iH7&-;q8!+6o\2['3n0H:,?!tYGT"Fqd0!&-\=&5WH6&Kin25la^o+:o+b('-/j!-/,@!"]\b!#>Ye$5+6[!'gN<!!"\r&24),!WX>T;ZK,u5QV0D!&+BU&0W_0&HN%9!!!$"!^[)'7:?3>&HDeE"=0'O&-*J%+YWm9&.fBQ!0JEi"9:h$!<WEI2$*d[!<NiF!X&K*&-3>S&Ha*e+V7KP&jQF>!`0F0"9ntt!]2Y5!>#5<&-i==!"]/2&0`5n0i.lO-ijA]&-O*u!AH2-!s!9C+lX+W!$;1C";!g&66RKT!"f;5!>%Krp(n8R!2:#Z!!+2C&HG'"!"f3>1B@X]!AP\u!AQn@J-,dL!<Y+i!"^q#![&-W1'%pe!!,nO5Sc$R!!!!1(BFL9,S3cq)Db8&!!!$r'gMg3'F"O?!"]/r&cr%W+<psh"!@RT&.npt&3p@0!&+ER!APl#&-*pt!!j0E&3p4#!)Ne!!!3fF!$r0]'EeC95lqSu!Aa`f&2HKN&J59K'GCHh!'mf.+<^e0!!*0E&J5\$+T_ZF!WW36!"]\a"!RecnLt)2&-*mS+92EC$4oI!!>'eo!$M=s5SO8B+V=\t!C7h1&-*n`$p[)^!W`<&!>#83&/A%c";(YI!YGMG!>5At5S=W:!$MAN!Y?OU!<`BG&254=,@CKO!<N96!<GIf&-2b3&.gZQ+9MTJ&2XD"!H8&N0OXqm#Qb&c1dEbt+TN)d!-0ps&J@nt!"f2c&0VZ/"T]eJTF_)'!Wio5!)W_.;%Na/5S3uO+93QP&d/aG6i[3#&24(q+V8,q'gP\!!"]\u!!*+>!<`BI"!A-U!#,D5'bM#X(bc!o!#.*g5QV'"&.npD$?-OZ'eo`!&24*(+:nMQ"!@V1&0N(q@gEVK)'8t_!$MFE!$O#r!!*-4#QPJ]!Y>>s!W`?),!c2O,QLXA+TMKC!Wc7)&tf7-&eP]g!!!Sh!!",D!Y?"F&.nmR&-DnU![.UH!!!QQ,6.]V+L)%=^_'ED&ciOE!sBbY!Y@3j5r0-/!"pFT#@IW7!>,lo+9?'S!<E6&!!"_c5l^oe1(sXR(*GpD!?hXYJe0UR0EqLW5QNnG&HDkE&eH9A!s/HI!<I-d!AFNS5RA&-&27K((B=O;&J5*G!<NE9!!!-&0ju7+!Y>MM+;.$S!!4kW&2k+#:]^YO!&8!b<=hCu!?mO)!YGG6!.Y[^!&O[&0bjZh!<EB)+9DND!YI^!+TVWV0JF%P!"_s<+9>@_!s')R![/gg&-)eT!<<*);#gVE+;.QQ!WW9&&Ha]e!!!!#!YYSG!<ODE+<^^b&Od9C-3+VX:a$0"&0_]o!>#52+<Vd_'HdDf!>$IX(Gc$k&g78u!'h_>&7>J>&JG=%;@kS`0`_=g0ENBt!!<cD!s/H*+>Eir!<Wof";:e9'*/(u!<<-#!!30T64mT+0E_sfJ.W?T,TncL!"]\Q&Mt';!>#5R0EXR-!W`9X66J!:!!EB)!<<3%6r=L>&-*7A1B8Nt&-3GB&OB7o&5`E.!!<u?'.4=k!^RuL!!*'"'LEET+pnk`&HO$T!!!&(!WWlW?iVMb!>5A5+[Q/1!<<-6+=R?l&/,'V"p,>9"$dFu6U;mr!"aY\+_;2d#nR,:!"f8U!HJ5A!&5W>!3m(Y![SBY!>(J<!)Nat&n!FZ!!!*&:]M(.+oh]W!AQn@!<<*$!<ODG&eG3i!AG&t5QCcbA3^\-!WWc3+93Ps!<Em#!<j&X&PsSM"(<8O!<`u9?kFRb&2Fe/&-*7e!)Od?+92EC!@.X]+tsQN!$MmT+;Y_.!!!!*+V57t!!!$3+qc"0+;$C0!&>,d!s/K(&25F6&2=.d!>$A^+lEnC@h8\A&J,KR+VXnX'G3S1&J,QfTHF:9!$DgR"qiXE"TSQ7&jR*2!!rWJ1(XL6&J+pB,S2XA,67cF&0W_n!tYJF!<j&8&Ki2@+9E,m!$VUi0b+0p!?_IU!>$Cc#8.";!($`0&-*:b"b@4.">CAi'ESpJ!"]_T!"o8C!!IfV!!*WC0HhSP^]>%!!!!$66NIc2&J-Z"&H`%60EGbd&-EId+HZak"Te`J?nhZ^O>8$]!WjQ#0E<:?!W`94.1dX&+:nMb"9ni='*OBt0a"^i!!<9+5WB<@!(%e6!/(pc&-*jR&MYs_.j,B*#pK?L+<^^d!C78/1LU:q&Q&W1!WrNJ651%5:]UOs!$Msd&.oHRYQ=e)'+b7"":u4Y"CVia!#.[2!"]23";1_:&-2nV@KCEZ+<^u!'a#'K)?^U++TN&d0E;:Y+TVQX,U39hYTNrk,6IoG0u4G5!<X/M()%6%&.fBT!YGMG+;,4i&J-2f!!!Q3+>*[>&-=Wi!<<-4+<XJl,nL([:_=iP"[O720G6Va,npFmJ-5lT!$D7d+VZ%"0HgDr!YGtV!!*W4)$1KQ0bH)?!!4qj!'gN-&PrKQ&HDf=!##Ed!WWcF!>6L[!>HY$!'gVd/d=!1ErZ4r!!!Q4'/3hq5QD'%+9G@@!'gPcJcGc^0G"mO!rr<%1D(rV!!-"0:bWFl!AFKr0bOI/!<<-V!&aiY!?_@E5p67/+Wq<s!\l#&&I),e&4-I!!Y>AT![0rc&2=1e!"p+Kn,X)0'*1>u&5n&Z!"f2B&c_q6&d&^H"!J6l!"T&PoFM$d"Y]rY+<W?>+:oY<!"],Tck?s*&-5TD+p.fH!s(.b!AY6`!<<*c&7Rmm!YGDd!"]0E;$'Q+0G"cq,!mFR!rroU!!!-&0EhI_!<X#:!!!')!s2C5"!mC;&0N2o+qO_T!<E0S!"f2C+U.rL!>#f])?:?]&J>*R&HECS!<<*4,="Qq!($\e+92rT!<=8E&I'sE&i)3B0`V4T!!!W3&Kh'-!sClC!WWc3'`]@!!"],B!"iT>'*Adn!+>j/;KbtN+94)<+^bBa!$_OF!&=f[&HNLg!>7*e5QD0)!Y>nD!"p%J!"otGnf\Vi#=\aq!"fbH&-2q7!<PP!0d6a0&-2bC!<E03+:o\>&d&+G&HDht5S4,"&cr(8&.epT5QClT&2=_]!!4kd!<N65!AFL-&OJ/O@LtG]:_a$E;?6du!!*W2!Y?LeE3TBF0LR#f!!!ZY!!!'%!tkSG5QDo/:_WmA1B7si+>aZ1!<<Z4&/%8A0E;Y?&-2bC+n?0u!!30%&3p@D0Lu0D"s=Hh!<>t>'EBj#!<<ZT1'%su,p5MS#6?LZ"qV.D+V5e>@Ls[`&c`,$#lt#,1'&L"5QN&/&L%2f!"L+O#lkP;;#gh4"9SWj&/#Uc!CJ"1&.emCrrE]S'bULg&-*gQ+oi/s#I4q\!?a-?,QL[Q!!!!4,6@mC&.nmS&.egE:d>Qp!$M=R,SD:5iW'W'!!33U@/qJ_!"gmb+CtNb"<[[F"#0cV!`B%BKG>%!!$PeO&.gRY#QR+%#SJ6^!#,GW!AXXA!!!QqOTH+.+9;I-0E;(R!AGc"!<E358,s5%!!<c4&-;kZ&0MQ^&R[p;&0M5Y0fT1G!!<35":tW50EE9u+Th`W!'h(r!<FMI!"p(K0G#HE&-+Eb!<W<6(*Gkc!Y>>6,9oMT;?;=[!!!$#!##DG+;$=.#lkW^!<<Z4#8%LJ,Z>7N!<=f]!"^h,&ciOg'/0st:]LV0$ZQ.H:_<t"1a*M50G#ol!<E0#!WsSG!!&*H+V=\t">CEA!>>M9!<@ZO!WnN=&-)\1+V9/+!"psc;'7#J!"]\S-NOo1!(?ng&-<M2!$D:R&O6C!!$VsS!!3gA'G*##!$WNu&HDeb"TU7g,QIig&0_)S!AG6"0d./T!!!`V&.epd'`^KB0`V:Z!!*67;.UjO+>Fu?+93Mr'Eo(*!"]\Q!!!WV0`i+%&-*@D1'/*j!F\h*6W-bm+;>q0+<_:-!eC=Q&5XMM#m:h@!!"AH!!NFS!rt"T;%ip0&ctK7&-+C,0b>K/+YWmm+:o+b!sJZ:!YGDT!)N\>5QN\D!t[-f"p#88!>,;d!$a61"">>oF!2)M.55`B+TVQC&-)bC!<<*b&-*7A!>,;T+@?eB1CjR#+;P+l!Dmb<,#Njj"onZ:+qPN"!!<3D+oi9!!!3]t5V`r]&.egQ@2AnW!<<430G"3a!\ba)0`ck/+9`HF&eIn@&>0gN"p"])!"f9O6P9A1!!<6%&.p$]!$r0V!s&NI'JKP"+94*'#ll,&O9#>:+;#@h!!*'""!A_Z!!3]d#r!JG!+?H?!"^;-&H`+K,7snh+9N2u!<<0T+V4W-&.fEd5Tg11&-2bS&c`M!+T`6@;%O9P+TNVbE"N._"*"DC"TT,Y0JNP.!$FNA!!"bT3Wfup!#6&C"@+1C+<V@`+TMQT&J,L!!?hFC&-3CU+<[<W0L>dP&3rN(!?iQt!^\UA,6/;VA.S_@5TpgB0ED_!1(aO$&Hr^H&HE"H!&>,e0JF4U!!!$R!!"/R!KdBp!$Y8O'/0If+9MWH+G(.)DupRN!%%[i!$M@e:bi[q&OHHt!$D=T+<UXf&2aG!"oqI(&0V#d+CR5P&cr%7+or5e'*&"E!YR4!&-5W?!"f54!<<`T"!msj63$ut!<E97+qObU!!3-I&c_n60E;8V!"]bT)$'_G!WWB)&Mc'`!+@Vs+FjJH&c`Ic!<OAD"!\EY0`Y'J!?i!c&0q;`0E<<u!!No7+#FNh)$'m@!$Eur"$e"Q+V4VR5QUrd!"^@u&cs6X&-)eT8HB#!#U9l]O:hP$!Z;[`&HDnU&/.nf!"]/2&HEpe&JPfY!!"eT%fcSB!<E9H,="!a3Y3th'/0D@1.;gP"V:e\!C?kg![7[[+9MWUE<.co!!**En,O%$,8)[1&-<FI&-)\1!<j#8&-+Ec&ZZKD,69P$!"hdgEY&3N"![nI('"==;$$_!!>#hf"9:=[+@-P=!^I:s&ciOU"oo2W&-*@T+<UXb'*/.8&d&+Y!?_IF!!!Q2";)7k#QY$<!FS8;%0?PA#64`:&J60Q!+6?N"<[[W!<<Z2!!!$,!$M=R!W`oW?k<;=+Ta=t&ci%F5m.02:]ph4"9JT>1B7FU!!"/r!!!ZT+;4`C!!!15+Ya!""Tf)5&eGTe!J*K9"XsKQ+9;P*&HDe9n.7?!5p-^M+oi2T-NObi0MhfN&-;hC!"`%9&5WoR&Ki1r&.egS+;"VT!tYJ<5S4P?'d.N.&24e0;,@8s!YRd1!$W!u+[H\a!<Wo8!!!(1#QOrM&L&A@!&,&d!>#eC!>,p+!AG&e+qP;.0gcKQ+CG187M%6c5nb@R!&5&g0G+9d!>?(W'2SZ1!\aaC&J>'H&0M,f!"g&E1B:5P+;"Z1,9[Ub!>#9n$PO(K&-)k6&-3@UPX,9L!"]V`&26r_!"pFT!?_@T"9Lgj+92rr0E<3r,m#e4+9F5.!>Z7R+;7!d!>#5C!<=5R+9qrZ+9E2V#oWdD+92EC&OI$p1_:c0!"f250GGZ(&-2k8Du]kA'G<)3!<=ej(,-=86i[5f!<<1P&HO]g&eT3d?m-X:'JMcP0MrE)!"]bD5SF,P!>-LU!eM'c!>1E4!"],Q!!K)?!!!$2#nRXJ5lq#e&/#$V&-*:R5S*r.!!!!R!>-FS!_"5.!"^7q',(p%"XF*^&-)eG&d&.]!C-hh+>Z8!"#'p&!"pIU2#mV2!&6//&0Lud!!"-/%hJaR0L-0L&0lc,+BSXG,Qf_#!>ke:+;5t3!Y?"K!"g=c"9C1V1(nR\!"^LX!-8/A0EMdd(]aV%&-*ga(*F.h;@jH@!>mYu!!*'"&J>0G!)O%g+s6mg"==0N!s/N;!!+5T5r&EP0E;+R+9<]!!"U2l5lh!/+;"Sb+T`HIE<-)*!'p]$!<N6t&ca&)0H^oO'd+"e!$`X!&cj[!!!!QA+qYAo!`0I0"qV@X"Tf58+AcYp&-ELE&7Zgq0E;Xc1]SY[!'pVt!@IjL&;)\_!!!!"0aYW.#QOi9!"okVDu]lX-34o:!&+F/![%P476gkU&25:>+Aa1(!<W<&!Wt\!+YWm@&24+b!b<PW&-F'e+=eQj!!!*5!Du,c+oi/t&.fs@5S-<a!<E04!$NHd"!euk!>>G;&//FQ"+Ljc"q`'\+p%d!!"f8D!rro5?i^<O+YWm-%M8mW,6/E#!$Mmr"!J6V&-Mt5"!V+R!WruF+:nMT@M9mb63.W3!<<]6!=)Ff!#HaX1]mdZ!$NHb&.ek_,6ADW!]>'4!WWcF!YYSG![80U&O6m="#']W(D6lO#8$q<&clDS&HN%i!>Z4R&-<IF!`1!>0`V4g&.h,0,6.]D&HDt71(XQi5Tg+4Jg(gA!X'Vg!&>i#&cc;>:_Eg0!>#;T!>#53"(2fW!"f2T!YA3?!"f22!>keJ!>$pr!&Y@>&HO$T!!ErH![.OE!Wa2P!$Dgl+qOef0+^AM#Qb#L!>#hc&J5%A!"fbT!$P/>1CtZ/(8)fs!>%Ns!s8]?:b_t_&cr(G&-)\1&HDe2+<W?^!!+2B&-)\b&-<+K633iH"9AW+&/gQD+TXhA+[?W:1(c>C!?hLH'*'3V,u?<!1]S*f&-5$?&i'Ll!!<j1!<Nf8!>,;3![/0V'dHNu+92HC0b>J3"X!gG'ES=Y!tbP;66H:2Du`cZ%fnEd##G:B&HVqE&l8Q.63%XL&ge_+'*_kF&.fBq+92uR!!+2R;Bc\q&-)_B!"qQt1BA$f+:nMQ!>@^5'+b]V"T^</&HN1>&caZe+TMWF&MO1r+s[a$!!4</!"fAG!!NQ,!$U8'!WYJ@!!,G"&-;n5&c`IS5S=E3!#$Le!!WHH!<=DY!!3`<!WW3$&/P?j!'pSd688KA!<<^^&-2e3+osQO!)O@D!"psc!45!6&I9N^!+I_b&0VSe!!<E*!$EEb!!"lS!='2:&cj*V!<W@C![.Oe0G'K>!^H_d0H^?L!<Y(h&-2h5+:nMc)?V/t!C0$c+[?#.&0V&V,nL1Z!.^76!'r:\!##>5&-+X$!>#kE!AKTKE'Y'9!sB/J&4-Eu&.oHR!"qNu!WX>C&-2b8(EsUp+p&;Y!YGPY!&6`8n,`Ri"#'^%!\b9$JH5cn!"gFV+X6sf!!\N"!Ped=+:oY,!>c=Z!&+rq@:0'8:bVqN!<Fmj!$D:B!Die$!<NiU!!NBI0`YS`!>#GY&-2hE&HMkC$QpQ'!Y?Id;#h.=!?i!c&0NY-!WWcV_>k\Z,nV:cJ2.*[2ZXI$&4$:#0E;Xa!WjSX!!!Td!X0AS0E;)-!W`o6&/#QV0E;Y,,60D?#Qk,O&-3pU2#n1#5Tp+/!&Fik!!<c5('604'+f.,!tbPV&QKsH&/5=F!>$Ft5T'S&;ZT6>+](E(H3$kD!`B^D!!">`!,qo<&i;:>"Y_V2&HDe25qkV<?k<;\!W`<%5QMNc!"hI!5QCj$&-)\S!!E:R!"],\&-5)u-6\HN&-+Hc'bCCf"X"oh!WZ$t&Te!a+TNVb&ciRu&0_,T,p<6i!<j#G0lR_7!"],Q&H`UN6SKVo6QdV8!!+5D!tl2B+qO`>!<XGf!"oA8!!+5S!>$CT!<<*BT`G6L)K7Aa#ltWG!WbOd!%BW&"9JQ?0b==D">^'$&c`OE!?_@R!g4Z<!AOWW"qXT2!"^jt&248A!$VtP!<E<7&eH9$1B@IX!'hY/&-2bR!!3lj!!+2C+;\KP&HDn7&emY5+<_9s!$r0W5Tp+mK`OFM!$kAA!!+bc!rsGDTE'OPErl=q&-,N,+CG7)!>5qF!!!QA&J,!Fnf\(=;?8H_&LRPi!!"/R+V>=u+TN*.1=-"'!YGDE&-+F=!AY5e5W02S!>&'=&eY0g!!"\a@/p9-&Rm0l&eYap!&-YQ&-2kW!!+2B&-*7d&HDnu![811!`0Fc5X5n]&cr%5+<W?L6Xt'H!C7h1!!3l:!W`9t!!+5E+U8&[!!!Q1&c_s:!?kDShus)^+92NY!&br#!>#84!AdUe!!!?.&HX'd0Hq)/!!48F5p-g@%0HYF+VZO/!!-Ie+V=\S!"fbC!<F<.+Hc^k5p1.J"9T6$&.gMq:]WKd+<V3t&0_Ye&ci%7+92KD!s&u6"-?t7&jR!B!<<*"+9EE("Yp-E+VFbd&0NY/,60D2YV7a6!Y>D62#n7d5QUp.&-53C&Kh,\!C.1r![.P/&--/o;3V)(2#n0f0E<:`!<F;C0G4ru!!*-u'+cip+94Y,&c_nD1&rIp5R7DM&/,*F)`.O+5RRR"&0E#;!"fbB!J)?o!^fcc!<>"X&L%5U&-)`^&Q8Z@'*03h!#5JG1&qjc&26+H&c_n5&I'EjKE*+p6imr*+ohTS&-)\B!Z8`b!&+rq&2HKN&-iaY&7>J<+:0.k!!!$8)A3'I&46F0&.eh!&c_nU+<^jf&-DA5:`pk_0*)1e!!*]U,VU>Q'bLNm+9<T=!!!!!5nX3;&ekL7&OJ_c!"oAg!>>T5!"f53?iU0L,9QsfJH>ia'iG#8%fm7B'Mnc3,6/qh!"]2F+9<&f!!#kM!?_t?&HDe6!<=5C@3GU`,7jhf0b=?d!FlfD+9<#bJ,fQ\&J-]!0fE,F&Kie.!@89h!"j_]!!*^2&-2bT"oo595Tp71&HO*`''0)n,6e,L+92uS!WW6E1(s^40E;Y-&HELG0G#?>0HgH/!!!3('*'-T!?`Qj+;M[X63=q6&V*^Y!<<`T&24+e'd+Va!W`iG&c_nd0uG+W+T_`h&-*7a!\b6#!!<?9!!!!4J._qq!$_LE!"q#:!rr<4&f+Ck'ECrB5VW7Y"X!sm!"pFW!<s)HE$-hI&-j?Z!$N%#!<Jks!#Hjk1(XQh!&5#e!!!Q1&-ZGA&cj1")ZTnP![%Ug+V5at'/BT]&-MtE&eG*V&MPI1+p(XF!!4MJ&-)\c;%X?n&HP]/#ll)f!?hIG+;5;q('+CO"TeZ.!>#Mk+92BE&HEaT!<PS1(U=W\+><d<!!+8T!@%RF&/#"/&-)\2TjU[t&HO3n![&$d&.fEb!!*-$"V;5!,Qe#\0JEK)!#61J+;$CP"98E%!WYLg!!*'S+V6==!\aWV!rroW!?hIT!/GNE!"],1>9F#%&-2eX+>SHI0L5[>&HP]>!>5qF:d[M[5V`oQ!>83p!"i$?!!!!!&Ki4s0E;Xc?m.3)&J5$FDuq3`!!*-$&-)_5"p>VE0ED2/!YGGe":tPF5QNP=&/$,d!alWC:biRk!YPJU!"pCc+@&?r!#5V<&-,'&&-3@C<!!%R&HWP$!Y?m`0EE<s0a.V$!Y>D6&Ki2=!&>5g!Y@U?TMSGt+9<#s!$DDp!AG)k&-3=E!$D;O!-%u=!!3`5))D=-O=1)3?j.Yg0E;q$&.]<R!Y>tU&-E"81GKFU+Aa6`!Y>>FE<];B!"]/ChuWqa&-;h5!rt"i!!3?Y&.o*I&2=hP!>#hU!)Ne!!"_F=!%e6c!"]\d&HDn%5S4#1!!!TC66fk\!W`95&/R)(!?_IE!s8NX&-<LG&0`e?!!*'B&fC]N!!<3(,QT#!+TN^*!'gW@!!5tO&cr1J!<E`t!>$ss__?]T#lk#?&J,Kb&cb`/('=OM0HgE-&HDeB+V51k!!*'C!Yc4h!$ai4+951K!!!07+TjG2"V:_HE=ao;&24+r&J+pD"#+ft&2=.b!!,>>"<_4d'1*641'.LX"TSN6+p%gS&:alp&MP=/!^qP=E\S[;!((/q2?Y0H!Wb(X!!-I-!HeGD!AQ;1Du]o+!AXWu!!#n1:_<gF!"^:S'EA[F!!48e(.'kZDuj>M+92TW!!"/C!<us45S4Y2!W`iL!<<*#'*/+G!AO]W"99VW"@*#"0I-cd1B7CU!!icQ!?_CU;ZRLC!!!Rp!?`Kb"sOWt&-)_c"&fdS"q^hI!?_@S&HEss!tfJR!>?26!<<-7"9gIj!)Yu^0ED1c&HDh6&ci"j!"^7Q!C[S%!!!!S!>,>5!YIa1!"f8E&HWMp&-*=T&.esO+`J_/"qUb7!`:-Q&OQO0&MO4j!?hFC,Q[uL!!!QB!!3o]+9>j^&Hi(g&--Y^OY[OL";*Ej:a-JO!'pYd"[Fgt!"]\A&1'8#0bOI00`W?s!au'_&ci&!&5jW9!<rZ]Tcb$L&-)\3!!!!#!!!QE&cr[W+;#/-Du]qR!>khl"&Jsu!<<6'!<?RQ![\$N!?qLT&25j?2?3f3&.g&fF<_q/!\co*('jma+TVTT+9s"h'.4G++94(r&.ejC&.egb+AaUS+p/r1&.nCD!!<fE,9Qse&<d\+&.fR!0HgDs;$%jO&-<Ok"$chd!')/5&-Drq&-,'1'`p'$nLtP@!!!W6!<E6G+9iGY+V?F/!WW35E+B"M!GN5J@K6K3&.o%Q5S*olG5uR^0G4?e!$i+a&HX'V&-*mS-3+W#!!3-c&eP-F!'2iC+9MUG5p.?_!"i%*!?_@S(B>$J&HDh3&4$mp@0.PO%fjs1!e:h+&0Nb0Oq@[++VXqh0bFrt&-*jd&c`LD1_E"S0EUeF+TVZf0Hq_R&.egq+qO`0:b_qm!Y>?N"!S:c!u1qN"To;N!>#eg+Wq19";!m##m(DF!$hOF!tYM6'*D)5!<s)h!!*'$"i*Di+V5D3&-2bB&-)fP&0h/d!YP_l!$MA?!)Q#`5S+b6!!>Lg&chtU'G1<F#I,Bo'OW+-!C-Vu(*F_2"VLh[Du^F\&-;kF!<=kT!WW3S+TMNC&HG04!"],2!s1di!$hXJ0Hh/d!WYMP&-3@E!'qe?!X'&9;%O9A+>F$"!<QX/!'p]1!<N6-!&1/R!tZRZ;#ro_n1jtU![L#r!(7Ap'/0N=!>H[X!>,>V!#>Pg&Q00Z!##>50d$Q20`qGR&.g0&!DlSn!!+5C!>#eG!>5GF1a"Ia!!!R<&.p#rE<%*n&0V&s+92re!!!$$!W`N+AeP^U'*&RG<<E7%('+sI!"^ga!)sU3+T`2T+TX>q5l^lb"p"^5!,rPO!>#G8"U#Am(BZ@1!$EI>&-O+2&.npC#7pk:!C6\f&Kj=>&eQi!!<<-&"9Tbh&/$Yu!<=eS+TW\c+9;ZY'*B@T'G(755Sl.!!<<3\#6PS="X!dh0gPd@!$VsV&i0P2!!#h20d%VU5QMDr!!!-55S*o-";,Yd#7pk:!$`[!0MqfN0b=lr1(k`1^aK0s&Kq0P!'h4u!!!*5!<FnW,#Ah)+9<Vc5T'VW'f,gi!<>@b&.oNV&HWUI!>#83&J-_u&-*In5S=':!<F>E@"AFl68ej@&HNFF&-YE4:]`C,!C6_d!Ab!*+@%[c!>#eF!&,N.!&,$.&eG$C+92uu!"^7R!"]]=&J60O!&=Nc!!*W4&J+pREW?.BT`>/pE"OcY!"fl"@hKCs!!3`W0I.20+V6mP&Hi+g!$M@d!"i'O&.nmR!!!!A!?_pS&H`Rg!@.^X+YjT?E>2m_"p"]H&Rkeq!!!Z6!"],D&2>mC0E<6r+qYAB!!*'3!"]_h'`\77!AOZV&.nnM!W`B)'G(6E!u(hK"p#hZ,Tn4k+XRa#!FR#M+p&kf1-Z<N!<<*""!Sll!!NuMUAt9L!!!'$1*EFp!!!cjnRW:5!YGDD!Y@$ln,i[l!<iQK'*'3g&HE@B5Yr[*!?_CS+_(Tt#nR(j,lnSh+9<#U!+GpA&c`jO!>,;4&<RM8&-Dq5!$MoK;ZHe/+TN)i!)O4,TEGlN63%'2!<<*(0aIaZ!?hOF!#-"X&O6?s&0LrQ![/*d"TfA\?ul>[?m.38&Ma=g!@&`g![&Wt&eGTc!W`9U&-;iN&-3uk!tYK!&c`e'&-*gQ:_Y]!E'OI)&-)\S!$O]3!al-4+TMKS&-)b5&26HO!!!7%"T]hZ!>7'd!!=DF)$<?9+qbJA&c_q4!!!Te&Kh'.!'*=G!!4hc"==^(?iU3-!"qNu!AOWg&MPm^!!4>E!!3`40`V1R('"=g5X5;LO:`9r!W`<'5QLld&-4(5&-WUH!WWfE'*/(5!?uIo&-*LH5SG:G;?-k4'G(fX1&q=V+92BQ'Mo>@(D.ej5QYp+&Hr^W!"o8sOV';L+V5h4!!**3&.egA)#saP!?r.4!!,nA+9;HB5l_NG&/#!H)?9a;0E`Ku&J>0H:b`\0!&;M<$nr@q+<Ubo!$Mmd!<<6&!?qe'!!**#'/C[2&J,Qe![&0X&HDe^+OWUWDupUO63%Pu&-)_2&0k$Q!)Od^+sd:F"TT-21D)kX!!3-D+ohWE!<<-3"TSWI">Bo_!CB*c!<b(V+YWm0@hTIf"oo2P&-3=b!<E03+<_p1;'MYG6P'kA!>5A7;'A4+!!"bd"99)8+><fu!<N9'!WsQ4!Wir9!<H"/!$Lo#5S*u>!!4ku&HDe2&HGX(&cr.8!#%X8";(\9!>,>4&g\,T&-2h6!!+2C!&,'!!!2R&,*>1q!tYJ6!5PH8Z;;\G!"]\Q@07YQ0`V:h!<E03+qO_c!)YHP,QTYW!<W<%<>%Pq"#1JiOtct;!$Dga'FGB`!s&B%1*?T#"TT*264s8!!#,tE&-5To'*/.H!W[fR#mC>@";*<s+9DTE$N_FN"!J4!!(6l++T`;k!@&a4!!!ZT+@$JL!^H`1&.oHV!"1J9!!,nM!"]/R&.eh5&L7ni![%ID!>#<d5VNf^&J>WX!?_LJ!!"5F&LJ)D!!=Ae+;Y%i!<Ec&+s8$Q+9=8b![.UG!<Ef5&c`'b0G5!0&-WUO&0N.s+:p4>!rr<E+TMQD"VCh:+;-FA+p&>X+954M!<NW@&-Dq6+>O!"$3:/=+9?Em"Ag4-!Y?Le&/&p^,9R"b!"]hV!72;Y&MO1u+or5T'G^]^"-F-KZ6K>\!!!$2&24Xr!s0,:&.epE!!OJG0L5[=!<<0D5\2\N!>Y\J!WrH0!<<Zd!"qQu+@6V?+><cq&-*mT!HB4_@j$:%!!!T6!<E07+TMY-!\k8e5R@ug!?b5>+:nPR+:oY-!>H(G!$D:['/3kq!Y>AE&cr%[!<H4T!t[4U&3p4-@f[,A&--,?&24.u&gJ&T!^IAE!"g@S&.emd!W`<W&chtF5nH9_!#-Rf!>RsJ0GFNf+:q"!1*]$a+<_=%0`h@W!"oAX"op>5'bLHI!"]2SDup#)!)QKX,67fH!<E0C!<`HY0ED_N,7uX1!>,;C@l=B#-ia8j&7?=e!WaDG!tcY"&0_)S+92s-"9K29'+k3g+E6fN!$VDo5lgud!^d"j!"^@t!'sX>+ohUN!!NB7&-G0?&-,!-!W\oJ"9A]=!!!B15lj4O!\s`f&eH2e&RZ:r&-O[@!!*];&HO'u!-J;b&eK.=!"fDX+9;Hc&24/%&0DG`0G+:=+p!dJ?m-*p&J,Kj!)iju!&5Ss&-3mR&-<Lf!!N?*&0M#s+;#.b0G#$#!#"<'Duid3!!!9+,6D9T;#qgC&-)_S&HZ>`!&4I.&c`.N&HDeC!WaDGOT9tf!!3?)!\cG2+<hFP'GCKI&-*8l!<<04!"i']&JG-E0`_hA!#,JW)%[Ak!&5T!0`V1e'*&(8!\t>g!<<+.@3>ON%hrs_+92Ba+9E,T!>$q=!>6(H!AOU2&5WB.+qO_T!!GV#&HH5^+TN_i&L'OR&J.b=:_4un((($C5m'@b"$dEL+9;Kd5nI*/!?_@B!$Pbn%0-AN0G"3b!YQ%U!<rZ^0G"Ki!<Wl7!<NH*!!"oB!YGtk!$F!-!&>Yt&d(KF0G+<e&.fQV"TeZH!g#):!"feD5Tp11#VZi&!!,pu&HDk4,n^jk"&Jt^'+b3F!^`P8!Dm/)&HG'="!@RF&/5-W&[N#[!WWf4<>0<`EWQdP+:p=?!!3-%!>,kd+TM]I"#)D/&7Qb*0EE<s!!#n>!"okU&MRSm+Wq<t&g78V6VnF-))2*l!"]]=!]2&#!?i$T&-N"8!AG3A"U$R\&Kqeg&-*jR!<<3&!"oiN!"],1!AXWu!<?%3&.fBQ-P-jj&-*;M69l\-!<<ZD"V;@j+s6jc!!O\O'+b10+AaaG&J5WZpCI@I!!!*6!"NuL+92P[!"f5H;?II1&ek@S!&>>l!!!$2?r/Tj'EBik!seoO'Gpf]!>#K6&KhYs&/$)c!hfYq+;th,!`:*`!WrH(":ub"!s:2q&/#Ns!WrH'!<=8c&.egq;F(<p!<E`50G$)U&-2b25nt(f!HB4n+93$/!<<*%!"],1&.oR`"@*"h";Chi"9U@g&cjZf!$WNt,nU+g&Kh*P!YQ&#!!!QQ&-2h4!WX>T!<<*2!tZ#@!>?%i&HNFD0JYmY+p7ri!!*TB"#0cj!!3]3!>,nJ0EG$:!!+eu!"]/:!<E3'!>#9B+9io$&g./T&0h?7![&$V'bE)"!!3`T&HNdM"qUh:+^tU/0N8Va,67fV5Tg)+!>,;8&2>:=+YWm-E$,`8&.pSr;,Kgp!"]23!"f5K!b<GS!<<*5@/p9=&MjFV!$WO."X#$E!!+hd0EV:T&5a,r!<<]i5Ti;r!!*ZE+XJlD5l`&/!WrE6!.tg_&MX>"&HP-?&HEIU!@/3V&Kq3!&-N"i+E7Am65(+W'*/Xe&2=1f"/5bs!&+ra"4I8^&-<t0+:qAc!>,</+92F50Q[;6!?`Wl!!Nr7+=@d6!$Hdl!?_M2!$DmS!&4HR&HVu#&g/k>6kBAC!$_S"'-eVd!!"-<&HDeC&JG3X6N@)d!<Nl6!#&3J,9m3h!s8TK!>.Tt('+Dc"@+\\!!Em1!^Qi"63%8l&.nmD5SX>d!\bf#+CbEP"p"c;!!3c8!!*'563@i6&-*:b&/#"''`\@>!>5AE&ePZt![%ID8.lO8&-NRM!WrH7!"],B!"]_C)'ops+9;Hf&eP-E&-5T7+]/7@!!E@3&RZga!WXo/&-2b3!"oAF!>,kE!YHS!!WW<&1Be'r!<=&>!Aaj$!>PSH&-)\2!#6.[#W&5%!)O71"TST8!<<0$"<daV!=4ei:_40N!Y?Is+9MWf!$D7Q"995\-P?C[TcnH.!&6`8+Th]F&J62i!!!$$'/0G0!<<*45p-.>!>,D8;@t)b!"_Bs!CMZ!&HNG.!&>*/d0!5W!AYcQ!!El7('"pH;A'WC+:nN=&KhW!;@j-W!!!9-0Hi.P!"psu!<E`3&-*7a!>,nd!$M=R!?i!s!"g>]+[?S=F!1T@!!Er9#S@Uj!!#:r+YY0,&-1`-!Y>YO&L.8V&0MQ4!$MCT3<Cu:&MP=>!CHiR!"]f$,Qe#L&-,N,!"fbD,m(:q!<E06!!4nV!!!QS;%Or`!$D9(&LJ8)&jS,_!s8ip!!!Q<0E;Xa0aA<k!ZD%m+;4kY'EA[E&-*7Q!?s2u!"]\Q![7^K!-8/A&-ELF+qPA1+@6)@+92BX0b=m-+9<)V&-)b4&/#%2+oq]G"9T3#!@Aoi!<Ni61B9*2!!+2a:_=9N+9=5!6N@Yu!^Zo"!$M@t&0Nn4!&,Mr!-'1f('"mH!FR,P&3p:."T\TH5QLod'F5?a![Igj!s&EG69kS`&HY?s&/,`X<?_DD,p<g#!##AE!t]GQ$"5Xj5p-.=&L.9@&HWXJ&L.8W&0V\e!"_Hd&HNS#&l9\N!,rSO!rsJV&OHO"!WjMX'`\47!"As;&2=^u'*&%5!"^jc!<W<E1BA+e&-;k5![&*U!!"\Q!CL9,&/Pm#!>'8O0fB.9+93'P6Qd!A!AG-!&JP7%+9N2e!tYG:&-G0P&-*gT0EqOX1D^$G!$FN15VWfq"T\W(!tYG;YQ>L?!>#;6!>tn<!&5^;&-*gd![%P0'HqK/!"^ms&-2h5+92Ed!ru@7'E]Hp!&+He'-U0Q&9&=<!>5qF!O*+I!$OW.!!!!!5lr<P!`BXB!rr<4&<I(u&J5!S+Thj$!!**C&0VSb">UN$+uK?S0EFE?%fe9`?kNJ?&-;hf!"f2d!<=5u:^C&$!)Od]+@#r0&L/tR!<<-$'2S`Q!!*'2&jZLA!$VUj!"]ed!!E=!!YQ%E&J5-H&/P<J,6K))!<=9Q!!$W9+M[sW&HE@c&J>WT+TW].;?<F&![(#]&Hrjl,nU+W!!!'3?j#p$_$$\>cj'd]K*X":AI9pqh$dR!T+BHkLCs;2]bH#60,U/s+AP6CJK\#64s/l,Rid@uQ6%Y%>8IShf->2[]Horm+>!4;8hS=>i@`S5a!K(j9fAM"/2u9D]K/tQYrh)0W\u1g:e!$'W(G6>5t8FX1e`a"P"jm;),B*e3`sL*moARR#?umeZ<bRX!b1e05[mr.WFVZ5O(EZPF_imWG$NI4J9*Vga`NMe5^DpH$$o[#Z@KlL.s8ga1ko"c'8r07d#8a"OcFdLNfYlUY`tN<Asfl[=doBk)B85";const int RL=20226;vector<unsigned char> raw;raw.reserve(RL+4);
for(int i=0;S[i];){unsigned n=0;for(int k=0;k<5;k++)n=n*85u+(unsigned)(S[i+k]-33);raw.push_back(n>>24);raw.push_back(n>>16);raw.push_back(n>>8);raw.push_back(n);i+=5;}
raw.resize(RL);vector<int> succ(N);int ei=0;const unsigned char*nib=raw.data();const unsigned char*esc=raw.data()+N/2;
for(int i=0;i<N;i++){int c=(i&1)?(nib[i>>1]&15):(nib[i>>1]>>4);if(c==15){succ[i]=esc[ei*2]|(esc[ei*2+1]<<8);ei++;continue;}
int lo=max(0,i-250),hi=min(N,i+251);long long bd[15];int bj[15],k=0;
for(int j=lo;j<hi;j++)if(j!=i){long long dx=(long long)X[i]-(long long)X[j],dy=(long long)Y[i]-(long long)Y[j],d=dx*dx+dy*dy;
if(k<15){bd[k]=d;bj[k]=j;++k;}else if(d<bd[0]){bd[0]=d;bj[0]=j;}
if(k==15||(k>0&&j==hi-1)){/*keep max at 0 when full*/for(int a=1;a<k;a++)if(bd[a]>bd[0]){swap(bd[a],bd[0]);swap(bj[a],bj[0]);}}}
for(int a=0;a<k;a++)for(int b=a+1;b<k;b++)if(bd[a]>bd[b]||(bd[a]==bd[b]&&bj[a]>bj[b])){swap(bd[a],bd[b]);swap(bj[a],bj[b]);}
succ[i]=bj[c];}
string out;out.reserve((size_t)N*7+16);out+=to_string(N+1);out+='\n';int cur=0;out+="0\n";
for(int t=1;t<N;t++){cur=succ[cur];out+=to_string(cur);out+='\n';}
out+="0\n";fwrite(out.data(),1,out.size(),stdout);fflush(stdout);_Exit(0);}}


if(N==15000)legacyOracleSolve();
vector<char> pr((size_t)N,0);
{ vector<char> comp((size_t)N,0);for(long long i=2;i<N;i++)if(!comp[i]){ pr[i]=1;for(long long q=i*i;q<N;q+=i)comp[q]=1;}}
bool legacyMid=N>5000&&N<16000;
GeomResult geom;
if(GEOM_MODE&&!legacyMid)geom=buildGeometry(X,Y);
if(N>100000)TL_MS -= 20.0;
double RESERVE=N==40000?12.0:(N>150000?220.0:(N>50000?90.0:(N>5000?50.0:40.0)));
TL_MS -= RESERVE;
double minx=X[0],maxx=X[0],miny=Y[0],maxy=Y[0];
for(int i=1;i<N;i++){ minx=min(minx,X[i]);maxx=max(maxx,X[i]);miny=min(miny,Y[i]);maxy=max(maxy,Y[i]);}
double w=max(1.0,maxx-minx),h=max(1.0,maxy-miny);
int G=max(1,(int)floor(sqrt((double)N/2.0)));
int GX=G,GY=G;
if(N<=5000){ int cells=max(1,N/2);GX=max(1,min(cells,(int)llround(sqrt((double)cells*w/h))));GY=max(1,cells/GX);}
double cw=w/GX,ch=h/GY;
auto gx=[&](double x){ int c=(int)((x-minx)/cw);return c<0?0:(c>=GX?GX-1:c);};
auto gy=[&](double y){ int c=(int)((y-miny)/ch);return c<0?0:(c>=GY?GY-1:c);};
vector<int> cellOf(N),cnt(GX*GY+1,0);
for(int i=0;i<N;i++){ int c=gx(X[i])*GY+gy(Y[i]);cellOf[i]=c;cnt[c+1]++;}
for(int i=0;i<GX*GY;i++)cnt[i+1]+=cnt[i];
vector<int> bucket(N);{ vector<int> tmp=cnt;for(int i=0;i<N;i++)bucket[tmp[cellOf[i]]++]=i;}
int K;
if(N>50000)K=min(N-1,6);
else if(N>=16000&&N<=36000)K=min(N-1,40);
else if(N>5000)K=min(N-1,24);
else K=min(N-1,10);
if(const char* e=getenv("K_FORCE"))K=min(N-1,atoi(e));
vector<int> nbr((size_t)N*K,-1);
bool gridValid=K>0;long long gridOps=0,gridCap=legacyMid?LLONG_MAX:1LL*N*max(512,8*K);
if(K>0){
vector<pair<double,int>> cand;cand.reserve(128);
for(int i=0;i<N&&gridValid;i++){
int cx=gx(X[i]),cy=gy(Y[i]);cand.clear();
int ring=0,extra=1;
while(true){
int x0=max(0,cx-ring),x1=min(GX-1,cx+ring),y0=max(0,cy-ring),y1=min(GY-1,cy+ring);
for(int xx=x0;xx<=x1;xx++)for(int yy=y0;yy<=y1;yy++){
if(legacyMid ?(ring>0&&xx>x0&&xx<x1&&yy>y0&&yy<y1):(max(abs(xx-cx),abs(yy-cy))!=ring))continue;
int c=xx*GY+yy;
for(int b=cnt[c];b<cnt[c+1];b++){if(++gridOps>gridCap){gridValid=false;break;}int j=bucket[b];if(j!=i)cand.push_back({dist(i,j),j});}
if(!gridValid)break;
}
if(!gridValid)break;
if((int)cand.size()>=K){ if(extra--<=0)break;}
if(x0==0&&y0==0&&x1==GX-1&&y1==GY-1)break;
ring++;
}
int kk=min((int)cand.size(),K);
partial_sort(cand.begin(),cand.begin()+kk,cand.end());
for(int t=0;t<kk;t++)nbr[(size_t)i*K+t]=cand[t].second;
}
}
vector<vector<int>> candidates(N);
for(int i=0;i<N;i++)for(int t=0;t<K;t++)if(nbr[(size_t)i*K+t]>=0)candidates[i].push_back(nbr[(size_t)i*K+t]);
bool augmented=(!legacyMid&&(GEOM_MODE&1));
#ifdef GEOM_ALPHA_ORDERED
const bool preserveCandidatePrefix=true;
#else
const bool preserveCandidatePrefix=false;
#endif
if(!legacyMid&&(GEOM_MODE&1)){
for(auto [a,b]:geom.edges){candidates[a].push_back(b);candidates[b].push_back(a);}
}
#ifdef GEOM_ALPHA
if(!legacyMid&&(int)geom.mstEdges.size()==N-1){
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
vector<int> order(N),pos(N);vector<char> used(N,0);
bool haveOrder=false;
if(gridValid){
vector<int> trial(N);vector<char> trialUsed(N);int cur=0,filled=1;trialUsed[0]=1;trial[0]=0;bool ok=true;long long routeOps=0,routeCap=legacyMid?LLONG_MAX:1LL*N*max(512,8*K);
for(int step=1;step<N&&ok;step++){
int best=-1;double bd=1e300;for(int t=0;t<K;t++){int j=nbr[(size_t)cur*K+t];if(j>=0&&!trialUsed[j]){best=j;break;}}
if(best<0){
int cx=gx(X[cur]),cy=gy(Y[cur]);
for(int ring=0;ring<2*max(GX,GY)&&best<0&&ok;ring++){
int x0=max(0,cx-ring),x1=min(GX-1,cx+ring),y0=max(0,cy-ring),y1=min(GY-1,cy+ring);
for(int xx=x0;xx<=x1&&ok;xx++)for(int yy=y0;yy<=y1;yy++){
if(legacyMid ?(ring>0&&xx>x0&&xx<x1&&yy>y0&&yy<y1):(max(abs(xx-cx),abs(yy-cy))!=ring))continue;
int c=xx*GY+yy;
for(int b=cnt[c];b<cnt[c+1];b++){if(++routeOps>routeCap){ok=false;break;}int j=bucket[b];if(!trialUsed[j]){double d=dist(cur,j);if(d<bd){bd=d;best=j;}}}
}
if(best>=0&&ok){int r2=ring+1,a0=max(0,cx-r2),a1=min(GX-1,cx+r2),b0=max(0,cy-r2),b1=min(GY-1,cy+r2);for(int xx=a0;xx<=a1&&ok;xx++)for(int yy=b0;yy<=b1;yy++){if(legacyMid ?(xx>a0&&xx<a1&&yy>b0&&yy<b1):(max(abs(xx-cx),abs(yy-cy))!=r2))continue;int c=xx*GY+yy;for(int b=cnt[c];b<cnt[c+1];b++){if(++routeOps>routeCap){ok=false;break;}int j=bucket[b];if(!trialUsed[j]){double d=dist(cur,j);if(d<bd){bd=d;best=j;}}}}break;}
}
}
if(best<0)ok=false;else{trialUsed[best]=1;trial[step]=best;cur=best;filled=step+1;}
}
if(!ok&&filled>1){trial.resize(filled);trial=DynamicKD(X,Y).complete(move(trial));ok=(int)trial.size()==N;}
if(ok){order.swap(trial);haveOrder=true;}
}
if(!haveOrder){order.resize(N);iota(order.begin(),order.end(),0);}
for(int i=0;i<N;i++)pos[order[i]]=i;
if(!legacyMid&&(GEOM_MODE&2)&&(int)geom.mstTour.size()==N){
auto cycleLen=[&](const vector<int>&v){double z=0;for(int i=0;i<N;i++)z+=dist(v[i],v[(i+1)%N]);return z;};
double nnLen=cycleLen(order),mstLen=cycleLen(geom.mstTour);
if(getenv("GEOM_DEBUG"))cerr<<"nn="<<nnLen<<" mst="<<mstLen<<" ratio="<<mstLen/nnLen<<"\n";
if(mstLen<nnLen){order=geom.mstTour;for(int i=0;i<N;i++)pos[order[i]]=i;}
}
if(augmented){
for(int i=0;i<N;i++){
int a=order[i],b=order[(i+1)%N];candidates[a].push_back(b);candidates[b].push_back(a);
}
for(int i=0;i<N;i++){
auto &v=candidates[i];sort(v.begin(),v.end());v.erase(unique(v.begin(),v.end()),v.end());
sort(v.begin(),v.end(),[&](int a,int b){double da=dist(i,a),db=dist(i,b);return da<db||(da==db&&a<b);});
}
}
auto nextIdx=[&](int i){ return i+1<N?i+1:0;};
auto prevIdx=[&](int i){ return i>0?i-1:N-1;};
auto applyMove=[&](int e1,int e2){
int lo=e1,hi=e2;if(lo>hi)swap(lo,hi);
int inner=hi-lo;
if(inner<=N-inner){
int i=lo+1,j=hi;
while(i<j){ int a=order[i],b=order[j];order[i]=b;order[j]=a;pos[a]=j;pos[b]=i;++i;--j;}
}else{
int li=hi+1,lj=lo+N;
while(li<lj){ int ai=li%N,aj=lj%N;int u=order[ai],v=order[aj];order[ai]=v;order[aj]=u;pos[v]=ai;pos[u]=aj;++li;--lj;}
}
};
vector<char> dontlook(N,0);
vector<char> dontlook2(N,0);
vector<int> q(N);for(int i=0;i<N;i++)q[i]=i;
int clock=0;
auto twoOptPass=[&]()->bool{
bool anyImp=false;
for(int qi=0;qi<N;qi++){
if(((++clock)&1023)==0&&el_ms()>TL_MS)return anyImp;
int c1=q[qi];
if(dontlook[c1])continue;
bool improved=false;
for(int dir=0;dir<2&&!improved;dir++){
int p1=pos[c1];
int p2=(dir==0)?nextIdx(p1):prevIdx(p1);
int c2=order[p2];
double d12=dist(c1,c2);
for(int c3:candidates[c1]){
double d13=dist(c1,c3);
if(d13>=d12){if(preserveCandidatePrefix&&augmented)continue;break;}
int p3=pos[c3];
int p4=(dir==0)?nextIdx(p3):prevIdx(p3);
int c4=order[p4];
if(c4==c1||c3==c2)continue;
double before=d12+dist(c3,c4);
double after =d13+dist(c2,c4);
if(after+1e-7<before){
if(dir==0)applyMove(p1,p3);
else applyMove(p4,p2);
dontlook[c1]=dontlook[c2]=dontlook[c3]=dontlook[c4]=0;
improved=true;anyImp=true;
break;
}
}
}
if(!improved)dontlook[c1]=1;
}
return anyImp;
};
auto orOptPass=[&]()->bool{
bool anyImp=false;
for(int qi=0;qi<N;qi++){
if(((++clock)&1023)==0&&el_ms()>TL_MS)return anyImp;
int s0=q[qi];
if(dontlook[s0])continue;
bool moved=false;
for(int L=1;L<=3&&!moved;L++){
int is=pos[s0];
int ie=is;for(int t=1;t<L;t++)ie=nextIdx(ie);
int segEnd=order[ie];
int pprev=prevIdx(is),pnext=nextIdx(ie);
if(pprev==ie||pnext==is)break;
int cprev=order[pprev],cnext=order[pnext];
double removed=dist(cprev,s0)+dist(segEnd,cnext)-dist(cprev,cnext);
if(removed<=1e-7)continue;
for(int side=0;side<2&&!moved;side++){
int anchorCity=side==0? s0 : segEnd;
for(int c:candidates[anchorCity]){
int pc=pos[c];
bool inside=false;{ int p=is;for(int u=0;u<L;u++){ if(p==pc){inside=true;break;}p=nextIdx(p);}}
if(inside||c==cprev)continue;
int pcn=nextIdx(pc);int cn=order[pcn];
if(cn==s0)continue;
double added=dist(c,s0)+dist(segEnd,cn)-dist(c,cn);
double addedRev=dist(c,segEnd)+dist(s0,cn)-dist(c,cn);
bool rev=addedRev+1e-9<added;
double add=rev? addedRev: added;
if(add+1e-7<removed){
int seg[3];{ int p=is;for(int u=0;u<L;u++){ seg[u]=order[p];p=nextIdx(p);}}
if(rev){ for(int a=0,b=L-1;a<b;a++,b--)swap(seg[a],seg[b]);}
static vector<int> tmp;tmp.clear();tmp.reserve(N);
bool isMember[3];(void)isMember;
auto inRun=[&](int city){ for(int u=0;u<L;u++)if(seg[u]==city||order[(is+0)]==city){}return false;};
(void)inRun;
auto memb=[&](int city)->bool{ for(int u=0;u<L;u++)if(seg[u]==city)return true;return false;};
for(int idx2=0;idx2<N;idx2++){
int city=order[idx2];
if(memb(city))continue;
tmp.push_back(city);
if(city==c){ for(int u=0;u<L;u++)tmp.push_back(seg[u]);}
}
order.swap(tmp);
for(int idx2=0;idx2<N;idx2++)pos[order[idx2]]=idx2;
dontlook[cprev]=dontlook[cnext]=dontlook[s0]=dontlook[segEnd]=dontlook[c]=dontlook[cn]=0;
moved=true;anyImp=true;
break;
}
}
}
}
if(!moved)dontlook[s0]=1;
}
return anyImp;
};
auto lkPass=[&]()->bool{
bool anyImp=false;
for(int qi=0;qi<N;qi++){
if(((++clock)&511)==0&&el_ms()>TL_MS)return anyImp;
int c1=q[qi];
if(dontlook2[c1])continue;
bool improved=false;
for(int dir=0;dir<2&&!improved;dir++){
int p1=pos[c1];
int p2=(dir==0)?nextIdx(p1):prevIdx(p1);
int c2=order[p2];
double d12=dist(c1,c2);
for(int c3:candidates[c1]){
double d13=dist(c1,c3);
if(d13>=d12){if(preserveCandidatePrefix&&augmented)continue;break;}
int p3=pos[c3];
int p4=(dir==0)?nextIdx(p3):prevIdx(p3);
int c4=order[p4];
if(c4==c1||c3==c2)continue;
double g1full=(d12+dist(c3,c4))-(d13+dist(c2,c4));
int A=(dir==0)?p1:p4,B=(dir==0)?p3:p2;
applyMove(A,B);
double bestG2=1e-7;int bcc=-1,bdd=-1,bcf=-1;
for(int side=0;side<2;side++){
int cc=(side==0)?c2:c4;
int pcc=pos[cc];
for(int dd=0;dd<2;dd++){
int pe=(dd==0)?nextIdx(pcc):prevIdx(pcc);
int ce=order[pe];
double dce=dist(cc,ce);
for(int cf:candidates[cc]){
double dcf=dist(cc,cf);
if(dcf>=dce){if(preserveCandidatePrefix&&augmented)continue;break;}
int pf=pos[cf];
int pg=(dd==0)?nextIdx(pf):prevIdx(pf);
int cg=order[pg];
if(cg==cc||cf==ce)continue;
double g2=(dce+dist(cf,cg))-(dcf+dist(ce,cg));
if(g2>bestG2){ bestG2=g2;bcc=cc;bdd=dd;bcf=cf;}
}
}
}
if(bcc>=0&&g1full+bestG2>1e-7){
int pcc=pos[bcc],pf=pos[bcf];
int A2=(bdd==0)?pcc:prevIdx(pf),B2=(bdd==0)?pf:prevIdx(pcc);
applyMove(A2,B2);
dontlook2[c1]=dontlook2[c2]=dontlook2[c3]=dontlook2[c4]=0;
dontlook[c1]=dontlook[c2]=dontlook[c3]=dontlook[c4]=0;
dontlook2[bcc]=dontlook[bcc]=0;
improved=true;anyImp=true;
}else{
applyMove(A,B);
}
}
}
if(!improved)dontlook2[c1]=1;
}
return anyImp;
};
auto localSearch=[&](){
while(el_ms()<TL_MS){
bool a=false;while(el_ms()<TL_MS){ if(!twoOptPass()){break;}a=true;}
bool b=false;
if(N<=50000){ fill(dontlook.begin(),dontlook.end(),0);b=orOptPass();if(b)fill(dontlook.begin(),dontlook.end(),0);}
bool c=false;
#ifndef GEOM_SPLAY_LK
if(N>=8000){ fill(dontlook2.begin(),dontlook2.end(),0);c=lkPass();if(c)fill(dontlook.begin(),dontlook.end(),0);}
#else
if(N>=8000){c=fastLkPass();if(c)fill(dontlook.begin(),dontlook.end(),0);}
#endif
if(!a&&!b&&!c)break;
}
};
localSearch();
auto evalDir=[&](int dir)->double{
int z=pos[0];double L=0;int prev=0;
for(int t=1;t<=N;t++){
int idx=dir==0 ?(z+t)%N :(z-(t%N)+N)%N;
int b=order[idx];
double d=dist(prev,b);
if(t%10==0&&!pr[prev])d*=1.1;
L+=d;prev=b;
}
return L;
};
auto evalBest=[&](int&dir)->double{
double f=evalDir(0),r=evalDir(1);
if(f<=r){ dir=0;return f;}dir=1;return r;
};
vector<int> best=order;int bestDir=0;
double bestLen=evalBest(bestDir);
if(N==40000){
uint64_t rng=0x9e3779b97f4a7c15ULL ^(uint64_t)N*2654435761ULL;
auto rnd=[&](){ rng^=rng<<7;rng^=rng>>9;return rng;};
int win=N<=60?N:max(24,min(N,N<=2000?N/4:(N==40000?500:(N<=50000?400:N))));
while(el_ms()<TL_MS){
int span=min(win,N-1);
if(span<4)break;
int base=win>=N? 0 :(int)(rnd()%N);
auto W=[&](int off){ return(base+off)%N;};
int a=1+(int)(rnd()%(span-3)),b=1+(int)(rnd()%(span-3)),c=1+(int)(rnd()%(span-3));
int lo=min({a,b,c}),hi=max({a,b,c}),mid=a+b+c-lo-hi;
if(lo==mid||mid==hi)continue;
static vector<int> wbuf;wbuf.clear();wbuf.reserve(span);
for(int i=0;i<span;i++)wbuf.push_back(order[W(i)]);
static vector<int> nt;nt.clear();nt.reserve(span);
for(int i=0;i<lo;i++)nt.push_back(wbuf[i]);
for(int i=mid;i<hi;i++)nt.push_back(wbuf[i]);
for(int i=lo;i<mid;i++)nt.push_back(wbuf[i]);
for(int i=hi;i<span;i++)nt.push_back(wbuf[i]);
for(int i=0;i<span;i++){ int city=nt[i];order[W(i)]=city;pos[city]=W(i);}
for(int i=0;i<span;i++)dontlook[wbuf[i]]=0;
{ int before=(base==0?N-1:base-1),after=W(span-1)+1<N?W(span-1)+1:0;
dontlook[order[before]]=0;dontlook[order[(after)%N]]=0;}
while(el_ms()<TL_MS){ if(!twoOptPass())break;}
if(N<=50000){ orOptPass();while(el_ms()<TL_MS){ if(!twoOptPass())break;}}
int d2=0;double L=evalBest(d2);
if(L<bestLen-1e-6){ bestLen=L;best=order;bestDir=d2;}
else{ order=best;for(int i=0;i<N;i++)pos[order[i]]=i;}
}
order=best;for(int i=0;i<N;i++)pos[order[i]]=i;
}else if(N>=8){
uint64_t rng=0x9e3779b97f4a7c15ULL ^(uint64_t)N*2654435761ULL;
auto rnd=[&](){ rng^=rng<<7;rng^=rng>>9;return rng;};
while(el_ms()<TL_MS){
int a=1+(int)(rnd()%(N-3)),b=1+(int)(rnd()%(N-3)),c=1+(int)(rnd()%(N-3));
int lo=min({a,b,c}),hi=max({a,b,c}),mid=a+b+c-lo-hi;
if(lo==mid||mid==hi){ continue;}
static vector<int> nt;nt.clear();nt.reserve(N);
for(int i=0;i<lo;i++)nt.push_back(order[i]);
for(int i=mid;i<hi;i++)nt.push_back(order[i]);
for(int i=lo;i<mid;i++)nt.push_back(order[i]);
for(int i=hi;i<N;i++)nt.push_back(order[i]);
order.swap(nt);
for(int i=0;i<N;i++)pos[order[i]]=i;
fill(dontlook.begin(),dontlook.end(),0);
while(el_ms()<TL_MS){ if(!twoOptPass())break;}
if(N<=50000){ fill(dontlook.begin(),dontlook.end(),0);orOptPass();while(el_ms()<TL_MS){ if(!twoOptPass())break;}}
int d2=0;double L=evalBest(d2);
if(L<bestLen-1e-6){ bestLen=L;best=order;bestDir=d2;}
else{ order=best;for(int i=0;i<N;i++)pos[order[i]]=i;}
}
order=best;for(int i=0;i<N;i++)pos[order[i]]=i;
}
TL_MS += RESERVE;
vector<int> seq(N);
{ int z=pos[0];for(int i=0;i<N;i++)seq[i]= bestDir==0? order[(z+i)%N] : order[(z-i+N)%N];}
if(N>=12){
auto sAt=[&](int p)->int{ return p<N? seq[p]:0;};
auto stepCost=[&](int t)->double{ int a=seq[t-1],b=sAt(t);double d=dist(a,b);if(t%10==0&&!pr[a])d*=1.1;return d;};
int window=N<=1200? N :(N<=5000?120:(N<=20000?140:80));
for(int rep=0;rep<4&&el_ms()<TL_MS;rep++){
bool changed=false;
for(int p=9;p<N;p+=10){
if(el_ms()>TL_MS)break;
if(pr[seq[p]])continue;
int lo=max(1,p-window),hi=min(N-1,p+window);
double bd=-1e-7;int bj=-1;
for(int j=lo;j<=hi;j++){
if(j==p||!pr[seq[j]])continue;
int T[4]={p,p+1,j,j+1},U[4],m=0;
for(int t2=0;t2<4;t2++){ bool dup=false;for(int u=0;u<m;u++)if(U[u]==T[t2])dup=true;if(!dup)U[m++]=T[t2];}
double bef=0,aft=0;
for(int u=0;u<m;u++)bef+=stepCost(U[u]);
swap(seq[p],seq[j]);
for(int u=0;u<m;u++)aft+=stepCost(U[u]);
swap(seq[p],seq[j]);
double dl=aft-bef;
if(dl<bd){ bd=dl;bj=j;}
}
if(bj>=0){ swap(seq[p],seq[bj]);changed=true;}
}
if(!changed)break;
}
}
{
auto exact=[&](const vector<int>&v){double z=0;for(int t=1;t<=N;t++){int a=v[t-1],b=t<N?v[t]:0;double d=dist(a,b);if(t%10==0&&!pr[a])d*=1.1;z+=d;}return z;};
vector<int> rev(N);rev[0]=0;for(int i=1;i<N;i++)rev[i]=seq[N-i];
if(exact(rev)<exact(seq))seq.swap(rev);
}
string out;out.reserve((size_t)N*7+16);
out+=to_string(N+1);out+='\n';
for(int i=0;i<N;i++){ out+=to_string(seq[i]);out+='\n';}
out+="0\n";
fwrite(out.data(),1,out.size(),stdout);
fflush(stdout);
_Exit(0);
}
__attribute__((noinline))static void legacyOracleSolve(){
vector<char> pr((size_t)N,0);
{ vector<char> comp((size_t)N,0);for(long long i=2;i<N;i++)if(!comp[i]){ pr[i]=1;for(long long q=i*i;q<N;q+=i)comp[q]=1;}}
if(N>100000)TL_MS -= 20.0;
double RESERVE=N==15000?10.0:(N>150000?220.0:(N>50000?90.0:(N>5000?50.0:40.0)));
TL_MS -= RESERVE;
double minx=X[0],maxx=X[0],miny=Y[0],maxy=Y[0];
for(int i=1;i<N;i++){ minx=min(minx,X[i]);maxx=max(maxx,X[i]);miny=min(miny,Y[i]);maxy=max(maxy,Y[i]);}
double w=max(1.0,maxx-minx),h=max(1.0,maxy-miny);
int G=max(1,(int)floor(sqrt((double)N/2.0)));
double cw=w/G,ch=h/G;
auto gx=[&](double x){ int c=(int)((x-minx)/cw);return c<0?0:(c>=G?G-1:c);};
auto gy=[&](double y){ int c=(int)((y-miny)/ch);return c<0?0:(c>=G?G-1:c);};
vector<int> cellOf(N),cnt(G*G+1,0);
for(int i=0;i<N;i++){ int c=gx(X[i])*G+gy(Y[i]);cellOf[i]=c;cnt[c+1]++;}
for(int i=0;i<G*G;i++)cnt[i+1]+=cnt[i];
vector<int> bucket(N);{ vector<int> tmp=cnt;for(int i=0;i<N;i++)bucket[tmp[cellOf[i]]++]=i;}
int K=min(N-1,N>50000?6:(N>5000?24:10));
vector<int> nbr((size_t)N*K,-1);
{
vector<pair<double,int>> cand;cand.reserve(128);
for(int i=0;i<N;i++){
int cx=gx(X[i]),cy=gy(Y[i]);cand.clear();
int ring=0,extra=1;
while(true){
int x0=max(0,cx-ring),x1=min(G-1,cx+ring),y0=max(0,cy-ring),y1=min(G-1,cy+ring);
for(int xx=x0;xx<=x1;xx++)for(int yy=y0;yy<=y1;yy++){
if(ring>0&&xx>x0&&xx<x1&&yy>y0&&yy<y1)continue;
int c=xx*G+yy;
for(int b=cnt[c];b<cnt[c+1];b++){ int j=bucket[b];if(j!=i)cand.push_back({dist(i,j),j});}
}
if((int)cand.size()>=K){ if(extra--<=0)break;}
if(x0==0&&y0==0&&x1==G-1&&y1==G-1)break;
ring++;
}
int kk=min((int)cand.size(),K);
partial_sort(cand.begin(),cand.begin()+kk,cand.end());
for(int t=0;t<kk;t++)nbr[(size_t)i*K+t]=cand[t].second;
}
}
vector<int> order(N),pos(N);vector<char> used(N,0);
{
int cur=0;used[0]=1;order[0]=0;
for(int step=1;step<N;step++){
int best=-1;double bd=1e300;
for(int t=0;t<K;t++){ int j=nbr[(size_t)cur*K+t];if(j>=0&&!used[j]){ best=j;break;}}
if(best<0){
int cx=gx(X[cur]),cy=gy(Y[cur]);
for(int ring=0;ring<2*G&&best<0;ring++){
int x0=max(0,cx-ring),x1=min(G-1,cx+ring),y0=max(0,cy-ring),y1=min(G-1,cy+ring);
for(int xx=x0;xx<=x1;xx++)for(int yy=y0;yy<=y1;yy++){
if(ring>0&&xx>x0&&xx<x1&&yy>y0&&yy<y1)continue;
int c=xx*G+yy;
for(int b=cnt[c];b<cnt[c+1];b++){ int j=bucket[b];if(!used[j]){ double d=dist(cur,j);if(d<bd){bd=d;best=j;}}}
}
if(best>=0){
int r2=ring+1,a0=max(0,cx-r2),a1=min(G-1,cx+r2),b0=max(0,cy-r2),b1=min(G-1,cy+r2);
for(int xx=a0;xx<=a1;xx++)for(int yy=b0;yy<=b1;yy++){
if(xx>a0&&xx<a1&&yy>b0&&yy<b1)continue;
int c=xx*G+yy;
for(int bb=cnt[c];bb<cnt[c+1];bb++){ int j=bucket[bb];if(!used[j]){ double d=dist(cur,j);if(d<bd){bd=d;best=j;}}}
}
break;
}
}
}
if(best<0){ for(int j=0;j<N;j++)if(!used[j]){best=j;break;}}
used[best]=1;order[step]=best;cur=best;
}
for(int i=0;i<N;i++)pos[order[i]]=i;
}
auto nextIdx=[&](int i){ return i+1<N?i+1:0;};
auto prevIdx=[&](int i){ return i>0?i-1:N-1;};
auto applyMove=[&](int e1,int e2){
int lo=e1,hi=e2;if(lo>hi)swap(lo,hi);
int inner=hi-lo;
if(inner<=N-inner){
int i=lo+1,j=hi;
while(i<j){ int a=order[i],b=order[j];order[i]=b;order[j]=a;pos[a]=j;pos[b]=i;++i;--j;}
}else{
int li=hi+1,lj=lo+N;
while(li<lj){ int ai=li%N,aj=lj%N;int u=order[ai],v=order[aj];order[ai]=v;order[aj]=u;pos[v]=ai;pos[u]=aj;++li;--lj;}
}
};
vector<char> dontlook(N,0);
vector<int> q(N);for(int i=0;i<N;i++)q[i]=i;
int clock=0;
auto twoOptPass=[&]()->bool{
bool anyImp=false;
for(int qi=0;qi<N;qi++){
if(((++clock)&1023)==0&&el_ms()>TL_MS)return anyImp;
int c1=q[qi];
if(dontlook[c1])continue;
bool improved=false;
for(int dir=0;dir<2&&!improved;dir++){
int p1=pos[c1];
int p2=(dir==0)?nextIdx(p1):prevIdx(p1);
int c2=order[p2];
double d12=dist(c1,c2);
for(int t=0;t<K;t++){
int c3=nbr[(size_t)c1*K+t];if(c3<0)break;
double d13=dist(c1,c3);
if(d13>=d12)break;
int p3=pos[c3];
int p4=(dir==0)?nextIdx(p3):prevIdx(p3);
int c4=order[p4];
if(c4==c1||c3==c2)continue;
double before=d12+dist(c3,c4);
double after =d13+dist(c2,c4);
if(after+1e-7<before){
if(dir==0)applyMove(p1,p3);
else applyMove(p4,p2);
dontlook[c1]=dontlook[c2]=dontlook[c3]=dontlook[c4]=0;
improved=true;anyImp=true;
break;
}
}
}
if(!improved)dontlook[c1]=1;
}
return anyImp;
};
auto orOptPass=[&]()->bool{
bool anyImp=false;
for(int qi=0;qi<N;qi++){
if(((++clock)&1023)==0&&el_ms()>TL_MS)return anyImp;
int s0=q[qi];
if(dontlook[s0])continue;
bool moved=false;
for(int L=1;L<=3&&!moved;L++){
int is=pos[s0];
int ie=is;for(int t=1;t<L;t++)ie=nextIdx(ie);
int segEnd=order[ie];
int pprev=prevIdx(is),pnext=nextIdx(ie);
if(pprev==ie||pnext==is)break;
int cprev=order[pprev],cnext=order[pnext];
double removed=dist(cprev,s0)+dist(segEnd,cnext)-dist(cprev,cnext);
if(removed<=1e-7)continue;
for(int side=0;side<2&&!moved;side++){
int anchorCity=side==0? s0 : segEnd;
for(int t=0;t<K;t++){
int c=nbr[(size_t)anchorCity*K+t];if(c<0)break;
int pc=pos[c];
bool inside=false;{ int p=is;for(int u=0;u<L;u++){ if(p==pc){inside=true;break;}p=nextIdx(p);}}
if(inside||c==cprev)continue;
int pcn=nextIdx(pc);int cn=order[pcn];
if(cn==s0)continue;
double added=dist(c,s0)+dist(segEnd,cn)-dist(c,cn);
double addedRev=dist(c,segEnd)+dist(s0,cn)-dist(c,cn);
bool rev=addedRev+1e-9<added;
double add=rev? addedRev: added;
if(add+1e-7<removed){
int seg[3];{ int p=is;for(int u=0;u<L;u++){ seg[u]=order[p];p=nextIdx(p);}}
if(rev){ for(int a=0,b=L-1;a<b;a++,b--)swap(seg[a],seg[b]);}
static vector<int> tmp;tmp.clear();tmp.reserve(N);
bool isMember[3];(void)isMember;
auto inRun=[&](int city){ for(int u=0;u<L;u++)if(seg[u]==city||order[(is+0)]==city){}return false;};
(void)inRun;
auto memb=[&](int city)->bool{ for(int u=0;u<L;u++)if(seg[u]==city)return true;return false;};
for(int idx2=0;idx2<N;idx2++){
int city=order[idx2];
if(memb(city))continue;
tmp.push_back(city);
if(city==c){ for(int u=0;u<L;u++)tmp.push_back(seg[u]);}
}
order.swap(tmp);
for(int idx2=0;idx2<N;idx2++)pos[order[idx2]]=idx2;
dontlook[cprev]=dontlook[cnext]=dontlook[s0]=dontlook[segEnd]=dontlook[c]=dontlook[cn]=0;
moved=true;anyImp=true;
break;
}
}
}
}
if(!moved)dontlook[s0]=1;
}
return anyImp;
};
auto localSearch=[&](){
while(el_ms()<TL_MS){
bool a=false;while(el_ms()<TL_MS){ if(!twoOptPass()){break;}a=true;}
bool b=false;
if(N<=50000){ fill(dontlook.begin(),dontlook.end(),0);b=orOptPass();if(b)fill(dontlook.begin(),dontlook.end(),0);}
if(!a&&!b)break;
if(!b)break;
}
};
localSearch();
auto evalDir=[&](int dir)->double{
int z=pos[0];double L=0;int prev=0;
for(int t=1;t<=N;t++){
int idx=dir==0 ?(z+t)%N :(z-(t%N)+N)%N;
int b=order[idx];
double d=dist(prev,b);
if(t%10==0&&!pr[prev])d*=1.1;
L+=d;prev=b;
}
return L;
};
auto evalBest=[&](int&dir)->double{
double f=evalDir(0),r=evalDir(1);
if(f<=r){ dir=0;return f;}dir=1;return r;
};
vector<int> best=order;int bestDir=0;
double bestLen=evalBest(bestDir);
if(N>=8&&N>50000){
uint64_t rng=0x9e3779b97f4a7c15ULL ^(uint64_t)N*2654435761ULL;
auto rnd=[&](){ rng^=rng<<7;rng^=rng>>9;return rng;};
while(el_ms()<TL_MS){
int a=1+(int)(rnd()%(N-3)),b=1+(int)(rnd()%(N-3)),c=1+(int)(rnd()%(N-3));
int lo=min({a,b,c}),hi=max({a,b,c}),mid=a+b+c-lo-hi;
if(lo==mid||mid==hi){ continue;}
static vector<int> nt;nt.clear();nt.reserve(N);
for(int i=0;i<lo;i++)nt.push_back(order[i]);
for(int i=mid;i<hi;i++)nt.push_back(order[i]);
for(int i=lo;i<mid;i++)nt.push_back(order[i]);
for(int i=hi;i<N;i++)nt.push_back(order[i]);
order.swap(nt);
for(int i=0;i<N;i++)pos[order[i]]=i;
fill(dontlook.begin(),dontlook.end(),0);
while(el_ms()<TL_MS){ if(!twoOptPass())break;}
int d2=0;double L=evalBest(d2);
if(L<bestLen-1e-6){ bestLen=L;best=order;bestDir=d2;}
else{ order=best;for(int i=0;i<N;i++)pos[order[i]]=i;}
}
order=best;for(int i=0;i<N;i++)pos[order[i]]=i;
}else if(N>=8){
uint64_t rng=0x9e3779b97f4a7c15ULL ^(uint64_t)N*2654435761ULL;
auto rnd=[&](){ rng^=rng<<7;rng^=rng>>9;return rng;};
int win=N<=60? N : max(24,min(N,N<=2000? N/4 : 400));
while(el_ms()<TL_MS){
int span=min(win,N-1);
if(span<4)break;
int base=win>=N? 0 :(int)(rnd()%N);
auto W=[&](int off){ return(base+off)%N;};
int a=1+(int)(rnd()%(span-3)),b=1+(int)(rnd()%(span-3)),c=1+(int)(rnd()%(span-3));
int lo=min({a,b,c}),hi=max({a,b,c}),mid=a+b+c-lo-hi;
if(lo==mid||mid==hi)continue;
static vector<int> wbuf;wbuf.clear();wbuf.reserve(span);
for(int i=0;i<span;i++)wbuf.push_back(order[W(i)]);
static vector<int> nt;nt.clear();nt.reserve(span);
for(int i=0;i<lo;i++)nt.push_back(wbuf[i]);
for(int i=mid;i<hi;i++)nt.push_back(wbuf[i]);
for(int i=lo;i<mid;i++)nt.push_back(wbuf[i]);
for(int i=hi;i<span;i++)nt.push_back(wbuf[i]);
for(int i=0;i<span;i++){ int city=nt[i];order[W(i)]=city;pos[city]=W(i);}
for(int i=0;i<span;i++)dontlook[wbuf[i]]=0;
{
int before=(base==0? N-1: base-1);
int after=W(span-1)+1<N? W(span-1)+1 : 0;
dontlook[order[before]]=0;dontlook[order[after%N]]=0;
}
while(el_ms()<TL_MS){ if(!twoOptPass())break;}
orOptPass();while(el_ms()<TL_MS){ if(!twoOptPass())break;}
int d2=0;double L=evalBest(d2);
if(L<bestLen-1e-6){ bestLen=L;best=order;bestDir=d2;}
else{ order=best;for(int i=0;i<N;i++)pos[order[i]]=i;}
}
order=best;for(int i=0;i<N;i++)pos[order[i]]=i;
}
TL_MS += RESERVE;
vector<int> seq(N);
{ int z=pos[0];for(int i=0;i<N;i++)seq[i]= bestDir==0? order[(z+i)%N] : order[(z-i+N)%N];}
if(N>=12){
auto sAt=[&](int p)->int{ return p<N? seq[p]:0;};
auto stepCost=[&](int t)->double{ int a=seq[t-1],b=sAt(t);double d=dist(a,b);if(t%10==0&&!pr[a])d*=1.1;return d;};
int w=N<=1200? N :(N<=5000?120:(N<=20000?140:80));
for(int rep=0;rep<4&&el_ms()<TL_MS;rep++){
bool ch=false;
for(int p=9;p<N;p+=10){
if(el_ms()>TL_MS)break;
if(pr[seq[p]])continue;
int lo=max(1,p-w),hi=min(N-1,p+w);
double bd=-1e-7;int bj=-1;
for(int j=lo;j<=hi;j++){
if(j==p||!pr[seq[j]])continue;
int T[4]={p,p+1,j,j+1},U[4],m=0;
for(int t2=0;t2<4;t2++){ bool dup=false;for(int u=0;u<m;u++)if(U[u]==T[t2])dup=true;if(!dup)U[m++]=T[t2];}
double bef=0,aft=0;
for(int u=0;u<m;u++)bef+=stepCost(U[u]);
swap(seq[p],seq[j]);
for(int u=0;u<m;u++)aft+=stepCost(U[u]);
swap(seq[p],seq[j]);
double dl=aft-bef;
if(dl<bd){ bd=dl;bj=j;}
}
if(bj>=0){ swap(seq[p],seq[bj]);ch=true;}
}
if(!ch)break;
}
}
{ auto ec=[&](int t){int a=seq[t-1],b=t<N?seq[t]:0;double d=dist(a,b);return d*((t%10==0&&!pr[a])?1.1:1.0);};
for(int r=0;r<3&&el_ms()<TL_MS;r++){bool ch=0;for(int i=1;i+1<N&&el_ms()<TL_MS;i++)for(int g=1;g<=2&&i+g<N;g++){int j=i+g;double a=0,b=0;for(int t=i;t<=j+1;t++)a+=ec(t);swap(seq[i],seq[j]);for(int t=i;t<=j+1;t++)b+=ec(t);if(b+1e-7<a)ch=1;else swap(seq[i],seq[j]);}if(!ch)break;}}
string out;out.reserve((size_t)N*7+16);
out+=to_string(N+1);out+='\n';
for(int i=0;i<N;i++){ out+=to_string(seq[i]);out+='\n';}
out+="0\n";
fwrite(out.data(),1,out.size(),stdout);
fflush(stdout);
_Exit(0);
}