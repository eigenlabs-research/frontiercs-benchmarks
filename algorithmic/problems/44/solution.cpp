// Model South Core
#include <bits/stdc++.h>
using namespace std;
static chrono::steady_clock::time_point T0;
static double TL_MS=2400.0;
static inline double el_ms(){ return chrono::duration<double,milli>(chrono::steady_clock::now()-T0).count(); }
static int N;
static vector<double> X,Y,XS,YS,XW,YW;
static double baseSc=1.0,curSc=1.0,acx=0,acy=0;
static inline double distOrig(int a,int b){ double dx=X[a]-X[b],dy=Y[a]-Y[b]; return sqrt(dx*dx+dy*dy); }
static inline double dist(int a,int b){ double dx=XW[a]-XW[b],dy=YW[a]-YW[b]; return sqrt(dx*dx+dy*dy); }
static vector<char> inbuf; static size_t ip=0;
static long long readLL(){
    while(ip<inbuf.size()&&(inbuf[ip]<'0'||inbuf[ip]>'9')&&inbuf[ip]!='-') ip++;
    bool neg=false; if(ip<inbuf.size()&&inbuf[ip]=='-'){neg=true;ip++;}
    long long v=0; while(ip<inbuf.size()&&inbuf[ip]>='0'&&inbuf[ip]<='9'){v=v*10+(inbuf[ip]-'0');ip++;}
    return neg?-v:v;
}
static uint64_t rngS=0x9e3779b97f4a7c15ULL;
static inline uint64_t rng(){ rngS^=rngS<<7; rngS^=rngS>>9; return rngS; }
static vector<char> rigid;
static void detectRigid(const vector<int>& order, const vector<int>& nbr, int K){
    rigid.assign(N,0); if(N<4)return;
    vector<int> nnOf(N,-1); vector<double> nnD(N,1e100);
    for(int i=0;i<N;i++) for(int t=0;t<K;t++){ int j=nbr[(size_t)i*K+t]; if(j<0)break; double d=dist(i,j); if(d<nnD[i]){nnD[i]=d;nnOf[i]=j;} }
    for(int p=0;p<N;p++){ int a=order[p],b=order[(p+1)%N]; if(nnOf[a]==b&&nnOf[b]==a) rigid[p]=1; }
}
static inline bool rigBrk(int a,int b,const vector<int>& o,const vector<int>& p,const vector<char>& r){
    int pa=p[a],pb=p[b],na=(pa+1)%N; if(o[na]==b&&r[pa])return 1;
    int pra=(pa-1+N)%N; if(o[pra]==b&&r[pra])return 1; return 0;
}
static vector<double> eHeat;
static void injHeat(const vector<int>& order, double frac){
    fill(eHeat.begin(),eHeat.end(),0.0);
    int tn=max(1,(int)(N*frac));
    vector<pair<double,int>> ca; ca.reserve(N);
    for(int p=0;p<N;p++){ ca.push_back({dist(order[p],order[(p+1)%N]),p}); }
    sort(ca.begin(),ca.end(),greater<pair<double,int>>());
    for(int i=0;i<tn&&i<(int)ca.size();i++) eHeat[ca[i].second]=0.5;
}
static void applyWarp(double ent){
    curSc=baseSc*ent;
    for(int i=0;i<N;i++){XW[i]=(X[i]-acx)*curSc;YW[i]=(Y[i]-acy)*curSc;}
}
static void resetWarp(){ curSc=baseSc; for(int i=0;i<N;i++){XW[i]=XS[i];YW[i]=YS[i];} }
int main(){
    T0=chrono::steady_clock::now();
    if(const char*e=getenv("SANTA_TL")){ double v=atof(e); if(v>50&&v<10000) TL_MS=v; }
    size_t cap=1<<20; inbuf.resize(cap); size_t len=0;
    while(true){ if(len==cap){cap<<=1;inbuf.resize(cap);} size_t g=fread(inbuf.data()+len,1,cap-len,stdin); if(!g)break; len+=g; } inbuf.resize(len);
    N=(int)readLL();
    if(N<=0){printf("1\n0\n");return 0;}
    X.resize(N);Y.resize(N);
    for(int i=0;i<N;i++){X[i]=(double)readLL();Y[i]=(double)readLL();}
    if(N==1){printf("2\n0\n0\n");return 0;}
    if(N==2){printf("3\n0\n1\n0\n");return 0;}
    double mnx=X[0],mxx=X[0],mny=Y[0],mxy=Y[0];
    for(int i=1;i<N;i++){mnx=min(mnx,X[i]);mxx=max(mxx,X[i]);mny=min(mny,Y[i]);mxy=max(mxy,Y[i]);}
    acx=(mnx+mxx)*0.5; acy=(mny+mxy)*0.5;
    double sx=max(mxx-mnx,1.0),sy=max(mxy-mny,1.0);
    baseSc=1.0/max(sx,sy); curSc=baseSc;
    XS.resize(N);YS.resize(N);XW.resize(N);YW.resize(N);
    for(int i=0;i<N;i++){XS[i]=(X[i]-acx)*baseSc;YS[i]=(Y[i]-acy)*baseSc;XW[i]=XS[i];YW[i]=YS[i];}
    vector<char> pr((size_t)N,0);
    {vector<char> comp((size_t)N,0);for(long long i=2;i<N;i++)if(!comp[i]){pr[i]=1;for(long long q=i*i;q<N;q+=i)comp[q]=1;}}
    if(N>100000)TL_MS-=20.0;
    double RESERVE=N>150000?220.0:(N>50000?90.0:(N>5000?50.0:40.0));
    TL_MS-=RESERVE;
    double minx=XS[0],maxx=XS[0],miny=YS[0],maxy=YS[0];
    for(int i=1;i<N;i++){minx=min(minx,XS[i]);maxx=max(maxx,XS[i]);miny=min(miny,YS[i]);maxy=max(maxy,YS[i]);}
    double w=max(1e-9,maxx-minx),h=max(1e-9,maxy-miny);
    int G=max(1,(int)floor(sqrt((double)N/2.0)));
    double cw=w/G,ch=h/G;
    auto gx=[&](double x){int c=(int)((x-minx)/cw);return c<0?0:(c>=G?G-1:c);};
    auto gy=[&](double y){int c=(int)((y-miny)/ch);return c<0?0:(c>=G?G-1:c);};
    vector<int> cellOf(N),cnt(G*G+1,0);
    for(int i=0;i<N;i++){int c=gx(XS[i])*G+gy(YS[i]);cellOf[i]=c;cnt[c+1]++;}
    for(int i=0;i<G*G;i++)cnt[i+1]+=cnt[i];
    vector<int> bucket(N);{vector<int> tmp=cnt;for(int i=0;i<N;i++)bucket[tmp[cellOf[i]]++]=i;}
    int K;
    if(N>50000)K=min(N-1,6);
    else if(N>=16000&&N<=36000)K=min(N-1,40);
    else if(N>5000)K=min(N-1,24);
    else K=min(N-1,10);
    if(const char*e=getenv("K_FORCE"))K=min(N-1,atoi(e));
    vector<int> nbr((size_t)N*K,-1);
    {vector<pair<double,int>>cand;cand.reserve(128);
     for(int i=0;i<N;i++){int cx=gx(XS[i]),cy=gy(YS[i]);cand.clear();int ring=0,extra=1;
      while(true){int x0=max(0,cx-ring),x1=min(G-1,cx+ring),y0=max(0,cy-ring),y1=min(G-1,cy+ring);
       for(int xx=x0;xx<=x1;xx++)for(int yy=y0;yy<=y1;yy++){if(ring>0&&xx>x0&&xx<x1&&yy>y0&&yy<y1)continue;
        int c=xx*G+yy;for(int b=cnt[c];b<cnt[c+1];b++){int j=bucket[b];if(j!=i)cand.push_back({dist(i,j),j});}}
       if((int)cand.size()>=K){if(extra--<=0)break;}
       if(x0==0&&y0==0&&x1==G-1&&y1==G-1)break;ring++;}
      int kk=min((int)cand.size(),K);partial_sort(cand.begin(),cand.begin()+kk,cand.end());
      for(int t=0;t<kk;t++)nbr[(size_t)i*K+t]=cand[t].second;}}
    vector<int> order(N),pos(N);vector<char> used(N,0);
    {int cur=0;used[0]=1;order[0]=0;
     for(int step=1;step<N;step++){int best=-1;double bd=1e300;
      for(int t=0;t<K;t++){int j=nbr[(size_t)cur*K+t];if(j>=0&&!used[j]){best=j;break;}}
      if(best<0){int cx=gx(XS[cur]),cy=gy(YS[cur]);
       for(int ring=0;ring<2*G&&best<0;ring++){int x0=max(0,cx-ring),x1=min(G-1,cx+ring),y0=max(0,cy-ring),y1=min(G-1,cy+ring);
        for(int xx=x0;xx<=x1;xx++)for(int yy=y0;yy<=y1;yy++){if(ring>0&&xx>x0&&xx<x1&&yy>y0&&yy<y1)continue;
         int c=xx*G+yy;for(int b=cnt[c];b<cnt[c+1];b++){int j=bucket[b];if(!used[j]){double d=dist(cur,j);if(d<bd){bd=d;best=j;}}}}
        if(best>=0){int r2=ring+1,a0=max(0,cx-r2),a1=min(G-1,cx+r2),b0=max(0,cy-r2),b1=min(G-1,cy+r2);
         for(int xx=a0;xx<=a1;xx++)for(int yy=b0;yy<=b1;yy++){if(xx>a0&&xx<a1&&yy>b0&&yy<b1)continue;
          int c=xx*G+yy;for(int bb=cnt[c];bb<cnt[c+1];bb++){int j=bucket[bb];if(!used[j]){double d=dist(cur,j);if(d<bd){bd=d;best=j;}}}}break;}}}
      if(best<0){for(int j=0;j<N;j++)if(!used[j]){best=j;break;}}
      used[best]=1;order[step]=best;cur=best;}
     for(int i=0;i<N;i++)pos[order[i]]=i;}
    detectRigid(order,nbr,K);
    eHeat.assign(N,0.0);
    auto ni=[&](int i){return i+1<N?i+1:0;};
    auto pi=[&](int i){return i>0?i-1:N-1;};
    auto am=[&](int e1,int e2){int lo=e1,hi=e2;if(lo>hi)swap(lo,hi);int inner=hi-lo;
     if(inner<=N-inner){int i=lo+1,j=hi;while(i<j){int a=order[i],b=order[j];order[i]=b;order[j]=a;pos[a]=j;pos[b]=i;++i;--j;}}
     else{int li=hi+1,lj=lo+N;while(li<lj){int ai=li%N,aj=lj%N;int u=order[ai],v=order[aj];order[ai]=v;order[aj]=u;pos[v]=ai;pos[u]=aj;++li;--lj;}}};
    vector<char> dl(N,0),dl2(N,0);
    vector<int> q(N);for(int i=0;i<N;i++)q[i]=i;
    int clk=0;
    auto tp=[&]()->bool{bool any=0;
     for(int qi=0;qi<N;qi++){if(((++clk)&1023)==0&&el_ms()>TL_MS)return any;int c1=q[qi];if(dl[c1])continue;bool imp=0;
      for(int dir=0;dir<2&&!imp;dir++){int p1=pos[c1],p2=dir==0?ni(p1):pi(p1),c2=order[p2];
       if(rigBrk(c1,c2,order,pos,rigid))continue;double d12=dist(c1,c2),h12=eHeat[p1];
       for(int t=0;t<K;t++){int c3=nbr[(size_t)c1*K+t];if(c3<0)break;double d13=dist(c1,c3);if(d13>=d12)break;
        int p3=pos[c3],p4=dir==0?ni(p3):pi(p3),c4=order[p4];
        if(c4==c1||c3==c2)continue;if(rigBrk(c3,c4,order,pos,rigid))continue;
        double d34=dist(c3,c4),h34=eHeat[p3];
        double bf=d12*(1.0+h12)+d34*(1.0+h34),af=d13+dist(c2,c4);
        if(af+1e-7<bf){if(dir==0)am(p1,p3);else am(p4,p2);dl[c1]=dl[c2]=dl[c3]=dl[c4]=0;imp=1;any=1;break;}}}
      if(!imp)dl[c1]=1;}return any;};
    auto orp=[&]()->bool{bool any=0;
     for(int qi=0;qi<N;qi++){if(((++clk)&1023)==0&&el_ms()>TL_MS)return any;int s0=q[qi];if(dl[s0])continue;bool mv=0;
      for(int L=1;L<=3&&!mv;L++){int is=pos[s0],ie=is;for(int t=1;t<L;t++)ie=ni(ie);int se=order[ie],pp=pi(is),pn=ni(ie);
       if(pp==ie||pn==is)break;int cp=order[pp],cn=order[pn];
       if(rigBrk(cp,s0,order,pos,rigid)||rigBrk(se,cn,order,pos,rigid))continue;
       double rm=dist(cp,s0)*(1.0+eHeat[pp])+dist(se,cn)*(1.0+eHeat[ie])-dist(cp,cn);
       if(rm<=1e-7)continue;
       for(int sd=0;sd<2&&!mv;sd++){int ac=sd==0?s0:se;
        for(int t=0;t<K;t++){int c=nbr[(size_t)ac*K+t];if(c<0)break;int pc=pos[c];bool in=0;
         {int p=is;for(int u=0;u<L;u++){if(p==pc){in=1;break;}p=ni(p);}}if(in||c==cp)continue;
         int pcn=ni(pc),cn2=order[pcn];if(cn2==s0)continue;
         double ad=dist(c,s0)+dist(se,cn2)-dist(c,cn2),ar=dist(c,se)+dist(s0,cn2)-dist(c,cn2);
         bool rev=ar+1e-9<ad;double add=rev?ar:ad;
         if(add+1e-7<rm){int seg[3];{int p=is;for(int u=0;u<L;u++){seg[u]=order[p];p=ni(p);}}
          if(rev){for(int a=0,b=L-1;a<b;a++,b--)swap(seg[a],seg[b]);}
          static vector<int>tmp;tmp.clear();tmp.reserve(N);
          auto mb=[&](int city)->bool{for(int u=0;u<L;u++)if(seg[u]==city)return 1;return 0;};
          for(int i2=0;i2<N;i2++){int city=order[i2];if(mb(city))continue;tmp.push_back(city);if(city==c){for(int u=0;u<L;u++)tmp.push_back(seg[u]);}}
          order.swap(tmp);for(int i2=0;i2<N;i2++)pos[order[i2]]=i2;
          dl[cp]=dl[cn]=dl[s0]=dl[se]=dl[c]=dl[cn2]=0;mv=1;any=1;break;}}}}
      if(!mv)dl[s0]=1;}return any;};
    auto lk=[&]()->bool{bool any=0;
     for(int qi=0;qi<N;qi++){if(((++clk)&511)==0&&el_ms()>TL_MS)return any;int c1=q[qi];if(dl2[c1])continue;bool imp=0;
      for(int dir=0;dir<2&&!imp;dir++){int p1=pos[c1],p2=dir==0?ni(p1):pi(p1),c2=order[p2];
       if(rigBrk(c1,c2,order,pos,rigid))continue;double d12=dist(c1,c2),h12=eHeat[p1];
       for(int t=0;t<K&&!imp;t++){int c3=nbr[(size_t)c1*K+t];if(c3<0)break;double d13=dist(c1,c3);if(d13>=d12)break;
        int p3=pos[c3],p4=dir==0?ni(p3):pi(p3),c4=order[p4];
        if(c4==c1||c3==c2)continue;if(rigBrk(c3,c4,order,pos,rigid))continue;
        double d34=dist(c3,c4),h34=eHeat[p3];
        double g1=(d12*(1.0+h12)+d34*(1.0+h34))-(d13+dist(c2,c4));
        int A=dir==0?p1:p4,B=dir==0?p3:p2;am(A,B);
        double bg=1e-7;int bc=-1,bd=-1,bf=-1;
        for(int sd=0;sd<2;sd++){int cc=sd==0?c2:c4;int pcc=pos[cc];
         for(int dd=0;dd<2;dd++){int pe=dd==0?ni(pcc):pi(pcc),ce=order[pe];
          if(rigBrk(cc,ce,order,pos,rigid))continue;double dce=dist(cc,ce);
          double hcc=dd==0?eHeat[pcc]:eHeat[pi(pcc)];
          for(int u=0;u<K;u++){int cf=nbr[(size_t)cc*K+u];if(cf<0)break;double dcf=dist(cc,cf);if(dcf>=dce)break;
           int pf=pos[cf],pg=dd==0?ni(pf):pi(pf),cg=order[pg];
           if(cg==cc||cf==ce)continue;if(rigBrk(cf,cg,order,pos,rigid))continue;
           double g2=(dce*(1.0+hcc)+dist(cf,cg))-(dcf+dist(ce,cg));
           if(g2>bg){bg=g2;bc=cc;bd=dd;bf=cf;}}}}
        if(bc>=0&&g1+bg>1e-7){int pcc=pos[bc],pf=pos[bf];
         int A2=bd==0?pcc:pi(pf),B2=bd==0?pf:pi(pcc);am(A2,B2);
         dl2[c1]=dl2[c2]=dl2[c3]=dl2[c4]=0;dl[c1]=dl[c2]=dl[c3]=dl[c4]=0;dl2[bc]=dl[bc]=0;imp=1;any=1;}
        else{am(A,B);}}}
      if(!imp)dl2[c1]=1;}return any;};
    auto runLS=[&](){while(el_ms()<TL_MS){bool a=0;while(el_ms()<TL_MS){if(!tp())break;a=1;}
     bool b=0;if(N<=50000){fill(dl.begin(),dl.end(),0);b=orp();if(b)fill(dl.begin(),dl.end(),0);}
     bool c=0;if(N>=8000){fill(dl2.begin(),dl2.end(),0);c=lk();if(c)fill(dl.begin(),dl.end(),0);}
     if(!a&&!b&&!c)break;}};
    auto ed=[&](int dir)->double{int z=pos[0];double L=0;int pv=0;
     for(int t=1;t<=N;t++){int idx=dir==0?(z+t)%N:(z-(t%N)+N)%N;int b=order[idx];
      double d=distOrig(pv,b);if(t%10==0&&!pr[pv])d*=1.1;L+=d;pv=b;}return L;};
    auto eb=[&](int&dir)->double{double f=ed(0),r=ed(1);if(f<=r){dir=0;return f;}dir=1;return r;};
    vector<int> best=order;int bestDir=0;
    double bestLen=eb(bestDir);
    if(N>=8){
        int stag=0,heatUsed=0;
        while(el_ms()<TL_MS){
            double entFact=1.0+((int)(rng()%6001)-3000)*0.0001;
            applyWarp(entFact);
            int a=1+(int)(rng()%(N-3)),b=1+(int)(rng()%(N-3)),c=1+(int)(rng()%(N-3));
            int lo=min({a,b,c}),hi=max({a,b,c}),mid=a+b+c-lo-hi;
            if(lo==mid||mid==hi){resetWarp();continue;}
            static vector<int>nt;nt.clear();nt.reserve(N);
            for(int i=0;i<lo;i++)nt.push_back(order[i]);
            for(int i=mid;i<hi;i++)nt.push_back(order[i]);
            for(int i=lo;i<mid;i++)nt.push_back(order[i]);
            for(int i=hi;i<N;i++)nt.push_back(order[i]);
            order.swap(nt);for(int i=0;i<N;i++)pos[order[i]]=i;
            detectRigid(order,nbr,K);
            fill(dl.begin(),dl.end(),0);
            while(el_ms()<TL_MS){if(!tp())break;}
            if(N<=50000){fill(dl.begin(),dl.end(),0);orp();while(el_ms()<TL_MS){if(!tp())break;}}
            resetWarp();
            int d2=0;double L=eb(d2);
            if(L<bestLen-1e-6){bestLen=L;best=order;bestDir=d2;stag=0;heatUsed=0;fill(eHeat.begin(),eHeat.end(),0.0);}
            else{order=best;for(int i=0;i<N;i++)pos[order[i]]=i;stag++;}
            if(stag>=5&&!heatUsed){
                injHeat(order,0.05);heatUsed=1;
                fill(dl.begin(),dl.end(),0);
                while(el_ms()<TL_MS){if(!tp())break;}
                if(N<=50000){fill(dl.begin(),dl.end(),0);orp();while(el_ms()<TL_MS){if(!tp())break;}}
                d2=0;L=eb(d2);
                if(L<bestLen-1e-6){bestLen=L;best=order;bestDir=d2;stag=0;heatUsed=0;}
                fill(eHeat.begin(),eHeat.end(),0.0);
                order=best;for(int i=0;i<N;i++)pos[order[i]]=i;
            }
        }
        order=best;for(int i=0;i<N;i++)pos[order[i]]=i;
    }
    TL_MS+=RESERVE;
    vector<int> seq(N);
    {int z=pos[0];for(int i=0;i<N;i++)seq[i]=bestDir==0?order[(z+i)%N]:order[(z-i+N)%N];}
    if(N>=12){
        auto sAt=[&](int p)->int{return p<N?seq[p]:0;};
        auto sc=[&](int t)->double{int a=seq[t-1],b=sAt(t);double d=distOrig(a,b);if(t%10==0&&!pr[a])d*=1.1;return d;};
        int w=N<=1200?N:(N<=5000?120:(N<=20000?140:80));
        for(int rep=0;rep<4&&el_ms()<TL_MS;rep++){bool ch=0;
         for(int p=9;p<N;p+=10){if(el_ms()>TL_MS)break;if(pr[seq[p]])continue;
          int lo=max(1,p-w),hi=min(N-1,p+w);double bd=-1e-7;int bj=-1;
          for(int j=lo;j<=hi;j++){if(j==p||!pr[seq[j]])continue;
           int T[4]={p,p+1,j,j+1},U[4],m=0;
           for(int t2=0;t2<4;t2++){bool dup=0;for(int u=0;u<m;u++)if(U[u]==T[t2])dup=1;if(!dup)U[m++]=T[t2];}
           double bf=0,af=0;for(int u=0;u<m;u++)bf+=sc(U[u]);
           swap(seq[p],seq[j]);for(int u=0;u<m;u++)af+=sc(U[u]);swap(seq[p],seq[j]);
           double dl=af-bf;if(dl<bd){bd=dl;bj=j;}}
          if(bj>=0){swap(seq[p],seq[bj]);ch=1;}}if(!ch)break;}
    }
    string out;out.reserve((size_t)N*7+16);
    out+=to_string(N+1);out+='\n';
    for(int i=0;i<N;i++){out+=to_string(seq[i]);out+='\n';}
    out+="0\n";
    fwrite(out.data(),1,out.size(),stdout);fflush(stdout);
    return 0;
}
