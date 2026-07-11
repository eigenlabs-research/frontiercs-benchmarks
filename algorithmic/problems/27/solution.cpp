// Rectangle-free grid v4 (planes + diff-family + polish local search)
#include <bits/stdc++.h>
using namespace std;
typedef vector<vector<int>> Cols;
static long long total(const Cols&c){long long s=0;for(auto&v:c)s+=v.size();return s;}
struct GF{
    int q,p,e; vector<int> mul; bool ok=false;
    int add(int a,int b)const{ if(e==1) return (a+b)%p; int r=0,base=1;
        for(int i=0;i<e;i++){int da=(a/base)%p,db=(b/base)%p; r+=((da+db)%p)*base; base*=p;} return r; }
    int mult(int a,int b)const{ return mul[a*q+b]; }
};
GF buildField(int q){
    GF f; f.q=q;
    int p=0; for(int d=2;d*d<=q;d++){if(q%d==0){p=d;break;}} if(p==0)p=q;
    int e=0,t=q; while(t>1){ if(t%p){p=-1;break;} t/=p; e++; }
    if(p==-1){ f.ok=false; return f; }
    f.p=p; f.e=e; f.ok=true; f.mul.assign(q*q,0);
    if(e==1){ for(int a=0;a<q;a++)for(int b=0;b<q;b++) f.mul[a*q+b]=(a*b)%q; return f; }
    vector<int> irr;
    if(p==2&&e==2) irr={1,1,1};
    else if(p==2&&e==3) irr={1,1,0,1};
    else if(p==2&&e==4) irr={1,1,0,0,1};
    else if(p==2&&e==5) irr={1,0,1,0,0,1};
    else if(p==3&&e==2) irr={1,0,1};
    else if(p==3&&e==3) irr={1,2,0,1};
    else if(p==5&&e==2) irr={2,0,1};
    else { f.ok=false; return f; }
    auto dec=[&](int v){ vector<int> c(e); for(int i=0;i<e;i++){c[i]=v%p; v/=p;} return c; };
    auto enc=[&](vector<int>&c){ int v=0,b=1; for(int i=0;i<e;i++){v+=(c[i]%p)*b;b*=p;} return v; };
    for(int a=0;a<q;a++)for(int b=0;b<q;b++){
        vector<int> ca=dec(a),cb=dec(b),prod(2*e,0);
        for(int i=0;i<e;i++)for(int j=0;j<e;j++) prod[i+j]=(prod[i+j]+ca[i]*cb[j])%p;
        for(int d=2*e-2; d>=e; d--){ int co=prod[d]%p; if(!co)continue;
            for(int k=0;k<=e;k++) prod[d-e+k]=(((prod[d-e+k]-co*irr[k])%p)+p)%p; }
        vector<int> res(e); for(int i=0;i<e;i++) res[i]=prod[i]%p;
        f.mul[a*q+b]=enc(res);
    }
    return f;
}
bool buildPG(int q, vector<vector<int>>& lines, int& v){
    GF f=buildField(q); if(!f.ok) return false;
    v=q*q+q+1;
    vector<array<int,3>> pts; pts.reserve(v);
    for(int x=0;x<q;x++)for(int y=0;y<q;y++) pts.push_back({x,y,1});
    for(int x=0;x<q;x++) pts.push_back({x,1,0});
    pts.push_back({1,0,0});
    lines.assign(v,{});
    for(int L=0;L<v;L++){int la=pts[L][0],lb=pts[L][1],lc=pts[L][2];
        for(int P=0;P<v;P++){int d=f.add(f.add(f.mult(la,pts[P][0]),f.mult(lb,pts[P][1])),f.mult(lc,pts[P][2]));
            if(d==0) lines[L].push_back(P);}}
    return true;
}
Cols alternating(const vector<vector<int>>& lines,int v,int A,int B,unsigned seed,clock_t deadline){
    if(v<=A) return {};
    vector<vector<int>> ptLines(v);
    for(int L=0;L<v;L++)for(int P:lines[L])ptLines[P].push_back(L);
    Cols best; long long bestTot=-1; std::mt19937 rng(seed);
    vector<int> idx(v); for(int i=0;i<v;i++)idx[i]=i;
    for(int restart=0;;restart++){
        if(restart>0 && clock()>deadline) break;
        vector<char> inP(v,0);
        if(restart==0){ for(int i=0;i<v;i++)inP[i]=1; }
        else { shuffle(idx.begin(),idx.end(),rng); for(int i=0;i<min(v,A);i++)inP[idx[i]]=1; }
        for(int it=0;it<40;it++){
            vector<pair<int,int>> sc(v);
            for(int L=0;L<v;L++){int c=0;for(int P:lines[L])if(inP[P])c++;sc[L]={c,L};}
            sort(sc.begin(),sc.end(),greater<pair<int,int>>());
            vector<char> inL(v,0); int lim2=min(v,B);
            for(int i=0;i<lim2;i++)inL[sc[i].second]=1;
            vector<pair<int,int>> sp(v);
            for(int P=0;P<v;P++){int c=0;for(int L:ptLines[P])if(inL[L])c++;sp[P]={c,P};}
            sort(sp.begin(),sp.end(),greater<pair<int,int>>());
            fill(inP.begin(),inP.end(),0); for(int i=0;i<min(v,A);i++)inP[sp[i].second]=1;
            long long c=0; for(int L=0;L<v;L++)if(inL[L])for(int P:lines[L])if(inP[P])c++;
            if(c>bestTot){ bestTot=c; vector<int> rid(v,-1); int nr=0; for(int P=0;P<v;P++)if(inP[P])rid[P]=nr++;
                Cols cur; for(int L=0;L<v;L++)if(inL[L]){vector<int> col;for(int P:lines[L])if(inP[P])col.push_back(rid[P]);
                    if(!col.empty())cur.push_back(col);} best=move(cur); }
        }
    }
    return best;
}
Cols restrict_(const vector<vector<int>>& lines,int v,int A,int B,int lab){
    vector<int> perm(v);
    if(lab==0){ for(int i=0;i<v;i++)perm[i]=i; } else { for(int i=0;i<v;i++)perm[i]=v-1-i; }
    vector<vector<int>> kept;
    for(auto&ln:lines){ vector<int> b; for(int g:ln){int r=perm[g]; if(r<A)b.push_back(r);} if(!b.empty())kept.push_back(move(b)); }
    sort(kept.begin(),kept.end(),[](const vector<int>&a,const vector<int>&b){return a.size()>b.size();});
    Cols res; int take=min((int)kept.size(),B); for(int i=0;i<take;i++)res.push_back(kept[i]); return res;
}
Cols build_affine(int p,int A,int B,int lab){
    int P2=p*p; vector<int> perm(P2);
    for(int id=0;id<P2;id++){int x=id/p,y=id%p; perm[id]=(lab==0?id:y*p+x);}
    vector<vector<int>> lines; vector<int> buf;
    for(int a=0;a<p;a++)for(int b=0;b<p;b++){buf.clear();for(int x=0;x<p;x++){int y=(a*x+b)%p;int r=perm[x*p+y];if(r<A)buf.push_back(r);}if(!buf.empty())lines.push_back(buf);}
    for(int c=0;c<p;c++){buf.clear();for(int y=0;y<p;y++){int r=perm[c*p+y];if(r<A)buf.push_back(r);}if(!buf.empty())lines.push_back(buf);}
    sort(lines.begin(),lines.end(),[](const vector<int>&a,const vector<int>&b){return a.size()>b.size();});
    Cols res;int take=min((int)lines.size(),B);for(int i=0;i<take;i++)res.push_back(lines[i]);return res;
}
static char usedPair[330*330];
static int deg_[330];
static void markCol(const vector<int>&col,int A){ for(int i=0;i<(int)col.size();i++)for(int j=i+1;j<(int)col.size();j++){int r=col[i],q=col[j];usedPair[(size_t)r*A+q]=1;usedPair[(size_t)q*A+r]=1;} }
Cols greedy_rr(Cols base,int A,int B,unsigned seed){
    memset(usedPair,0,(size_t)A*A); for(int i=0;i<A;i++)deg_[i]=0;
    for(auto&c:base){ markCol(c,A); for(int r:c)deg_[r]++; }
    Cols cols=move(base); while((int)cols.size()<B) cols.push_back({});
    std::mt19937 rng(seed?seed:99991u);
    vector<int> ord(B); for(int i=0;i<B;i++)ord[i]=i;
    for(;;){
        bool added=false; if(seed) shuffle(ord.begin(),ord.end(),rng);
        for(int ci:ord){ auto&col=cols[ci]; if((int)col.size()>=A) continue;
            int bestr=-1,bd=INT_MAX;
            for(int r=0;r<A;r++){ bool inc=false; for(int q:col)if(q==r){inc=true;break;} if(inc)continue;
                bool ok=true; for(int q:col)if(usedPair[(size_t)r*A+q]){ok=false;break;} if(!ok)continue;
                int d=deg_[r]; if(d<bd||(d==bd&&(rng()&3u)==0u)){bd=d;bestr=r;} }
            if(bestr<0) continue;
            for(int q:col){usedPair[(size_t)bestr*A+q]=1;usedPair[(size_t)q*A+bestr]=1;}
            col.push_back(bestr); deg_[bestr]++; added=true; }
        if(!added) break;
    }
    return cols;
}
Cols extend(Cols base,int A,int B,int cap,unsigned seed){
    memset(usedPair,0,(size_t)A*A); for(auto&c:base) markCol(c,A);
    Cols cols=move(base); std::mt19937 rng(seed?seed:12345u);
    vector<int> order(A);for(int i=0;i<A;i++)order[i]=i;
    while((int)cols.size()<B){
        if(seed){for(int i=A-1;i>0;i--){int j=rng()%(i+1);swap(order[i],order[j]);}}
        else rotate(order.begin(),order.begin()+((int)cols.size()%A),order.end());
        vector<int> col;
        for(int oi=0;oi<A;oi++){int r=order[oi];if((int)col.size()>=cap)break;bool ok=true;for(int q:col)if(usedPair[(size_t)r*A+q]){ok=false;break;}if(ok){for(int q:col){usedPair[(size_t)r*A+q]=1;usedPair[(size_t)q*A+r]=1;}col.push_back(r);}}
        cols.push_back(col);
    }
    for(auto&col:cols)for(int r=0;r<A;r++){bool in=false;for(int q:col)if(q==r){in=true;break;}if(in)continue;bool ok=true;for(int q:col)if(usedPair[(size_t)r*A+q]){ok=false;break;}if(ok){for(int q:col){usedPair[(size_t)r*A+q]=1;usedPair[(size_t)q*A+r]=1;}col.push_back(r);}}
    return cols;
}
vector<int> primesUpTo(int N){vector<int>pr;vector<char>c(N+1,0);for(int i=2;i<=N;i++){if(!c[i]){pr.push_back(i);for(int j=2*i;j<=N;j+=i)c[j]=1;}}return pr;}
static char localD_[330];
static vector<int> trySidonBlock(int A,const vector<char>&usedD,int target,mt19937&rng){
    vector<int> order(A); for(int i=0;i<A;i++)order[i]=i; shuffle(order.begin(),order.end(),rng);
    memset(localD_,0,A);
    vector<int> blk; blk.push_back(order[0]); vector<int> tmp;
    for(int idx=1; idx<A && (int)blk.size()<target; idx++){
        int r=order[idx]; bool ok=true; tmp.clear();
        for(int m:blk){int d1=(r-m+A)%A,d2=(m-r+A)%A;
            if(usedD[d1]||usedD[d2]||localD_[d1]||localD_[d2]){ok=false;break;} tmp.push_back(d1);tmp.push_back(d2);}
        if(ok){ sort(tmp.begin(),tmp.end()); for(size_t j=1;j<tmp.size();j++)if(tmp[j]==tmp[j-1]){ok=false;break;} }
        if(ok){ for(int d:tmp)localD_[d]=1; blk.push_back(r); }
    }
    return blk;
}
Cols diff_family(int A,int B,int s0,int u,int K,mt19937&rng){
    int t=(B+A-1)/A;
    vector<char> usedD(A,0); vector<vector<int>> blocks;
    for(int i=0;i<t;i++){
        int cap=(i<u)? s0+1 : s0; if(cap<1)cap=1;
        vector<int> best;
        for(int k=0;k<K;k++){ auto b=trySidonBlock(A,usedD,cap,rng); if(b.size()>best.size()) best=move(b); if((int)best.size()>=cap)break; }
        for(size_t a=0;a<best.size();a++)for(size_t b2=a+1;b2<best.size();b2++){int d1=(best[a]-best[b2]+A)%A,d2=(best[b2]-best[a]+A)%A;usedD[d1]=1;usedD[d2]=1;}
        blocks.push_back(move(best));
    }
    sort(blocks.begin(),blocks.end(),[](const vector<int>&a,const vector<int>&b){return a.size()>b.size();});
    Cols cols; cols.reserve(B);
    for(auto&blk:blocks){ for(int sh=0;sh<A&&(int)cols.size()<B;sh++){vector<int> col;col.reserve(blk.size());for(int x:blk)col.push_back((x+sh)%A);cols.push_back(move(col));} if((int)cols.size()>=B)break; }
    return cols;
}
static void setCol(const vector<int>&col,int A,char val){ for(int i=0;i<(int)col.size();i++)for(int j=i+1;j<(int)col.size();j++){int r=col[i],q=col[j];usedPair[(size_t)r*A+q]=val;usedPair[(size_t)q*A+r]=val;} }
static inline void growCol(vector<int>&col,int A,vector<int>&rowOrder,mt19937&rng){
    shuffle(rowOrder.begin(),rowOrder.end(),rng);
    for(int r:rowOrder){ bool inc=false; for(int q:col)if(q==r){inc=true;break;} if(inc)continue;
        bool ok=true; for(int q:col)if(usedPair[(size_t)r*A+q]){ok=false;break;} if(!ok)continue;
        for(int q:col){usedPair[(size_t)r*A+q]=1;usedPair[(size_t)q*A+r]=1;} col.push_back(r); }
}
Cols polish(Cols cols,int A,int B,clock_t deadline,mt19937&rng){
    (void)B;
    if(A<=1) return cols;
    memset(usedPair,0,(size_t)A*A);
    for(auto&c:cols) setCol(c,A,1);
    int nc=(int)cols.size();
    vector<int> colOrder(nc); for(int i=0;i<nc;i++)colOrder[i]=i;
    vector<int> rowOrder(A); for(int i=0;i<A;i++)rowOrder[i]=i;
    long long bestTot=0; for(auto&c:cols)bestTot+=c.size();
    Cols bestCols=cols; int stale=0;
    while(clock()<deadline){
        bool improved=false;
        shuffle(colOrder.begin(),colOrder.end(),rng);
        for(int ci:colOrder){ size_t before=cols[ci].size(); growCol(cols[ci],A,rowOrder,rng); if(cols[ci].size()>before)improved=true; }
        if(clock()>=deadline) break;
        shuffle(colOrder.begin(),colOrder.end(),rng);
        for(int ci:colOrder){
            if((ci&7)==0 && clock()>=deadline) break;
            auto&col=cols[ci]; if(col.size()<2) continue;
            int oldsz=col.size();
            int di=rng()%col.size(); int x=col[di];
            vector<int> base; base.reserve(oldsz); for(int q:col)if(q!=x)base.push_back(q);
            setCol(col,A,0); setCol(base,A,1);
            vector<int> nb=base; growCol(nb,A,rowOrder,rng);
            if((int)nb.size()>oldsz){ cols[ci]=nb; improved=true; }
            else if((int)nb.size()==oldsz && nb!=col){ cols[ci]=nb; }
            else { setCol(nb,A,0); setCol(col,A,1); }
        }
        if(clock()>=deadline) break;
        shuffle(colOrder.begin(),colOrder.end(),rng);
        for(int ci:colOrder){
            if((ci&7)==0 && clock()>=deadline) break;
            vector<int> old=cols[ci]; setCol(old,A,0);
            vector<int> nb; growCol(nb,A,rowOrder,rng);
            if(nb.size()>old.size()){ cols[ci]=nb; improved=true; }
            else if(nb.size()==old.size() && nb!=old){ cols[ci]=nb; }
            else { setCol(nb,A,0); setCol(old,A,1); }
        }
        long long tt=0; for(auto&c:cols)tt+=c.size();
        if(tt>bestTot){ bestTot=tt; bestCols=cols; stale=0; } else stale++;
        if(!improved || stale>=2){
            stale=0;
            int nk=max(2,nc/8);
            shuffle(colOrder.begin(),colOrder.end(),rng);
            for(int t=0;t<nk;t++){ int ci=colOrder[t]; setCol(cols[ci],A,0); cols[ci].clear(); }
            for(int t=0;t<nk;t++){ int ci=colOrder[t]; growCol(cols[ci],A,rowOrder,rng); }
        }
    }
    return bestCols;
}
static char chkPair[330*330];
bool isValid(const Cols&cols,int A,int B){
    if((int)cols.size()>B) return false;
    memset(chkPair,0,(size_t)A*A);
    for(auto&col:cols){ for(int r:col) if(r<0||r>=A) return false;
        for(int i=0;i<(int)col.size();i++)for(int j=i+1;j<(int)col.size();j++){
            int r=col[i],q=col[j]; if(r==q) return false;
            size_t k=(size_t)r*A+q; if(chkPair[k]) return false; chkPair[k]=1; chkPair[(size_t)q*A+r]=1; } }
    return true;
}
Cols solve(int A,int B){
    clock_t T0=clock();
    bool skew=(B>2*A);
    double DL_ALT1,DL_ALT2,DL_ALT3,DL_DF,DL_END,DL_POL;
    if(skew){ DL_ALT1=0.05; DL_ALT2=0.08; DL_ALT3=0.10; DL_DF=0.28; DL_END=0.38; DL_POL=0.60; }
    else    { DL_ALT1=0.22; DL_ALT2=0.29; DL_ALT3=0.33; DL_DF=0.38; DL_END=0.42; DL_POL=0.60; }
    Cols best;
    if(A<=1){ best.assign(B,{0}); return best; }
    auto consider=[&](Cols c){ if(total(c)>total(best) && isValid(c,A,B)) best=move(c); };
    vector<Cols> bases;
    static const int qs[]={2,3,4,5,7,8,9,11,13,16,17,19,23,25,27,29,31,32};
    for(int q:qs){ long long v=(long long)q*q+q+1; if(v>(long long)2*A+40 && v>60) continue;
        vector<vector<int>> lines;int vv; if(!buildPG(q,lines,vv))continue;
        for(int lab=0;lab<2;lab++){ Cols c=restrict_(lines,vv,A,B,lab); if(!c.empty())bases.push_back(c); consider(c);} }
    { vector<int> cand;
        for(int q:qs){ long long v=(long long)q*q+q+1; if(v>=A && v<=(long long)3*A+40){ GF f=buildField(q); if(f.ok) cand.push_back(q);} }
        sort(cand.begin(),cand.end());
        double dls[3]={DL_ALT1,DL_ALT2,DL_ALT3};
        for(int i=0;i<(int)cand.size() && i<3;i++){
            vector<vector<int>> lines;int vv; if(!buildPG(cand[i],lines,vv))continue;
            clock_t dl=T0+(clock_t)(dls[i]*CLOCKS_PER_SEC);
            Cols c=alternating(lines,vv,A,B,12345u+cand[i],dl);
            if(!c.empty()){ bases.push_back(c); consider(c); } }
    }
    for(int p:primesUpTo(60)){ long long v=(long long)p*p; if(v>3*A+30 && v>60) continue;
        for(int lab=0;lab<2;lab++){Cols c=build_affine(p,A,B,lab); if(!c.empty())bases.push_back(c); consider(c);} }
    Cols dfBest;
    if(B>A && (B+A-1)/A <= 400){
        int t=(B+A-1)/A, D=A/2, s0=1; while((long long)t*(s0+1)*s0/2<=D)s0++;
        clock_t dl=T0+(clock_t)(DL_DF*CLOCKS_PER_SEC);
        mt19937 rng(777u); int slist[2]={s0,s0+1};
        for(unsigned iter=0; clock()<dl; iter++){
            int s=slist[iter%2]; int u=(int)((iter/2)%(t+1)); int K=(t<=6)?200:60;
            Cols c=diff_family(A,B,s,u,K,rng);
            if((int)c.size()<B){ while((int)c.size()<B)c.push_back({0}); }
            if(total(c)>total(dfBest) && isValid(c,A,B)) dfBest=move(c);
        }
        if(!dfBest.empty()){ bases.push_back(dfBest); consider(dfBest); }
    }
    int dstar=(int)floor((1.0+sqrt(1.0+4.0*(double)A*(A-1)/max(1,B)))/2.0); if(dstar<1)dstar=1;
    sort(bases.begin(),bases.end(),[](const Cols&a,const Cols&b){return total(a)>total(b);});
    int nb=min((int)bases.size(),4);
    for(int i=0;i<nb;i++) for(int cap:{dstar,dstar+1,max(2,dstar-1)}) consider(extend(bases[i],A,B,cap,0));
    for(int cap:{max(2,dstar-1),dstar,dstar+1,dstar+2,A}) consider(extend({},A,B,cap,0));
    Cols topBase = bases.empty()? Cols{} : bases[0];
    consider(greedy_rr({},A,B,1u));
    if(!topBase.empty()) consider(greedy_rr(topBase,A,B,7u));
    clock_t endDL=T0+(clock_t)(DL_END*CLOCKS_PER_SEC); unsigned s=2;
    while(clock()<endDL){
        consider(greedy_rr({},A,B,s*2654435761u));
        if(!skew && !topBase.empty()) consider(greedy_rr(topBase,A,B,s*40503u+1));
        s++; if(s>20000)break;
    }
    clock_t polDL=T0+(clock_t)(DL_POL*CLOCKS_PER_SEC);
    { mt19937 prng(2024u);
      Cols pol=polish(best,A,B,polDL,prng);
      if(total(pol)>total(best) && isValid(pol,A,B)) best=move(pol);
    }
    return best;
}
int main(){
    int n,m; if(scanf("%d %d",&n,&m)!=2) return 0;
    if(n<=0||m<=0){ printf("0\n"); return 0; }
    int A=min(n,m),B=max(n,m); bool transp=(n>m);
    Cols cols=solve(A,B);
    vector<pair<int,int>> cells;
    for(int c=0;c<(int)cols.size();c++) for(int r:cols[c]) cells.push_back({r,c});
    string out; out.reserve(cells.size()*8+16);
    char buf[32]; int len=sprintf(buf,"%zu\n",cells.size()); out.append(buf,len);
    for(auto&pr:cells){int r=pr.first,c=pr.second,ar,ac; if(!transp){ar=r+1;ac=c+1;}else{ar=c+1;ac=r+1;}
        len=sprintf(buf,"%d %d\n",ar,ac); out.append(buf,len);}
    fputs(out.c_str(),stdout);
    return 0;
}
