#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>
#include <cstdint>
#include <random>
#include <chrono>
using namespace std;
typedef long long ll;
static chrono::steady_clock::time_point _t0;
static double _el(){return chrono::duration<double>(chrono::steady_clock::now()-_t0).count();}
static double _tl=0.95;
struct IT{string n;int w,h;ll v;int lmt;double d;};
struct Pl{int ti,x,y,r;};
static int _W,_H;static bool _ar;static vector<IT> _I;
struct Strip{
    int x0,x1,H;
    struct S{int y;int h;};
    vector<S> s;
    void init(int x0_,int x1_){x0=x0_;x1=x1_;H=_H;s={{0,0}};}
    bool fb(int rw,int rh,int&ox,int&oy){
        if(rw<=0||rh<=0||x0+rw>x1)return false;
        int by=1<<30,bx=x0;
        for(int i=0;i<(int)s.size();++i){int mx=s[i].h;int y=mx;
        if(y+rh<=_H&&(y<by||(y==by&&x0<bx))){by=y;bx=x0;}}
        if(by==1<<30)return false;ox=bx;oy=by;return true;}
    void pl(int x,int y,int rw,int rh){
        int nH=y+rh;vector<S> o;int i=0;
        while(i<(int)s.size()&&s[i].y<nH){o.push_back(s[i]);++i;}
        o.push_back({nH,0});s=o;}
};
vector<Pl> split_plan(vector<IT>&items,int W,int H,bool ar){
    vector<Pl> all;int sw=W/3;
    for(int strip=0;strip<3;++strip){
        int x0=strip*sw,x1=(strip==2)?W:(strip+1)*sw;
        Strip st;st.init(x0,x1);vector<int>ct(items.size(),0);
        vector<int>ord(items.size());iota(ord.begin(),ord.end(),0);
        sort(ord.begin(),ord.end(),[&](int a,int b){return items[a].d>items[b].d;});
        for(int idx:ord){IT&it=items[idx];while(ct[idx]<it.lmt){
            bool placed=false;int ori[2][2]={{it.w,it.h},{it.h,it.w}};
            int nO=ar&&it.w!=it.h?2:1;
            for(int o=0;o<nO&&!placed;++o){int rw=ori[o][0],rh=ori[o][1];int px,py;
            if(st.fb(rw,rh,px,py)){st.pl(px,py,rw,rh);all.push_back({idx,px,py,o});ct[idx]++;placed=true;}}
            if(!placed)break;}}}
    return all;}
vector<Pl> geo_rep(vector<Pl>init,vector<IT>&items,int iters){
    vector<Pl>best=init;ll bp=0;for(auto&p:best)bp+=items[p.ti].v;
    mt19937 rng(42);
    for(int it=0;it<iters;++it){if(_el()>_tl)break;if(best.empty())break;
    int i=rng()%best.size(),j=rng()%best.size();if(i==j)continue;
    swap(best[i],best[j]);ll np=0;for(auto&p:best)np+=items[p.ti].v;
    if(np>bp)bp=np;else swap(best[i],best[j]);}
    return best;}
struct JP{
    const string&s;size_t i=0;
    explicit JP(const string&x):s(x){}
    void sw(){while(i<s.size()&&(s[i]==32||s[i]==9||s[i]==10))++i;}
    char pk(){sw();return i<s.size()?s[i]:0;}
    char gc(){sw();return i<s.size()?s[i++]:0;}
    string rs(){
        sw();
        if(i<s.size()&&s[i]=='\"')++i;
        string o;
        while(i<s.size()&&s[i]!='\"'){
            if(s[i]=='\\'){
                ++i;
                if(i<s.size()){char e=s[i];++i;
                switch(e){
                    case'\"':o+='\"';break;
                    case'\\':o+='\\';break;
                    default:o+=e;break;
                }}
            }else{
                o+=s[i];++i;
            }
        }
        if(i<s.size())++i;
        return o;
    }
    ll rl(){
        sw();
        bool neg=false;
        if(i<s.size()&&(s[i]=='+'||s[i]=='-')){neg=(s[i]=='-');++i;}
        ll v=0;
        while(i<s.size()&&s[i]>='0'&&s[i]<='9'){v=v*10+(s[i]-'0');++i;}
        return neg?-v:v;
    }
    bool rb(){
        sw();
        if(i<s.size()&&s[i]=='t'){i+=4;return true;}
        if(i<s.size()&&s[i]=='f'){i+=5;return false;}
        return false;
    }
};
int main(){
    _t0=chrono::steady_clock::now();
    string inp;char buf[8192];while(fgets(buf,sizeof(buf),stdin))inp+=buf;
    JP jp(inp);jp.gc();jp.gc();string k=jp.rs();jp.gc();jp.gc();
    _W=0;_H=0;_ar=false;
    while(jp.pk()!='}'){jp.gc();k=jp.rs();jp.gc();
    if(k=="W")_W=jp.rl();else if(k=="H")_H=jp.rl();
    else if(k=="allow_rotate")_ar=jp.rb();if(jp.pk()==',')jp.gc();}
    jp.gc();jp.gc();k=jp.rs();jp.gc();jp.gc();
    while(jp.pk()!=']'){jp.gc();IT it;
    while(jp.pk()!='}'){jp.gc();k=jp.rs();jp.gc();
    if(k=="type")it.n=jp.rs();else if(k=="w")it.w=jp.rl();
    else if(k=="h")it.h=jp.rl();else if(k=="v")it.v=jp.rl();
    else if(k=="limit")it.lmt=jp.rl();
    if(jp.pk()==',')jp.gc();}
    jp.gc();it.d=(double)it.v/(it.w*it.h);_I.push_back(it);
    if(jp.pk()==',')jp.gc();}
    jp.gc();
    vector<Pl>best=split_plan(_I,_W,_H,_ar);
    best=geo_rep(best,_I,10000);
    printf("{\"placements\":[");
    for(int i=0;i<(int)best.size();++i){
        if(i)printf(",");
        Pl&p=best[i];
        printf("{\"type\":\"%s\",\"x\":%d,\"y\":%d,\"rot\":%d}",_I[p.ti].n.c_str(),p.x,p.y,p.r);}
    printf("]}");
    return 0;
}