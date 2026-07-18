#include <bits/stdc++.h>
using namespace std;

struct JsonParser {
    string s; size_t p=0;
    JsonParser(const string& in):s(in){}
    void ws(){ while(p<s.size() && isspace((unsigned char)s[p])) p++; }
    void expect(char c){ ws(); if(p>=s.size()||s[p]!=c) exit(0); p++; }
    string str(){ ws(); expect('"'); string r; while(p<s.size()){ char c=s[p++]; if(c=='"') break; if(c=='\\' && p<s.size()){ char e=s[p++]; if(e=='n') r+='\n'; else if(e=='t') r+='\t'; else if(e=='r') r+='\r'; else r+=e; } else r+=c; } return r; }
    long long num(){ ws(); long long sign=1, v=0; if(s[p]=='-'){sign=-1;p++;} while(p<s.size() && isdigit((unsigned char)s[p])) v=v*10+(s[p++]-'0'); return sign*v; }
    bool boolean(){ ws(); if(s.compare(p,4,"true")==0){p+=4; return true;} p+=5; return false; }
};

struct Item { string id; int w,h; long long v; int limit; };
struct Bin { int W,H; bool rot; };
struct Rect { int x,y,w,h; };
struct Placement { int t,x,y,rot,w,h; long long v; };

static bool contains(const Rect&a,const Rect&b){return a.x<=b.x&&a.y<=b.y&&a.x+a.w>=b.x+b.w&&a.y+a.h>=b.y+b.h;}
static bool intersect(const Rect&a,const Rect&b){return a.x<b.x+b.w&&b.x<a.x+a.w&&a.y<b.y+b.h&&b.y<a.y+a.h;}

struct Packer {
    Bin bin; vector<Item> it; vector<Rect> freeR; vector<int> used; vector<Placement> out;
    Packer(const Bin& b,const vector<Item>& items):bin(b),it(items){ freeR.push_back({0,0,b.W,b.H}); used.assign(items.size(),0); }
    bool place(int ti,int rot, bool bssf){
        int w = rot ? it[ti].h : it[ti].w, h = rot ? it[ti].w : it[ti].h;
        int best=-1, bx=0, by=0; long long bs1=LLONG_MAX, bs2=LLONG_MAX;
        for(int i=0;i<(int)freeR.size();++i){ Rect r=freeR[i]; if(w<=r.w && h<=r.h){
            long long dw=r.w-w, dh=r.h-h;
            long long s1 = bssf ? min(dw,dh) : (long long)r.y*1000000LL + r.x;
            long long s2 = bssf ? max(dw,dh) : min(dw,dh);
            if(s1<bs1 || (s1==bs1 && s2<bs2)){ best=i; bx=r.x; by=r.y; bs1=s1; bs2=s2; }
        }}
        if(best<0) return false;
        Rect pr{bx,by,w,h}; out.push_back({ti,bx,by,rot,w,h,it[ti].v}); used[ti]++;
        vector<Rect> nf;
        for(Rect r: freeR){
            if(!intersect(r,pr)){ nf.push_back(r); continue; }
            if(pr.x>r.x) nf.push_back({r.x,r.y,pr.x-r.x,r.h});
            if(pr.x+pr.w<r.x+r.w) nf.push_back({pr.x+pr.w,r.y,r.x+r.w-(pr.x+pr.w),r.h});
            if(pr.y>r.y) nf.push_back({r.x,r.y,r.w,pr.y-r.y});
            if(pr.y+pr.h<r.y+r.h) nf.push_back({r.x,pr.y+pr.h,r.w,r.y+r.h-(pr.y+pr.h)});
        }
        freeR.clear();
        for(Rect r:nf) if(r.w>0&&r.h>0) freeR.push_back(r);
        vector<char> rem(freeR.size(),0);
        for(int i=0;i<(int)freeR.size();++i) for(int j=0;j<(int)freeR.size();++j) if(i!=j&&!rem[i]&&contains(freeR[j],freeR[i])) rem[i]=1;
        vector<Rect> pruned; for(int i=0;i<(int)freeR.size();++i) if(!rem[i]) pruned.push_back(freeR[i]); freeR.swap(pruned);
        return true;
    }
    long long val() const { long long s=0; for(auto&p:out) s+=p.v; return s; }
};

struct Ori { int t,rot,w,h; long long v; double den; };
vector<Ori> orientations(const Bin& bin,const vector<Item>& items){
    vector<Ori> o;
    for(int i=0;i<(int)items.size();++i){
        auto &a=items[i]; if(a.limit<=0||a.v<=0) continue;
        if(a.w<=bin.W&&a.h<=bin.H) o.push_back({i,0,a.w,a.h,a.v,(double)a.v/(a.w*a.h)});
        if(bin.rot && a.w!=a.h && a.h<=bin.W&&a.w<=bin.H) o.push_back({i,1,a.h,a.w,a.v,(double)a.v/(a.w*a.h)});
    }
    return o;
}

vector<Placement> run_order(const Bin& bin,const vector<Item>& items,int mode,bool bssf){
    auto os=orientations(bin,items);
    sort(os.begin(), os.end(), [&](const Ori&a,const Ori&b){
        if(mode==0) return tuple<double,long long,int>(-a.den,-a.v,a.h) < tuple<double,long long,int>(-b.den,-b.v,b.h);
        if(mode==1) return tuple<long long,double>(-a.v,-a.den) < tuple<long long,double>(-b.v,-b.den);
        if(mode==2) return tuple<int,double>(-max(a.w,a.h),-a.den) < tuple<int,double>(-max(b.w,b.h),-b.den);
        if(mode==3) return tuple<int,double>(a.h,-a.den) < tuple<int,double>(b.h,-b.den);
        if(mode==4) return tuple<int,double>(-a.h,-a.den) < tuple<int,double>(-b.h,-b.den);
        return tuple<double,int>(-(a.den*sqrt((double)a.w*a.h)), a.w) < tuple<double,int>(-(b.den*sqrt((double)b.w*b.h)), b.w);
    });
    Packer pk(bin,items);
    bool any=true;
    while(any){ any=false; for(auto &o: os) while(pk.used[o.t]<items[o.t].limit && pk.place(o.t,o.rot,bssf)) any=true; }
    return pk.out;
}

vector<Placement> run_bestfit(const Bin& bin,const vector<Item>& items,int mode){
    Packer pk(bin,items); auto os=orientations(bin,items);
    while(true){
        int bestT=-1,bestR=0,bestF=-1,bx=0,by=0,bw=0,bh=0; double best=-1e100; long long tie=LLONG_MAX;
        for(auto&o:os) if(pk.used[o.t]<items[o.t].limit){
            for(int fi=0;fi<(int)pk.freeR.size();++fi){ Rect r=pk.freeR[fi]; if(o.w<=r.w&&o.h<=r.h){
                int dw=r.w-o.w, dh=r.h-o.h; double sc;
                if(mode==0) sc=o.den*1000000.0 + o.v*1e-3 - (dw+dh)*0.01;
                else if(mode==1) sc=o.v - (double)(dw*dh)*o.den*0.05;
                else sc=o.den*1000000.0 + (double)o.v/(1+min(dw,dh));
                long long ti=(long long)r.y*1000000LL+r.x+max(dw,dh);
                if(sc>best || (fabs(sc-best)<1e-9 && ti<tie)){ best=sc; tie=ti; bestT=o.t; bestR=o.rot; bestF=fi; bx=r.x; by=r.y; bw=o.w; bh=o.h; }
            }}
        }
        if(bestT<0) break;
        pk.place(bestT,bestR,true);
    }
    return pk.out;
}

vector<Placement> run_shelves(const Bin& bin,const vector<Item>& items,int mode){
    auto os=orientations(bin,items); vector<int> hs; for(auto&o:os) hs.push_back(o.h); sort(hs.begin(),hs.end()); hs.erase(unique(hs.begin(),hs.end()),hs.end());
    vector<Placement> best;
    for(int baseH: hs){
        vector<int> used(items.size(),0); vector<Placement> out; int y=0;
        while(y<bin.H){
            int sh=baseH; if(mode==1){ // choose best remaining shelf height by potential density
                double bd=-1; for(int h:hs) if(y+h<=bin.H){ double d=0; for(auto&o:os) if(o.h<=h&&used[o.t]<items[o.t].limit) d=max(d,o.den); if(d>bd){bd=d; sh=h;} }
            }
            if(y+sh>bin.H) break;
            vector<long long> dp(bin.W+1,-1); vector<int> pre(bin.W+1,-1), prew(bin.W+1); dp[0]=0;
            for(int x=0;x<=bin.W;x++) if(dp[x]>=0){
                for(int k=0;k<(int)os.size();++k){ auto&o=os[k]; if(o.h<=sh && used[o.t]<items[o.t].limit && x+o.w<=bin.W){ long long nv=dp[x]+o.v; if(nv>dp[x+o.w]){dp[x+o.w]=nv; pre[x+o.w]=k; prew[x+o.w]=x;} } }
            }
            int bx=max_element(dp.begin(),dp.end())-dp.begin(); if(dp[bx]<=0) break;
            vector<int> seq; int cur=bx; while(cur>0&&pre[cur]>=0){seq.push_back(pre[cur]); cur=prew[cur];}
            reverse(seq.begin(),seq.end()); int x=0; bool placed=false;
            for(int k:seq){ auto&o=os[k]; if(used[o.t]>=items[o.t].limit) continue; out.push_back({o.t,x,y,o.rot,o.w,o.h,o.v}); used[o.t]++; x+=o.w; placed=true; }
            if(!placed) break; y+=sh;
        }
        long long vb=0,vo=0; for(auto&p:best) vb+=p.v; for(auto&p:out) vo+=p.v; if(vo>vb) best=out;
    }
    return best;
}

string jsonEscape(const string& s){
    string r;
    for(char c:s){
        if(c=='"'||c=='\\'){ r.push_back('\\'); r.push_back(c); }
        else if(c=='\n') r += "\\n";
        else if(c=='\r') r += "\\r";
        else if(c=='\t') r += "\\t";
        else r.push_back(c);
    }
    return r;
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    string input((istreambuf_iterator<char>(cin)), istreambuf_iterator<char>());
    JsonParser jp(input); Bin bin{}; vector<Item> items;
    jp.expect('{');
    for(int top=0;top<2;top++){
        if(top) jp.expect(','); string key=jp.str(); jp.expect(':');
        if(key=="bin"){
            jp.expect('{'); for(int i=0;i<3;i++){ if(i) jp.expect(','); string k=jp.str(); jp.expect(':'); if(k=="W") bin.W=jp.num(); else if(k=="H") bin.H=jp.num(); else bin.rot=jp.boolean(); } jp.expect('}');
        }else{
            jp.expect('['); bool first=true; while(true){ jp.ws(); if(jp.p<input.size()&&input[jp.p]==']'){jp.p++; break;} if(!first) jp.expect(','); first=false; Item it; jp.expect('{'); for(int f=0;f<5;f++){ if(f) jp.expect(','); string k=jp.str(); jp.expect(':'); if(k=="type") it.id=jp.str(); else if(k=="w") it.w=jp.num(); else if(k=="h") it.h=jp.num(); else if(k=="v") it.v=jp.num(); else it.limit=jp.num(); } jp.expect('}'); items.push_back(it); }
        }
    }
    vector<Placement> best; auto upd=[&](vector<Placement> v){ long long a=0,b=0; for(auto&p:best)a+=p.v; for(auto&p:v)b+=p.v; if(b>a) best=move(v); };
    for(int m=0;m<6;m++){ upd(run_order(bin,items,m,true)); upd(run_order(bin,items,m,false)); }
    for(int m=0;m<3;m++) upd(run_bestfit(bin,items,m));
    upd(run_shelves(bin,items,0)); upd(run_shelves(bin,items,1));
    cout << "{\"placements\":[";
    for(size_t i=0;i<best.size();++i){ auto&p=best[i]; if(i) cout<<','; cout << "{\"type\":\"" << jsonEscape(items[p.t].id) << "\",\"x\":" << p.x << ",\"y\":" << p.y << ",\"rot\":" << p.rot << "}"; }
    cout << "]}\n";
    return 0;
}
