#include <bits/stdc++.h>
using namespace std;

struct JsonParser{
    string s; size_t p=0;
    JsonParser(const string& in):s(in){}
    void ws(){ while(p<s.size() && isspace((unsigned char)s[p])) p++; }
    void expect(char c){ ws(); if(p>=s.size()||s[p]!=c) exit(0); p++; }
    string str(){ ws(); expect('"'); string r; while(p<s.size()){ char c=s[p++]; if(c=='"') break; if(c=='\\'&&p<s.size()){ char e=s[p++]; if(e=='"'||e=='\\'||e=='/') r+=e; else if(e=='n') r+='\n'; else if(e=='t') r+='\t'; else if(e=='r') r+='\r'; else r+=e; } else r+=c; } return r; }
    long long num(){ ws(); int sign=1; if(p<s.size()&&(s[p]=='-'||s[p]=='+')){ if(s[p]=='-') sign=-1; p++; } long long v=0; while(p<s.size()&&isdigit((unsigned char)s[p])) v=v*10+(s[p++]-'0'); return sign*v; }
    bool boolean(){ ws(); if(s.compare(p,4,"true")==0){p+=4; return true;} if(s.compare(p,5,"false")==0){p+=5; return false;} return false; }
};

struct Item{ string id; int w,h; long long v; int lim; };
struct Placement{ int idx,x,y,rot,w,h; };
struct Seg{ int x,h; };

int W,H; bool allowRot; vector<Item> items;

struct Skyline{
    vector<Seg> seg;
    Skyline(){ seg.push_back({0,0}); seg.push_back({W,0}); }
    int heightAt(int x) const { int i=int(upper_bound(seg.begin(), seg.end(), Seg{x, INT_MAX}, [](const Seg&a,const Seg&b){return a.x<b.x;})-seg.begin())-1; return seg[i].h; }
    bool eval(int x,int w,int h,int &y) const{
        if(x<0 || x+w>W) return false;
        y=0;
        int i=int(upper_bound(seg.begin(), seg.end(), Seg{x, INT_MAX}, [](const Seg&a,const Seg&b){return a.x<b.x;})-seg.begin())-1;
        int pos=x;
        while(pos<x+w){
            y=max(y, seg[i].h);
            int nx=seg[i+1].x;
            if(nx<=pos) return false;
            pos=min(nx,x+w); i++;
        }
        return y+h<=H;
    }
    void place(int x,int w,int newh){
        int r=x+w;
        vector<Seg> pts;
        pts.reserve(seg.size()+2);
        bool ax=false, ar=false;
        for(size_t i=0;i<seg.size();++i){
            if(!ax && seg[i].x>x){ pts.push_back({x,heightAt(x)}); ax=true; }
            if(!ar && seg[i].x>r){ pts.push_back({r,heightAt(r)}); ar=true; }
            if(seg[i].x==x) ax=true;
            if(seg[i].x==r) ar=true;
            pts.push_back(seg[i]);
        }
        if(!ax) pts.push_back({x,heightAt(x)});
        if(!ar) pts.push_back({r,heightAt(r)});
        sort(pts.begin(), pts.end(), [](const Seg&a,const Seg&b){return a.x<b.x;});
        vector<Seg> out;
        for(auto sg: pts){
            if(!out.empty() && out.back().x==sg.x) out.back().h=sg.h; else out.push_back(sg);
        }
        vector<Seg> upd;
        for(auto sg: out){
            if(sg.x>=x && sg.x<r) sg.h=newh;
            if(upd.empty() || upd.back().h!=sg.h) upd.push_back(sg);
            else if(upd.back().x==sg.x) upd.back()=sg;
        }
        if(upd.empty() || upd.back().x!=W) upd.push_back({W,0});
        seg.swap(upd);
    }
};

string esc(const string& a){ string r; for(char c:a){ if(c=='"'||c=='\\') {r+='\\'; r+=c;} else r+=c; } return r; }

vector<Placement> runGreedy(int mode, int pass){
    vector<int> rem(items.size()); for(size_t i=0;i<items.size();++i) rem[i]=items[i].lim;
    Skyline sky; vector<Placement> res;
    while(true){
        double best=-1e300; int bi=-1,bx=0,by=0,br=0,bw=0,bh=0;
        for(int i=0;i<(int)items.size();++i) if(rem[i]>0 && items[i].v>0){
            for(int rot=0; rot<=(allowRot?1:0); ++rot){
                int w=rot?items[i].h:items[i].w, h=rot?items[i].w:items[i].h;
                if(w>W||h>H) continue;
                vector<int> xs; xs.reserve(sky.seg.size()*2);
                for(size_t si=0; si+1<sky.seg.size(); ++si){
                    xs.push_back(sky.seg[si].x);
                    xs.push_back(sky.seg[si].x - w); // right edge flush with a skyline break
                }
                sort(xs.begin(), xs.end()); xs.erase(unique(xs.begin(), xs.end()), xs.end());
                for(int x: xs){
                    int y;
                    if(x<0 || x+w>W || !sky.eval(x,w,h,y)) continue;
                    long long area=1LL*w*h;
                    double dens=(double)items[i].v/(double)area;
                    double score;
                    if(mode==0) score=dens*1e9 + items[i].v*0.02 - y*20 - x*0.01;
                    else if(mode==1) score=items[i].v + dens*200000.0 - y*50;
                    else if(mode==2) score=dens*1e9 - (double)(W-(x+w))*200.0 - y*5 + items[i].v*0.005;
                    else if(mode==3) score=(double)items[i].v/(w+pass*h) * 1000000.0 - y*30 - abs(w-h)*0.1;
                    else score=dens*8e8 + items[i].v*0.05 - (double)(y+h)*100 - sky.seg.size()*10;
                    // Prefer lower/left placements on nearly equal scores.
                    score -= (y*0.001 + x*0.000001);
                    if(score>best){ best=score; bi=i; bx=x; by=y; br=rot; bw=w; bh=h; }
                }
            }
        }
        if(bi<0) break;
        sky.place(bx,bw,by+bh);
        rem[bi]--;
        res.push_back({bi,bx,by,br,bw,bh});
        if(res.size()>20000) break;
    }
    return res;
}

long long valueOf(const vector<Placement>& v){ long long s=0; for(auto &p:v) s+=items[p.idx].v; return s; }

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    string input((istreambuf_iterator<char>(cin)), istreambuf_iterator<char>());
    JsonParser jp(input);
    jp.expect('{');
    for(int top=0; top<2; ++top){
        if(top) jp.expect(',');
        string k=jp.str(); jp.expect(':');
        if(k=="bin"){
            jp.expect('{');
            for(int f=0; f<3; ++f){ if(f) jp.expect(','); string bk=jp.str(); jp.expect(':'); if(bk=="W") W=jp.num(); else if(bk=="H") H=jp.num(); else if(bk=="allow_rotate") allowRot=jp.boolean(); }
            jp.expect('}');
        }else if(k=="items"){
            jp.expect('['); jp.ws(); bool first=true;
            while(jp.p<jp.s.size() && jp.s[jp.p]!=']'){
                if(!first) jp.expect(','); first=false;
                Item it; jp.expect('{');
                for(int f=0; f<5; ++f){ if(f) jp.expect(','); string ik=jp.str(); jp.expect(':'); if(ik=="type") it.id=jp.str(); else if(ik=="w") it.w=jp.num(); else if(ik=="h") it.h=jp.num(); else if(ik=="v") it.v=jp.num(); else if(ik=="limit") it.lim=jp.num(); }
                jp.expect('}'); items.push_back(it); jp.ws();
            }
            jp.expect(']');
        }
    }
    vector<Placement> best; long long bestv=-1;
    for(int mode=0; mode<5; ++mode){
        for(int pass=1; pass<=3; ++pass){
            auto cand=runGreedy(mode,pass);
            long long val=valueOf(cand);
            if(val>bestv){ bestv=val; best.swap(cand); }
        }
    }
    cout << "{\"placements\":[";
    for(size_t i=0;i<best.size();++i){
        if(i) cout << ',';
        auto &p=best[i];
        cout << "{\"type\":\"" << esc(items[p.idx].id) << "\",\"x\":" << p.x << ",\"y\":" << p.y << ",\"rot\":" << p.rot << "}";
    }
    cout << "]}\n";
    return 0;
}
