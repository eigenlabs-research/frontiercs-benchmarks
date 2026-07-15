#include <bits/stdc++.h>
using namespace std;

struct Item { string id; long long w,h,v,lim; };
struct Pl { int t; long long x,y; int r; };
struct Seg { long long x, w, y; };

// Small JSON reader: input only contains objects, arrays, strings, integers and booleans.
struct Json {
    string s; size_t p=0;
    Json(){ s.assign(istreambuf_iterator<char>(cin), {}); }
    void ws(){ while(p<s.size() && isspace((unsigned char)s[p])) ++p; }
    void ch(char c){ ws(); if(p<s.size() && s[p]==c) ++p; }
    string str(){ ws(); string r; if(p<s.size() && s[p]=='\"') ++p; while(p<s.size() && s[p]!='\"'){ if(s[p]=='\\' && p+1<s.size()) ++p; r+=s[p++]; } if(p<s.size()) ++p; return r; }
    long long num(){ ws(); long long z=0, sg=1; if(s[p]=='-') sg=-1,++p; while(p<s.size()&&isdigit((unsigned char)s[p])) z=z*10+s[p++]-'0'; return z*sg; }
    bool boolean(){ ws(); bool q=s.compare(p,4,"true")==0; p+=q?4:5; return q; }
};

static bool findSpot(const vector<Seg>& a, long long W, long long H, long long rw, long long rh, long long &bx, long long &by) {
    bool ok=false; bx=by=0;
    for(int i=0;i<(int)a.size();++i){
        long long x=a[i].x;
        if(x+rw>W) continue;
        long long end=x+rw, y=0;
        for(int j=i;j<(int)a.size() && a[j].x<end;++j) y=max(y,a[j].y);
        if(y+rh>H) continue;
        // Bottom-left placement, with a small preference for a tight right edge.
        if(!ok || y<by || (y==by && x<bx)){ ok=true; bx=x; by=y; }
    }
    return ok;
}
static void raiseSkyline(vector<Seg>& a,long long x,long long w,long long y){
    long long e=x+w; vector<Seg> b;
    for(auto q:a){
        long long qe=q.x+q.w;
        if(qe<=x || q.x>=e) b.push_back(q);
        else {
            if(q.x<x) b.push_back({q.x,x-q.x,q.y});
            if(qe>e) b.push_back({e,qe-e,q.y});
        }
    }
    b.push_back({x,w,y});
    sort(b.begin(),b.end(),[](const Seg&A,const Seg&B){return A.x<B.x;});
    vector<Seg> c;
    for(auto q:b){
        if(q.w==0) continue;
        if(!c.empty() && c.back().x+c.back().w==q.x && c.back().y==q.y) c.back().w+=q.w;
        else c.push_back(q);
    }
    a.swap(c);
}

static pair<long long,vector<Pl>> pack(const vector<Item>& it,long long W,long long H,bool rot, vector<int> ord, int mode){
    vector<Seg> sky={{0,W,0}}; vector<int> used(it.size()); vector<Pl> ans; long long val=0;
    // Take the highest ranked still-placeable type.  This lets lower ranked types consume
    // awkward residual skyline gaps once a dense type no longer fits.
    while(true){
        int choose=-1, cr=0; long long cx=0,cy=0;
        for(int z:ord) if(used[z]<it[z].lim){
            bool have=false; long long ax,ay,bx,by; int rr=0;
            if(findSpot(sky,W,H,it[z].w,it[z].h,ax,ay)){ have=true; bx=ax;by=ay; }
            if(rot && it[z].w!=it[z].h && findSpot(sky,W,H,it[z].h,it[z].w,ax,ay)) {
                // The two orientations are selected by bottom-left, then either low skyline
                // profile or (in alternate trials) the orientation with the smaller height.
                if(!have || ay<by || (ay==by && ((mode%4==2) || (mode%4==1 && it[z].w<it[z].h) || (mode%4==0 && ax<bx)))){
                    have=true; bx=ax;by=ay;rr=1;
                }
            }
            if(have){ choose=z; cx=bx;cy=by;cr=rr; break; }
        }
        if(choose<0) break;
        long long rw=cr?it[choose].h:it[choose].w, rh=cr?it[choose].w:it[choose].h;
        ans.push_back({choose,cx,cy,cr}); val+=it[choose].v; ++used[choose];
        raiseSkyline(sky,cx,rw,cy+rh);
    }
    return {val,ans};
}
int main(){
    Json j; vector<Item> a; long long W=0,H=0; bool rotate=false;
    j.ch('{');
    while(true){ j.ws(); if(j.p>=j.s.size()||j.s[j.p]=='}'){j.ch('}');break;} string key=j.str(); j.ch(':');
        if(key=="bin") { j.ch('{'); for(int k=0;k<3;k++){ string q=j.str();j.ch(':'); if(q=="W")W=j.num(); else if(q=="H")H=j.num(); else rotate=j.boolean(); j.ws();if(j.p<j.s.size()&&j.s[j.p]==',')j.p++; }j.ch('}'); }
        else { j.ch('['); while(true){j.ws();if(j.p<j.s.size()&&j.s[j.p]==']'){j.p++;break;} j.ch('{'); Item x; for(int k=0;k<5;k++){string q=j.str();j.ch(':');if(q=="type")x.id=j.str();else if(q=="w")x.w=j.num();else if(q=="h")x.h=j.num();else if(q=="v")x.v=j.num();else x.lim=j.num();j.ws();if(j.p<j.s.size()&&j.s[j.p]==',')j.p++;}j.ch('}');a.push_back(x);j.ws();if(j.p<j.s.size()&&j.s[j.p]==',')j.p++; } }
        j.ws(); if(j.p<j.s.size()&&j.s[j.p]==',')j.p++;
    }
    vector<Pl> best; long long bestv=-1; int n=a.size();
    // Several deterministic value/area tradeoffs retain one packing mechanism but avoid a
    // single arbitrary ordering deciding every residual gap.
    for(int trial=0;trial<10;trial++){
        double exponent = trial==0 ? 1.0 : (trial==1 ? 0.0 : (trial==2 ? .55 : (trial==3 ? 1.45 : .85 + .06*trial)));
        vector<int> o(n); iota(o.begin(),o.end(),0);
        sort(o.begin(),o.end(),[&](int x,int y){
            long double sx=(long double)a[x].v/powl((long double)a[x].w*a[x].h,exponent);
            long double sy=(long double)a[y].v/powl((long double)a[y].w*a[y].h,exponent);
            if(fabsl(sx-sy)>1e-12L*max((long double)1,max(fabsl(sx),fabsl(sy)))) return sx>sy;
            return a[x].v>a[y].v;
        });
        auto got=pack(a,W,H,rotate,o,trial); if(got.first>bestv) bestv=got.first,best=move(got.second);
    }
    cout << "{\"placements\":[";
    for(size_t k=0;k<best.size();k++){ if(k) cout<<','; auto &p=best[k];
        cout<<"{\"type\":\""<<a[p.t].id<<"\",\"x\":"<<p.x<<",\"y\":"<<p.y<<",\"rot\":"<<p.r<<'}'; }
    cout << "]}\n";
}
