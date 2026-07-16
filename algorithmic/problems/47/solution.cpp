#include <bits/stdc++.h>
using namespace std;

struct Item { string id; int w,h; long long v; int lim; };
struct Rect { int x,y,w,h; };
struct Place { int t,x,y,rot,w,h; };

// Tiny JSON reader: input is the restricted JSON format specified by the task.
struct Json {
    string s; size_t p=0;
    Json(){ s.assign(istreambuf_iterator<char>(cin), {}); }
    void ws(){ while(p<s.size() && isspace((unsigned char)s[p])) ++p; }
    char ch(){ ws(); return p<s.size()?s[p++]:0; }
    string str(){
        ws(); if(p>=s.size() || s[p++]!='"') return "";
        string r;
        while(p<s.size() && s[p]!='"') {
            if(s[p]=='\\' && p+1<s.size()) { ++p; r+=s[p++]; }
            else r+=s[p++];
        }
        if(p<s.size()) ++p;
        return r;
    }
    long long num(){ ws(); long long sign=1,x=0; if(s[p]=='-') sign=-1,++p; while(p<s.size() && isdigit((unsigned char)s[p])) x=x*10+s[p++]-'0'; return sign*x; }
    bool boolean(){ ws(); if(s.compare(p,4,"true")==0) { p+=4; return true; } p+=5; return false; }
};

static bool intersects(const Rect& a,const Rect& b) {
    return a.x < b.x+b.w && b.x < a.x+a.w && a.y < b.y+b.h && b.y < a.y+a.h;
}
static bool contains(const Rect& a,const Rect& b) {
    return a.x<=b.x && a.y<=b.y && a.x+a.w>=b.x+b.w && a.y+a.h>=b.y+b.h;
}

struct Result { long long value=0; vector<Place> p; };

// Maximal-free-rectangle packing.  Every selected rectangle lies in an explicitly empty
// rectangle; splitting every intersected free rectangle maintains that invariant.
Result pack(int W,int H,bool rotok,const vector<Item>& it,double exponent,double fitweight,int mode) {
    vector<Rect> freeR(1,{0,0,W,H});
    vector<int> used(it.size());
    Result out;
    while(!freeR.empty()) {
        int bi=-1,bf=-1,br=0,bw=0,bh=0;
        double best=-1e300, bestwaste=1e300;
        for(int f=0;f<(int)freeR.size();++f) {
            const Rect &q=freeR[f];
            for(int t=0;t<(int)it.size();++t) if(used[t]<it[t].lim) {
                for(int r=0;r<=(rotok && it[t].w!=it[t].h);++r) {
                    int w=r?it[t].h:it[t].w, h=r?it[t].w:it[t].h;
                    if(w>q.w || h>q.h) continue;
                    double area=(double)w*h;
                    double val=(double)it[t].v / pow(area,exponent);
                    double shortfit=(double)min(q.w-w,q.h-h)/(max(q.w,q.h)+1.0);
                    double cover=area/((double)q.w*q.h);
                    // mode changes the harmless tie geometry preference, giving distinct layouts.
                    double shape = mode==0 ? cover : (mode==1 ? 1.0-shortfit : (double)min(w,h)/max(w,h));
                    double key=val*(1.0+fitweight*shape);
                    double waste=(double)q.w*q.h-area;
                    if(key>best+1e-12 || (fabs(key-best)<1e-12 && waste<bestwaste)) {
                        best=key; bestwaste=waste; bi=t; bf=f; br=r; bw=w; bh=h;
                    }
                }
            }
        }
        if(bi<0) break;
        Rect placed{freeR[bf].x,freeR[bf].y,bw,bh};
        out.p.push_back({bi,placed.x,placed.y,br,bw,bh});
        out.value+=it[bi].v; ++used[bi];
        vector<Rect> next;
        next.reserve(freeR.size()+4);
        for(const Rect& q:freeR) {
            if(!intersects(q,placed)) { next.push_back(q); continue; }
            if(placed.x>q.x) next.push_back({q.x,q.y,placed.x-q.x,q.h});
            if(placed.x+placed.w<q.x+q.w) next.push_back({placed.x+placed.w,q.y,q.x+q.w-(placed.x+placed.w),q.h});
            if(placed.y>q.y) next.push_back({q.x,q.y,q.w,placed.y-q.y});
            if(placed.y+placed.h<q.y+q.h) next.push_back({q.x,placed.y+placed.h,q.w,q.y+q.h-(placed.y+placed.h)});
        }
        // Remove contained free rectangles. Overlapping free rectangles are intentional in MaxRects.
        vector<char> gone(next.size());
        for(int i=0;i<(int)next.size();++i) if(!gone[i])
            for(int j=0;j<(int)next.size();++j) if(i!=j && !gone[i] && contains(next[j],next[i])) gone[i]=1;
        freeR.clear();
        for(int i=0;i<(int)next.size();++i) if(!gone[i] && next[i].w>0 && next[i].h>0) freeR.push_back(next[i]);
    }
    return out;
}

int main(){
    Json j;
    if(j.s.empty()) return 0;
    int W=0,H=0; bool allow=false; vector<Item> a;
    j.ch(); // root {
    while(true) {
        string key=j.str(); j.ch();
        if(key=="bin") {
            j.ch();
            for(int k=0;k<3;k++) { string q=j.str(); j.ch(); if(q=="W") W=j.num(); else if(q=="H") H=j.num(); else allow=j.boolean(); j.ws(); if(j.p<j.s.size()&&j.s[j.p]==',') ++j.p; }
            j.ch();
        } else if(key=="items") {
            j.ch(); j.ws();
            while(j.p<j.s.size() && j.s[j.p]!=']') {
                j.ch(); Item z;
                for(int k=0;k<5;k++) { string q=j.str(); j.ch(); if(q=="type") z.id=j.str(); else if(q=="w") z.w=j.num(); else if(q=="h") z.h=j.num(); else if(q=="v") z.v=j.num(); else z.lim=j.num(); j.ws(); if(j.p<j.s.size()&&j.s[j.p]==',') ++j.p; }
                j.ch(); a.push_back(z); j.ws(); if(j.p<j.s.size()&&j.s[j.p]==',') ++j.p; j.ws();
            }
            j.ch();
        }
        j.ws(); if(j.p>=j.s.size() || j.s[j.p]=='}') break; if(j.s[j.p]==',') ++j.p;
    }
    Result ans;
    const double exps[]={0.55,0.75,0.95,1.12,1.30};
    const double fits[]={0.0,0.22};
    for(int k=0;k<5;k++) for(double f:fits) {
        Result r=pack(W,H,allow,a,exps[k],f,k&1);
        if(r.value>ans.value) ans=move(r);
    }
    cout << "{\"placements\":[";
    for(size_t i=0;i<ans.p.size();++i) {
        if(i) cout<<',';
        const Place& q=ans.p[i];
        cout << "{\"type\":\"" << a[q.t].id << "\",\"x\":" << q.x << ",\"y\":" << q.y << ",\"rot\":" << q.rot << '}';
    }
    cout << "]}\n";
}
