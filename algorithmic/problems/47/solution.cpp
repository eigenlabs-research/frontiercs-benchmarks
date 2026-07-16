#include <bits/stdc++.h>
using namespace std;

// Small JSON reader sufficient for the fixed input schema.
struct Tok {
    string s; size_t p=0;
    Tok(string z):s(move(z)) {}
    void ws(){ while(p<s.size() && isspace((unsigned char)s[p])) ++p; }
    char ch(){ ws(); return p<s.size()?s[p]:0; }
    void eat(char c){ ws(); if(p<s.size() && s[p]==c) ++p; }
    string str(){
        ws(); string r; if(p==s.size() || s[p++]!='"') return r;
        while(p<s.size() && s[p]!='"') { if(s[p]=='\\' && p+1<s.size()) ++p; r+=s[p++]; }
        if(p<s.size()) ++p; return r;
    }
    long long num(){ ws(); long long z=0, sign=1; if(p<s.size()&&s[p]=='-') sign=-1,++p; while(p<s.size()&&isdigit((unsigned char)s[p])) z=z*10+s[p++]-'0'; return sign*z; }
    bool boolean(){ ws(); bool r=s.compare(p,4,"true")==0; while(p<s.size()&&isalpha((unsigned char)s[p])) ++p; return r; }
};
struct Item { string id; int w,h; long long v,lim; };
struct Pl { int id,x,y,rot; };
struct Rect { int x,y,w,h; };
struct Cand { long long val=0; vector<Pl> p; };
int W,H; bool canrot; vector<Item> a;

static bool contains(const Rect& u,const Rect& v) {
    return u.x<=v.x && u.y<=v.y && u.x+u.w>=v.x+v.w && u.y+u.h>=v.y+v.h;
}

// Maximal rectangles: unlike a skyline, every remaining rectangular cavity stays available.
Cand pack(int mode) {
    vector<Rect> freeR(1,{0,0,W,H});
    vector<long long> used(a.size()); vector<Pl> out;
    long long value=0;
    const int CAP=5000;
    while(!freeR.empty() && (int)out.size()<CAP) {
        int bi=-1, bf=-1, br=0, bw=0,bh=0;
        long double best=-1e100L;
        for(int f=0; f<(int)freeR.size(); ++f) {
            const Rect &q=freeR[f];
            for(int i=0;i<(int)a.size();++i) if(used[i]<a[i].lim) {
                for(int r=0;r<=(canrot?1:0);++r) {
                    if(r && a[i].w==a[i].h) continue;
                    int w=r?a[i].h:a[i].w, h=r?a[i].w:a[i].h;
                    if(w>q.w || h>q.h) continue;
                    long double den=(long double)a[i].v/(w*(long long)h);
                    long double fill=(long double)(w*(long long)h)/(q.w*(long long)q.h);
                    int shortLeft=min(q.w-w,q.h-h), longLeft=max(q.w-w,q.h-h);
                    long double sc;
                    // The modes differ only in free-space fit preference; density remains dominant.
                    if(mode==0) sc=den + 1e-12L*(1000000-shortLeft);
                    else if(mode==1) sc=den*(1.0L+0.45L*fill);
                    else if(mode==2) sc=den*(1.0L+0.20L*(long double)min(w,h)/max(w,h));
                    else if(mode==3) sc=den*(1.0L+0.30L*(long double)(q.w-w+q.h-h==0));
                    else sc=den + 1e-13L*(1000000-longLeft);
                    // Stable ties give the variants genuinely different geometric choices.
                    sc += (mode&1 ? w : h)*1e-18L;
                    if(sc>best) best=sc,bi=i,bf=f,br=r,bw=w,bh=h;
                }
            }
        }
        if(bi<0) break;
        Rect put{freeR[bf].x,freeR[bf].y,bw,bh};
        out.push_back({bi,put.x,put.y,br}); value+=a[bi].v; ++used[bi];
        vector<Rect> next; next.reserve(freeR.size()+4);
        for(const Rect &q: freeR) {
            int ix=max(q.x,put.x), iy=max(q.y,put.y);
            int ax=min(q.x+q.w,put.x+put.w), ay=min(q.y+q.h,put.y+put.h);
            if(ix>=ax || iy>=ay) { next.push_back(q); continue; }
            // Split every free rectangle touched by the placed item. Each result is known empty.
            if(put.x>q.x) next.push_back({q.x,q.y,put.x-q.x,q.h});
            if(put.x+put.w<q.x+q.w) next.push_back({put.x+put.w,q.y,q.x+q.w-(put.x+put.w),q.h});
            if(put.y>q.y) next.push_back({q.x,q.y,q.w,put.y-q.y});
            if(put.y+put.h<q.y+q.h) next.push_back({q.x,put.y+put.h,q.w,q.y+q.h-(put.y+put.h)});
        }
        // MaxRects free rectangles can overlap, but contained ones add no possible placement.
        vector<char> dead(next.size());
        for(int i=0;i<(int)next.size();++i) for(int j=0;j<(int)next.size();++j)
            if(i!=j && contains(next[j],next[i]) && (next[j].x!=next[i].x || next[j].y!=next[i].y || next[j].w!=next[i].w || next[j].h!=next[i].h || j<i)) { dead[i]=1; break; }
        freeR.clear();
        for(int i=0;i<(int)next.size();++i) if(!dead[i] && next[i].w>0 && next[i].h>0) freeR.push_back(next[i]);
        // Guard the quadratic pruning cost on pathological tiny-rectangle inputs.
        if(freeR.size()>9000) break;
    }
    return {value,move(out)};
}
static string esc(const string& s) {
    string r; for(char c:s) { if(c=='"'||c=='\\') r+='\\'; r+=c; } return r;
}
int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    string in((istreambuf_iterator<char>(cin)),{}); Tok t(in); t.eat('{');
    while(t.ch()!='}' && t.ch()) {
        string k=t.str(); t.eat(':');
        if(k=="bin") { t.eat('{'); while(t.ch()!='}') { string z=t.str(); t.eat(':'); if(z=="W")W=t.num(); else if(z=="H")H=t.num(); else if(z=="allow_rotate")canrot=t.boolean(); if(t.ch()==',')t.eat(','); } t.eat('}'); }
        else if(k=="items") { t.eat('['); while(t.ch()!=']') { t.eat('{'); Item it; while(t.ch()!='}') { string z=t.str(); t.eat(':'); if(z=="type")it.id=t.str(); else if(z=="w")it.w=t.num(); else if(z=="h")it.h=t.num(); else if(z=="v")it.v=t.num(); else if(z=="limit")it.lim=t.num(); if(t.ch()==',')t.eat(','); } t.eat('}'); a.push_back(it); if(t.ch()==',')t.eat(','); } t.eat(']'); }
        if(t.ch()==',')t.eat(',');
    }
    Cand ans;
    for(int mode=0;mode<5;++mode) { Cand c=pack(mode); if(c.val>ans.val) ans=move(c); }
    cout << "{\"placements\":[";
    for(size_t i=0;i<ans.p.size();++i) { const Pl &p=ans.p[i]; if(i) cout<<','; cout<<"{\"type\":\""<<esc(a[p.id].id)<<"\",\"x\":"<<p.x<<",\"y\":"<<p.y<<",\"rot\":"<<p.rot<<'}'; }
    cout << "]}\n";
}
