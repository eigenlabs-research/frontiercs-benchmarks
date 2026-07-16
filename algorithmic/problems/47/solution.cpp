#include <bits/stdc++.h>
using namespace std;

struct Item { string id; int w,h; long long v; int lim; };
struct Rect { int x,y,w,h; };
struct Pl { int t,x,y,rot,w,h; };

// Small JSON reader: the instance schema contains only objects, arrays, strings, numbers and bools.
struct Reader {
    string s; size_t p=0;
    Reader(){ s.assign(istreambuf_iterator<char>(cin), {}); }
    void ws(){ while(p<s.size() && isspace((unsigned char)s[p])) ++p; }
    void ch(char c){ ws(); if(p<s.size() && s[p]==c) ++p; }
    string str(){ ws(); string r; if(p<s.size() && s[p]=='"') ++p; while(p<s.size() && s[p]!='"'){ char c=s[p++]; if(c=='\\' && p<s.size()){ char e=s[p++]; if(e=='n') r+='\n'; else if(e=='t') r+='\t'; else r+=e; } else r+=c; } if(p<s.size()) ++p; return r; }
    long long num(){ ws(); int z=1; if(s[p]=='-'){z=-1;++p;} long long r=0; while(p<s.size()&&isdigit((unsigned char)s[p])) r=r*10+s[p++]-'0'; return z*r; }
    bool boolean(){ ws(); bool r=p<s.size()&&s[p]=='t'; while(p<s.size()&&isalpha((unsigned char)s[p]))++p; return r; }
};
static bool contains(const Rect&a,const Rect&b){ return b.x>=a.x&&b.y>=a.y&&b.x+b.w<=a.x+a.w&&b.y+b.h<=a.y+a.h; }
static bool hit(const Rect&a,const Rect&b){ return a.x<b.x+b.w&&b.x<a.x+a.w&&a.y<b.y+b.h&&b.y<a.y+a.h; }

struct Solver {
    int W,H; bool rotate; vector<Item> a;
    vector<Pl> run(int mode) {
        vector<Rect> fr(1,{0,0,W,H}); vector<int> used(a.size()); vector<Pl> out;
        const int FREE_CAP=700;
        while(!fr.empty()) {
            int bi=-1,bf=-1,br=0,bw=0,bh=0; long double best=-1;
            for(int i=0;i<(int)a.size();++i) if(used[i]<a[i].lim && a[i].v>0) {
                for(int ro=0;ro<=(rotate && a[i].w!=a[i].h);++ro) {
                    int iw=ro?a[i].h:a[i].w, ih=ro?a[i].w:a[i].h;
                    for(int j=0;j<(int)fr.size();++j) if(iw<=fr[j].w&&ih<=fr[j].h) {
                        long double den=(long double)a[i].v/(iw*(long double)ih);
                        long double fill=(long double)iw*ih/(fr[j].w*(long double)fr[j].h);
                        long double side=min((long double)iw/fr[j].w,(long double)ih/fr[j].h);
                        long double q;
                        if(mode==0) q=den;
                        else if(mode==1) q=den*(1.0L+0.32L*fill);
                        else if(mode==2) q=den*(1.0L+0.50L*side);
                        else if(mode==3) q=den*(1.0L+0.20L*fill+0.25L*side);
                        else q=den + (long double)a[i].v*1e-7L*fill;
                        // Best short-side fit is a stable tie breaker and avoids gratuitous slivers.
                        int shortleft=min(fr[j].w-iw,fr[j].h-ih);
                        int oldleft=(bf<0?INT_MAX:min(fr[bf].w-bw,fr[bf].h-bh));
                        if(q>best+1e-15L || (fabsl(q-best)<1e-15L && shortleft<oldleft))
                            best=q,bi=i,bf=j,br=ro,bw=iw,bh=ih;
                    }
                }
            }
            if(bi<0) break;
            Rect take{fr[bf].x,fr[bf].y,bw,bh};
            out.push_back({bi,take.x,take.y,br,bw,bh}); ++used[bi];
            vector<Rect> nf; nf.reserve(fr.size()*2+4);
            for(const Rect& f:fr) {
                if(!hit(f,take)) { nf.push_back(f); continue; }
                if(take.x>f.x) nf.push_back({f.x,f.y,take.x-f.x,f.h});
                if(take.x+take.w<f.x+f.w) nf.push_back({take.x+take.w,f.y,f.x+f.w-(take.x+take.w),f.h});
                if(take.y>f.y) nf.push_back({f.x,f.y,f.w,take.y-f.y});
                if(take.y+take.h<f.y+f.h) nf.push_back({f.x,take.y+take.h,f.w,f.y+f.h-(take.y+take.h)});
            }
            vector<char> dead(nf.size());
            for(int i=0;i<(int)nf.size();++i) if(!dead[i]) for(int j=0;j<(int)nf.size();++j)
                if(i!=j&&!dead[j]&&contains(nf[i],nf[j])) dead[j]=1;
            fr.clear(); for(int i=0;i<(int)nf.size();++i) if(!dead[i]&&nf[i].w>0&&nf[i].h>0) fr.push_back(nf[i]);
            // This only discards legal empty regions; it is a runtime bound, never a feasibility shortcut.
            if((int)fr.size()>FREE_CAP) {
                nth_element(fr.begin(),fr.begin()+FREE_CAP,fr.end(),[](const Rect&x,const Rect&y){return (long long)x.w*x.h>(long long)y.w*y.h;});
                fr.resize(FREE_CAP);
            }
        }
        return out;
    }
};
static string esc(const string&s){ string r; for(char c:s){ if(c=='"'||c=='\\')r+='\\'; r+=c; } return r; }
int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    Reader r; int W=0,H=0; bool rot=false; vector<Item> it;
    r.ch('{');
    while(true){ r.ws(); if(r.p>=r.s.size()||r.s[r.p]=='}'){r.ch('}');break;} string k=r.str(); r.ch(':');
        if(k=="bin") { r.ch('{'); while(true){ r.ws(); if(r.s[r.p]=='}'){r.ch('}');break;} string q=r.str(); r.ch(':'); if(q=="W")W=(int)r.num(); else if(q=="H")H=(int)r.num(); else if(q=="allow_rotate")rot=r.boolean(); r.ws(); if(r.s[r.p]==',')r.ch(','); } }
        else if(k=="items") { r.ch('['); while(true){ r.ws(); if(r.s[r.p]==']'){r.ch(']');break;} r.ch('{'); Item z; while(true){ r.ws(); if(r.s[r.p]=='}'){r.ch('}');break;} string q=r.str(); r.ch(':'); if(q=="type")z.id=r.str(); else if(q=="w")z.w=(int)r.num(); else if(q=="h")z.h=(int)r.num(); else if(q=="v")z.v=r.num(); else if(q=="limit")z.lim=(int)r.num(); r.ws(); if(r.s[r.p]==',')r.ch(','); } it.push_back(z); r.ws(); if(r.s[r.p]==',')r.ch(','); } }
        r.ws(); if(r.p<r.s.size()&&r.s[r.p]==',')r.ch(',');
    }
    Solver s{W,H,rot,it}; vector<Pl> ans; long long val=-1;
    // Policies differ only in the max-rectangle candidate ranking, making this a bounded ablation of that mechanism.
    for(int k=0;k<5;++k){ vector<Pl> q=s.run(k); long long v=0; for(auto&p:q)v+=it[p.t].v; if(v>val){val=v;ans.swap(q);} }
    cout<<"{\"placements\":[";
    for(size_t i=0;i<ans.size();++i){ if(i)cout<<','; auto&p=ans[i]; cout<<"{\"type\":\""<<esc(it[p.t].id)<<"\",\"x\":"<<p.x<<",\"y\":"<<p.y<<",\"rot\":"<<p.rot<<'}'; }
    cout<<"]}\n";
}
