#include <bits/stdc++.h>
using namespace std;

struct Item { string id; int w,h; long long v; int lim; };
struct Seg { int x,y,w; };
struct Put { int t,x,y,r; };

static string esc(const string& s) { string r; for(char c:s) { if(c=='"'||c=='\\') r+='\\'; r+=c; } return r; }

static vector<Seg> raised(const vector<Seg>& a, int x, int wid, int ny) {
    int z=x+wid; vector<Seg> b; b.reserve(a.size()+2);
    for (auto q:a) {
        int l=q.x, r=q.x+q.w;
        if (r<=x || l>=z) b.push_back(q);
        else {
            if(l<x) b.push_back({l,q.y,x-l});
            if(r>z) b.push_back({z,q.y,r-z});
        }
    }
    b.push_back({x,ny,wid});
    sort(b.begin(),b.end(),[](const Seg&A,const Seg&B){return A.x<B.x;});
    vector<Seg> c;
    for(auto q:b) {
        if(q.w<=0) continue;
        if(!c.empty() && c.back().x+c.back().w==q.x && c.back().y==q.y) c.back().w+=q.w;
        else c.push_back(q);
    }
    return c;
}

struct Result { long long val=-1; vector<Put> p; };

static Result run(const vector<Item>& it, int W, int H, bool rotok, long double denExp, long double wasteK, long double topK, int tieMode) {
    vector<Seg> sky(1,{0,0,W}); vector<int> used(it.size()); Result out; out.val=0;
    // A skyline candidate is always supported by the current contour.  Unlike free-rectangle
    // splitting, this representation has no overlapping residual-region bookkeeping.
    for(int step=0; step<6000 && sky.size()<3500; ++step) {
        long double best=-1; int bi=-1,bx=0,by=0,br=0,bw=0,bh=0;
        for(int i=0;i<(int)it.size();++i) if(used[i]<it[i].lim && it[i].v>0) {
            for(int rr=0;rr<=(rotok && it[i].w!=it[i].h);++rr) {
                int ww=rr?it[i].h:it[i].w, hh=rr?it[i].w:it[i].h;
                if(ww>W || hh>H) continue;
                for(auto s:sky) {
                    int x=s.x; if(x+ww>W) continue;
                    int top=0; long long delta=0; int z=x+ww;
                    for(auto q:sky) if(q.x<z && q.x+q.w>x) top=max(top,q.y);
                    if(top+hh>H) continue;
                    for(auto q:sky) {
                        int l=max(x,q.x), r=min(z,q.x+q.w);
                        if(l<r) delta += 1LL*(top+hh-q.y)*(r-l);
                    }
                    long double area=(long double)ww*hh;
                    long double density=(long double)it[i].v/area;
                    long double waste=max((long double)0,(long double)delta-area);
                    long double key=pow(density,denExp)*pow((long double)it[i].v,1.0L-denExp);
                    key/=1.0L+wasteK*waste/(area+1);
                    key/=1.0L+topK*(long double)(top+hh)/H;
                    // deterministic but deliberately varied final ties make profile families diverge
                    key += 1e-12L*(tieMode==0 ? -top : tieMode==1 ? -ww : tieMode==2 ? hh : it[i].v%97);
                    if(key>best) { best=key; bi=i; bx=x; by=top; br=rr; bw=ww; bh=hh; }
                }
            }
        }
        if(bi<0) break;
        out.p.push_back({bi,bx,by,br}); out.val+=it[bi].v; ++used[bi];
        sky=raised(sky,bx,bw,by+bh);
    }
    return out;
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    string in((istreambuf_iterator<char>(cin)),{});
    smatch m; int W=0,H=0; bool allow=false;
    // JSON object member order is not significant; accept every ordering allowed by the contract.
    regex binR("\\\"bin\\\"\\s*:\\s*\\{([^{}]*)\\}");
    if(!regex_search(in,m,binR)) { cout << "{\"placements\":[]}"; return 0; }
    string bs=m[1]; smatch f;
    regex wr("\\\"W\\\"\\s*:\\s*([0-9]+)"), hr("\\\"H\\\"\\s*:\\s*([0-9]+)"), ar("\\\"allow_rotate\\\"\\s*:\\s*(true|false)");
    if(!regex_search(bs,f,wr)) { cout << "{\"placements\":[]}"; return 0; } W=stoi(f[1]);
    if(!regex_search(bs,f,hr)) { cout << "{\"placements\":[]}"; return 0; } H=stoi(f[1]);
    if(!regex_search(bs,f,ar)) { cout << "{\"placements\":[]}"; return 0; } allow=f[1]=="true";
    vector<Item> a;
    regex objR("\\{([^{}]*)\\}"), tr("\\\"type\\\"\\s*:\\s*\\\"([^\\\"]*)\\\""),
          irw("\\\"w\\\"\\s*:\\s*([0-9]+)"), irh("\\\"h\\\"\\s*:\\s*([0-9]+)"),
          irv("\\\"v\\\"\\s*:\\s*([0-9]+)"), irl("\\\"limit\\\"\\s*:\\s*([0-9]+)");
    for(sregex_iterator q(in.begin(),in.end(),objR),e;q!=e;++q) {
        string o=(*q)[1]; smatch z;
        if(regex_search(o,z,tr)) { string id=z[1];
            if(!regex_search(o,z,irw)) continue; int w=stoi(z[1]);
            if(!regex_search(o,z,irh)) continue; int h=stoi(z[1]);
            if(!regex_search(o,z,irv)) continue; long long v=stoll(z[1]);
            if(!regex_search(o,z,irl)) continue; int lim=stoi(z[1]);
            a.push_back({id,w,h,v,lim});
        }
    }
    Result ans;
    // The passes are alternate contour objectives, rather than random perturbations of one layout.
    const long double par[][3]={{1,0,0},{1,.08L,.02L},{1,.35L,.05L},{.72L,.10L,.08L},{.45L,.16L,.12L},{1.25L,.12L,.02L},{.82L,.65L,.16L},{.2L,.25L,.2L},{1.08L,1.2L,.1L},{.62L,1.8L,.25L}};
    for(int k=0;k<10;k++) { Result r=run(a,W,H,allow,par[k][0],par[k][1],par[k][2],k%4); if(r.val>ans.val) ans=move(r); }
    cout << "{\"placements\":[";
    for(size_t i=0;i<ans.p.size();++i) { auto p=ans.p[i]; if(i) cout<<','; cout<<"{\"type\":\""<<esc(a[p.t].id)<<"\",\"x\":"<<p.x<<",\"y\":"<<p.y<<",\"rot\":"<<p.r<<'}'; }
    cout << "]}";
}
