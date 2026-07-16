#include <bits/stdc++.h>
using namespace std;

struct Item { string id; int w,h,lim; long long v; };
struct Rect { int x,y,w,h; };
struct Pl { int t,x,y,r; };
struct Sol { long long value=0; vector<Pl> p; };

static bool intersects(const Rect&a,const Rect&b){
    return a.x < b.x+b.w && b.x < a.x+a.w && a.y < b.y+b.h && b.y < a.y+a.h;
}

/* Maximal free rectangles.  Free rectangles are allowed to overlap each other,
   but each individual one is empty; therefore any placement chosen in one is safe. */
static void cut(vector<Rect>& fr, const Rect& q){
    vector<Rect> a; a.reserve(fr.size()*2+4);
    for(const Rect& f:fr){
        if(!intersects(f,q)){ a.push_back(f); continue; }
        if(q.x>f.x) a.push_back({f.x,f.y,q.x-f.x,f.h});
        if(q.x+q.w<f.x+f.w) a.push_back({q.x+q.w,f.y,f.x+f.w-q.x-q.w,f.h});
        if(q.y>f.y) a.push_back({f.x,f.y,f.w,q.y-f.y});
        if(q.y+q.h<f.y+f.h) a.push_back({f.x,q.y+q.h,f.w,f.y+f.h-q.y-q.h});
    }
    vector<Rect> b; b.reserve(a.size());
    for(int i=0;i<(int)a.size();++i){
        if(a[i].w<=0 || a[i].h<=0) continue;
        bool inside=false;
        long long ai=1LL*a[i].w*a[i].h;
        for(int j=0;j<(int)a.size();++j) if(i!=j){
            const Rect& z=a[j]; long long aj=1LL*z.w*z.h;
            if(z.x<=a[i].x && z.y<=a[i].y && z.x+z.w>=a[i].x+a[i].w && z.y+z.h>=a[i].y+a[i].h &&
               (aj>ai || (aj==ai && j<i))){ inside=true; break; }
        }
        if(!inside) b.push_back(a[i]);
    }
    fr.swap(b);
}

static Sol build(const vector<Item>& it,int W,int H,bool rot,int policy){
    vector<Rect> fr={{0,0,W,H}};
    vector<int> used(it.size(),0);
    Sol out;
    // Different exponents deliberately trade density against individual profit.
    static const long double alpha[]={0.35L,0.55L,0.70L,0.85L,1.0L,1.12L,1.28L,0.90L,0.65L,1.05L,0.48L,1.20L};
    static const long double beta []={0.00L,0.08L,0.22L,0.42L,0.12L,0.32L,0.05L,0.55L,0.35L,0.18L,0.65L,0.48L};
    int k=policy%12, corner=policy%4;
    const long double binarea=(long double)W*H;
    while(true){
        long double best=-1e100L; int bi=-1,bf=-1,br=0,bw=0,bh=0;
        for(int i=0;i<(int)it.size();++i){
            if(used[i]>=it[i].lim || it[i].v<=0) continue;
            for(int r=0;r<=(rot?1:0);++r){
                int w=r?it[i].h:it[i].w, h=r?it[i].w:it[i].h;
                for(int f=0;f<(int)fr.size();++f) if(w<=fr[f].w && h<=fr[f].h){
                    long double area=(long double)w*h;
                    long double residual=((long double)fr[f].w*fr[f].h-area)/binarea;
                    long double key=log((long double)it[i].v+1.0L)-alpha[k]*log(area)-beta[k]*residual;
                    // Prefer a tight short-side fit only as a deterministic tie-scale term.
                    key-=1e-7L*(min(fr[f].w-w,fr[f].h-h));
                    if(key>best+1e-15L || (fabsl(key-best)<=1e-15L && it[i].v>(bi<0?-1:it[bi].v))){
                        best=key; bi=i; bf=f; br=r; bw=w; bh=h;
                    }
                }
            }
        }
        if(bi<0) break;
        Rect f=fr[bf];
        int x=(corner&1)?f.x+f.w-bw:f.x;
        int y=(corner&2)?f.y+f.h-bh:f.y;
        Rect q={x,y,bw,bh};
        out.p.push_back({bi,x,y,br}); out.value+=it[bi].v; ++used[bi];
        cut(fr,q);
    }
    return out;
}

static string esc(const string&s){
    string r; for(char c:s){ if(c=='"'||c=='\\') r+='\\'; r+=c; } return r;
}
int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    string in((istreambuf_iterator<char>(cin)),{});
    smatch m;
    regex rb("\\\"W\\\"\\s*:\\s*([0-9]+)\\s*,\\s*\\\"H\\\"\\s*:\\s*([0-9]+)\\s*,\\s*\\\"allow_rotate\\\"\\s*:\\s*(true|false)");
    if(!regex_search(in,m,rb)){ cout << "{\"placements\":[]}"; return 0; }
    int W=stoi(m[1]),H=stoi(m[2]); bool allow=m[3]=="true";
    regex ri("\\{\\s*\\\"type\\\"\\s*:\\s*\\\"([^\\\"]*)\\\"\\s*,\\s*\\\"w\\\"\\s*:\\s*([0-9]+)\\s*,\\s*\\\"h\\\"\\s*:\\s*([0-9]+)\\s*,\\s*\\\"v\\\"\\s*:\\s*([0-9]+)\\s*,\\s*\\\"limit\\\"\\s*:\\s*([0-9]+)\\s*\\}");
    vector<Item> it;
    for(sregex_iterator q(in.begin(),in.end(),ri),e;q!=e;++q){
        const smatch& z=*q; it.push_back({z[1],stoi(z[2]),stoi(z[3]),stoi(z[5]),stoll(z[4])});
    }
    Sol ans;
    // The twelve deterministic policies differ only in priorities and corner anchoring.
    for(int p=0;p<12;++p){ Sol s=build(it,W,H,allow,p); if(s.value>ans.value) ans=move(s); }
    cout << "{\"placements\":[";
    for(size_t i=0;i<ans.p.size();++i){
        if(i) cout<<','; const Pl& q=ans.p[i];
        cout << "{\"type\":\"" << esc(it[q.t].id) << "\",\"x\":" << q.x << ",\"y\":" << q.y << ",\"rot\":" << q.r << '}';
    }
    cout << "]}";
}
