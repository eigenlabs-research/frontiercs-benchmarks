#include <bits/stdc++.h>
using namespace std;

struct Item { string id; long long w,h,v,lim; };
struct Seg { long long l,r,y; };
struct Put { int t, rot; long long x,y; };
// Range assignment/max tree makes contour queries logarithmic even after many skyline splits.
struct Tree {
    int n; vector<long long> mx, tag;
    Tree(int N):n(N),mx(4*N+8),tag(4*N+8,-1) {}
    void push(int p){ if(tag[p]>=0) { mx[p*2]=mx[p*2+1]=tag[p]; tag[p*2]=tag[p*2+1]=tag[p]; tag[p]=-1; } }
    void set(int p,int l,int r,int a,int b,long long v){ if(b<=l||r<=a)return; if(a<=l&&r<=b){mx[p]=v;tag[p]=v;return;} push(p);int m=(l+r)/2;set(p*2,l,m,a,b,v);set(p*2+1,m,r,a,b,v);mx[p]=max(mx[p*2],mx[p*2+1]); }
    long long get(int p,int l,int r,int a,int b){ if(b<=l||r<=a)return 0; if(a<=l&&r<=b)return mx[p];push(p);int m=(l+r)/2;return max(get(p*2,l,m,a,b),get(p*2+1,m,r,a,b)); }
    void set(long long a,long long b,long long v){set(1,0,n,(int)a,(int)b,v);} long long get(long long a,long long b){return get(1,0,n,(int)a,(int)b);}
};

static string unescape(const string& s) {
    string r; bool e=false;
    for(char c:s) { if(e) { if(c=='n') r+='\n'; else if(c=='t') r+='\t'; else r+=c; e=false; } else if(c=='\\') e=true; else r+=c; }
    return r;
}
static bool fieldNum(const string& s,const string& k,long long& z) {
    regex re("\\\""+k+"\\\"\\s*:\\s*(-?[0-9]+)"); smatch m;
    if(!regex_search(s,m,re)) return false; z=stoll(m[1]); return true;
}
static bool fieldStr(const string& s,const string& k,string& z) {
    regex re("\\\""+k+"\\\"\\s*:\\s*\\\"((?:\\\\.|[^\\\"])*)\\\""); smatch m;
    if(!regex_search(s,m,re)) return false; z=unescape(m[1]); return true;
}
static long long topAt(const vector<Seg>& a,long long x,long long w) {
    long long z=0, e=x+w;
    for(const auto& q:a) { if(q.r<=x) continue; if(q.l>=e) break; z=max(z,q.y); }
    return z;
}
static void raiseSky(vector<Seg>& a,long long l,long long r,long long y) {
    vector<Seg> b; b.reserve(a.size()+2);
    for(const auto& q:a) {
        if(q.r<=l || q.l>=r) b.push_back(q);
        else {
            if(q.l<l) b.push_back({q.l,l,q.y});
            if(q.r>r) b.push_back({r,q.r,q.y});
        }
    }
    b.push_back({l,r,y});
    sort(b.begin(),b.end(),[](const Seg&A,const Seg&B){return A.l<B.l;});
    vector<Seg> c; c.reserve(b.size());
    for(auto q:b) {
        if(!c.empty() && c.back().r==q.l && c.back().y==q.y) c.back().r=q.r;
        else c.push_back(q);
    }
    a.swap(c);
}
static unsigned mix(unsigned x) { x^=x>>16; x*=0x7feb352dU; x^=x>>15; x*=0x846ca68bU; return x^(x>>16); }

static pair<long long,vector<Put>> solve(long long W,long long H,bool rot,const vector<Item>& it,int mode) {
    int n=it.size(); vector<long long> used(n); vector<Seg> sky(1,{0,W,0}); Tree contour((int)W); vector<Put> ans;
    long long value=0, binArea=max(1LL,W*H);
    for(int step=0; step<50000; ++step) {
        int bt=-1, br=0; long long bx=0,by=0,bw=0,bh=0;
        long double best=-1e100L; long long bestTop=LLONG_MAX, bestX=LLONG_MAX;
        for(int t=0;t<n;t++) if(used[t]<it[t].lim) {
            for(int rr=0;rr<=(rot?1:0);rr++) {
                long long w=rr?it[t].h:it[t].w, h=rr?it[t].w:it[t].h;
                if(w>W||h>H) continue;
                long long py=LLONG_MAX, px=LLONG_MAX;
                for(const auto& s:sky) {
                    long long x=s.l; if(x+w>W) continue;
                    long long y=contour.get(x,x+w); if(y+h>H) continue;
                    if(py==LLONG_MAX || y+h<py+h || (y+h==py+h && x<px)) py=y,px=x;
                }
                if(px==LLONG_MAX) continue;
                long double area=(long double)w*h, den=(long double)it[t].v/area;
                long double q;
                if(mode==0) q=den;
                else if(mode==1) q=(long double)it[t].v;
                else if(mode==2) q=den*sqrt(area);
                else if(mode==3) q=den/sqrt(area);
                else if(mode==4) q=(long double)it[t].v/(max(w,h));
                else if(mode==5) q=den*(1.0L+0.35L*(long double)min(w,h)/max(w,h));
                else {
                    unsigned z=mix((unsigned)(t*11939+mode*7919));
                    q=den*(0.78L+0.44L*(z%10000)/9999.0L);
                }
                // A tiny compactness tie break only; objective rank remains the main decision.
                if(q>best+1e-14L || (fabsl(q-best)<1e-14L && (py+h<bestTop || (py+h==bestTop && px<bestX)))) {
                    best=q; bt=t;br=rr;bx=px;by=py;bw=w;bh=h;bestTop=py+h;bestX=px;
                }
            }
        }
        if(bt<0) break;
        ans.push_back({bt,br,bx,by}); used[bt]++; value+=it[bt].v;
        raiseSky(sky,bx,bx+bw,by+bh); contour.set(bx,bx+bw,by+bh);
    }
    return {value,ans};
}
static string esc(const string& s) { string r; for(char c:s) { if(c=='"'||c=='\\') r+='\\'; r+=c; } return r; }
int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    string in((istreambuf_iterator<char>(cin)),{}); long long W=0,H=0;
    fieldNum(in,"W",W); fieldNum(in,"H",H);
    bool allow = regex_search(in,regex("\\\"allow_rotate\\\"\\s*:\\s*true"));
    vector<Item> a;
    regex obj("\\{[^{}]*\\}");
    for(sregex_iterator p(in.begin(),in.end(),obj),e;p!=e;++p) {
        string z=(*p).str(), id; long long w,h,v,l;
        if(fieldStr(z,"type",id)&&fieldNum(z,"w",w)&&fieldNum(z,"h",h)&&fieldNum(z,"v",v)&&fieldNum(z,"limit",l)) a.push_back({id,w,h,v,l});
    }
    pair<long long,vector<Put>> best={0,{}};
    for(int mode=0;mode<7;mode++) { auto r=solve(W,H,allow,a,mode); if(r.first>best.first) best=move(r); }
    cout << "{\"placements\":[";
    for(size_t i=0;i<best.second.size();i++) { auto&p=best.second[i]; if(i) cout<<','; cout<<"{\"type\":\""<<esc(a[p.t].id)<<"\",\"x\":"<<p.x<<",\"y\":"<<p.y<<",\"rot\":"<<p.rot<<'}'; }
    cout << "]}\n";
}
