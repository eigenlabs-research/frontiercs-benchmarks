#include <bits/stdc++.h>
using namespace std;

struct Item { string id; int w,h,lim; long long v; };
struct R { int x,y,w,h; };
struct P { int t,x,y,rot; };
struct Answer { long long value=0; vector<P> p; };
int W,H; bool canrot; vector<Item> it;

static bool intersects(const R& a, int x,int y,int w,int h) {
    return a.x < x+w && x < a.x+a.w && a.y < y+h && y < a.y+a.h;
}
static bool contains(const R& a,const R& b) {
    return a.x<=b.x && a.y<=b.y && a.x+a.w>=b.x+b.w && a.y+a.h>=b.y+b.h;
}

// Remove the occupied part of every maximal empty rectangle.  Free rectangles may
// overlap each other, but none intersects an already placed item.
static void occupy(vector<R>& fr,int x,int y,int w,int h) {
    vector<R> q; q.reserve(fr.size()*2+4);
    for (const R& a: fr) {
        if (!intersects(a,x,y,w,h)) { q.push_back(a); continue; }
        int ax2=a.x+a.w, ay2=a.y+a.h, bx2=x+w, by2=y+h;
        if (x>a.x) q.push_back({a.x,a.y,x-a.x,a.h});
        if (bx2<ax2) q.push_back({bx2,a.y,ax2-bx2,a.h});
        if (y>a.y) q.push_back({a.x,a.y,a.w,y-a.y});
        if (by2<ay2) q.push_back({a.x,by2,a.w,ay2-by2});
    }
    vector<char> dead(q.size());
    for (int i=0;i<(int)q.size();++i) if(q[i].w>0 && q[i].h>0)
        for (int j=0;j<(int)q.size();++j) if(i!=j && !dead[i] && contains(q[j],q[i])) {
            // identical rectangles are removed only once
            if (q[j].x!=q[i].x || q[j].y!=q[i].y || q[j].w!=q[i].w || q[j].h!=q[i].h || j<i)
                dead[i]=1;
        }
    fr.clear(); fr.reserve(q.size());
    for(int i=0;i<(int)q.size();++i) if(!dead[i] && q[i].w>0 && q[i].h>0) fr.push_back(q[i]);
}

struct Choice { int t=-1,f=-1,rot=0,w=0,h=0; long double key=-1e100L; int s1=INT_MAX,s2=INT_MAX; };

static long double priority(const Item& a, int mode) {
    long double area=(long double)a.w*a.h;
    long double d=a.v/area;
    long double v=a.v;
    switch(mode) {
        case 0: return d;
        case 1: return v;
        case 2: return d + v/((long double)W*H)*0.45L;
        case 3: return v/sqrtl(area);
        case 4: return v/powl(area,0.72L);
        case 5: return d*powl((long double)min(a.w,a.h)/max(a.w,a.h),0.25L);
        default: return d*powl(area,0.12L);
    }
}

static Answer build(int mode,int first) {
    vector<R> fr(1,{0,0,W,H});
    vector<int> used(it.size());
    Answer ans;
    bool must = first>=0;
    // The cap is only a runtime guard. All provided item limits are normally much smaller.
    for(int step=0; step<1600 && !fr.empty(); ++step) {
        Choice best;
        for(int t=0;t<(int)it.size();++t) {
            if(used[t]>=it[t].lim || (must && t!=first)) continue;
            for(int ro=0;ro<=(canrot?1:0);++ro) {
                if(ro && it[t].w==it[t].h) continue;
                int w=ro?it[t].h:it[t].w, h=ro?it[t].w:it[t].h;
                long double k=priority(it[t],mode);
                for(int f=0;f<(int)fr.size();++f) if(w<=fr[f].w && h<=fr[f].h) {
                    int a=min(fr[f].w-w,fr[f].h-h), b=max(fr[f].w-w,fr[f].h-h);
                    bool take=false;
                    if(k>best.key+1e-18L) take=true;
                    else if(fabsl(k-best.key)<=1e-18L) {
                        if(a<best.s1 || (a==best.s1 && (b<best.s2 || (b==best.s2 && (fr[f].y<fr[best.f].y || (fr[f].y==fr[best.f].y && fr[f].x<fr[best.f].x)))))) take=true;
                    }
                    if(take) best={t,f,ro,w,h,k,a,b};
                }
            }
        }
        if(best.t<0) break;
        const R r=fr[best.f];
        ans.p.push_back({best.t,r.x,r.y,best.rot});
        ans.value+=it[best.t].v; ++used[best.t];
        occupy(fr,r.x,r.y,best.w,best.h);
        must=false;
    }
    return ans;
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    string s((istreambuf_iterator<char>(cin)),{});
    smatch m;
    regex bw("\\\"W\\\"\\s*:\\s*([0-9]+)");
    regex bh("\\\"H\\\"\\s*:\\s*([0-9]+)");
    regex br("\\\"allow_rotate\\\"\\s*:\\s*(true|false)");
    if(!regex_search(s,m,bw)) return 0; W=stoi(m[1]);
    regex_search(s,m,bh); H=stoi(m[1]);
    regex_search(s,m,br); canrot=m[1]=="true";
    regex ri("\\{\\s*\\\"type\\\"\\s*:\\s*\\\"([^\\\"]*)\\\"\\s*,\\s*\\\"w\\\"\\s*:\\s*([0-9]+)\\s*,\\s*\\\"h\\\"\\s*:\\s*([0-9]+)\\s*,\\s*\\\"v\\\"\\s*:\\s*([0-9]+)\\s*,\\s*\\\"limit\\\"\\s*:\\s*([0-9]+)\\s*\\}");
    for(sregex_iterator q(s.begin(),s.end(),ri), e; q!=e; ++q) {
        auto z=*q; it.push_back({z[1],stoi(z[2]),stoi(z[3]),stoi(z[5]),stoll(z[4])});
    }
    Answer best;
    // Different exponents express one mechanism: deciding how much area should be
    // reserved for scarce high-profit pieces before ordinary density filling.
    for(int mode=0;mode<5;++mode) {
        Answer a=build(mode,-1); if(a.value>best.value) best=move(a);
    }
    // A forced seed exposes layouts which greedy density cannot recover after it has
    // fragmented the initial empty bin.  Density filling after that one commitment
    // keeps this a focused test of the scarce-large-item hypothesis.
    for(int t=0;t<(int)it.size();++t) if(it[t].lim>0) {
        Answer a=build(0,t); if(a.value>best.value) best=move(a);
    }
    cout << "{\"placements\":[";
    for(int i=0;i<(int)best.p.size();++i) {
        if(i) cout<<',';
        const P& p=best.p[i];
        cout << "{\"type\":\"" << it[p.t].id << "\",\"x\":" << p.x
             << ",\"y\":" << p.y << ",\"rot\":" << p.rot << '}';
    }
    cout << "]}\n";
}
