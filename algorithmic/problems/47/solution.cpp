#include <bits/stdc++.h>
using namespace std;
struct Item { string id; int w,h; long long v; int lim; };
struct Put { int t,x,y,r; };
struct Choice { int t,w,h,r; };
struct Row { vector<Put> p; vector<int> used; long long val=0; int h=0, width=0; };
int W,H; bool allowRot; vector<Item> a;

static string esc(const string& s) { string r; for(char c:s) { if(c=='"'||c=='\\') r+='\\'; r+=c; } return r; }

// In this representation a solution is a sequence of complete horizontal shelves.
// A shelf is seeded by each possible height, then its remaining width is a small
// one-dimensional greedy knapsack; this is deliberately a guillotine/shelf family.
Row makeRow(const vector<Choice>& c, int seed, const vector<int>& left, int fillMode) {
    Row z; z.h=c[seed].h; z.used.assign(a.size(),0);
    int x=0;
    auto add=[&](int k) { const Choice &q=c[k]; z.p.push_back({q.t,x,0,q.r}); z.used[q.t]++; z.val+=a[q.t].v; x+=q.w; };
    add(seed);
    while (true) {
        int best=-1; long double bs=-1;
        for(int k=0;k<(int)c.size();k++) { const Choice&q=c[k];
            if(q.h>z.h || q.w>W-x || z.used[q.t]>=left[q.t]) continue;
            long double s;
            if(fillMode==0) s=(long double)a[q.t].v/(q.w*q.h);            // area density
            else if(fillMode==1) s=(long double)a[q.t].v/q.w;              // valuable row width
            else if(fillMode==2) s=(long double)a[q.t].v/(q.w*z.h);        // shelf density
            else if(fillMode==3) s=(long double)a[q.t].v/(q.w*q.h)*1.0L + (long double)q.w/W*1e-8L;
            else s=(long double)a[q.t].v/(q.w*q.h) + (long double)q.h/H*1e-7L;
            if(s>bs+1e-18L || (fabsl(s-bs)<1e-18L && (best<0 || q.w>c[best].w))) bs=s,best=k;
        }
        if(best<0) break;
        add(best);
    }
    z.width=x; return z;
}

vector<Put> build(bool transpose, int fillMode, int rowMode) {
    int ow=W, oh=H; if(transpose) swap(W,H);
    vector<Choice> c;
    for(int i=0;i<(int)a.size();i++) for(int r=0;r<= (allowRot?1:0);r++) {
        int ww=r?a[i].h:a[i].w, hh=r?a[i].w:a[i].h;
        if(transpose) swap(ww,hh);
        if(ww<=W && hh<=H) c.push_back({i,ww,hh,r});
    }
    vector<int> left(a.size()); for(int i=0;i<(int)a.size();i++) left[i]=a[i].lim;
    vector<Put> ans; long long total=0; int y=0;
    while(y<H) {
        int best=-1; Row br; long double bs=-1;
        for(int s=0;s<(int)c.size();s++) {
            if(left[c[s].t]<=0 || c[s].h>H-y) continue;
            Row r=makeRow(c,s,left,fillMode);
            if(r.width==0) continue;
            long double score;
            if(rowMode==0) score=(long double)r.val/(r.h*W);              // filled shelf density
            else if(rowMode==1) score=(long double)r.val/r.h;              // value per consumed height
            else if(rowMode==2) score=(long double)r.val/(r.h*max(1,r.width));
            else score=(long double)r.val/(r.h*W) + (long double)r.width/W*1e-9L;
            if(score>bs+1e-18L || (fabsl(score-bs)<1e-18L && r.val>br.val)) bs=score,best=s,br=move(r);
        }
        if(best<0) break;
        for(auto q:br.p) { q.y=y; ans.push_back(q); left[q.t]--; total+=a[q.t].v; }
        y+=br.h;
    }
    // restore global dimensions before mapping output coordinates
    if(transpose) swap(W,H);
    if(transpose) for(auto &q:ans) { int nx=q.y, ny=q.x; q.x=nx; q.y=ny; }
    return ans;
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    string in((istreambuf_iterator<char>(cin)),{});
    smatch m;
    regex binR("\\\"W\\\"\\s*:\\s*([0-9]+)[\\s\\S]*?\\\"H\\\"\\s*:\\s*([0-9]+)[\\s\\S]*?\\\"allow_rotate\\\"\\s*:\\s*(true|false)");
    if(!regex_search(in,m,binR)) { cout<<"{\"placements\":[]}"; return 0; }
    W=stoi(m[1]); H=stoi(m[2]); allowRot=m[3]=="true";
    regex itemR("\\{\\s*\\\"type\\\"\\s*:\\s*\\\"([^\\\"]*)\\\"\\s*,\\s*\\\"w\\\"\\s*:\\s*([0-9]+)\\s*,\\s*\\\"h\\\"\\s*:\\s*([0-9]+)\\s*,\\s*\\\"v\\\"\\s*:\\s*([0-9]+)\\s*,\\s*\\\"limit\\\"\\s*:\\s*([0-9]+)\\s*\\}");
    for(sregex_iterator it(in.begin(),in.end(),itemR),e;it!=e;++it) { auto q=*it; a.push_back({q[1],stoi(q[2]),stoi(q[3]),stoll(q[4]),stoi(q[5])}); }
    vector<Put> best; long long bv=-1;
    // Multiple scoring lenses are cheap, but all candidates use the same shelf decomposition.
    for(int tr=0;tr<2;tr++) for(int f=0;f<5;f++) for(int rm=0;rm<4;rm++) {
        auto z=build(tr,f,rm); long long v=0; for(auto&p:z)v+=a[p.t].v;
        if(v>bv) bv=v,best=move(z);
    }
    cout<<"{\"placements\":[";
    for(size_t i=0;i<best.size();i++) { if(i) cout<<','; auto&p=best[i];
        cout<<"{\"type\":\""<<esc(a[p.t].id)<<"\",\"x\":"<<p.x<<",\"y\":"<<p.y<<",\"rot\":"<<p.r<<'}'; }
    cout<<"]}";
}
