#include <bits/stdc++.h>
using namespace std;
struct It { string id; int w,h; long long v; int lim; };
struct Pl { int t,x,y,r,w,h; };
struct Node { int x,h; };
int W,H; bool rotok; vector<It> a;

// A contour (skyline) is a compact representation of all already occupied area:
// each new rectangle sits above the maximum contour under its horizontal span.
struct Sol {
    vector<Node> s; vector<int> used; vector<Pl> p; long long val=0;
    Sol(){ s={{0,0},{W,INT_MAX/4}}; used.assign(a.size(),0); }
    int split(int x) {
        int k=upper_bound(s.begin(),s.end(),x,[](int z,const Node& n){return z<n.x;})-s.begin()-1;
        if(s[k].x==x) return k;
        s.insert(s.begin()+k+1,{x,s[k].h}); return k+1;
    }
    int top(int x,int w) const {
        int e=x+w, z=0;
        int k=upper_bound(s.begin(),s.end(),x,[](int q,const Node& n){return q<n.x;})-s.begin()-1;
        for(; s[k].x<e; ++k) z=max(z,s[k].h);
        return z;
    }
    void put(int t,int x,int y,int r,int w,int h) {
        int l=split(x), rr=split(x+w);
        s[l].h=y+h;
        s.erase(s.begin()+l+1,s.begin()+rr);
        for(int i=1;i+1<(int)s.size();) {
            if(s[i-1].h==s[i].h) s.erase(s.begin()+i);
            else ++i;
        }
        p.push_back({t,x,y,r,w,h}); used[t]++; val+=a[t].v;
    }
};

Sol run(int mode) {
    Sol z; int cap=1200; // keeps contour search bounded on instances with huge limits
    for(int step=0; step<cap; ++step) {
        int bt=-1,bx=0,by=0,br=0,bw=0,bh=0; long double bq=-1e100L;
        for(int k=0;k+1<(int)z.s.size();++k) {
            int x=z.s[k].x; if(x>=W) continue;
            for(int t=0;t<(int)a.size();++t) if(z.used[t]<a[t].lim) {
                for(int r=0;r<(rotok?2:1);++r) {
                    if(r && a[t].w==a[t].h) continue;
                    int w=r?a[t].h:a[t].w, h=r?a[t].w:a[t].h;
                    if(x+w>W) continue;
                    int y=z.top(x,w); if(y+h>H) continue;
                    long double den=(long double)a[t].v/(w*(long double)h);
                    long double q;
                    if(mode==0) q=den;
                    else if(mode==1) q=a[t].v;
                    else if(mode==2) q=(long double)a[t].v/(w+h);
                    else if(mode==3) q=den*(1.0L+0.18L*(long double)min(w,h)/max(w,h));
                    else if(mode==4) q=den*(1.0L+0.12L*(long double)w/W);
                    else q=den*(1.0L+0.12L*(long double)h/H);
                    // Two modes deliberately use strict bottom-left placement first; the
                    // others are value-led skyline variants.
                    bool take=false;
                    if(bt<0) take=true;
                    else if(mode>=4) {
                        if(y<by || (y==by && (x<bx || (x==bx && q>bq)))) take=true;
                    } else if(q>bq+1e-18L || (fabsl(q-bq)<1e-18L && (y<by || (y==by&&x<bx)))) take=true;
                    if(take) bt=t,bx=x,by=y,br=r,bw=w,bh=h,bq=q;
                }
            }
        }
        if(bt<0) break;
        z.put(bt,bx,by,br,bw,bh);
    }
    return z;
}
int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    string in((istreambuf_iterator<char>(cin)),{});
    auto num=[&](const string& key)->long long { regex re("\\\""+key+"\\\"\\s*:\\s*([0-9]+)"); smatch m; return regex_search(in,m,re)?stoll(m[1]):0; };
    W=num("W"); H=num("H");
    regex rotateRe("\\\"allow_rotate\\\"\\s*:\\s*true");
    rotok=regex_search(in,rotateRe);
    regex re("\\{\\s*\\\"type\\\"\\s*:\\s*\\\"([^\\\"]*)\\\"\\s*,\\s*\\\"w\\\"\\s*:\\s*([0-9]+)\\s*,\\s*\\\"h\\\"\\s*:\\s*([0-9]+)\\s*,\\s*\\\"v\\\"\\s*:\\s*([0-9]+)\\s*,\\s*\\\"limit\\\"\\s*:\\s*([0-9]+)\\s*\\}");
    for(sregex_iterator it(in.begin(),in.end(),re), e;it!=e;++it) a.push_back({(*it)[1],stoi((*it)[2]),stoi((*it)[3]),stoll((*it)[4]),stoi((*it)[5])});
    Sol best; for(int m=0;m<6;m++){ Sol q=run(m); if(q.val>best.val) best=move(q); }
    cout<<"{\"placements\":[";
    for(size_t i=0;i<best.p.size();++i){ auto &q=best.p[i]; if(i) cout<<','; cout<<"{\"type\":\""<<a[q.t].id<<"\",\"x\":"<<q.x<<",\"y\":"<<q.y<<",\"rot\":"<<q.r<<'}'; }
    cout<<"]}\n";
}
