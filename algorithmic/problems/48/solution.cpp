#include <bits/stdc++.h>
using namespace std;

struct P { int x, y, z; };
struct Candidate {
    vector<P> p;
    double d = 0, value = -1;
};

static Candidate scaledCandidate(vector<P> p, double d) {
    int lx = p[0].x, hx = p[0].x, ly = p[0].y, hy = p[0].y, lz = p[0].z, hz = p[0].z;
    for (const P &q : p) {
        lx=min(lx,q.x); hx=max(hx,q.x);
        ly=min(ly,q.y); hy=max(hy,q.y);
        lz=min(lz,q.z); hz=max(hz,q.z);
    }
    int w = max({hx-lx, hy-ly, hz-lz});
    Candidate c;
    c.p = move(p); c.d = d;
    // With lattice spacing h, the closest pair is d*h and the widest
    // coordinate span is w*h.  Thus h=1/(w+d) makes both constraints tight.
    c.value = d / (2.0 * (w + d));
    return c;
}

static vector<P> compactFCC(int n) {
    int k=0;
    while (true) {
        vector<P> v;
        for (int x=-k;x<=k;x++) for (int y=-k;y<=k;y++) for (int z=-k;z<=k;z++)
            if (((x+y+z)&1)==0) v.push_back({x,y,z});
        if ((int)v.size() >= n) {
            sort(v.begin(),v.end(),[](const P&a,const P&b) {
                int ca=max({abs(a.x),abs(a.y),abs(a.z)}), cb=max({abs(b.x),abs(b.y),abs(b.z)});
                if(ca!=cb) return ca<cb;
                int na=a.x*a.x+a.y*a.y+a.z*a.z, nb=b.x*b.x+b.y*b.y+b.z*b.z;
                if(na!=nb) return na<nb;
                if(a.x!=b.x) return a.x<b.x; if(a.y!=b.y) return a.y<b.y; return a.z<b.z;
            });
            v.resize(n); return v;
        }
        ++k;
    }
}

static Candidate boxCandidate(int a,int b,int c,int n) {
    vector<P> v;
    for(int x=0;x<a;x++) for(int y=0;y<b;y++) for(int z=0;z<c;z++) v.push_back({x,y,z});
    // Take a centered, nearly cubical prefix rather than a long lexicographic slab.
    sort(v.begin(),v.end(),[&](const P&u,const P&v) {
        int ku=max({abs(2*u.x-(a-1)),abs(2*u.y-(b-1)),abs(2*u.z-(c-1))});
        int kv=max({abs(2*v.x-(a-1)),abs(2*v.y-(b-1)),abs(2*v.z-(c-1))});
        if(ku!=kv) return ku<kv;
        if(u.x!=v.x) return u.x<v.x; if(u.y!=v.y) return u.y<v.y; return u.z<v.z;
    });
    v.resize(n);
    return scaledCandidate(move(v),1.0);
}

int main() {
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int n;
    if (!(cin >> n) || n < 2) return 0;
    vector<Candidate> all;
    all.push_back(scaledCandidate(compactFCC(n), sqrt(2.0)));

    int t=1; while (t*t*t<n) ++t;
    for(int a : {max(1,t-1),t}) for(int b : {max(1,t-1),t}) for(int c : {max(1,t-1),t})
        if (a*b*c>=n) all.push_back(boxCandidate(a,b,c,n));

    // Binary cube vertices are especially effective for the small cases.
    vector<P> corners={{0,0,0},{1,1,1},{0,1,1},{1,0,1},{0,0,1},{0,1,0},{1,0,0},{1,1,0}};
    if(n<=8) {
        vector<P> q(corners.begin(),corners.begin()+n);
        int md=3;
        for(int i=0;i<n;i++) for(int j=0;j<i;j++)
            md=min(md, abs(q[i].x-q[j].x)+abs(q[i].y-q[j].y)+abs(q[i].z-q[j].z));
        all.push_back(scaledCandidate(move(q),sqrt((double)md)));
    }
    Candidate best=all[0];
    for(const auto &c:all) if(c.value>best.value) best=c;

    int lx=best.p[0].x, ly=best.p[0].y, lz=best.p[0].z;
    int hx=lx,hy=ly,hz=lz;
    for(auto q:best.p) { lx=min(lx,q.x);hx=max(hx,q.x);ly=min(ly,q.y);hy=max(hy,q.y);lz=min(lz,q.z);hz=max(hz,q.z); }
    int w=max({hx-lx,hy-ly,hz-lz});
    double h=1.0/(w+best.d), r=best.d*h/2;
    cout<<setprecision(17);
    for(auto q:best.p)
        cout << r+(q.x-lx)*h << ' ' << r+(q.y-ly)*h << ' ' << r+(q.z-lz)*h << '\n';
}
