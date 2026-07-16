#include <bits/stdc++.h>
using namespace std;

struct P { double x, y, z; };

static vector<P> cubicGrid(int n) {
    int q=1; while (1LL*q*q*q<n) ++q;
    int a=q, b=q, c=(n + q*q-1)/(q*q);
    vector<P> v; v.reserve(n);
    for (int i=0; i<a && (int)v.size()<n; ++i)
        for (int j=0; j<b && (int)v.size()<n; ++j)
            for (int k=0; k<c && (int)v.size()<n; ++k)
                v.push_back({(i+.5)/a,(j+.5)/b,(k+.5)/c});
    return v;
}

// Points with a fixed parity of x+y+z form an FCC lattice with nearest distance sqrt(2).
static vector<P> fcc(int n, int parity, double& achieved) {
    int m=0;
    vector<array<int,3>> all;
    for (;;) {
        ++m;
        all.clear();
        for (int x=0;x<=m;++x) for (int y=0;y<=m;++y) for (int z=0;z<=m;++z)
            if (((x+y+z)&1)==parity) all.push_back({x,y,z});
        if ((int)all.size()>=n) break;
    }
    // Taking the points nearest the middle keeps partially filled lattice boxes compact.
    sort(all.begin(), all.end(), [m](const auto& a, const auto& b) {
        int da=(2*a[0]-m)*(2*a[0]-m)+(2*a[1]-m)*(2*a[1]-m)+(2*a[2]-m)*(2*a[2]-m);
        int db=(2*b[0]-m)*(2*b[0]-m)+(2*b[1]-m)*(2*b[1]-m)+(2*b[2]-m)*(2*b[2]-m);
        if (da!=db) return da<db;
        return a<b;
    });
    all.resize(n);
    int lx=all[0][0], ly=all[0][1], lz=all[0][2];
    int hx=lx, hy=ly, hz=lz;
    for (auto p: all) {
        lx=min(lx,p[0]); ly=min(ly,p[1]); lz=min(lz,p[2]);
        hx=max(hx,p[0]); hy=max(hy,p[1]); hz=max(hz,p[2]);
    }
    int w=max({hx-lx,hy-ly,hz-lz});
    int mind2=INT_MAX;
    for (int i=0;i<n;++i) for (int j=0;j<i;++j) {
        int dx=all[i][0]-all[j][0], dy=all[i][1]-all[j][1], dz=all[i][2]-all[j][2];
        mind2=min(mind2,dx*dx+dy*dy+dz*dz);
    }
    double d=sqrt((double)mind2);
    double r=d/(2.0*(w+d));
    achieved=r;
    double s=(1.0-2.0*r)/w;
    vector<P> out; out.reserve(n);
    for (auto p: all) out.push_back({r+s*(p[0]-lx), r+s*(p[1]-ly), r+s*(p[2]-lz)});
    return out;
}

int main() {
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int n;
    if (!(cin>>n)) return 0;
    int q=1; while (1LL*q*q*q<n) ++q;
    vector<P> best=cubicGrid(n);
    double bestR=.5/q;
    for (int parity=0; parity<2; ++parity) {
        double r;
        vector<P> v=fcc(n, parity, r);
        if (r>bestR) bestR=r, best.swap(v);
    }
    // The two-sphere case has a better body-diagonal realization than a finite FCC cell.
    if (n==2) {
        double r=sqrt(3.0)/(2.0*(1.0+sqrt(3.0)));
        vector<P> v={{r,r,r},{1-r,1-r,1-r}};
        if (r>bestR) best.swap(v);
    }
    cout<<setprecision(17);
    for (const P& p: best) cout<<p.x<<' '<<p.y<<' '<<p.z<<'\n';
}
