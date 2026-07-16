#include <bits/stdc++.h>
using namespace std;

using P = array<double,3>;

static long long cubeCap(long long a, long long lim) {
    if (a <= 0) return 0;
    if (a > lim / a) return lim;
    long long v = a * a;
    if (v > lim / a) return lim;
    return v * a;
}

// Number of usable coordinates r+a*(i+s), i >= 0, in one direction.
static long long levels(double r, double a, double s) {
    double t = (1.0 - 2.0*r) / a - s;
    if (t < -1e-11) return 0;
    return max(0LL, (long long)floor(t + 1e-11) + 1);
}

static long long squareCap(long long a, long long lim) {
    if (a <= 0) return 0;
    return a > lim/a ? lim : a*a;
}

static long long bccCap(double r, long long lim) {
    double a = 4.0*r / sqrt(3.0);
    long long u = levels(r,a,0.0), v = levels(r,a,0.5);
    long long x = cubeCap(u,lim), y = cubeCap(v,lim);
    return min(lim, x > lim-y ? lim : x+y);
}
static long long fccCap(double r, long long lim) {
    double a = 2.0*sqrt(2.0)*r;
    long long u = levels(r,a,0.0), v = levels(r,a,0.5);
    long long x = cubeCap(u,lim);
    long long vv = squareCap(v,lim);
    long long uv = min(lim, u > 0 && vv > lim/u ? lim : u*vv);
    long long add = min(lim, 3 > 0 && uv > lim/3 ? lim : 3*uv);
    return min(lim, x > lim-add ? lim : x+add);
}

template<class F> static double bestRadius(int n, F cap) {
    double lo=0.0, hi=0.5;
    for (int it=0; it<80; ++it) {
        double mid=(lo+hi)*0.5;
        if (cap(mid,n) >= n) lo=mid; else hi=mid;
    }
    // Leave a small numerical gap from a lattice/boundary contact.
    return lo*(1.0-1e-12);
}

static void appendLattice(vector<P>& out, int n, double r, bool bcc) {
    double a = bcc ? 4.0*r/sqrt(3.0) : 2.0*sqrt(2.0)*r;
    vector<P> shifts;
    if (bcc) shifts={{0,0,0},{.5,.5,.5}};
    else shifts={{0,0,0},{0,.5,.5},{.5,0,.5},{.5,.5,0}};
    for (P s: shifts) {
        long long nx=levels(r,a,s[0]), ny=levels(r,a,s[1]), nz=levels(r,a,s[2]);
        for (long long i=0;i<nx && (int)out.size()<n;++i)
            for (long long j=0;j<ny && (int)out.size()<n;++j)
                for (long long k=0;k<nz && (int)out.size()<n;++k)
                    out.push_back({r+a*(i+s[0]),r+a*(j+s[1]),r+a*(k+s[2])});
    }
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n;
    if (!(cin>>n)) return 0;

    double rb=bestRadius(n,bccCap);
    double rf=bestRadius(n,fccCap);
    int q=1;
    while (1LL*q*q*q<n) ++q;
    double rg=1.0/(2.0*q);

    vector<P> ans;
    if (rb >= rf && rb >= rg) appendLattice(ans,n,rb,true);
    else if (rf >= rg) appendLattice(ans,n,rf,false);
    else {
        ans.reserve(n);
        for (int i=0;i<q && (int)ans.size()<n;++i)
            for (int j=0;j<q && (int)ans.size()<n;++j)
                for (int k=0;k<q && (int)ans.size()<n;++k)
                    ans.push_back({(i+.5)/q,(j+.5)/q,(k+.5)/q});
    }
    // The capacity formula and generator are intentionally identical; this is a
    // defensive fallback only for an exotic floating-point boundary case.
    if ((int)ans.size()!=n) {
        ans.clear();
        for (int i=0;i<q && (int)ans.size()<n;++i)
            for (int j=0;j<q && (int)ans.size()<n;++j)
                for (int k=0;k<q && (int)ans.size()<n;++k)
                    ans.push_back({(i+.5)/q,(j+.5)/q,(k+.5)/q});
    }
    cout<<setprecision(17);
    for (const P& p: ans) cout<<p[0]<<' '<<p[1]<<' '<<p[2]<<'\n';
}
