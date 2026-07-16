#include <bits/stdc++.h>
using namespace std;

struct P { double x,y,z; };
static const double rt2=sqrt(2.0);
static const double B[4][3]={{0,0,0},{0,.5,.5},{.5,0,.5},{.5,.5,0}};

// Number of integer i for which r <= a*(i+b+s) <= 1-r.
static long long oneCount(double r, double b, double s) {
    double a=2.0*rt2*r;
    double lo=r/a-b-s, hi=(1.0-r)/a-b-s;
    long long L=(long long)ceil(lo-1e-12);
    long long R=(long long)floor(hi+1e-12);
    return max(0LL,R-L+1);
}
static long long capacity(double r, double sx, double sy, double sz) {
    long long ans=0;
    for(int q=0;q<4;q++) {
        long long x=oneCount(r,B[q][0],sx), y=oneCount(r,B[q][1],sy), z=oneCount(r,B[q][2],sz);
        // The official sizes make this small; saturating also makes the predicate safe.
        if(x && y && z) ans += min(1000000007LL,x*y*z);
    }
    return ans;
}
static vector<P> makeFCC(int n,double r,double sx,double sy,double sz) {
    vector<P> v; v.reserve(n);
    double a=2.0*rt2*r;
    for(int q=0;q<4 && (int)v.size()<n;q++) {
        double bb[3]={B[q][0]+sx,B[q][1]+sy,B[q][2]+sz};
        long long L[3],R[3];
        for(int d=0;d<3;d++) {
            L[d]=(long long)ceil(r/a-bb[d]-1e-12);
            R[d]=(long long)floor((1-r)/a-bb[d]+1e-12);
        }
        for(long long i=L[0];i<=R[0] && (int)v.size()<n;i++)
        for(long long j=L[1];j<=R[1] && (int)v.size()<n;j++)
        for(long long k=L[2];k<=R[2] && (int)v.size()<n;k++) {
            P p{a*(i+bb[0]),a*(j+bb[1]),a*(k+bb[2])};
            // The predicate above is deliberately conservative, but roundoff should never leak outside.
            p.x=min(1.0,max(0.0,p.x)); p.y=min(1.0,max(0.0,p.y)); p.z=min(1.0,max(0.0,p.z));
            v.push_back(p);
        }
    }
    return v;
}
int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int n; if(!(cin>>n)) return 0;
    vector<P> out;
    if(n==2) {
        double r=sqrt(3.0)/(2.0*(1.0+sqrt(3.0)));
        out={{r,r,r},{1-r,1-r,1-r}};
    } else if(n==3 || n==4) {
        double r=sqrt(2.0)/(2.0*(1.0+sqrt(2.0)));
        out={{r,r,r},{r,1-r,1-r},{1-r,r,1-r},{1-r,1-r,r}};
        out.resize(n);
    } else {
        // Certified grid fallback: its radius is precisely the published simple baseline.
        int q=1; while(1LL*q*q*q<n) ++q;
        double bestR=1.0/(2.0*q), bsx=0,bsy=0,bsz=0;
        // Translating the FCC lattice is the only experimental degree of freedom here.
        // Eighth-cell phases cover distinct ways its boundary shells meet the cube.
        for(int ix=0;ix<8;ix++) for(int iy=0;iy<8;iy++) for(int iz=0;iz<8;iz++) {
            double sx=ix/8.0, sy=iy/8.0, sz=iz/8.0;
            double lo=0.0, hi=0.5;
            for(int it=0;it<58;it++) {
                double mid=(lo+hi)*.5;
                if(capacity(mid,sx,sy,sz)>=n) lo=mid; else hi=mid;
            }
            if(lo>bestR*(1.0+1e-10)) bestR=lo,bsx=sx,bsy=sy,bsz=sz;
        }
        if(bestR>1.0/(2.0*q)*(1.0+1e-10)) out=makeFCC(n,bestR*(1.0-2e-14),bsx,bsy,bsz);
        else {
            out.reserve(n);
            for(int i=0;i<q && (int)out.size()<n;i++) for(int j=0;j<q && (int)out.size()<n;j++) for(int k=0;k<q && (int)out.size()<n;k++)
                out.push_back({(i+.5)/q,(j+.5)/q,(k+.5)/q});
        }
    }
    cout<<setprecision(17);
    for(auto p:out) cout<<p.x<<' '<<p.y<<' '<<p.z<<'\n';
}
