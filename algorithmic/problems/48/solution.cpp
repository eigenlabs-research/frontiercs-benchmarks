#include <bits/stdc++.h>
using namespace std;
struct P { double x,y,z; };
static unsigned long long rngstate;
static unsigned long long nextrnd(){ rngstate+=0x9e3779b97f4a7c15ULL; unsigned long long z=rngstate; z=(z^(z>>30))*0xbf58476d1ce4e5b9ULL; z=(z^(z>>27))*0x94d049bb133111ebULL; return z^(z>>31); }
static double urand(){ return (nextrnd()>>11)*(1.0/9007199254740992.0); }
static double radiusOf(const vector<P>& a){
    int n=a.size(); double ans=.5;
    for(auto p:a) ans=min(ans,min(min(p.x,1-p.x),min(p.y,1-p.y)));
    for(int i=0;i<n;i++) for(int j=0;j<i;j++) { double x=a[i].x-a[j].x,y=a[i].y-a[j].y,z=a[i].z-a[j].z; ans=min(ans,.5*sqrt(x*x+y*y+z*z)); }
    return ans;
}
static vector<P> relax(int n, int seed) {
    rngstate=0x123456789abcdefULL+seed*0x9e3779b97f4a7c15ULL;
    vector<P>a(n);
    // A scrambled uniform cloud deliberately has no lattice/cell occupancy structure.
    for(auto &p:a) p={urand(),urand(),urand()};
    double base=pow((double)n,-1.0/3.0);
    for(int stage=0;stage<15;stage++) {
        double want=base*(.16+.027*stage); // progressively test larger equal-sphere margins
        double reach=2*want;
        // Cell width is at least the interaction diameter, so 27 neighboring cells are complete.
        int C=max(1,(int)floor(1.0/reach));
        vector<vector<int>> box(C*C*C);
        for(int it=0;it<42;it++) {
            for(auto &v:box) v.clear();
            for(int i=0;i<n;i++) { int x=min(C-1,(int)(a[i].x*C)),y=min(C-1,(int)(a[i].y*C)),z=min(C-1,(int)(a[i].z*C)); box[(x*C+y)*C+z].push_back(i); }
            vector<P> f(n,{0,0,0});
            for(int ix=0;ix<C;ix++) for(int iy=0;iy<C;iy++) for(int iz=0;iz<C;iz++) {
                auto &v=box[(ix*C+iy)*C+iz];
                for(int dx=-1;dx<=1;dx++) for(int dy=-1;dy<=1;dy++) for(int dz=-1;dz<=1;dz++) { int X=ix+dx,Y=iy+dy,Z=iz+dz; if(X<0||Y<0||Z<0||X>=C||Y>=C||Z>=C) continue; auto&w=box[(X*C+Y)*C+Z];
                    for(int ii:v) for(int jj:w) if(ii>jj) { double x=a[ii].x-a[jj].x,y=a[ii].y-a[jj].y,z=a[ii].z-a[jj].z; double d=sqrt(x*x+y*y+z*z); if(d<reach) { if(d<1e-12){ x=urand()-.5;y=urand()-.5;z=urand()-.5;d=sqrt(x*x+y*y+z*z); } double q=(reach-d)/d*.42; f[ii].x+=x*q;f[ii].y+=y*q;f[ii].z+=z*q; f[jj].x-=x*q;f[jj].y-=y*q;f[jj].z-=z*q; }
                    }
                }
            }
            for(int i=0;i<n;i++) { // hard wall violations encode the same radius as pair constraints
                if(a[i].x<want) f[i].x+=(want-a[i].x)*.85; if(a[i].x>1-want) f[i].x-=(a[i].x-(1-want))*.85;
                if(a[i].y<want) f[i].y+=(want-a[i].y)*.85; if(a[i].y>1-want) f[i].y-=(a[i].y-(1-want))*.85;
                if(a[i].z<want) f[i].z+=(want-a[i].z)*.85; if(a[i].z>1-want) f[i].z-=(a[i].z-(1-want))*.85;
                a[i].x=min(1.0,max(0.0,a[i].x+f[i].x)); a[i].y=min(1.0,max(0.0,a[i].y+f[i].y)); a[i].z=min(1.0,max(0.0,a[i].z+f[i].z));
            }
        }
    }
    return a;
}
int main(){ ios::sync_with_stdio(false);cin.tie(nullptr); int n;if(!(cin>>n))return 0; vector<P> best; double br=-1; int tries=(n<=128?4:1); for(int s=0;s<tries;s++){auto q=relax(n,s);double r=radiusOf(q);if(r>br)br=r,best.swap(q);} cout<<setprecision(17); for(auto p:best) cout<<p.x<<' '<<p.y<<' '<<p.z<<'\n'; }
