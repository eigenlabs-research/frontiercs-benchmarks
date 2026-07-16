#include <bits/stdc++.h>
using namespace std;
struct P { double x,y,z; };

static double halton(int v,int b) {
    double f=1, r=0;
    while(v) { f/=b; r+=f*(v%b); v/=b; }
    return r;
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int n; if(!(cin>>n)) return 0;
    vector<P> a(n);
    // For large instances this is a safe, inexpensive volumetric seed.  The
    // continuous optimizer below is deliberately used on the sizes where its
    // all-pairs interaction is affordable.
    if(n>320){
        int nx=ceil(cbrt((double)n)), ny=nx, nz=nx;
        while((long long)nx*ny*nz<n) ++nz;
        int k=0;
        for(int z=0;z<nz && k<n;z++) for(int y=0;y<ny && k<n;y++) for(int x=0;x<nx && k<n;x++)
            a[k++]={(x+.5)/nx,(y+.5)/ny,(z+.5)/nz};
    } else {
        // Radical-inverse seeds avoid the axis-aligned symmetry of a lattice.
        for(int i=0;i<n;i++) a[i]={.04+.92*halton(i+1,2),.04+.92*halton(i+1,3),.04+.92*halton(i+1,5)};
        vector<P> f(n), best=a;
        double bestR=-1;
        auto radius = [&](const vector<P>& v) {
            double r=.5;
            for(auto p:v) r=min(r,min(min(p.x,1-p.x),min(p.y,1-p.y)));
            for(int u=0;u<n;u++) for(int w=u+1;w<n;w++) {
                double x=v[u].x-v[w].x, y=v[u].y-v[w].y, z=v[u].z-v[w].z;
                r=min(r,.5*sqrt(x*x+y*y+z*z));
            }
            return r;
        };
        const int rounds = n<20 ? 900 : (n<100 ? 520 : 260);
        for(int it=0;it<rounds;it++){
            double q = 3.0 + 7.0*(double)it/(rounds-1);
            double wall = pow(2.0,-q-1.0);
            fill(f.begin(),f.end(),P{0,0,0});
            for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){
                double dx=a[i].x-a[j].x, dy=a[i].y-a[j].y, dz=a[i].z-a[j].z;
                double d2=dx*dx+dy*dy+dz*dz+1e-18;
                double c=q*pow(d2,-.5*q-1.0);
                f[i].x+=c*dx; f[i].y+=c*dy; f[i].z+=c*dz;
                f[j].x-=c*dx; f[j].y-=c*dy; f[j].z-=c*dz;
            }
            for(int i=0;i<n;i++){
                auto addwall=[&](double v, double &g){ g += wall*q*(pow(v,-q-1)-pow(1-v,-q-1)); };
                addwall(a[i].x,f[i].x); addwall(a[i].y,f[i].y); addwall(a[i].z,f[i].z);
            }
            // Unit-length projected descent makes the method stable despite
            // the changing inverse-power scale.
            double step=.028*(1.0-(double)it/rounds)+.0015;
            for(int i=0;i<n;i++){
                double s=sqrt(f[i].x*f[i].x+f[i].y*f[i].y+f[i].z*f[i].z);
                if(s>0){
                    a[i].x += step*f[i].x/s; a[i].y += step*f[i].y/s; a[i].z += step*f[i].z/s;
                }
                a[i].x=min(.999999999,max(.000000001,a[i].x));
                a[i].y=min(.999999999,max(.000000001,a[i].y));
                a[i].z=min(.999999999,max(.000000001,a[i].z));
            }
            if(it%8==7 || it+1==rounds) {
                double r=radius(a);
                if(r>bestR) bestR=r, best=a;
            }
        }
        a=best;
    }
    cout<<setprecision(17);
    for(auto p:a) cout<<p.x<<' '<<p.y<<' '<<p.z<<'\n';
}
