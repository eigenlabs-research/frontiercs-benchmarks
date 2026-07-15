#include <bits/stdc++.h>
using namespace std;

struct P { double x, y, z; };

static double nearestDistance(const vector<P>& a, double guaranteed) {
    if (a.size() > 1024) return guaranteed;
    double d2 = numeric_limits<double>::infinity();
    for (size_t i = 0; i < a.size(); ++i)
        for (size_t j = 0; j < i; ++j) {
            double dx=a[i].x-a[j].x, dy=a[i].y-a[j].y, dz=a[i].z-a[j].z;
            d2=min(d2, dx*dx+dy*dy+dz*dz);
        }
    return sqrt(d2);
}

static vector<P> fit(vector<P> a, double guaranteed) {
    double mn[3]={1e100,1e100,1e100}, mx[3]={-1e100,-1e100,-1e100};
    for (auto p:a) {
        mn[0]=min(mn[0],p.x); mn[1]=min(mn[1],p.y); mn[2]=min(mn[2],p.z);
        mx[0]=max(mx[0],p.x); mx[1]=max(mx[1],p.y); mx[2]=max(mx[2],p.z);
    }
    double span[3]={mx[0]-mn[0],mx[1]-mn[1],mx[2]-mn[2]};
    double wide=max(span[0],max(span[1],span[2]));
    double d=nearestDistance(a, guaranteed), s=1.0/(wide+d);
    for (auto &p:a) {
        double v[3]={p.x,p.y,p.z};
        for (int q=0;q<3;q++) v[q]=(v[q]-mn[q]+(wide-span[q])*0.5+d*0.5)*s;
        p={v[0],v[1],v[2]};
    }
    return a;
}

static double value(const vector<P>& a) {
    double face=1;
    for(auto p:a) face=min(face,min(p.x,min(1-p.x,min(p.y,min(1-p.y,min(p.z,1-p.z))))));
    double d2=1e100;
    for(size_t i=0;i<a.size();++i) for(size_t j=0;j<i;++j) {
        double x=a[i].x-a[j].x,y=a[i].y-a[j].y,z=a[i].z-a[j].z;
        d2=min(d2,x*x+y*y+z*z);
    }
    return min(face, .5*sqrt(d2));
}

int main() {
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int n; if (!(cin>>n)) return 0;
    if (n==2) { // The diameter joins opposite cube corners; balance it against face clearance.
        double r=sqrt(3.0)/(2.0*(1.0+sqrt(3.0)));
        cout<<setprecision(17)<<r<<' '<<r<<' '<<r<<'\n'<<1-r<<' '<<1-r<<' '<<1-r<<'\n';
        return 0;
    }
    int R=1;
    while ((2*R+1LL)*(2*R+1)*(2*R+1) < 5LL*n) ++R;
    vector<P> hcp, cubic, fcc, bcc;
    const double sy=sqrt(3.0)/2.0, hz=sqrt(2.0/3.0);
    for(int k=-R;k<=R;k++) for(int j=-R;j<=R;j++) for(int i=-R;i<=R;i++) {
        // ABAB close-packed triangular layers, with unit nearest-neighbour distance.
        double ox=(k&1)?0.5:0.0, oy=(k&1)?sqrt(3.0)/6.0:0.0;
        hcp.push_back({i+0.5*j+ox, sy*j+oy, hz*k});
        cubic.push_back({(double)i,(double)j,(double)k});
        if (((i+j+k)&1)==0) fcc.push_back({(double)i,(double)j,(double)k});
        if (((i&1)==(j&1)) && ((j&1)==(k&1))) bcc.push_back({(double)i,(double)j,(double)k});
    }
    auto shell=[](vector<P>& v, int n) {
        stable_sort(v.begin(),v.end(),[](const P&a,const P&b) {
            return a.x*a.x+a.y*a.y+a.z*a.z < b.x*b.x+b.y*b.y+b.z*b.z;
        });
        v.resize(n);
    };
    shell(hcp,n); shell(cubic,n); shell(fcc,n); shell(bcc,n);
    vector<vector<P>> candidates;
    candidates.push_back(fit(hcp,1.0));
    candidates.push_back(fit(cubic,1.0));
    candidates.push_back(fit(fcc,sqrt(2.0)));
    candidates.push_back(fit(bcc,sqrt(3.0)));
    int best=0;
    for(int i=1;i<(int)candidates.size();++i)
        if(value(candidates[i])>value(candidates[best])) best=i;
    cout<<setprecision(17);
    for(auto p:candidates[best]) cout<<p.x<<' '<<p.y<<' '<<p.z<<'\n';
}
