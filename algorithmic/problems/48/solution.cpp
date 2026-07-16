#include <bits/stdc++.h>
using namespace std;

struct P{ double x,y,z; };

static vector<P> small_pack(int n){
    vector<array<int,3>> bits;
    if(n==2) bits={{0,0,0},{1,1,1}};
    else if(n==3) bits={{0,0,0},{0,1,1},{1,0,1}};
    else if(n==4) bits={{0,0,0},{0,1,1},{1,0,1},{1,1,0}};
    else bits={{0,0,0},{0,0,1},{0,1,0},{0,1,1},{1,0,0},{1,0,1},{1,1,0},{1,1,1}};
    double d = 1.0;
    if(n==2) d=sqrt(3.0);
    else if(n<=4) d=sqrt(2.0);
    double r = d/(2.0+2.0*d);
    vector<P> v;
    for(int i=0;i<n;i++){
        auto b=bits[i];
        v.push_back({b[0]?1-r:r, b[1]?1-r:r, b[2]?1-r:r});
    }
    return v;
}

static vector<P> gen_fcc(double r){
    vector<P> v;
    if(r<=0 || r>=0.5) return v;
    const double a = 2.0*sqrt(2.0)*r;
    const double hi = 1.0-r+1e-12;
    const double bas[4][3]={{0,0,0},{0,0.5,0.5},{0.5,0,0.5},{0.5,0.5,0}};
    int M = (int)floor((1.0-2.0*r)/a)+2;
    for(int i=0;i<=M;i++) for(int j=0;j<=M;j++) for(int k=0;k<=M;k++){
        for(auto &b: bas){
            double x=r+a*(i+b[0]), y=r+a*(j+b[1]), z=r+a*(k+b[2]);
            if(x<=hi && y<=hi && z<=hi) v.push_back({x,y,z});
        }
    }
    return v;
}

static vector<P> gen_bcc(double r){
    vector<P> v;
    if(r<=0 || r>=0.5) return v;
    const double a = 4.0*r/sqrt(3.0);
    const double hi = 1.0-r+1e-12;
    const double bas[2][3]={{0,0,0},{0.5,0.5,0.5}};
    int M = (int)floor((1.0-2.0*r)/a)+2;
    for(int i=0;i<=M;i++) for(int j=0;j<=M;j++) for(int k=0;k<=M;k++){
        for(auto &b: bas){
            double x=r+a*(i+b[0]), y=r+a*(j+b[1]), z=r+a*(k+b[2]);
            if(x<=hi && y<=hi && z<=hi) v.push_back({x,y,z});
        }
    }
    return v;
}

static vector<P> gen_grid(int n){
    int bestM=INT_MAX, ba=1,bb=1,bc=n;
    int lim = (int)ceil(cbrt((double)n))+100;
    for(int a=1; a<=lim; ++a) for(int b=a; b<=lim; ++b){
        int c=(n+a*b-1)/(a*b); if(c<b) c=b;
        int m=max({a,b,c});
        if(1LL*a*b*c>=n && m<bestM){bestM=m; ba=a; bb=b; bc=c;}
    }
    vector<P> v; v.reserve(n);
    for(int i=0;i<ba && (int)v.size()<n;i++) for(int j=0;j<bb && (int)v.size()<n;j++) for(int k=0;k<bc && (int)v.size()<n;k++)
        v.push_back({(i+0.5)/ba,(j+0.5)/bb,(k+0.5)/bc});
    return v;
}

static double radius_of(const vector<P>& v){
    int n=v.size(); double r=1e100;
    for(auto&p:v) r=min(r, min({p.x,1-p.x,p.y,1-p.y,p.z,1-p.z}));
    for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){
        double dx=v[i].x-v[j].x, dy=v[i].y-v[j].y, dz=v[i].z-v[j].z;
        r=min(r, 0.5*sqrt(dx*dx+dy*dy+dz*dz));
    }
    return r;
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int n; if(!(cin>>n)) return 0;
    vector<P> best;
    if(n<=8) best=small_pack(n); else best=gen_grid(n);
    double bestR = radius_of(best);

    auto try_lattice = [&](auto gen){
        double lo=0, hi=0.5;
        vector<P> pts;
        for(int it=0; it<55; ++it){
            double mid=(lo+hi)/2;
            auto v=gen(mid);
            if((int)v.size()>=n){ lo=mid; pts.swap(v); }
            else hi=mid;
        }
        auto v=gen(lo*0.999999999);
        if((int)v.size()>=n){
            v.resize(n);
            double rr=radius_of(v);
            if(rr>bestR){ bestR=rr; best.swap(v); }
        }
    };
    try_lattice(gen_fcc);
    try_lattice(gen_bcc);

    cout.setf(ios::fixed); cout<<setprecision(12);
    for(auto &p: best) cout<<p.x<<' '<<p.y<<' '<<p.z<<'\n';
    return 0;
}
