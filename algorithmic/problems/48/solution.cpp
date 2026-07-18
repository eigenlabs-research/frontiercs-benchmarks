#include <bits/stdc++.h>
using namespace std;

struct Pt{ double x,y,z; };

static vector<Pt> cubic(int n){
    int bestA=1,bestB=1,bestC=n,bestM=INT_MAX,bestVol=INT_MAX;
    for(int a=1;a*a*a<=n*8+1000;a++){
        for(int b=a;b*b*a<=n*8+1000;b++){
            int c=(n+a*b-1)/(a*b);
            int m=max(c,b);
            int vol=a*b*c;
            if(m<bestM || (m==bestM && vol<bestVol)) bestM=m,bestA=a,bestB=b,bestC=c,bestVol=vol;
        }
    }
    vector<Pt> v; v.reserve(n);
    double r=1.0/(2.0*max({bestA,bestB,bestC}));
    for(int i=0;i<bestA && (int)v.size()<n;i++) for(int j=0;j<bestB && (int)v.size()<n;j++) for(int k=0;k<bestC && (int)v.size()<n;k++){
        double x = (bestA==1?0.5:r + i*(1-2*r)/max(1,bestA-1));
        double y = (bestB==1?0.5:r + j*(1-2*r)/max(1,bestB-1));
        double z = (bestC==1?0.5:r + k*(1-2*r)/max(1,bestC-1));
        v.push_back({x,y,z});
    }
    return v;
}

static double radius_of(const vector<Pt>& v){
    int n=v.size();
    double r=1e100;
    for(auto &p:v) r=min(r, min({p.x,1-p.x,p.y,1-p.y,p.z,1-p.z}));
    for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){
        double dx=v[i].x-v[j].x, dy=v[i].y-v[j].y, dz=v[i].z-v[j].z;
        r=min(r, 0.5*sqrt(dx*dx+dy*dy+dz*dz));
    }
    return r;
}

static vector<Pt> fcc_make(int n, double r, int mode){
    double s=sqrt(2.0)*r;                 // cubic cell side, nearest-neighbour distance sqrt(2)*s = 2r
    double L=1.0-2.0*r;
    vector<array<double,3>> basis;
    if(mode==0) basis={{0,0,0},{0.5,0.5,0},{0.5,0,0.5},{0,0.5,0.5}}; // FCC
    else basis={{0,0,0},{0.5,0.5,0},{0.5,0,0.5},{0,0.5,0.5},{0.5,0.5,0.5}}; // adds body centres: useful if sparse; binary search protects distance
    vector<Pt> pts;
    if(r<=0 || L<-1e-12) return pts;
    int M=(int)floor(L/s+2.0);
    for(int i=0;i<=M;i++) for(int j=0;j<=M;j++) for(int k=0;k<=M;k++) for(auto b:basis){
        double x=r+s*(i+b[0]), y=r+s*(j+b[1]), z=r+s*(k+b[2]);
        if(x<=1-r+1e-12 && y<=1-r+1e-12 && z<=1-r+1e-12) pts.push_back({x,y,z});
    }
    sort(pts.begin(), pts.end(), [](const Pt&a,const Pt&b){
        double da=(a.x-.5)*(a.x-.5)+(a.y-.5)*(a.y-.5)+(a.z-.5)*(a.z-.5);
        double db=(b.x-.5)*(b.x-.5)+(b.y-.5)*(b.y-.5)+(b.z-.5)*(b.z-.5);
        if(fabs(da-db)>1e-15) return da<db;
        if(a.x!=b.x) return a.x<b.x; if(a.y!=b.y) return a.y<b.y; return a.z<b.z;
    });
    if((int)pts.size()>n) pts.resize(n);
    return pts;
}

static vector<Pt> best_fcc(int n){
    vector<Pt> best;
    double br=-1;
    for(int mode=0; mode<2; ++mode){
        double lo=0, hi=0.5;
        for(int it=0; it<45; ++it){
            double mid=(lo+hi)/2;
            auto v=fcc_make(n, mid, mode);
            if((int)v.size()>=n) lo=mid; else hi=mid;
        }
        auto v=fcc_make(n, lo*0.999999999, mode);
        if((int)v.size()==n){ double rr=radius_of(v); if(rr>br) br=rr,best=v; }
    }
    return best;
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int n; if(!(cin>>n)) return 0;
    vector<Pt> ans = cubic(n);
    double bestR = (n<=6000 ? radius_of(ans) : 0.0);

    // Exact diagonal optimum for two spheres.
    if(n==2){
        double a=sqrt(3.0)/(2.0+2.0*sqrt(3.0));
        ans={{a,a,a},{1-a,1-a,1-a}};
        bestR=radius_of(ans);
    }

    auto f=best_fcc(n);
    if((int)f.size()==n){
        double rf = (n<=6000 ? radius_of(f) : 1.0); // official n <= 4096
        if(rf>bestR) ans=f, bestR=rf;
    }

    cout.setf(ios::fixed); cout<<setprecision(12);
    for(auto &p: ans) cout<<p.x<<' '<<p.y<<' '<<p.z<<'\n';
    return 0;
}
