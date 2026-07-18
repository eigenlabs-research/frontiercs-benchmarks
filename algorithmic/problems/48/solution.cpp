#include <bits/stdc++.h>
using namespace std;

struct P { double x,y,z; };

static double radius_of(const vector<P>& v){
    double r = 1e100;
    int n = (int)v.size();
    for (auto &p: v) {
        r = min(r, min(min(p.x, 1.0-p.x), min(p.y, min(1.0-p.y, min(p.z, 1.0-p.z)))));
    }
    double md2 = 1e100;
    for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){
        double dx=v[i].x-v[j].x, dy=v[i].y-v[j].y, dz=v[i].z-v[j].z;
        md2 = min(md2, dx*dx+dy*dy+dz*dz);
    }
    if(n>=2) r = min(r, 0.5*sqrt(md2));
    return r;
}

static vector<P> best_cubic(int n){
    int ba=1,bb=1,bc=n, bestM=n, bestProd=INT_MAX;
    int lim = (int)ceil(cbrt((double)n))+100;
    for(int a=1;a<=lim;a++) for(int b=a;b<=lim;b++){
        int c = (n + a*b - 1)/(a*b);
        if(c < b) c = b;
        int M = c, prod = a*b*c;
        if(M < bestM || (M==bestM && prod < bestProd)){
            bestM=M; bestProd=prod; ba=a; bb=b; bc=c;
        }
    }
    vector<int> dims = {ba,bb,bc};
    sort(dims.begin(), dims.end(), greater<int>());
    int M = dims[0];
    double r = 1.0/(2.0*M);
    vector<vector<double>> xs(3);
    for(int d=0; d<3; d++){
        int c=dims[d]; xs[d].reserve(c);
        if(c==1) xs[d].push_back(0.5);
        else for(int i=0;i<c;i++) xs[d].push_back(r + (1.0-2.0*r)*i/(c-1));
    }
    vector<P> v; v.reserve(n);
    for(double x: xs[0]) for(double y: xs[1]) for(double z: xs[2]){
        if((int)v.size()<n) v.push_back({x,y,z});
    }
    return v;
}

static vector<P> cube_vertices(int n){
    vector<P> verts;
    for(int mask=0; mask<8; mask++) verts.push_back({double(mask&1), double((mask>>1)&1), double((mask>>2)&1)});
    double bestd=-1; int bestmask=0;
    for(int mask=0; mask<256; mask++) if(__builtin_popcount((unsigned)mask)==n){
        double md=1e9;
        for(int i=0;i<8;i++) if(mask>>i&1) for(int j=i+1;j<8;j++) if(mask>>j&1){
            double dx=verts[i].x-verts[j].x, dy=verts[i].y-verts[j].y, dz=verts[i].z-verts[j].z;
            md=min(md, sqrt(dx*dx+dy*dy+dz*dz));
        }
        if(md>bestd){ bestd=md; bestmask=mask; }
    }
    double r = bestd / (2.0*(1.0+bestd));
    vector<P> v;
    for(int i=0;i<8;i++) if(bestmask>>i&1){
        v.push_back({r + (1-2*r)*verts[i].x, r + (1-2*r)*verts[i].y, r + (1-2*r)*verts[i].z});
    }
    return v;
}

static long long fcc_count(double r){
    if(r<=0) return (long long)4e18;
    double a = sqrt(2.0)*r;
    int L = (int)floor((1.0-2.0*r)/a + 1e-12);
    if(L<0) return 0;
    long long m = (long long)L+1, total=m*m*m;
    return (total+1)/2; // larger parity class
}

static vector<P> fcc_lattice(int n){
    double lo=0, hi=0.5;
    for(int it=0; it<70; it++){
        double mid=(lo+hi)/2;
        if(fcc_count(mid)>=n) lo=mid; else hi=mid;
    }
    double r=lo*0.999999999; // tiny safety against floor/rounding at thresholds
    double a=sqrt(2.0)*r;
    int L=(int)floor((1.0-2.0*r)/a + 1e-12);
    long long cnt[2]={0,0};
    for(int i=0;i<=L;i++) for(int j=0;j<=L;j++) for(int k=0;k<=L;k++) cnt[(i+j+k)&1]++;
    int par = cnt[1] > cnt[0];
    double used = a*L;
    double off = (1.0-used)/2.0;
    vector<P> all; all.reserve((size_t)min<long long>(cnt[par], n+1000));
    for(int i=0;i<=L;i++) for(int j=0;j<=L;j++) for(int k=0;k<=L;k++) if(((i+j+k)&1)==par){
        all.push_back({off+a*i, off+a*j, off+a*k});
    }
    sort(all.begin(), all.end(), [](const P& A, const P& B){
        double da=(A.x-.5)*(A.x-.5)+(A.y-.5)*(A.y-.5)+(A.z-.5)*(A.z-.5);
        double db=(B.x-.5)*(B.x-.5)+(B.y-.5)*(B.y-.5)+(B.z-.5)*(B.z-.5);
        return da < db;
    });
    if((int)all.size()>n) all.resize(n);
    return all;
}

static long long bcc_count(double r){
    if(r<=0) return (long long)4e18;
    double a = 4.0*r/sqrt(3.0);
    double t = (1.0-2.0*r)/a;
    if(t < -1e-14) return 0;
    long long L0 = (long long)floor(t + 1e-12) + 1;
    long long L1 = (long long)floor(t - 0.5 + 1e-12) + 1;
    if(L1<0) L1=0;
    return L0*L0*L0 + L1*L1*L1;
}

static vector<P> bcc_lattice(int n){
    double lo=0, hi=0.5;
    for(int it=0; it<70; it++){
        double mid=(lo+hi)/2;
        if(bcc_count(mid)>=n) lo=mid; else hi=mid;
    }
    double r=lo*0.999999999, a=4.0*r/sqrt(3.0);
    double t=(1.0-2.0*r)/a;
    int L0=(int)floor(t+1e-12), L1=(int)floor(t-0.5+1e-12);
    double used = a*max((double)L0, (L1>=0? L1+0.5 : 0.0));
    double off=(1.0-used)/2.0;
    vector<P> all;
    for(int i=0;i<=L0;i++) for(int j=0;j<=L0;j++) for(int k=0;k<=L0;k++) all.push_back({off+a*i,off+a*j,off+a*k});
    for(int i=0;i<=L1;i++) for(int j=0;j<=L1;j++) for(int k=0;k<=L1;k++) all.push_back({off+a*(i+.5),off+a*(j+.5),off+a*(k+.5)});
    sort(all.begin(), all.end(), [](const P& A, const P& B){
        double da=(A.x-.5)*(A.x-.5)+(A.y-.5)*(A.y-.5)+(A.z-.5)*(A.z-.5);
        double db=(B.x-.5)*(B.x-.5)+(B.y-.5)*(B.y-.5)+(B.z-.5)*(B.z-.5);
        return da < db;
    });
    if((int)all.size()>n) all.resize(n);
    return all;
}

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n; if(!(cin>>n)) return 0;
    vector<vector<P>> cand;
    cand.push_back(best_cubic(n));
    if(n<=8) cand.push_back(cube_vertices(n));
    cand.push_back(fcc_lattice(n));
    cand.push_back(bcc_lattice(n));
    int bi=0; double br=-1;
    for(int i=0;i<(int)cand.size();i++) if((int)cand[i].size()==n){
        double rr=radius_of(cand[i]);
        if(rr>br){ br=rr; bi=i; }
    }
    cout.setf(ios::fixed); cout<<setprecision(12);
    for(auto &p: cand[bi]) cout << p.x << ' ' << p.y << ' ' << p.z << '\n';
    return 0;
}
