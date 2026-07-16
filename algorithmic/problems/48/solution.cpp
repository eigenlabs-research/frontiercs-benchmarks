#include <bits/stdc++.h>
using namespace std;

struct Pt { double x,y,z; };

static double radius_of(const vector<Pt>& p){
    double r = 1e100;
    int n = (int)p.size();
    for (auto &a: p) {
        r = min(r, min(min(a.x, 1.0-a.x), min(min(a.y, 1.0-a.y), min(a.z, 1.0-a.z))));
    }
    for (int i=0;i<n;i++) for (int j=i+1;j<n;j++) {
        double dx=p[i].x-p[j].x, dy=p[i].y-p[j].y, dz=p[i].z-p[j].z;
        r = min(r, 0.5*sqrt(dx*dx+dy*dy+dz*dz));
    }
    return r;
}

static vector<Pt> cube_corners(int n, double r){
    vector<Pt> v;
    for (int mask=0; mask<8 && (int)v.size()<n; ++mask) {
        v.push_back({(mask&1)?1.0-r:r, (mask&2)?1.0-r:r, (mask&4)?1.0-r:r});
    }
    return v;
}

static vector<Pt> cubic_grid(int n){
    int ba=1,bb=1,bc=n, bestM=n, bestCap=n;
    int lim = 1;
    while (lim*lim*lim < n) ++lim;
    lim += 3;
    for (int a=1;a<=lim;a++) for (int b=1;b<=lim;b++) {
        int c = (n + a*b - 1)/(a*b);
        int M = max(a, max(b,c));
        int cap = a*b*c;
        if (M < bestM || (M==bestM && cap < bestCap)) {
            bestM=M; bestCap=cap; ba=a; bb=b; bc=c;
        }
    }
    vector<Pt> v; v.reserve(n);
    double r = 1.0/(2.0*bestM);
    auto coord = [&](int i, int cnt)->double{
        if (cnt==1) return 0.5;
        return r + (1.0-2.0*r) * (double)i / (double)(cnt-1);
    };
    for (int i=0;i<ba && (int)v.size()<n;i++)
        for (int j=0;j<bb && (int)v.size()<n;j++)
            for (int k=0;k<bc && (int)v.size()<n;k++)
                v.push_back({coord(i,ba), coord(j,bb), coord(k,bc)});
    return v;
}

static int parity_count(int nx,int ny,int nz,int par){
    int total = nx*ny*nz;
    if ((nx%2)==0 || (ny%2)==0 || (nz%2)==0) return total/2;
    int even = (total + 1) / 2; // (0,0,0) has even parity
    return par==0 ? even : total-even;
}

static vector<Pt> fcc_lattice(int n){
    int bestM = INT_MAX, bx=1, by=1, bz=1, bp=0, bestCap=INT_MAX;
    for (int M=0; M<=80; ++M) {
        bool any=false;
        for (int nx=1; nx<=M+1; ++nx) for (int ny=1; ny<=M+1; ++ny) for (int nz=1; nz<=M+1; ++nz) {
            if (max({nx,ny,nz}) != M+1) continue;
            for (int par=0; par<2; ++par) {
                int cap = parity_count(nx,ny,nz,par);
                if (cap >= n) {
                    any=true;
                    if (M < bestM || (M==bestM && cap < bestCap)) {
                        bestM=M; bestCap=cap; bx=nx; by=ny; bz=nz; bp=par;
                    }
                }
            }
        }
        if (any) break;
    }
    vector<Pt> v; v.reserve(n);
    double s = 1.0 / (bestM + sqrt(2.0));
    double r = s / sqrt(2.0);
    for(int i=0;i<bx && (int)v.size()<n;i++) for(int j=0;j<by && (int)v.size()<n;j++) for(int k=0;k<bz && (int)v.size()<n;k++) {
        if (((i+j+k)&1) == bp) v.push_back({r + i*s, r + j*s, r + k*s});
    }
    return v;
}

static vector<Pt> bcc_lattice(int n){
    int bestM=INT_MAX;
    for(int M=1; M<=80; ++M){
        long long cap = 1LL*(M+1)*(M+1)*(M+1) + 1LL*M*M*M;
        if(cap >= n){ bestM=M; break; }
    }
    int M=bestM;
    double a = 1.0 / (M + sqrt(3.0)/2.0);
    double r = sqrt(3.0)*a/4.0;
    vector<Pt> v; v.reserve(n);
    for(int i=0;i<=M && (int)v.size()<n;i++) for(int j=0;j<=M && (int)v.size()<n;j++) for(int k=0;k<=M && (int)v.size()<n;k++)
        v.push_back({r+i*a, r+j*a, r+k*a});
    for(int i=0;i<M && (int)v.size()<n;i++) for(int j=0;j<M && (int)v.size()<n;j++) for(int k=0;k<M && (int)v.size()<n;k++)
        v.push_back({r+(i+0.5)*a, r+(j+0.5)*a, r+(k+0.5)*a});
    return v;
}

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n; if(!(cin>>n)) return 0;
    vector<vector<Pt>> cand;
    if (n==1) cand.push_back({{0.5,0.5,0.5}});
    if (n==2) {
        double r = sqrt(3.0)/(2.0+2.0*sqrt(3.0));
        cand.push_back({{r,r,r},{1.0-r,1.0-r,1.0-r}});
    }
    if (n>=3 && n<=4) {
        double r = sqrt(2.0)/(2.0+2.0*sqrt(2.0));
        vector<Pt> t={{r,r,r},{r,1-r,1-r},{1-r,r,1-r},{1-r,1-r,r}};
        t.resize(n); cand.push_back(t);
    }
    if (n>=5 && n<=8) cand.push_back(cube_corners(n, 0.25));
    cand.push_back(cubic_grid(n));
    cand.push_back(fcc_lattice(n));
    cand.push_back(bcc_lattice(n));

    int bi=0; double br=-1;
    for(int i=0;i<(int)cand.size();++i) if((int)cand[i].size()==n) {
        double rr = radius_of(cand[i]);
        if(rr > br) { br=rr; bi=i; }
    }
    cout.setf(ios::fixed); cout << setprecision(17);
    for(auto &p: cand[bi]) {
        double x=min(1.0,max(0.0,p.x)), y=min(1.0,max(0.0,p.y)), z=min(1.0,max(0.0,p.z));
        cout << x << ' ' << y << ' ' << z << '\n';
    }
    return 0;
}
