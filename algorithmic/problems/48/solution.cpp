#include <bits/stdc++.h>
using namespace std;

struct Pt { double x,y,z; };

static vector<Pt> make_cubic(int n) {
    int bestA=1,bestB=1,bestC=n;
    double bestR=-1;
    int lim = (int)ceil(cbrt((double)n))*8 + 20;
    for (int a=1; a<=lim; ++a) for (int b=1; b<=lim; ++b) {
        int c = (n + a*b - 1)/(a*b);
        double r = 0.5 / max({a,b,c});
        if (r > bestR) bestR=r, bestA=a, bestB=b, bestC=c;
    }
    vector<Pt> p; p.reserve(n);
    double r = bestR;
    for (int i=0; i<bestA && (int)p.size()<n; ++i)
        for (int j=0; j<bestB && (int)p.size()<n; ++j)
            for (int k=0; k<bestC && (int)p.size()<n; ++k) {
                double x = (bestA==1)?0.5:(r + i*(1-2*r)/(bestA-1));
                double y = (bestB==1)?0.5:(r + j*(1-2*r)/(bestB-1));
                double z = (bestC==1)?0.5:(r + k*(1-2*r)/(bestC-1));
                p.push_back({x,y,z});
            }
    return p;
}

// Finite face-centred-cubic lattice.  Cubic cell side a gives nearest distance a/sqrt(2).
static vector<Pt> make_fcc_for(int n, double r, int mode) {
    const double a = 2.0 * sqrt(2.0) * r;
    const double lo = r, hi = 1.0 - r;
    const double eps = 1e-10;
    vector<array<double,3>> off = {{ {0,0,0}, {0,0.5,0.5}, {0.5,0,0.5}, {0.5,0.5,0} }};
    // A few harmless orderings/phase choices improve truncated cells near the boundary.
    if (mode&1) reverse(off.begin(), off.end());
    double phase = (mode>=2 ? 0.5*a : 0.0);
    vector<Pt> p; p.reserve(n);
    int M = (int)ceil((1.0-2.0*r + a)/a) + 3;
    for (int i=-1; i<=M && (int)p.size()<n; ++i) {
        for (int j=-1; j<=M && (int)p.size()<n; ++j) {
            for (int k=-1; k<=M && (int)p.size()<n; ++k) {
                for (auto o: off) {
                    double x = lo + phase + (i + o[0]) * a;
                    double y = lo + phase + (j + o[1]) * a;
                    double z = lo + phase + (k + o[2]) * a;
                    if (x >= lo-eps && x <= hi+eps && y >= lo-eps && y <= hi+eps && z >= lo-eps && z <= hi+eps) {
                        p.push_back({min(1.0,max(0.0,x)), min(1.0,max(0.0,y)), min(1.0,max(0.0,z))});
                        if ((int)p.size()==n) break;
                    }
                }
            }
        }
    }
    return p;
}

static vector<Pt> best_fcc(int n) {
    vector<Pt> best;
    double bestR = -1;
    for (int mode=0; mode<4; ++mode) {
        double L=0, U=0.5;
        for (int it=0; it<55; ++it) {
            double mid=(L+U)/2;
            if ((int)make_fcc_for(n, mid, mode).size() >= n) L=mid; else U=mid;
        }
        vector<Pt> p = make_fcc_for(n, L*(1-2e-12), mode);
        if ((int)p.size()>=n && L>bestR) bestR=L, best=p;
    }
    return best;
}

static double radius_of(const vector<Pt>& p) {
    double r = 1e100;
    int n=p.size();
    for (auto &q:p) r = min(r, min({q.x,1-q.x,q.y,1-q.y,q.z,1-q.z}));
    double md2 = 1e100;
    for (int i=0;i<n;i++) for (int j=i+1;j<n;j++) {
        double dx=p[i].x-p[j].x, dy=p[i].y-p[j].y, dz=p[i].z-p[j].z;
        md2 = min(md2, dx*dx+dy*dy+dz*dz);
    }
    return min(r, 0.5*sqrt(md2));
}

static vector<Pt> small_n(int n) {
    vector<Pt> p;
    if (n==2) {
        double s=sqrt(3.0), r=s/(2.0*(1.0+s));
        p={{r,r,r},{1-r,1-r,1-r}};
    } else if (n==3 || n==4) {
        // Three/four vertices of a regular tetrahedron on alternating cube corners, optimally inset.
        double s=sqrt(2.0), r=s/(2.0*(1.0+s));
        p={{r,r,r},{r,1-r,1-r},{1-r,r,1-r}};
        if (n==4) p.push_back({1-r,1-r,r});
    } else if (n==8) {
        double r=0.25; for(int i: {0,1}) for(int j:{0,1}) for(int k:{0,1}) p.push_back({r+i*0.5,r+j*0.5,r+k*0.5});
    }
    return p;
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int n; if(!(cin>>n)) return 0;
    vector<Pt> best = make_cubic(n);
    vector<Pt> fcc = best_fcc(n);
    if ((int)fcc.size()==n && radius_of(fcc) > radius_of(best)) best = fcc;
    vector<Pt> sm = small_n(n);
    if ((int)sm.size()==n && radius_of(sm) > radius_of(best)) best = sm;
    cout.setf(ios::fixed); cout<<setprecision(12);
    for (int i=0;i<n;i++) cout << best[i].x << ' ' << best[i].y << ' ' << best[i].z << '\n';
    return 0;
}
