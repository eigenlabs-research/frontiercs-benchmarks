#include <bits/stdc++.h>
using namespace std;

struct P { double x,y,z; };
struct Mat { double a[3][3]; };
static P rot(const Mat& m, P p) {
    return {m.a[0][0]*p.x+m.a[0][1]*p.y+m.a[0][2]*p.z,
            m.a[1][0]*p.x+m.a[1][1]*p.y+m.a[1][2]*p.z,
            m.a[2][0]*p.x+m.a[2][1]*p.y+m.a[2][2]*p.z};
}
static Mat axisRot(int axis, double t) {
    Mat m{}; double c=cos(t), s=sin(t);
    for(int i=0;i<3;i++) m.a[i][i]=1;
    if(axis==0) { m.a[1][1]=c;m.a[1][2]=-s;m.a[2][1]=s;m.a[2][2]=c; }
    if(axis==1) { m.a[0][0]=c;m.a[0][2]=s;m.a[2][0]=-s;m.a[2][2]=c; }
    if(axis==2) { m.a[0][0]=c;m.a[0][1]=-s;m.a[1][0]=s;m.a[1][1]=c; }
    return m;
}
int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int n; if(!(cin>>n)) return 0;
    const double sy=sqrt(3.0)/2.0, sz=sqrt(2.0/3.0);
    vector<Mat> rots;
    rots.push_back(axisRot(0,0));
    for(int ax=0;ax<3;ax++) for(double d: {15.0,30.0,45.0}) rots.push_back(axisRot(ax,d*acos(-1.0)/180.0));
    vector<P> best;
    double bestScale=-1;
    int q=(int)ceil(cbrt((double)n));
    // The useful dimensions are close to q, but the extended range also covers
    // thin cases where an incomplete final HCP layer is advantageous.
    int lim=max(3,2*q+4);
    for(int A=1;A<=lim;A++) for(int B=1;B<=lim;B++) {
        int C=(n+A*B-1)/(A*B);
        if(C<1 || C>lim*2) continue;
        int total=A*B*C;
        // Do not spend the time budget on extremely anisotropic boxes: for
        // non-tiny n they cannot win a max-side-length cube scaling.
        if(n>=30) {
            double ex=A-.5, ey=(B-1)*sy+sy/3.0, ez=(C-1)*sz;
            double mn=min({ex,ey,ez}), mx=max({ex,ey,ez});
            if(mn<=0.0 || mx/mn>1.75) continue;
        }
        vector<pair<double,P>> v; v.reserve(total);
        double cx=(A-1)*.5, cy=(B-1)*sy*.5, cz=(C-1)*sz*.5;
        for(int k=0;k<C;k++) for(int j=0;j<B;j++) for(int i=0;i<A;i++) {
            double off=(k&1)?0.5:0.0;
            P p={(double)i+off, (double)j*sy+(k&1?sy/3.0:0.0), (double)k*sz};
            // Remove central points first.  The retained shell preserves every
            // extremal lattice point of this crop while allowing exact cardinality.
            double bd=min({(double)i,(double)(A-1-i),(double)j,(double)(B-1-j),(double)k,(double)(C-1-k)});
            double rr=(p.x-cx)*(p.x-cx)+(p.y-cy)*(p.y-cy)+(p.z-cz)*(p.z-cz);
            v.push_back({bd*1000000.0-rr,p});
        }
        if(total>n) {
            nth_element(v.begin(),v.begin()+n,v.end(),[](const auto& u,const auto& w){return u.first<w.first;});
            v.resize(n);
        }
        for(const Mat& m: rots) {
            double lo[3]={1e100,1e100,1e100}, hi[3]={-1e100,-1e100,-1e100};
            vector<P> w; w.reserve(n);
            for(auto &e:v) {
                P p=rot(m,e.second); w.push_back(p);
                lo[0]=min(lo[0],p.x); lo[1]=min(lo[1],p.y); lo[2]=min(lo[2],p.z);
                hi[0]=max(hi[0],p.x); hi[1]=max(hi[1],p.y); hi[2]=max(hi[2],p.z);
            }
            double span=max({hi[0]-lo[0],hi[1]-lo[1],hi[2]-lo[2]});
            double scale=1.0/(span+1.0);
            if(scale>bestScale) {
                bestScale=scale; best.clear(); best.reserve(n);
                for(P p:w) best.push_back({scale*(p.x-lo[0])+scale*.5,
                                             scale*(p.y-lo[1])+scale*.5,
                                             scale*(p.z-lo[2])+scale*.5});
            }
        }
    }
    // A finite HCP crop can lose to the checker’s cubic-grid baseline,
    // notably for small or strongly incomplete crops.  Keep that baseline
    // available rather than allowing a negative normalized contribution.
    double gridScale=1.0/q;
    bool useGrid=(gridScale>bestScale);
    cout<<setprecision(17);
    if(useGrid) {
        int put=0;
        for(int i=0;i<q && put<n;i++) for(int j=0;j<q && put<n;j++)
            for(int k=0;k<q && put<n;k++,put++)
                cout<<(i+.5)/q<<' '<<(j+.5)/q<<' '<<(k+.5)/q<<'\n';
    } else {
        for(P p:best) cout<<p.x<<' '<<p.y<<' '<<p.z<<'\n';
    }
}
