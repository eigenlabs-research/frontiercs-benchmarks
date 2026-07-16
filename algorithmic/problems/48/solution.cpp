#include <bits/stdc++.h>
using namespace std;

struct P { double x,y,z; };
struct Candidate { vector<P> p; double r = -1; };

static vector<P> takeClosest(vector<P> a, int n) {
    auto key=[](const P& p) { return p.x*p.x+p.y*p.y+p.z*p.z; };
    nth_element(a.begin(), a.begin()+n, a.end(), [&](const P& u,const P& v) {
        double ku=key(u), kv=key(v);
        if (ku != kv) return ku < kv;
        if (u.x != v.x) return u.x < v.x;
        if (u.y != v.y) return u.y < v.y;
        return u.z < v.z;
    });
    a.resize(n);
    return a;
}

// Fit a point set with known nearest-neighbour separation d into the cube.
static Candidate fit(const vector<P>& a, double d) {
    double lo[3]={1e100,1e100,1e100}, hi[3]={-1e100,-1e100,-1e100};
    for (auto q:a) {
        lo[0]=min(lo[0],q.x); hi[0]=max(hi[0],q.x);
        lo[1]=min(lo[1],q.y); hi[1]=max(hi[1],q.y);
        lo[2]=min(lo[2],q.z); hi[2]=max(hi[2],q.z);
    }
    double mid[3]={(lo[0]+hi[0])*.5,(lo[1]+hi[1])*.5,(lo[2]+hi[2])*.5};
    double R=max({(hi[0]-lo[0])*.5,(hi[1]-lo[1])*.5,(hi[2]-lo[2])*.5});
    double scale=1.0/(2.0*R+d);
    Candidate c; c.r=scale*d*.5; c.p.reserve(a.size());
    for(auto q:a) c.p.push_back({.5+scale*(q.x-mid[0]),.5+scale*(q.y-mid[1]),.5+scale*(q.z-mid[2])});
    return c;
}

static vector<P> rotateSet(const vector<P>& a, double az, double tilt) {
    // Rz(az) Ry(tilt): a small deterministic orientation net for the finite cluster.
    double ca=cos(az), sa=sin(az), ct=cos(tilt), st=sin(tilt);
    vector<P> b; b.reserve(a.size());
    for(auto p:a) {
        double X=ct*p.x+st*p.z, Y=p.y, Z=-st*p.x+ct*p.z;
        b.push_back({ca*X-sa*Y, sa*X+ca*Y, Z});
    }
    return b;
}

static Candidate fcc(int n) {
    int L=(int)ceil(cbrt((double)n))+4;
    vector<P> a; a.reserve((2*L+1)*(2*L+1)*(2*L+1)/2);
    for(int i=-L;i<=L;i++) for(int j=-L;j<=L;j++) for(int k=-L;k<=L;k++)
        if(((i+j+k)&1)==0) a.push_back({(double)i,(double)j,(double)k});
    a=takeClosest(move(a),n);
    Candidate best;
    // Include a body-diagonal view as well as cube-aligned FCC termination.
    const double PI=acos(-1.0);
    for (double t: {0.0, acos(1.0/sqrt(3.0)), PI/4})
        for (double a0: {0.0, PI/8, PI/4}) {
            Candidate c=fit(rotateSet(a,a0,t),sqrt(2.0));
            if(c.r>best.r) best=move(c);
        }
    return best;
}

static Candidate hcp(int n) {
    // AB-stacked triangular layers.  Its bulk density equals FCC but its finite shells differ.
    int L=(int)ceil(cbrt((double)n))+4;
    const double sy=sqrt(3.0)/2.0, sz=sqrt(2.0/3.0);
    vector<P> a;
    a.reserve((2*L+1)*(2*L+1)*(2*L+1));
    for(int k=-L;k<=L;k++) for(int j=-L;j<=L;j++) for(int i=-L;i<=L;i++) {
        int odd=(k&1); // works for negative layers too
        double x=i + .5*(j&1) + (odd?.5:0.0);
        double y=j*sy + (odd?sqrt(3.0)/6.0:0.0);
        a.push_back({x,y,k*sz});
    }
    a=takeClosest(move(a),n);
    Candidate best;
    const double PI=acos(-1.0);
    // Orientations deliberately do not alter pair distances, only layer termination at faces.
    const double tilts[]={0.0, PI/12, PI/6, acos(1.0/sqrt(3.0)), 5*PI/12};
    for(double t:tilts) for(int q=0;q<6;q++) {
        Candidate c=fit(rotateSet(a,q*PI/12,t),1.0);
        if(c.r>best.r) best=move(c);
    }
    return best;
}

static Candidate grid(int n) {
    int bestM=INT_MAX, ba=1,bb=1,bc=n;
    for(int a=1;a*a*a<=n*2+2;a++) for(int b=a;;b++) {
        int c=(n+a*b-1)/(a*b);
        if(c<b) c=b;
        if(a*b*c>=n && c<bestM) bestM=c,ba=a,bb=b,bc=c;
        if(b>=c) break;
    }
    Candidate out; out.r=1.0/(2.0*bestM); out.p.reserve(n);
    for(int i=0;i<ba && (int)out.p.size()<n;i++) for(int j=0;j<bb && (int)out.p.size()<n;j++)
        for(int k=0;k<bc && (int)out.p.size()<n;k++)
            out.p.push_back({(i+.5)/ba,(j+.5)/bb,(k+.5)/bc});
    return out;
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int n; if(!(cin>>n)) return 0;
    Candidate best=grid(n);
    Candidate a=fcc(n); if(a.r>best.r) best=move(a);
    Candidate b=hcp(n); if(b.r>best.r) best=move(b);
    cout<<setprecision(17);
    for(const P& p:best.p) cout<<p.x<<' '<<p.y<<' '<<p.z<<'\n';
}
