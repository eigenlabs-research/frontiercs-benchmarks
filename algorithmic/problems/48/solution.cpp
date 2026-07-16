#include <bits/stdc++.h>
using namespace std;

struct Candidate {
    vector<array<double,3>> p;
    double d, value;
};

static Candidate makeCandidate(vector<array<double,3>> p, double d) {
    double lo[3]={1e100,1e100,1e100}, hi[3]={-1e100,-1e100,-1e100};
    for (auto q:p) for(int k=0;k<3;k++) lo[k]=min(lo[k],q[k]), hi[k]=max(hi[k],q[k]);
    double e=0;
    for(int k=0;k<3;k++) e=max(e,hi[k]-lo[k]);
    // A uniform scale leaves a margin scale*d/2 at both opposite faces.
    double s=1.0/(e+d);
    for(auto &q:p) for(int k=0;k<3;k++) q[k]=0.5+s*(q[k]-(lo[k]+hi[k])*0.5);
    return {move(p),d,d*s*0.5};
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int n; if(!(cin>>n)) return 0;
    Candidate best; best.value=-1;
    auto consider=[&](vector<array<double,3>> p,double d){
        if((int)p.size()!=n) return;
        Candidate q=makeCandidate(move(p),d);
        if(q.value>best.value) best=move(q);
    };

    // The small binary cube codes remove the especially costly boundary shell of a lattice.
    if(n<=8){
        vector<int> code;
        if(n==2) code={0,7};
        else if(n==3) code={0,3,5};
        else if(n==4) code={0,3,5,6};
        else { code={0,3,5,6,1,2,4,7}; code.resize(n); }
        vector<array<double,3>> p;
        for(int v:code) p.push_back({double(v&1),double((v>>1)&1),double((v>>2)&1)});
        int mh=3;
        for(int i=0;i<n;i++) for(int j=0;j<i;j++) mh=min(mh,__builtin_popcount(code[i]^code[j]));
        consider(move(p),sqrt((double)mh));
    }

    int m=1; while((long long)m*m*m<n) ++m;
    // A centered order makes partially filled cells substantially less elongated than row order.
    {
        vector<array<double,3>> a;
        for(int x=0;x<m;x++) for(int y=0;y<m;y++) for(int z=0;z<m;z++) a.push_back({double(x),double(y),double(z)});
        double c=(m-1)*.5;
        sort(a.begin(),a.end(),[&](auto u,auto v){
            double du=(u[0]-c)*(u[0]-c)+(u[1]-c)*(u[1]-c)+(u[2]-c)*(u[2]-c);
            double dv=(v[0]-c)*(v[0]-c)+(v[1]-c)*(v[1]-c)+(v[2]-c)*(v[2]-c);
            return du<dv;
        });
        a.resize(n); consider(move(a),1.0);
    }

    // FCC: integer sites of even coordinate sum, whose nearest separation is sqrt(2).
    int mf=1; while(((long long)mf*mf*mf+1)/2<n) ++mf;
    for(int M=mf;M<=mf+2;M++) {
        vector<array<double,3>> a;
        for(int x=0;x<M;x++) for(int y=0;y<M;y++) for(int z=0;z<M;z++)
            if(((x+y+z)&1)==0) a.push_back({double(x),double(y),double(z)});
        if((int)a.size()<n) continue;
        double c=(M-1)*.5;
        sort(a.begin(),a.end(),[&](auto u,auto v){
            double du=(u[0]-c)*(u[0]-c)+(u[1]-c)*(u[1]-c)+(u[2]-c)*(u[2]-c);
            double dv=(v[0]-c)*(v[0]-c)+(v[1]-c)*(v[1]-c)+(v[2]-c)*(v[2]-c);
            return du<dv;
        });
        a.resize(n); consider(move(a),sqrt(2.0));
    }

    // Discriminating ablation: BCC has worse bulk density but a different finite boundary shell.
    auto bccCount=[](int M) { long long e=(M+1)/2, o=M/2; return e*e*e+o*o*o; };
    int mb=1; while(bccCount(mb)<n) ++mb;
    for(int M=mb;M<=mb+2;M++) {
        vector<array<double,3>> a;
        for(int x=0;x<M;x++) for(int y=0;y<M;y++) for(int z=0;z<M;z++)
            if((x&1)==(y&1) && (y&1)==(z&1)) a.push_back({double(x),double(y),double(z)});
        if((int)a.size()<n) continue;
        double c=(M-1)*.5;
        sort(a.begin(),a.end(),[&](auto u,auto v){
            double du=(u[0]-c)*(u[0]-c)+(u[1]-c)*(u[1]-c)+(u[2]-c)*(u[2]-c);
            double dv=(v[0]-c)*(v[0]-c)+(v[1]-c)*(v[1]-c)+(v[2]-c)*(v[2]-c);
            return du<dv;
        });
        a.resize(n); consider(move(a),sqrt(3.0));
    }
    cout<<setprecision(17);
    for(auto q:best.p) cout<<q[0]<<' '<<q[1]<<' '<<q[2]<<'\n';
}
