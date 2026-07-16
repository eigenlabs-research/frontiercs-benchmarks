#include <bits/stdc++.h>
using namespace std;
using P=array<double,3>;
struct Cand { vector<P> p; double r=-1; };

static Cand fit(vector<P> a, double mind) {
    Cand z; if(a.empty()) return z;
    double lo[3]={1e100,1e100,1e100}, hi[3]={-1e100,-1e100,-1e100};
    for(auto q:a) for(int k=0;k<3;k++) lo[k]=min(lo[k],q[k]),hi[k]=max(hi[k],q[k]);
    double span=max({hi[0]-lo[0],hi[1]-lo[1],hi[2]-lo[2]});
    double den=span+mind;
    for(auto &q:a) for(int k=0;k<3;k++) q[k]=(q[k]-lo[k]+mind*.5)/den;
    z.p=move(a); z.r=mind/(2*den); return z;
}
static double radius(const vector<P>& a){
    double r=1;
    for(auto p:a) for(int k=0;k<3;k++) r=min(r,min(p[k],1-p[k]));
    for(int i=0;i<(int)a.size();i++) for(int j=0;j<i;j++){
        double s=0; for(int k=0;k<3;k++){double d=a[i][k]-a[j][k];s+=d*d;}
        r=min(r,.5*sqrt(s));
    }
    return r;
}
static double halton(int x,int b){ double f=1,r=0; while(x){f/=b;r+=f*(x%b);x/=b;}return r; }

// A feasibility-style overlap relaxation: unlike a lattice prefix it explicitly
// tests the finite-population explanation, and is used only where it is cheap.
static Cand relaxPacking(int n){
    Cand best;
    if(n>72) return best;
    for(int seed=0;seed<4;seed++){
        vector<P> start(n);
        for(int i=0;i<n;i++) start[i]={halton(i+1+37*seed,2),halton(i+1+53*seed,3),halton(i+1+71*seed,5)};
        double low=0, high=.5; vector<P> answer;
        for(int step=0;step<17;step++){
            double rr=(low+high)/2; vector<P> a=start;
            // Repeated projected repulsion is a certificate experiment: after it
            // finishes we measure all distances rather than trusting rr.
            for(int it=0;it<420;it++){
                vector<P> d(n,{0,0,0});
                for(int i=0;i<n;i++) for(int k=0;k<3;k++){
                    if(a[i][k]<rr) d[i][k]+=(rr-a[i][k])*.55;
                    if(a[i][k]>1-rr) d[i][k]-=(a[i][k]-(1-rr))*.55;
                }
                for(int i=0;i<n;i++) for(int j=0;j<i;j++){
                    double v[3],s=0; for(int k=0;k<3;k++){v[k]=a[i][k]-a[j][k];s+=v[k]*v[k];}
                    double q=sqrt(s)+1e-12, need=2*rr-q;
                    if(need>0) for(int k=0;k<3;k++){ double w=.32*need*v[k]/q; d[i][k]+=w;d[j][k]-=w; }
                }
                for(int i=0;i<n;i++) for(int k=0;k<3;k++) a[i][k]=min(1.0,max(0.0,a[i][k]+d[i][k]));
            }
            if(radius(a)+2e-5>=rr){low=rr;answer=a;start=a;} else high=rr;
        }
        if(!answer.empty()){ Cand c{answer,radius(answer)}; if(c.r>best.r) best=move(c); }
    }
    return best;
}
int main(){
    ios::sync_with_stdio(false);cin.tie(nullptr);
    int n;if(!(cin>>n)) return 0; Cand best;
    auto take=[&](vector<P> v,double d){
        sort(v.begin(),v.end(),[](const P&a,const P&b){return a[0]*a[0]+a[1]*a[1]+a[2]*a[2]<b[0]*b[0]+b[1]*b[1]+b[2]*b[2];});
        v.resize(n); Cand c=fit(move(v),d); if(c.r>best.r) best=move(c);
    };
    int q=1;while(q*q*q<n)q++;
    vector<P> g; for(int x=0;x<q;x++)for(int y=0;y<q;y++)for(int z=0;z<q;z++)g.push_back({(double)x,(double)y,(double)z});
    take(g,1);
    // Centered finite windows of the two cubic Bravais lattices.
    int L=1; while((2*L+1)*(2*L+1)*(2*L+1)<5*n) L++;
    vector<P> f,b;
    for(int x=-L;x<=L;x++)for(int y=-L;y<=L;y++)for(int z=-L;z<=L;z++){
        if(((x+y+z)&1)==0) f.push_back({(double)x,(double)y,(double)z});
        if(((x-y)&1)==0 && ((y-z)&1)==0) b.push_back({(double)x,(double)y,(double)z});
    }
    take(f,sqrt(2.0)); take(b,sqrt(3.0));
    Cand local=relaxPacking(n); if(local.r>best.r) best=move(local);
    cout<<setprecision(17); for(auto p:best.p) cout<<p[0]<<' '<<p[1]<<' '<<p[2]<<'\n';
}
