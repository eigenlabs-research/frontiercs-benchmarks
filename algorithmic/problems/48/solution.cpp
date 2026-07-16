#include <bits/stdc++.h>
using namespace std;
using P = array<double,3>;

double radius(const vector<P>& a) {
    double r=1;
    for (auto p:a) for(int q=0;q<3;q++) r=min(r,min(p[q],1-p[q]));
    for (int i=0;i<(int)a.size();i++) for(int j=0;j<i;j++) {
        double s=0; for(int q=0;q<3;q++){ double d=a[i][q]-a[j][q]; s+=d*d; }
        r=min(r,.5*sqrt(s));
    }
    return r;
}
vector<P> grid(int n) {
    int lim=1; while(lim*lim*lim<n) ++lim;
    int bx=1,by=1,bz=n, best=INT_MAX;
    for(int x=1;x<=lim+1;x++) for(int y=1;y<=lim+1;y++) {
        int z=(n+x*y-1)/(x*y), mx=max(x,max(y,z));
        if(mx<best || (mx==best && x*y*z<bx*by*bz)) bx=x,by=y,bz=z,best=mx;
    }
    vector<P> a; a.reserve(n);
    for(int i=0;i<bx && (int)a.size()<n;i++) for(int j=0;j<by && (int)a.size()<n;j++) for(int k=0;k<bz && (int)a.size()<n;k++)
        a.push_back({(i+.5)/bx,(j+.5)/by,(k+.5)/bz});
    return a;
}
// A finite piece of the face-centred cubic lattice.  Farthest-point sampling
// makes its truncated boundary much less biased than simply taking a prefix.
vector<P> fcc(int n) {
    int K=0; vector<array<int,3>> all;
    do {
        ++K; all.clear();
        for(int x=0;x<=K;x++) for(int y=0;y<=K;y++) for(int z=0;z<=K;z++)
            if(((x+y+z)&1)==0) all.push_back({x,y,z});
    } while((int)all.size()<n);
    int M=all.size(); vector<int> take; take.reserve(n); vector<double> near(M,1e100); vector<char> used(M);
    // Start at a corner; subsequent choices are the points farthest from those kept.
    take.push_back(0); used[0]=1;
    for(int it=1;it<n;it++) {
        int last=take.back();
        for(int j=0;j<M;j++) { double s=0; for(int q=0;q<3;q++){double d=all[j][q]-all[last][q];s+=d*d;} near[j]=min(near[j],s); }
        int b=-1; for(int j=0;j<M;j++) if(!used[j] && (b<0 || near[j]>near[b])) b=j;
        take.push_back(b); used[b]=1;
    }
    double md=1e100; int lo[3]={INT_MAX,INT_MAX,INT_MAX}, hi[3]={INT_MIN,INT_MIN,INT_MIN};
    for(int u:take) for(int q=0;q<3;q++) lo[q]=min(lo[q],all[u][q]),hi[q]=max(hi[q],all[u][q]);
    for(int i=0;i<n;i++) for(int j=0;j<i;j++) { double s=0;for(int q=0;q<3;q++){double d=all[take[i]][q]-all[take[j]][q];s+=d*d;}md=min(md,sqrt(s)); }
    double span=0; for(int q=0;q<3;q++) span=max(span,double(hi[q]-lo[q]));
    double scale=1.0/(span+md); vector<P> out;
    for(int u:take) { P p; for(int q=0;q<3;q++) p[q]=.5+scale*(all[u][q]-(lo[q]+hi[q])*.5); out.push_back(p); }
    return out;
}
int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int n; if(!(cin>>n)) return 0;
    vector<P> best=grid(n); double br=radius(best);
    vector<P> b=fcc(n); double rr=radius(b); if(rr>br){best=b;br=rr;}
    // For small instances, polish both lattice constructions with a deterministic
    // annealing search.  This is intentionally only a supplement: large cases use
    // the O(n^2)-free lattice construction.
    if(n<=20) {
        mt19937_64 gen(0x9e3779b97f4a7c15ULL+n);
        vector<P> cur=best; double cr=br;
        int rounds=70000;
        for(int it=0;it<rounds;it++) {
            int v=gen()%n, q=gen()%3;
            P old=cur[v];
            double t=0.012*(1.0-double(it)/rounds)+0.00015;
            double step=(0.16*(1.0-double(it)/rounds)+0.003);
            cur[v][q]=min(0.999999999999,max(0.000000000001,old[q]+(double(gen()%20001)/10000.0-1.0)*step));
            double nr=radius(cur);
            double u=double(gen()%1000000)/1000000.0;
            if(nr>=cr || u<exp((nr-cr)/t)) { cr=nr; if(cr>br){br=cr;best=cur;} }
            else cur[v]=old;
        }
    }
    cout<<setprecision(17);
    for(auto p:best) cout<<p[0]<<' '<<p[1]<<' '<<p[2]<<'\n';
}
