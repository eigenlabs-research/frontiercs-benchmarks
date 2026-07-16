#include <bits/stdc++.h>
using namespace std;

struct Candidate {
    vector<array<int,3>> p;
    double d = 0;
    double value = -1;
};

// A box of points from one parity class of Z^3 is an FCC lattice:
// its shortest possible nonzero vector has length sqrt(2).
static Candidate fccBox(int a, int b, int c, int parity, int n) {
    Candidate q;
    q.d = sqrt(2.0);
    for (int x=0; x<a && (int)q.p.size()<n; ++x)
        for (int y=0; y<b && (int)q.p.size()<n; ++y)
            for (int z=0; z<c && (int)q.p.size()<n; ++z)
                if (((x+y+z)&1) == parity) q.p.push_back({x,y,z});
    if ((int)q.p.size() != n) { q.value = -1; return q; }
    int lo[3]={INT_MAX,INT_MAX,INT_MAX}, hi[3]={INT_MIN,INT_MIN,INT_MIN};
    for (auto v:q.p) for(int k=0;k<3;k++) lo[k]=min(lo[k],v[k]),hi[k]=max(hi[k],v[k]);
    int span=max({hi[0]-lo[0],hi[1]-lo[1],hi[2]-lo[2]});
    q.value=q.d/(2.0*(span+q.d));
    return q;
}

static Candidate gridBox(int a, int b, int c, int n) {
    Candidate q;
    q.d=1.0;
    for (int x=0;x<a && (int)q.p.size()<n;x++)
        for (int y=0;y<b && (int)q.p.size()<n;y++)
            for (int z=0;z<c && (int)q.p.size()<n;z++) q.p.push_back({x,y,z});
    if ((int)q.p.size()!=n) return q;
    int lo[3]={INT_MAX,INT_MAX,INT_MAX},hi[3]={INT_MIN,INT_MIN,INT_MIN};
    for(auto v:q.p) for(int k=0;k<3;k++) lo[k]=min(lo[k],v[k]),hi[k]=max(hi[k],v[k]);
    int span=max({hi[0]-lo[0],hi[1]-lo[1],hi[2]-lo[2]});
    q.value=1.0/(2.0*(span+1));
    return q;
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int n; if(!(cin>>n)) return 0;
    Candidate best;
    int dg=ceil(cbrt((double)n))+3;
    // The ordinary lattice is retained whenever a partial FCC box is unfavorable.
    for(int a=1;a<=dg;a++) for(int b=1;b<=dg;b++) for(int c=1;c<=dg;c++)
        if(1LL*a*b*c>=n) {
            Candidate q=gridBox(a,b,c,n);
            if(q.value>best.value) best=move(q);
        }
    int df=ceil(cbrt(2.0*n))+3;
    for(int a=1;a<=df;a++) for(int b=1;b<=df;b++) for(int c=1;c<=df;c++) {
        long long total=1LL*a*b*c;
        for(int par=0;par<2;par++) {
            long long cap=(total + ((par==0)?1:0))/2;
            // Formula above is exact only if the first cell has even parity.
            if(cap<n) continue;
            Candidate q=fccBox(a,b,c,par,n);
            if(q.value>best.value) best=move(q);
        }
    }
    int lo[3]={INT_MAX,INT_MAX,INT_MAX},hi[3]={INT_MIN,INT_MIN,INT_MIN};
    for(auto v:best.p) for(int k=0;k<3;k++) lo[k]=min(lo[k],v[k]),hi[k]=max(hi[k],v[k]);
    int span=max({hi[0]-lo[0],hi[1]-lo[1],hi[2]-lo[2]});
    double scale=1.0/(span+best.d), margin=scale*best.d/2.0;
    cout<<setprecision(17);
    for(auto v:best.p)
        cout << margin+scale*(v[0]-lo[0]) << ' '
             << margin+scale*(v[1]-lo[1]) << ' '
             << margin+scale*(v[2]-lo[2]) << '\n';
}
