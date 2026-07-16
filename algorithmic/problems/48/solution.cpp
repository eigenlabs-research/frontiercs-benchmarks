#include <bits/stdc++.h>
using namespace std;

struct Mat { double a[3][3]; };
struct Pick { double r=0, h=0, mn[3]; int A=0,B=0,C=0, rot=0; bool fcc=false; };

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n;
    if (!(cin >> n)) return 0;
    const double rt2=sqrt(2.0), rt3=sqrt(3.0), rt6=sqrt(6.0);
    vector<Mat> R;
    R.push_back({{{1,0,0},{0,1,0},{0,0,1}}});
    double q=1.0/rt2;
    R.push_back({{{q,-q,0},{q,q,0},{0,0,1}}});
    R.push_back({{{q,0,q},{0,1,0},{-q,0,q}}});
    R.push_back({{{1,0,0},{0,q,-q},{0,q,q}}});
    // One lattice direction is diagonal to the cube; useful for thin boxes.
    R.push_back({{{1/rt3,1/rt2,1/rt6},{1/rt3,-1/rt2,1/rt6},{1/rt3,0,-2/rt6}}});

    Pick best;
    // The ordinary grid is both a guaranteed baseline and preferable for small n.
    int limA=1;
    while ((long long)limA*limA*limA < n) ++limA;
    for (int A=1; A<=limA; ++A) for (int B=A; (long long)A*B*B<=n || B==A; ++B) {
        int C=(n + A*B-1)/(A*B);
        if (C<B) C=B;
        double r=1.0/(2.0*C);
        if (r>best.r) { best={r,0,{0,0,0},A,B,C,0,false}; }
        if (B>=C) break;
    }

    // A conventional FCC cell has four sites.  Search rectangular crops and
    // several orientations, then center the crop in the cube with its radius margin.
    int need=(n+3)/4;
    int la=1; while ((long long)la*la*la<need) ++la;
    for (int A=1; A<=la; ++A) {
        int lb=1; while ((long long)A*lb*lb<need) ++lb;
        for (int B=A; B<=lb; ++B) {
            int C=(need + A*B-1)/(A*B);
            if (C<B) C=B;
            // The output is a prefix of this crop.  Its bounding box can be
            // substantially smaller than that of all 4*A*B*C FCC sites.
            for (int q=0;q<(int)R.size();++q) {
                double lo[3]={1e100,1e100,1e100}, hi[3]={-1e100,-1e100,-1e100};
                int used=0;
                for (int i=0;i<A && used<n;++i) for (int j=0;j<B && used<n;++j)
                  for (int k=0;k<C && used<n;++k) for (int s=0;s<4 && used<n;++s,++used) {
                    double v[3]={i+(s>=2 ? .5 : 0), j+(s==1 || s==3 ? .5 : 0),
                                 k+(s==1 || s==2 ? .5 : 0)};
                    for (int u=0;u<3;++u) {
                        double x=0;
                        for (int z=0;z<3;++z) x+=R[q].a[u][z]*v[z];
                        lo[u]=min(lo[u],x); hi[u]=max(hi[u],x);
                    }
                }
                double W=max({hi[0]-lo[0],hi[1]-lo[1],hi[2]-lo[2]});
                double h=(1.0-1e-13)/(W+1.0/rt2);
                double r=h/(2.0*rt2);
                if (r>best.r) best={r,h,{lo[0],lo[1],lo[2]},A,B,C,q,true};
            }
        }
    }

    cout<<setprecision(17);
    if (!best.fcc) {
        int put=0;
        for(int i=0;i<best.A && put<n;++i)
          for(int j=0;j<best.B && put<n;++j)
            for(int k=0;k<best.C && put<n;++k,++put)
              cout<<(i+.5)/best.A<<' '<<(j+.5)/best.B<<' '<<(k+.5)/best.C<<'\n';
        return 0;
    }
    static const double off[4][3]={{0,0,0},{0,.5,.5},{.5,0,.5},{.5,.5,0}};
    int put=0;
    const Mat &M=R[best.rot];
    for(int i=0;i<best.A && put<n;++i) for(int j=0;j<best.B && put<n;++j)
      for(int k=0;k<best.C && put<n;++k) for(int s=0;s<4 && put<n;++s) {
        double v[3]={i+off[s][0],j+off[s][1],k+off[s][2]}, p[3];
        for(int u=0;u<3;++u) {
            p[u]=best.r-best.h*best.mn[u];
            for(int z=0;z<3;++z) p[u]+=best.h*M.a[u][z]*v[z];
        }
        cout<<p[0]<<' '<<p[1]<<' '<<p[2]<<'\n'; ++put;
      }
}
