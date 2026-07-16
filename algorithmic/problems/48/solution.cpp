#include <bits/stdc++.h>
using namespace std;

struct Cand { double r = -1; vector<array<double,3>> p; };

static long long fccCount(int a,int b,int c){
    long long ae=(a+1)/2, ao=a/2, be=(b+1)/2, bo=b/2, ce=(c+1)/2, co=c/2;
    return ae*be*ce + ae*bo*co + ao*be*co + ao*bo*ce;
}
static long long bccCount(int a,int b,int c){
    long long ae=(a+1)/2, ao=a/2, be=(b+1)/2, bo=b/2, ce=(c+1)/2, co=c/2;
    return ae*be*ce + ao*bo*co;
}
static double step(double r,int d) { return d==1 ? 0.0 : (1.0-2.0*r)/(d-1); }
static double minFCC(double r,int a,int b,int c) {
    double s[3]={step(r,a),step(r,b),step(r,c)}; int d[3]={a,b,c};
    double q=1e100;
    for(int i=0;i<3;i++) if(d[i]>=3) q=min(q,2*s[i]);
    for(int i=0;i<3;i++) for(int j=i+1;j<3;j++) if(d[i]>=2&&d[j]>=2) q=min(q,hypot(s[i],s[j]));
    return q;
}
static double minBCC(double r,int a,int b,int c) {
    double s[3]={step(r,a),step(r,b),step(r,c)}; int d[3]={a,b,c};
    double q=1e100;
    for(int i=0;i<3;i++) if(d[i]>=3) q=min(q,2*s[i]);
    if(a>=2&&b>=2&&c>=2) q=min(q,sqrt(s[0]*s[0]+s[1]*s[1]+s[2]*s[2]));
    return q;
}
static vector<array<double,3>> makePoints(int n,int a,int b,int c,double r,bool bcc){
    vector<array<double,3>> v; v.reserve(n);
    double sx=step(r,a),sy=step(r,b),sz=step(r,c);
    for(int i=0;i<a && (int)v.size()<n;i++) for(int j=0;j<b && (int)v.size()<n;j++) for(int k=0;k<c && (int)v.size()<n;k++) {
        bool take=bcc ? ((i&1)==(j&1) && (j&1)==(k&1)) : (((i+j+k)&1)==0);
        if(take) v.push_back({a==1?.5:r+i*sx,b==1?.5:r+j*sy,c==1?.5:r+k*sz});
    }
    return v;
}
int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int n; if(!(cin>>n)) return 0;
    Cand best;
    auto put=[&](double r, vector<array<double,3>> p){ if((int)p.size()==n && r>best.r){best.r=r;best.p=move(p);} };
    // These are geometric constructions, not instance-specific output tables.
    if(n==2){
        double r=sqrt(3.0)/(2.0*(1.0+sqrt(3.0))), t=1-r;
        put(r,{{{r,r,r},{t,t,t}}});
    }
    if(n==3){
        double r=sqrt(2.0)/(2.0+2.0*sqrt(2.0)), t=1-r;
        put(r,{{{r,r,r},{t,t,r},{t,r,t}}});
    }
    if(n==4){
        double r=sqrt(2.0)/(2.0+2.0*sqrt(2.0)), t=1-r;
        put(r,{{{r,r,r},{r,t,t},{t,r,t},{t,t,r}}});
    }
    // Cubic grid is an unconditional lower-bound safety candidate.
    int m=1; while(1LL*m*m*m<n) ++m;
    double gr=1.0/(2.0*m); vector<array<double,3>> gp; gp.reserve(n);
    for(int i=0;i<m&&(int)gp.size()<n;i++) for(int j=0;j<m&&(int)gp.size()<n;j++) for(int k=0;k<m&&(int)gp.size()<n;k++)
        gp.push_back({gr+i*2*gr,gr+j*2*gr,gr+k*2*gr});
    put(gr,move(gp));
    int lim=(int)ceil(cbrt((double)n))*2+6;
    for(int a=1;a<=lim;a++) for(int b=1;b<=lim;b++) for(int c=1;c<=lim;c++) for(int typ=0;typ<2;typ++) {
        if((typ?bccCount(a,b,c):fccCount(a,b,c))<n) continue;
        double lo=0,hi=.5;
        for(int it=0;it<52;it++) { double x=(lo+hi)/2; double md=typ?minBCC(x,a,b,c):minFCC(x,a,b,c); if(md>=2*x) lo=x; else hi=x; }
        if(lo>best.r+1e-14) put(lo,makePoints(n,a,b,c,lo,typ));
    }
    cout<<setprecision(17);
    for(auto q:best.p) cout<<q[0]<<' '<<q[1]<<' '<<q[2]<<'\n';
}
