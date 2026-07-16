#include <bits/stdc++.h>
using namespace std;
struct P { double x,y,z; };

static double radiusOf(const vector<P>& a) {
    double r=1.0;
    for (auto p:a) r=min(r,min({p.x,1-p.x,p.y,1-p.y,p.z,1-p.z}));
    for (int i=0;i<(int)a.size();++i) for(int j=0;j<i;++j) {
        double x=a[i].x-a[j].x, y=a[i].y-a[j].y, z=a[i].z-a[j].z;
        r=min(r,0.5*sqrt(x*x+y*y+z*z));
    }
    return r;
}
static unsigned long long rngstate;
static double rnd() { // reproducible, deliberately no dependence on clock
    rngstate=rngstate*6364136223846793005ULL+1442695040888963407ULL;
    return ((rngstate>>11)*(1.0/9007199254740992.0));
}
static void relax(vector<P>& a, double want, int rounds, double shake) {
    int n=a.size();
    double lo=want, hi=1.0-want, d=2*want;
    if (shake) for(P &p:a) {
        p.x=min(hi,max(lo,p.x+(rnd()-.5)*shake));
        p.y=min(hi,max(lo,p.y+(rnd()-.5)*shake));
        p.z=min(hi,max(lo,p.z+(rnd()-.5)*shake));
    }
    vector<P> f(n);
    for(int it=0;it<rounds;++it) {
        fill(f.begin(),f.end(),P{0,0,0});
        for(int i=0;i<n;++i) for(int j=0;j<i;++j) {
            double x=a[i].x-a[j].x, y=a[i].y-a[j].y, z=a[i].z-a[j].z;
            double q=sqrt(x*x+y*y+z*z);
            if(q<d) {
                if(q<1e-14) { x=rnd()-.5; y=rnd()-.5; z=rnd()-.5; q=sqrt(x*x+y*y+z*z); }
                double t=0.56*(d-q)/q;
                x*=t; y*=t; z*=t;
                f[i].x+=x; f[i].y+=y; f[i].z+=z;
                f[j].x-=x; f[j].y-=y; f[j].z-=z;
            }
        }
        for(int i=0;i<n;++i) {
            a[i].x=min(hi,max(lo,a[i].x+f[i].x));
            a[i].y=min(hi,max(lo,a[i].y+f[i].y));
            a[i].z=min(hi,max(lo,a[i].z+f[i].z));
        }
    }
}
int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int n; if(!(cin>>n)) return 0;
    int m=1; while(1LL*m*m*m<n) ++m;
    vector<P> best;
    for(int i=0;(int)best.size()<n && i<m;++i) for(int j=0;(int)best.size()<n && j<m;++j)
        for(int k=0;(int)best.size()<n && k<m;++k) best.push_back({(i+.5)/m,(j+.5)/m,(k+.5)/m});

    // Axis-aligned FCC remains the inexpensive, high-density large-n fallback.
    const double rt2=sqrt(2.0);
    auto info=[&](double s) {
        double r=s/rt2, span=1-2*r;
        if(span<0) return pair<long long,int>{0,0};
        int q=(int)floor(span/s)+1; if(q<=0) return pair<long long,int>{0,0};
        long long e=(q+1)/2,o=q/2;
        long long a=e*e*e+3*e*o*o, b=o*o*o+3*o*e*e;
        return pair<long long,int>{max(a,b),a>=b?0:1};
    };
    double l=0,h=1;
    for(int z=0;z<70;++z) { double s=(l+h)/2; if(info(s).first>=n) l=s; else h=s; }
    double fr=l/rt2, gr=radiusOf(best);
    if(fr>gr*(1+1e-13)) {
        vector<P> f; auto qinfo=info(l); int par=qinfo.second;
        int q=(int)floor((1-2*fr)/l)+1;
        for(int i=0;(int)f.size()<n && i<q;++i) for(int j=0;(int)f.size()<n && j<q;++j)
            for(int k=0;(int)f.size()<n && k<q;++k) if(((i+j+k)&1)==par)
                f.push_back({fr+i*l,fr+j*l,fr+k*l});
        if((int)f.size()==n && radiusOf(f)>gr) best=f;
    }

    // Finite cubes have large boundary effects.  This is a separate search
    // family: grow a requested radius and repeatedly project pair overlaps
    // outwards, retaining only arrangements whose measured radius improves.
    if(n<=32) {
        double start=radiusOf(best);
        const double delta=acos(-1.0)/sqrt(18.0);
        double upper=min(.499999, pow(delta*3.0/(4.0*acos(-1.0)*n),1.0/3.0));
        for(int trial=0;trial<6;++trial) {
            vector<P> cur=best;
            rngstate=0x9e3779b97f4a7c15ULL+(unsigned long long)n*1315423911ULL+trial;
            for(int step=1;step<=55;++step) {
                double want=start+(upper-start)*step/55.0;
                relax(cur,want,100, step==1 ? .045*(trial+1)/6.0 : 0.0);
                double got=radiusOf(cur);
                if(got>radiusOf(best)*(1+1e-12)) best=cur;
            }
        }
    }
    cout<<setprecision(17);
    for(P p:best) cout<<min(1.0,max(0.0,p.x))<<' '<<min(1.0,max(0.0,p.y))<<' '<<min(1.0,max(0.0,p.z))<<'\n';
}
