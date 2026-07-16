#include <bits/stdc++.h>
using namespace std;

struct Pt { long long x, y; };
int N;
vector<Pt> pts;
vector<char> isPrimeId;

static inline double d2(int a, int b){
    double dx = (double)pts[a].x - (double)pts[b].x;
    double dy = (double)pts[a].y - (double)pts[b].y;
    return hypot(dx, dy);
}
static inline double multForStep(int t, int src){
    return (t % 10 == 0 && !isPrimeId[src]) ? 1.1 : 1.0;
}
static inline double edgeCost(const vector<int>& r, int t){ // edge r[t-1] -> r[t], 1 <= t <= N
    return multForStep(t, r[t-1]) * d2(r[t-1], r[t]);
}
double routeCostOrder(const vector<int>& ord){
    double s = 0.0;
    int prev = 0;
    for(int i=0;i<(int)ord.size();++i){
        int t = i + 1;
        s += multForStep(t, prev) * d2(prev, ord[i]);
        prev = ord[i];
    }
    s += multForStep(N, prev) * d2(prev, 0);
    return s;
}
double routeCost(const vector<int>& r){
    double s=0; for(int t=1;t<=N;t++) s += edgeCost(r,t); return s;
}

vector<char> sievePrime(int n){
    vector<char> p(max(2,n), true);
    p[0]=false; if(n>1) p[1]=false;
    for(long long i=2;i*i<n;i++) if(p[i]) for(long long j=i*i;j<n;j+=i) p[(int)j]=false;
    return p;
}

uint64_t hilbertOrder(uint32_t x, uint32_t y){
    const int B = 21;
    uint64_t d = 0;
    for(int s=B-1; s>=0; --s){
        uint32_t rx = (x >> s) & 1u;
        uint32_t ry = (y >> s) & 1u;
        d += (uint64_t(1) << (2*s)) * ((3u * rx) ^ ry);
        if(ry == 0){
            if(rx == 1){
                uint32_t mask = (1u << B) - 1u;
                x = mask - x; y = mask - y;
            }
            swap(x,y);
        }
    }
    return d;
}

vector<int> idsAll(){ vector<int> v; v.reserve(max(0,N-1)); for(int i=1;i<N;i++) v.push_back(i); return v; }

vector<int> candidateHilbert(bool flipX, bool flipY, bool rev){
    long long minx=pts[0].x, maxx=pts[0].x, miny=pts[0].y, maxy=pts[0].y;
    for(auto &p: pts){ minx=min(minx,p.x); maxx=max(maxx,p.x); miny=min(miny,p.y); maxy=max(maxy,p.y); }
    const uint32_t M = (1u<<21) - 1u;
    vector<pair<uint64_t,int>> a; a.reserve(N-1);
    long double rx = max(1LL, maxx-minx), ry = max(1LL, maxy-miny);
    for(int i=1;i<N;i++){
        uint32_t xx = (uint32_t) llround(((long double)(pts[i].x-minx) / rx) * M);
        uint32_t yy = (uint32_t) llround(((long double)(pts[i].y-miny) / ry) * M);
        if(flipX) xx = M - xx;
        if(flipY) yy = M - yy;
        a.push_back({hilbertOrder(xx,yy), i});
    }
    sort(a.begin(), a.end());
    vector<int> v; v.reserve(N-1);
    for(auto &e:a) v.push_back(e.second);
    if(rev) reverse(v.begin(), v.end());
    return v;
}

vector<int> candidateSerpentine(bool primaryY, int bucket){
    vector<int> v = idsAll();
    if(primaryY) sort(v.begin(), v.end(), [&](int a,int b){ if(pts[a].y!=pts[b].y) return pts[a].y<pts[b].y; return pts[a].x<pts[b].x; });
    // else x order is already by id.
    for(int l=0, block=0; l<(int)v.size(); l+=bucket, block++){
        int r=min((int)v.size(), l+bucket);
        if(primaryY){
            sort(v.begin()+l, v.begin()+r, [&](int a,int b){ return block%2 ? pts[a].x>pts[b].x : pts[a].x<pts[b].x; });
        }else{
            sort(v.begin()+l, v.begin()+r, [&](int a,int b){ return block%2 ? pts[a].y>pts[b].y : pts[a].y<pts[b].y; });
        }
    }
    return v;
}

vector<int> candidateAngle(bool rev){
    long double cx = pts[0].x, cy = pts[0].y;
    vector<pair<long double,int>> a; a.reserve(N-1);
    for(int i=1;i<N;i++) a.push_back({atan2((long double)pts[i].y-cy, (long double)pts[i].x-cx), i});
    sort(a.begin(), a.end());
    vector<int> v; v.reserve(N-1); for(auto &e:a) v.push_back(e.second);
    if(rev) reverse(v.begin(), v.end());
    return v;
}

void trySwap(vector<int>& r, int i, int j){
    if(i==j || i<=0 || j<=0 || i>=N || j>=N) return;
    if(i>j) swap(i,j);
    int ts[4] = {i, i+1, j, j+1};
    double old=0, neu=0; vector<int> used;
    for(int t: ts) if(t>=1 && t<=N && find(used.begin(), used.end(), t)==used.end()) { used.push_back(t); old += edgeCost(r,t); }
    swap(r[i], r[j]);
    for(int t: used) neu += edgeCost(r,t);
    if(neu + 1e-7 < old) return;
    swap(r[i], r[j]);
}

void localPolish(vector<int>& r){
    // Cheap adjacent pass for geometric-order boundary mistakes.
    for(int pass=0; pass<2; ++pass){
        bool any=false;
        for(int i=1;i+1<N;i++){
            double before = 0, after = 0;
            vector<int> ts = {i, i+1, i+2};
            for(int t: ts) if(t>=1 && t<=N) before += edgeCost(r,t);
            swap(r[i], r[i+1]);
            for(int t: ts) if(t>=1 && t<=N) after += edgeCost(r,t);
            if(after + 1e-7 < before) any=true; else swap(r[i], r[i+1]);
        }
        if(!any) break;
    }
    // Move nearby prime IDs onto penalized source positions only when the full local cost improves.
    for(int pass=0; pass<2; ++pass){
        for(int p=9; p<N; p+=10){
            if(isPrimeId[r[p]]) continue;
            int lo=max(1,p-7), hi=min(N-1,p+7), best=-1; double bestDelta=0.0;
            for(int q=lo;q<=hi;q++) if(isPrimeId[r[q]]){
                int i=min(p,q), j=max(p,q);
                int ts[4]={i,i+1,j,j+1}; vector<int> used; double old=0, neu=0;
                for(int t: ts) if(t>=1 && t<=N && find(used.begin(),used.end(),t)==used.end()){ used.push_back(t); old+=edgeCost(r,t); }
                swap(r[p], r[q]);
                for(int t: used) neu += edgeCost(r,t);
                swap(r[p], r[q]);
                double delta=neu-old;
                if(delta < bestDelta){ bestDelta=delta; best=q; }
            }
            if(best!=-1) swap(r[p], r[best]);
        }
    }
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    if(!(cin>>N)) return 0;
    pts.resize(N);
    for(int i=0;i<N;i++) cin>>pts[i].x>>pts[i].y;
    isPrimeId = sievePrime(N);

    if(N<=10){
        vector<int> perm=idsAll(), best=perm; double bc=1e300;
        do{ double c=routeCostOrder(perm); if(c<bc){bc=c; best=perm;} }while(next_permutation(perm.begin(), perm.end()));
        cout << N+1 << '\n' << 0 << '\n';
        for(int v: best) cout << v << '\n';
        cout << 0 << '\n'; return 0;
    }

    vector<vector<int>> cand;
    cand.reserve(16);
    cand.push_back(idsAll());
    { auto v=idsAll(); reverse(v.begin(), v.end()); cand.push_back(move(v)); }
    { auto v=idsAll(); sort(v.begin(), v.end(), [&](int a,int b){ if(pts[a].y!=pts[b].y) return pts[a].y<pts[b].y; return pts[a].x<pts[b].x; }); cand.push_back(v); reverse(v.begin(), v.end()); cand.push_back(move(v)); }
    int b1=max(8,(int)sqrt((double)N));
    int b2=max(16, (int)(2*sqrt((double)N)));
    cand.push_back(candidateSerpentine(false,b1));
    cand.push_back(candidateSerpentine(false,b2));
    cand.push_back(candidateSerpentine(true,b1));
    cand.push_back(candidateSerpentine(true,b2));
    cand.push_back(candidateAngle(false));
    cand.push_back(candidateAngle(true));
    cand.push_back(candidateHilbert(false,false,false));
    cand.push_back(candidateHilbert(false,false,true));
    cand.push_back(candidateHilbert(true,false,false));
    cand.push_back(candidateHilbert(false,true,false));

    double bestCost=1e300; vector<int> bestOrd;
    for(auto &v: cand){
        double c=routeCostOrder(v);
        if(c<bestCost){ bestCost=c; bestOrd=v; }
    }
    vector<int> route; route.reserve(N+1); route.push_back(0); for(int v: bestOrd) route.push_back(v); route.push_back(0);
    localPolish(route);

    cout << N+1 << '\n';
    for(int v: route) cout << v << '\n';
    return 0;
}
