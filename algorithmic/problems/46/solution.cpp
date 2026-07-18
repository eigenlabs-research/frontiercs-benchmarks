#include <bits/stdc++.h>
using namespace std;

struct Instance {
    int J, M, N;
    vector<vector<int>> mach, pt, posOfMach, opId;
    vector<int> opJob, opPos, opMach, opPt;
    vector<long long> totalJob, totalMach;
};

static Instance in;

struct Sol {
    vector<vector<int>> ord;
    long long val = (1LL<<62);
};

long long evalOrder(const vector<vector<int>>& ord) {
    int J=in.J, M=in.M, N=in.N;
    vector<vector<int>> adj(N);
    vector<int> indeg(N,0);
    for (int j=0;j<J;j++) for (int k=0;k+1<M;k++) {
        int a=in.opId[j][k], b=in.opId[j][k+1]; adj[a].push_back(b); indeg[b]++;
    }
    for (int m=0;m<M;m++) {
        if ((int)ord[m].size()!=J) return (1LL<<61);
        vector<int> seen(J,0);
        int prev=-1;
        for (int x: ord[m]) {
            if (x<0 || x>=J || seen[x]++) return (1LL<<61);
            int k=in.posOfMach[x][m];
            if (k<0) return (1LL<<61);
            int id=in.opId[x][k];
            if (prev!=-1) { adj[prev].push_back(id); indeg[id]++; }
            prev=id;
        }
    }
    queue<int> q;
    vector<long long> dist(N,0);
    for (int i=0;i<N;i++) if (!indeg[i]) q.push(i);
    int cnt=0; long long best=0;
    while(!q.empty()) {
        int u=q.front(); q.pop(); cnt++;
        long long fin=dist[u]+in.opPt[u];
        if (fin>best) best=fin;
        for (int v: adj[u]) {
            if (dist[v]<fin) dist[v]=fin;
            if (--indeg[v]==0) q.push(v);
        }
    }
    if (cnt!=N) return (1LL<<61);
    return best;
}

Sol serialSchedule(int rule, uint64_t seed=1) {
    int J=in.J, M=in.M;
    vector<int> nxt(J,0);
    vector<long long> jr(J,0), mr(M,0);
    vector<vector<int>> ord(M); for (int m=0;m<M;m++) ord[m].reserve(J);
    int bottleneck = max_element(in.totalMach.begin(), in.totalMach.end()) - in.totalMach.begin();
    long long maxMach = max(1LL, in.totalMach[bottleneck]);
    mt19937_64 rng(seed);
    for (int step=0; step<J*M; step++) {
        int bj=-1; tuple<long long,long long,long long,long long> bk(LLONG_MAX,LLONG_MAX,LLONG_MAX,LLONG_MAX);
        for (int j=0;j<J;j++) if (nxt[j]<M) {
            int k=nxt[j], m=in.mach[j][k], p=in.pt[j][k];
            long long est=max(jr[j], mr[m]);
            long long rem=0, tail=0;
            for (int t=k;t<M;t++) rem += in.pt[j][t];
            for (int t=k+1;t<M;t++) tail += in.pt[j][t];
            long long before=0; for (int t=0;t<k;t++) before += in.pt[j][t];
            long long a=est, b=0, c=0, d=j;
            switch(rule) {
                case 0: b=-(long long)p; c=-rem; break;                 // long current work
                case 1: b=-rem; c=est+p; break;                         // most remaining work
                case 2: b=p; c=rem; break;                              // short current work
                case 3: b=-tail; c=est+p; break;                        // large delivery tail
                case 4: a=est+p; b=-rem; break;                         // earliest completion
                case 5: b=(m==bottleneck?0:1); c=(m==bottleneck?before:est); d=-tail; break;
                case 6: b=(m==bottleneck?0:1); c=-(long long)p; d=-rem; break;
                case 7: b=mr[m]; c=-in.totalMach[m]; d=-rem; break;     // feed busy/big machines
                case 8: b=-(long long)p*1000000/maxMach - rem/10; c=tail; break;
                case 9: b=(long long)(rng()%1000000) - rem/1000 - (long long)p/1000; break;
                case 10: b=(long long)(rng()%1000000) + p/1000 - tail/1000; break;
                default: b=(long long)(rng()%1000000) - (m==bottleneck?200000:0) - rem/2000; break;
            }
            auto key=make_tuple(a,b,c,d);
            if (key<bk) { bk=key; bj=j; }
        }
        int j=bj, k=nxt[j], m=in.mach[j][k], p=in.pt[j][k];
        long long st=max(jr[j], mr[m]);
        jr[j]=mr[m]=st+p;
        ord[m].push_back(j);
        nxt[j]++;
    }
    Sol s; s.ord=move(ord); s.val=evalOrder(s.ord); return s;
}

Sol sameOrderCandidate(vector<int> perm) {
    Sol s; s.ord.assign(in.M, perm); s.val=evalOrder(s.ord); return s;
}

void consider(Sol &best, const Sol &s) {
    if (s.val < best.val) best = s;
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int J,M; if(!(cin>>J>>M)) return 0;
    in.J=J; in.M=M; in.N=J*M;
    in.mach.assign(J, vector<int>(M)); in.pt.assign(J, vector<int>(M));
    in.posOfMach.assign(J, vector<int>(M,-1)); in.opId.assign(J, vector<int>(M));
    in.totalJob.assign(J,0); in.totalMach.assign(M,0);
    in.opJob.reserve(in.N); in.opPos.reserve(in.N); in.opMach.reserve(in.N); in.opPt.reserve(in.N);
    int id=0;
    for(int j=0;j<J;j++) for(int k=0;k<M;k++) {
        int m,p; cin>>m>>p; in.mach[j][k]=m; in.pt[j][k]=p; in.posOfMach[j][m]=k; in.opId[j][k]=id++;
        in.opJob.push_back(j); in.opPos.push_back(k); in.opMach.push_back(m); in.opPt.push_back(p);
        in.totalJob[j]+=p; in.totalMach[m]+=p;
    }
    Sol best;
    for(int r=0;r<=8;r++) consider(best, serialSchedule(r, 12345+r));
    int iters = (J*M<=1300 ? 60 : 30);
    for(int t=0;t<iters;t++) consider(best, serialSchedule(9 + (t%3), 88172645463325252ull + t*1000003ull));
    vector<int> perm(J); iota(perm.begin(), perm.end(), 0);
    sort(perm.begin(), perm.end(), [&](int a,int b){ return in.totalJob[a]>in.totalJob[b]; });
    consider(best, sameOrderCandidate(perm));
    int bottleneck = max_element(in.totalMach.begin(), in.totalMach.end()) - in.totalMach.begin();
    sort(perm.begin(), perm.end(), [&](int a,int b){
        int ka=in.posOfMach[a][bottleneck], kb=in.posOfMach[b][bottleneck];
        long long ba=0, bb=0, ta=0, tb=0;
        for(int t=0;t<ka;t++) ba+=in.pt[a][t]; for(int t=ka+1;t<M;t++) ta+=in.pt[a][t];
        for(int t=0;t<kb;t++) bb+=in.pt[b][t]; for(int t=kb+1;t<M;t++) tb+=in.pt[b][t];
        if (ba!=bb) return ba<bb; return ta>tb;
    });
    consider(best, sameOrderCandidate(perm));

    // A bounded adjacent-swap cleanup: accept any valid improving neighboring machine order.
    bool improved=true; int passes=0;
    while(improved && passes<2) {
        improved=false; passes++;
        for(int m=0;m<M;m++) for(int i=0;i+1<J;i++) {
            auto cand=best.ord; swap(cand[m][i], cand[m][i+1]);
            long long v=evalOrder(cand);
            if (v < best.val) { best.ord.swap(cand); best.val=v; improved=true; }
        }
    }

    if (best.ord.empty()) { best.ord.assign(M, vector<int>(J)); for(int m=0;m<M;m++) iota(best.ord[m].begin(), best.ord[m].end(), 0); }
    for(int m=0;m<M;m++) {
        for(int i=0;i<J;i++) { if(i) cout << ' '; cout << best.ord[m][i]; }
        cout << '\n';
    }
    return 0;
}
