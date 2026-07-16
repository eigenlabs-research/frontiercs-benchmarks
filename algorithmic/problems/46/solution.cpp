#include <bits/stdc++.h>
using namespace std;

struct Op { int m; long long p; };

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int J, M;
    if(!(cin >> J >> M)) return 0;
    vector<vector<Op>> ops(J, vector<Op>(M));
    for(int j=0;j<J;j++){
        for(int k=0;k<M;k++) cin >> ops[j][k].m >> ops[j][k].p;
    }

    vector<vector<long long>> rem(J, vector<long long>(M+1,0));
    vector<long long> totalJob(J,0), loadMach(M,0);
    for(int j=0;j<J;j++){
        for(int k=M-1;k>=0;k--) rem[j][k]=rem[j][k+1]+ops[j][k].p;
        totalJob[j]=rem[j][0];
        for(int k=0;k<M;k++) loadMach[ops[j][k].m]+=ops[j][k].p;
    }
    long long maxLoad = 1;
    for(long long x: loadMach) maxLoad=max(maxLoad,x);

    struct Result { long long cmax; vector<vector<int>> order; };

    auto run = [&](int rule, array<double,6> w, int seed)->Result{
        vector<int> nxt(J,0);
        vector<long long> jr(J,0), mr(M,0);
        vector<vector<int>> ord(M); for(int m=0;m<M;m++) ord[m].reserve(J);
        uint64_t rng = 1469598103934665603ULL ^ (uint64_t)(seed+1009*rule);
        auto rnd = [&](){ rng ^= rng<<7; rng ^= rng>>9; return (double)(rng & ((1ULL<<24)-1)) / (double)(1ULL<<24); };

        auto better = [&](int a, int b)->bool{
            if(b<0) return true;
            int ka=nxt[a], kb=nxt[b];
            const Op &A=ops[a][ka], &B=ops[b][kb];
            long long tailA=rem[a][ka+1], tailB=rem[b][kb+1];
            double scoreA=0, scoreB=0;
            switch(rule){
                case 0: scoreA=-A.p; scoreB=-B.p; break;                         // shortest processing time
                case 1: scoreA=A.p; scoreB=B.p; break;                           // longest processing time
                case 2: scoreA=rem[a][ka]; scoreB=rem[b][kb]; break;             // most work remaining
                case 3: scoreA=(M-ka)*1000000.0 + rem[a][ka]; scoreB=(M-kb)*1000000.0 + rem[b][kb]; break;
                case 4: scoreA=tailA; scoreB=tailB; break;                       // delivery time
                case 5: scoreA=-jr[a]; scoreB=-jr[b]; break;                     // earliest release
                case 6: scoreA=loadMach[A.m]*1.0/maxLoad*rem[a][ka] - 0.2*A.p; scoreB=loadMach[B.m]*1.0/maxLoad*rem[b][kb] - 0.2*B.p; break;
                default:
                    scoreA = w[0]*A.p + w[1]*rem[a][ka] + w[2]*tailA + w[3]*jr[a] + w[4]*ka + w[5]*(loadMach[A.m]*1.0/maxLoad) + 1e-6*rnd();
                    scoreB = w[0]*B.p + w[1]*rem[b][kb] + w[2]*tailB + w[3]*jr[b] + w[4]*kb + w[5]*(loadMach[B.m]*1.0/maxLoad) + 1e-6*rnd();
            }
            if(fabs(scoreA-scoreB)>1e-9) return scoreA>scoreB;
            if(jr[a]!=jr[b]) return jr[a]<jr[b];
            if(A.p!=B.p) return A.p<B.p;
            return a<b;
        };

        int N=J*M;
        for(int done=0; done<N; done++){
            long long bestC = LLONG_MAX; int crit=-1;
            for(int j=0;j<J;j++) if(nxt[j]<M){
                const Op &o=ops[j][nxt[j]];
                long long c=max(jr[j], mr[o.m]) + o.p;
                if(c<bestC || (c==bestC && better(j,crit))){ bestC=c; crit=j; }
            }
            int mach = ops[crit][nxt[crit]].m;
            int choose=-1;
            for(int j=0;j<J;j++) if(nxt[j]<M && ops[j][nxt[j]].m==mach && jr[j] < bestC){
                if(better(j, choose)) choose=j;
            }
            if(choose<0) choose=crit;
            const Op &o=ops[choose][nxt[choose]];
            long long st=max(jr[choose], mr[o.m]);
            long long en=st+o.p;
            jr[choose]=en; mr[o.m]=en; ord[o.m].push_back(choose); nxt[choose]++;
        }
        long long cmax=0; for(long long x: jr) cmax=max(cmax,x);
        return {cmax, ord};
    };

    vector<array<double,6>> weights;
    // A compact deterministic portfolio of dispatch priorities (one constructive mechanism).
    double vals[] = {-1.0,-0.4,0.0,0.4,1.0};
    for(double a: vals) for(double b: vals) {
        weights.push_back({a,b,0.7*b,-0.25,0.05,0.2});
        if((int)weights.size()>=80) break;
    }
    for(int i=0;i<40;i++){
        uint64_t x = 88172645463393265ULL + i*11995408973635179863ULL;
        array<double,6> w;
        for(int t=0;t<6;t++){ x ^= x<<7; x ^= x>>9; w[t] = ((int)(x%2001)-1000)/500.0; }
        weights.push_back(w);
    }

    Result best{LLONG_MAX, {}};
    array<double,6> z{0,0,0,0,0,0};
    for(int r=0;r<=6;r++){
        Result cur=run(r,z,r);
        if(cur.cmax<best.cmax) best=move(cur);
    }
    int seed=100;
    for(auto &w: weights){
        Result cur=run(7,w,seed++);
        if(cur.cmax<best.cmax) best=move(cur);
    }

    // Safety fallback should never be used, but preserves the required format.
    if(best.order.empty()){
        best.order.assign(M, {});
        for(int m=0;m<M;m++) for(int j=0;j<J;j++) best.order[m].push_back(j);
    }
    for(int m=0;m<M;m++){
        for(int i=0;i<J;i++){
            if(i) cout << ' ';
            cout << best.order[m][i];
        }
        cout << '\n';
    }
    return 0;
}
