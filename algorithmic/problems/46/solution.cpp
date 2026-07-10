// General Job-Shop Scheduling heuristic (FrontierCS open-optimization track).
// Giffler-Thompson active-schedule construction (multiple priority rules +
// randomized restarts) followed by iterated local search over the critical-path
// adjacent-swap neighbourhood, under a strict wall-clock budget. Outputs, per
// machine, the job processing order (a permutation of 0..J-1). Fully general.
#include <bits/stdc++.h>
using namespace std;

static chrono::steady_clock::time_point T0;
static double TL_MS = 880.0;
static inline double ms(){ return chrono::duration<double,milli>(chrono::steady_clock::now()-T0).count(); }

int J,M,NO;
vector<vector<int>> jm;        // jm[j][k] machine of op k
vector<vector<long long>> jp;  // jp[j][k] proc time
vector<vector<int>> kOf;       // kOf[j][m]
inline int OP(int j,int k){ return j*M+k; }

// --- fast allocation-free makespan evaluator ---
static vector<long long> dur, st;
static vector<int> jobPred, jobNext;   // fixed job-chain arcs (op ids or -1)
static vector<int> macPred, macNext;   // machine arcs (rebuilt per eval)
static vector<int> indeg, q, critPred; // topo helpers; critPred = pred that set start

void initEval(){
    dur.assign(NO,0); st.assign(NO,0); jobPred.assign(NO,-1); jobNext.assign(NO,-1);
    macPred.assign(NO,-1); macNext.assign(NO,-1); indeg.assign(NO,0); q.reserve(NO); critPred.assign(NO,-1);
    for(int j=0;j<J;j++) for(int k=0;k<M;k++){ int o=OP(j,k); dur[o]=jp[j][k];
        if(k>0) jobPred[o]=OP(j,k-1); if(k<M-1) jobNext[o]=OP(j,k+1); }
}
// build machine arcs from seq, compute longest path; return makespan or LLONG_MAX.
long long evalMakespan(const vector<vector<int>>& seq, bool wantCrit=false){
    for(int o=0;o<NO;o++){ macPred[o]=-1; macNext[o]=-1; indeg[o]=(jobPred[o]>=0?1:0); }
    for(int m=0;m<M;m++){ const auto&s=seq[m]; for(int i=1;i<(int)s.size();i++){
        int oa=OP(s[i-1],kOf[s[i-1]][m]), ob=OP(s[i],kOf[s[i]][m]);
        macNext[oa]=ob; macPred[ob]=oa; indeg[ob]++;
    }}
    q.clear();
    for(int o=0;o<NO;o++){ st[o]=0; if(wantCrit) critPred[o]=-1; if(indeg[o]==0) q.push_back(o); }
    int processed=0; long long mk=0; int endOp=-1;
    for(int qi=0; qi<(int)q.size(); qi++){
        int o=q[qi]; processed++; long long fin=st[o]+dur[o];
        if(fin>mk){ mk=fin; endOp=o; }
        int nx=jobNext[o]; if(nx>=0){ if(fin>st[nx]){ st[nx]=fin; if(wantCrit)critPred[nx]=o; } if(--indeg[nx]==0)q.push_back(nx); }
        nx=macNext[o];    if(nx>=0){ if(fin>st[nx]){ st[nx]=fin; if(wantCrit)critPred[nx]=o; } if(--indeg[nx]==0)q.push_back(nx); }
    }
    if(processed!=NO) return LLONG_MAX;
    return mk;
}

vector<vector<int>> giffler(int rule, mt19937_64& rng){
    vector<vector<int>> seq(M);
    vector<int> nextK(J,0);
    vector<long long> jobReady(J,0),macReady(M,0),rem(J,0);
    for(int j=0;j<J;j++) for(int k=0;k<M;k++) rem[j]+=jp[j][k];
    int total=J*M,done=0;
    while(done<total){
        long long bestC=LLONG_MAX; int bestM=-1;
        for(int j=0;j<J;j++){ int k=nextK[j]; if(k>=M)continue; int m=jm[j][k];
            long long est=max(jobReady[j],macReady[m]),c=est+jp[j][k]; if(c<bestC){bestC=c;bestM=m;} }
        int cj=-1; double bk=-1e18;
        for(int j=0;j<J;j++){ int k=nextK[j]; if(k>=M||jm[j][k]!=bestM)continue;
            long long est=max(jobReady[j],macReady[bestM]); if(est>=bestC)continue;
            double key; switch(rule){ case 0:key=-(double)jp[j][k];break; case 1:key=(double)jp[j][k];break;
                case 2:key=(double)rem[j];break; case 3:key=-(double)rem[j];break; default:key=(double)(rng()%100000); }
            key+=(double)(rng()%17)*1e-6; if(key>bk){bk=key;cj=j;} }
        int j=cj,k=nextK[j],m=bestM; long long est=max(jobReady[j],macReady[m]),fin=est+jp[j][k];
        seq[m].push_back(j); macReady[m]=fin; jobReady[j]=fin; rem[j]-=jp[j][k]; nextK[j]++; done++;
    }
    return seq;
}

// position of each job on each machine within seq (for swaps)
static vector<vector<int>> pos;
void buildPos(const vector<vector<int>>& seq){ pos.assign(M, vector<int>(J,-1)); for(int m=0;m<M;m++) for(int i=0;i<(int)seq[m].size();i++) pos[m][seq[m][i]]=i; }

// Iterated local search: hill-climb the critical-path machine-arc swaps, then
// perturb from the incumbent to escape local optima. Simple and effective.
long long localSearch(vector<vector<int>>& seq, long long mk, mt19937_64& rng, double budget){
    long long best=mk; vector<vector<int>> bestSeq=seq;
    while(ms()<budget){
        bool improved=false;
        long long cur=evalMakespan(seq,true);
        int endOp=-1; long long fmax=-1;
        for(int o=0;o<NO;o++){ long long f=st[o]+dur[o]; if(f>fmax){fmax=f;endOp=o;} }
        vector<pair<int,int>> arcs; int o=endOp;
        while(o>=0){ int p=critPred[o]; if(p>=0 && macNext[p]==o) arcs.push_back({p,o}); o=p; }
        for(auto&ar: arcs){
            if(ms()>=budget) break;
            int u=ar.first, v=ar.second; int m=jm[u/M][u%M];
            buildPos(seq); int iu=pos[m][u/M], iv=pos[m][v/M];
            if(iu<0||iv<0||abs(iu-iv)!=1) continue;
            swap(seq[m][iu], seq[m][iv]);
            long long nm=evalMakespan(seq);
            if(nm<cur){ cur=nm; improved=true; if(nm<best){best=nm;bestSeq=seq;} }
            else swap(seq[m][iu], seq[m][iv]);
        }
        if(!improved){
            seq=bestSeq;
            for(int t=0;t<3;t++){ int m=rng()%M; if(seq[m].size()<2)continue; int i=rng()%(seq[m].size()-1); swap(seq[m][i],seq[m][i+1]); }
            long long nm=evalMakespan(seq); if(nm<best){best=nm;bestSeq=seq;}
        }
    }
    seq=bestSeq; return best;
}

int main(){
    T0=chrono::steady_clock::now();
    ios::sync_with_stdio(false); cin.tie(nullptr);
    if(!(cin>>J>>M)) return 0;
    NO=J*M;
    jm.assign(J,vector<int>(M)); jp.assign(J,vector<long long>(M));
    for(int j=0;j<J;j++) for(int k=0;k<M;k++) cin>>jm[j][k]>>jp[j][k];
    kOf.assign(J,vector<int>(M,-1)); for(int j=0;j<J;j++) for(int k=0;k<M;k++) kOf[j][jm[j][k]]=k;
    initEval();

    mt19937_64 rng(0x9E3779B97F4A7C15ULL ^ ((uint64_t)J<<32) ^ (uint64_t)M);
    vector<vector<int>> best; long long bestMk=LLONG_MAX;
    for(int rule=0;rule<4;rule++){ auto s=giffler(rule,rng); long long mk=evalMakespan(s); if(mk<bestMk){bestMk=mk;best=s;} }
    while(ms()<TL_MS*0.35){ int rule=rng()%5; auto s=giffler(rule,rng); long long mk=evalMakespan(s); if(mk<bestMk){bestMk=mk;best=s;} }
    if(!best.empty()) bestMk=localSearch(best,bestMk,rng,TL_MS*0.97);

    for(int m=0;m<M;m++) for(int i=0;i<J;i++){ cout<<best[m][i]; cout<<(i+1<J?' ':'\n'); }
    return 0;
}
