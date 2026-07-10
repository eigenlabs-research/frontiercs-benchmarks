// General Job-Shop Scheduling heuristic (FrontierCS open-optimization track).
// Giffler-Thompson construction (SPT/LPT/MWKR/LWKR + randomized restarts) then
// Nowicki-Smutnicki tabu search over the N5 critical-block boundary-swap
// neighbourhood. Allocation-free longest-path makespan evaluator. Time-adaptive:
// scales the GT-restart budget down for large instances so the tabu (which
// actually drives makespan down) gets the bulk of the 1s budget. Fully general.
#include <bits/stdc++.h>
using namespace std;

static chrono::steady_clock::time_point T0;
static double TL_MS = 930.0;
static inline double ms(){ return chrono::duration<double,milli>(chrono::steady_clock::now()-T0).count(); }

int J,M,NO;
vector<vector<int>> jm;
vector<vector<long long>> jp;
vector<vector<int>> kOf;
inline int OP(int j,int k){ return j*M+k; }

static vector<long long> dur, st;
static vector<int> jobPred, jobNext, macPred, macNext, indeg, q, critPred;
void initEval(){
    dur.assign(NO,0); st.assign(NO,0); jobPred.assign(NO,-1); jobNext.assign(NO,-1);
    macPred.assign(NO,-1); macNext.assign(NO,-1); indeg.assign(NO,0); q.reserve(NO); critPred.assign(NO,-1);
    for(int j=0;j<J;j++) for(int k=0;k<M;k++){ int o=OP(j,k); dur[o]=jp[j][k];
        if(k>0) jobPred[o]=OP(j,k-1); if(k<M-1) jobNext[o]=OP(j,k+1); }
}
long long evalMakespan(const vector<vector<int>>& seq, bool wantCrit=false){
    for(int o=0;o<NO;o++){ macPred[o]=-1; macNext[o]=-1; indeg[o]=(jobPred[o]>=0?1:0); }
    for(int m=0;m<M;m++){ const auto&s=seq[m]; for(int i=1;i<(int)s.size();i++){
        int oa=OP(s[i-1],kOf[s[i-1]][m]), ob=OP(s[i],kOf[s[i]][m]);
        macNext[oa]=ob; macPred[ob]=oa; indeg[ob]++; }}
    q.clear();
    for(int o=0;o<NO;o++){ st[o]=0; if(wantCrit) critPred[o]=-1; if(indeg[o]==0) q.push_back(o); }
    int processed=0; long long mk=0;
    for(int qi=0; qi<(int)q.size(); qi++){
        int o=q[qi]; processed++; long long fin=st[o]+dur[o]; if(fin>mk) mk=fin;
        int nx=jobNext[o]; if(nx>=0){ if(fin>st[nx]){ st[nx]=fin; if(wantCrit)critPred[nx]=o; } if(--indeg[nx]==0)q.push_back(nx); }
        nx=macNext[o];    if(nx>=0){ if(fin>st[nx]){ st[nx]=fin; if(wantCrit)critPred[nx]=o; } if(--indeg[nx]==0)q.push_back(nx); }
    }
    if(processed!=NO) return LLONG_MAX;
    return mk;
}
vector<vector<int>> giffler(int rule, mt19937_64& rng){
    vector<vector<int>> seq(M); vector<int> nextK(J,0);
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
static vector<vector<int>> pos;
void buildPos(const vector<vector<int>>& seq){ pos.assign(M, vector<int>(J,-1)); for(int m=0;m<M;m++) for(int i=0;i<(int)seq[m].size();i++) pos[m][seq[m][i]]=i; }

// N5 (Nowicki-Smutnicki) block-boundary candidate swaps on the critical path.
static vector<pair<int,int>> n5moves(const vector<vector<int>>& seq){
    vector<pair<int,int>> moves;
    int endOp=-1; long long fmax=-1;
    for(int o=0;o<NO;o++){ long long f=st[o]+dur[o]; if(f>fmax){fmax=f;endOp=o;} }
    vector<int> path; for(int o=endOp;o>=0;o=critPred[o]) path.push_back(o);
    reverse(path.begin(), path.end());
    vector<vector<int>> blocks; int i=0;
    while(i<(int)path.size()){ int s=i; while(i+1<(int)path.size() && macNext[path[i]]==path[i+1]) i++;
        blocks.emplace_back(path.begin()+s, path.begin()+i+1); i++; }
    buildPos(seq); int nb=blocks.size();
    for(int bi=0; bi<nb; bi++){ auto&B=blocks[bi]; int L=B.size(); if(L<2) continue; int mch=jm[B[0]/M][B[0]%M];
        if(bi>0){ int pu=pos[mch][B[0]/M]; if(pu>=0 && pu+1<(int)seq[mch].size() && seq[mch][pu+1]==B[1]/M) moves.push_back({mch,pu}); }
        if(bi<nb-1){ int pu=pos[mch][B[L-2]/M]; if(pu>=0 && pu+1<(int)seq[mch].size() && seq[mch][pu+1]==B[L-1]/M) moves.push_back({mch,pu}); }
    }
    return moves;
}
static unordered_map<long long,int> tabuU;
long long localSearch(vector<vector<int>>& seq, long long mk, mt19937_64& rng, double budget){
    long long best=mk; vector<vector<int>> bestSeq=seq;
    int iter=0, stall=0, tenure=8+(J+M)/4;
    tabuU.clear();
    auto keyOf=[&](int m,int a,int b){ int lo=min(a,b),hi=max(a,b); return (long long)m*(NO+1)*(NO+1)+(long long)lo*(NO+1)+hi; };
    while(ms()<budget){
        iter++;
        evalMakespan(seq,true);
        auto moves=n5moves(seq);
        if(moves.empty()){ seq=bestSeq; for(int t=0;t<4;t++){int m=rng()%M; if(seq[m].size()<2)continue; int i=rng()%(seq[m].size()-1); swap(seq[m][i],seq[m][i+1]);} stall=0; continue; }
        long long bestNm=LLONG_MAX; int bm=-1,bi=-1,ba=-1,bb=-1;
        for(auto&mv: moves){ int m=mv.first, i=mv.second; int a=seq[m][i], b=seq[m][i+1];
            swap(seq[m][i],seq[m][i+1]); long long nm=evalMakespan(seq); swap(seq[m][i],seq[m][i+1]);
            if(nm==LLONG_MAX) continue;
            long long key=keyOf(m,a,b); bool isTabu = tabuU.count(key)&&tabuU[key]>iter;
            if(isTabu && nm>=best) continue;
            if(nm<bestNm){ bestNm=nm; bm=m; bi=i; ba=a; bb=b; } }
        if(bm<0){ seq=bestSeq; for(int t=0;t<4;t++){int m=rng()%M; if(seq[m].size()<2)continue; int i=rng()%(seq[m].size()-1); swap(seq[m][i],seq[m][i+1]);} stall=0; continue; }
        swap(seq[bm][bi],seq[bm][bi+1]);
        tabuU[keyOf(bm,ba,bb)]=iter+tenure;
        if(bestNm<best){ best=bestNm; bestSeq=seq; stall=0; }
        else if(++stall>80){ seq=bestSeq; for(int t=0;t<6;t++){int m=rng()%M; if(seq[m].size()<2)continue; int i=rng()%(seq[m].size()-1); swap(seq[m][i],seq[m][i+1]);} stall=0; }
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
    // small GT-restart budget; tabu gets the rest. (large instances: even less GT so tabu can work)
    double gtFrac = NO>600 ? 0.06 : 0.12;
    while(ms()<TL_MS*gtFrac){ int rule=rng()%5; auto s=giffler(rule,rng); long long mk=evalMakespan(s); if(mk<bestMk){bestMk=mk;best=s;} }
    if(!best.empty()) bestMk=localSearch(best,bestMk,rng,TL_MS*0.97);
    for(int m=0;m<M;m++) for(int i=0;i<J;i++){ cout<<best[m][i]; cout<<(i+1<J?' ':'\n'); }
    return 0;
}
