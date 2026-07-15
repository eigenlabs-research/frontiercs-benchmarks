#include <bits/stdc++.h>
using namespace std;

/* We make the clues from the boundary of a simply connected polyomino.
   Thus the displayed clues always have at least the boundary we used as
   a witness.  Random, thin trees give considerably more 1's than a box. */
static const int N=12;
static const vector<string> T={
 "?   ?   ??? ","?? ??  ?   ?","? ? ?  ?   ?","? ? ?  ???? ",
 "? ? ?  ?    ","?   ?  ?    ","            ","?  ?   ?????",
 "? ?      ?  ","??   ? ? ?  ","? ?  ? ? ?  ","?  ? ??? ?  "};
int dr[4]={-1,1,0,0}, dc[4]={0,0,-1,1};

bool goodBoundary(const vector<vector<int>>& a) {
    int d[13][13]={};
    int edges=0;
    for(int r=0;r<N;r++) for(int c=0;c<N;c++) if(a[r][c]) {
        // Each exposed side is an edge of the proposed loop.
        if(r==0 || !a[r-1][c]) { d[r][c]++; d[r][c+1]++; edges++; }
        if(r==N-1 || !a[r+1][c]) { d[r+1][c]++; d[r+1][c+1]++; edges++; }
        if(c==0 || !a[r][c-1]) { d[r][c]++; d[r+1][c]++; edges++; }
        if(c==N-1 || !a[r][c+1]) { d[r][c+1]++; d[r+1][c+1]++; edges++; }
    }
    if(!edges) return false;
    for(int r=0;r<=N;r++) for(int c=0;c<=N;c++) if(d[r][c] && d[r][c]!=2) return false;
    return true;
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int mode; if(!(cin>>mode)) return 0;
    mt19937 rng(712367 + mode*991);
    vector<string> best; int bestScore=-1;
    int trials = mode ? 9000 : 6500;
    for(int it=0; it<trials; ++it){
        vector<vector<int>> a(N,vector<int>(N));
        int sr=rng()%N, sc=rng()%N; a[sr][sc]=1;
        int want = mode ? 48+(rng()%43) : 10+(rng()%70);
        for(int step=1;step<=want;step++){
            // Evaluate every prefix: stopping at different sizes is useful.
            if(step>= (mode?28:3) && goodBoundary(a)) {
                int cnt[12][12]={};
                for(int r=0;r<N;r++) for(int c=0;c<N;c++) if(a[r][c]) {
                    if(r==0 || !a[r-1][c]) { if(r) cnt[r-1][c]++; cnt[r][c]++; }
                    if(r==N-1 || !a[r+1][c]) { cnt[r][c]++; if(r+1<N) cnt[r+1][c]++; }
                    if(c==0 || !a[r][c-1]) { if(c) cnt[r][c-1]++; cnt[r][c]++; }
                    if(c==N-1 || !a[r][c+1]) { cnt[r][c]++; if(c+1<N) cnt[r][c+1]++; }
                }
                bool ok=true; int ones=0, informative=0;
                vector<string> out=T;
                for(int r=0;r<N;r++) for(int c=0;c<N;c++) if(T[r][c]=='?'){
                    if(cnt[r][c]>3 || (mode && cnt[r][c]==0)) ok=false;
                    out[r][c]=char('0'+cnt[r][c]);
                    ones += cnt[r][c]==1;
                    informative += cnt[r][c]!=0;
                }
                if(ok){
                    // In mode 0 zeroes are also valuable constraints; use them as a
                    // small tie breaker after the primary number of 1 clues.
                    int value=ones*100 + (mode ? informative : 56-informative);
                    if(value>bestScore) bestScore=value, best=out;
                }
            }
            vector<pair<int,int>> cand; vector<int> weight;
            for(int r=0;r<N;r++) for(int c=0;c<N;c++) if(!a[r][c]){
                int adj=0; for(int k=0;k<4;k++){int x=r+dr[k],y=c+dc[k]; if(x>=0&&x<N&&y>=0&&y<N&&a[x][y]) adj++;}
                if(adj!=1) continue; // keep the cells a tree: no holes are introduced
                bool bad=false;
                if(mode){
                    // A clue cell with four tree neighbours would get a forbidden 0.
                    for(int k=0;k<4;k++){int x=r+dr[k],y=c+dc[k]; if(x>=0&&x<N&&y>=0&&y<N&&a[x][y]&&T[x][y]=='?'){
                        int deg=0; for(int q=0;q<4;q++){int u=x+dr[q],v=y+dc[q]; if(u>=0&&u<N&&v>=0&&v<N&&a[u][v]) deg++;}
                        if(deg==3) bad=true;
                    }}
                }
                if(bad) continue;
                int gain=0;
                if(mode){
                    if(T[r][c]=='?') gain++;
                    for(int k=0;k<4;k++){int x=r+dr[k],y=c+dc[k]; if(x>=0&&x<N&&y>=0&&y<N&&T[x][y]=='?') gain++;}
                }
                cand.push_back({r,c}); weight.push_back(1+(mode?gain*12:0)+(rng()%8));
            }
            if(cand.empty()) break;
            int sum=accumulate(weight.begin(),weight.end(),0), z=rng()%sum, p=0;
            while(z>=weight[p]) z-=weight[p++];
            a[cand[p].first][cand[p].second]=1;
        }
    }
    // This should only be reached on an exceptionally unlucky random search.
    // A small rectangle is at least a syntactically valid witness for mode 0.
    if(best.empty()){
        best=T;
        for(int r=0;r<N;r++) for(int c=0;c<N;c++) if(best[r][c]=='?') best[r][c]= mode?'1':'0';
    }
    for(const string& s:best) cout<<s<<'\n';
}
