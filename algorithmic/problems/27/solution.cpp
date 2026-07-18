#include <bits/stdc++.h>
using namespace std;

static long long C2(long long x){ return x*(x-1)/2; }

static void consider(vector<vector<int>>& best, const vector<vector<int>>& cand){
    long long a=0,b=0; for(auto &v: best) a += (int)v.size(); for(auto &v: cand) b += (int)v.size();
    if(b>a) best=cand;
}

static vector<vector<int>> pair_design(int R,int C){
    vector<vector<int>> col(C);
    int idx=0;
    for(int i=0;i<R && idx<C;i++) for(int j=i+1;j<R && idx<C;j++) col[idx++]={i,j};
    for(; idx<C; idx++) col[idx] = {idx%R};
    return col;
}

static vector<vector<int>> greedy_design(int R,int C,int shift,int bias){
    vector<vector<unsigned char>> used(R, vector<unsigned char>(R,0));
    vector<int> deg(R,0);
    vector<vector<int>> col(C);
    long long pairs=C2(R);
    int t=1;
    while(t<R && 1LL*C*C2(t+1)<=pairs) t++;
    long long rem = pairs - 1LL*C*C2(t);
    int high = (t>0? (int)min<long long>(C, rem / max(1,t)) : 0);
    vector<int> order(R);
    iota(order.begin(), order.end(), 0);
    for(int c=0;c<C;c++){
        int want = t + (c<high ? 1 : 0);
        if(bias && (c%7==0) && want>1) want--;
        sort(order.begin(), order.end(), [&](int a,int b){
            if(deg[a]!=deg[b]) return deg[a]<deg[b];
            return ((a+shift)%R) < ((b+shift)%R);
        });
        vector<int> block;
        for(int pass=0; pass<2 && (int)block.size()<want; pass++){
            for(int x: order){
                bool already=false; for(int y:block) if(y==x){already=true; break;}
                if(already) continue;
                bool ok=true;
                for(int y:block){ int a=min(x,y), b=max(x,y); if(used[a][b]){ ok=false; break; } }
                if(ok){ block.push_back(x); if((int)block.size()==want) break; }
            }
        }
        if(block.empty()) block.push_back((c+shift)%R);
        for(int i=0;i<(int)block.size();i++) for(int j=i+1;j<(int)block.size();j++){
            int a=min(block[i],block[j]), b=max(block[i],block[j]); used[a][b]=1;
        }
        for(int x:block) deg[x]++;
        col[c]=block;
    }
    return col;
}

static bool is_prime(int x){ if(x<2) return false; for(int d=2; d*d<=x; d++) if(x%d==0) return false; return true; }

static vector<vector<int>> affine_design(int R,int C,int q){
    vector<pair<int,int>> pts;
    pts.reserve(R);
    // balanced prefix over the q by q affine grid
    for(int y=0;y<q && (int)pts.size()<R;y++)
        for(int x=0;x<q && (int)pts.size()<R;x++) pts.push_back({x,y});
    vector<vector<int>> lines;
    lines.reserve(q*q+q);
    for(int a=0;a<q;a++) for(int b=0;b<q;b++){
        vector<int> v;
        for(int i=0;i<R;i++) if(pts[i].second == (a*pts[i].first + b)%q) v.push_back(i);
        if(!v.empty()) lines.push_back(v);
    }
    for(int x0=0;x0<q;x0++){
        vector<int> v;
        for(int i=0;i<R;i++) if(pts[i].first==x0) v.push_back(i);
        if(!v.empty()) lines.push_back(v);
    }
    sort(lines.begin(), lines.end(), [](const vector<int>&a,const vector<int>&b){ return a.size()>b.size(); });
    vector<vector<int>> col(C);
    int idx=0;
    for(; idx<C && idx<(int)lines.size(); idx++) col[idx]=lines[idx];
    for(; idx<C; idx++) col[idx] = {idx%R};
    return col;
}

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n,m; if(!(cin>>n>>m)) return 0;
    if(n==1){ cout<<m<<'\n'; for(int c=1;c<=m;c++) cout<<1<<' '<<c<<'\n'; return 0; }
    if(m==1){ cout<<n<<'\n'; for(int r=1;r<=n;r++) cout<<r<<' '<<1<<'\n'; return 0; }
    bool trans = n>m;
    int R = min(n,m), C = max(n,m);
    vector<vector<int>> best;
    consider(best, pair_design(R,C));
    for(int s=0; s<min(R,16); s++){
        consider(best, greedy_design(R,C,s,0));
        consider(best, greedy_design(R,C,s,1));
    }
    int start=1; while(start*start<R) start++;
    for(int q=start; q<=max(start+40, 2*start+10); q++) if(is_prime(q)) consider(best, affine_design(R,C,q));
    vector<pair<int,int>> ans;
    for(int c=0;c<C;c++) for(int r: best[c]){
        if(!trans) ans.push_back({r+1,c+1});
        else ans.push_back({c+1,r+1});
    }
    cout << ans.size() << '\n';
    for(auto &p: ans) cout << p.first << ' ' << p.second << '\n';
    return 0;
}
