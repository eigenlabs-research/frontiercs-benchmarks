#include <bits/stdc++.h>
using namespace std;

static bool isPrime(int x){
    if(x < 2) return false;
    for(int d=2; d*d<=x; ++d) if(x%d==0) return false;
    return true;
}

static vector<int> candidatePrimes(int s){
    int lo = 2;
    while(1LL*lo*lo < s) ++lo;
    vector<int> ps;
    int limit;
    if(s <= 1000) limit = max(lo, min(450, 3*lo + 50)); // small side: try varied embeddings
    else if(s <= 10000) limit = max(lo, min(260, 2*lo + 40));
    else limit = lo + 30;                        // large side: only nearby primes
    for(int q=lo; q<=limit; ++q) if(isPrime(q)) ps.push_back(q);
    if(ps.empty()){
        int q=lo;
        while(!isPrime(q)) ++q;
        ps.push_back(q);
    }
    return ps;
}

// Treat one grid dimension as points of an affine plane and the other as lines/blocks.
// Any two points lie on at most one selected line, so no two grid rows share two columns.
static vector<pair<int,int>> buildWithQ(int pointCount, int blockCount, bool transposed, int q){
    int lineCount = q*q + q; // y=a*x+b plus vertical x=c
    vector<vector<int>> lines(lineCount);

    for(int id=0; id<pointCount; ++id){
        int x = id / q;
        int y = id % q;
        lines[q*q + x].push_back(id + 1); // vertical line
        for(int a=0; a<q; ++a){
            int b = (y - (int)(1LL*a*x % q) + q) % q;
            lines[a*q + b].push_back(id + 1);
        }
    }

    vector<int> ord(lineCount);
    iota(ord.begin(), ord.end(), 0);
    stable_sort(ord.begin(), ord.end(), [&](int u, int v){
        return lines[u].size() > lines[v].size();
    });

    vector<pair<int,int>> ans;
    ans.reserve((size_t)pointCount * (size_t)blockCount);
    int usedBlocks = 0;
    for(int id: ord){
        if(usedBlocks >= blockCount) break;
        if(lines[id].empty()) break;
        ++usedBlocks;
        for(int p: lines[id]){
            if(!transposed) ans.push_back({p, usedBlocks});
            else ans.push_back({usedBlocks, p});
        }
    }

    // Extra block positions are filled with singletons; singletons introduce no row-pair.
    while(usedBlocks < blockCount){
        ++usedBlocks;
        if(!transposed) ans.push_back({1, usedBlocks});
        else ans.push_back({usedBlocks, 1});
    }
    return ans;
}

static vector<pair<int,int>> bestAffine(int points, int blocks, bool transposed){
    vector<pair<int,int>> best;
    for(int q: candidatePrimes(points)){
        if(1LL*q*q < points) continue;
        vector<pair<int,int>> cur = buildWithQ(points, blocks, transposed, q);
        if(cur.size() > best.size()) best.swap(cur);
    }
    return best;
}

// A complementary packing for very rectangular cases.  Use each unordered pair of
// rows in at most one column; all remaining columns get a singleton.  This is
// often better when one side is small and the other side has many spare columns.
static vector<pair<int,int>> buildPairColumns(int rows, int cols, bool transposed){
    vector<pair<int,int>> ans;
    ans.reserve((size_t)cols + (size_t)min<long long>(cols, 1LL*rows*(rows-1)/2));
    int used = 0;
    for(int i=1; i<=rows && used<cols; ++i){
        for(int j=i+1; j<=rows && used<cols; ++j){
            ++used;
            if(!transposed){
                ans.push_back({i, used});
                ans.push_back({j, used});
            }else{
                ans.push_back({used, i});
                ans.push_back({used, j});
            }
        }
    }
    while(used < cols){
        ++used;
        if(!transposed) ans.push_back({1, used});
        else ans.push_back({used, 1});
    }
    return ans;
}

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n, m;
    if(!(cin >> n >> m)) return 0;

    vector<pair<int,int>> ans;
    if(n == 1){
        for(int c=1;c<=m;++c) ans.push_back({1,c});
    } else if(m == 1){
        for(int r=1;r<=n;++r) ans.push_back({r,1});
    } else {
        vector<pair<int,int>> a = bestAffine(n, m, false); // rows are affine-plane points
        vector<pair<int,int>> b = bestAffine(m, n, true);  // columns are affine-plane points
        vector<pair<int,int>> c = buildPairColumns(n, m, false);
        vector<pair<int,int>> d = buildPairColumns(m, n, true);
        ans = std::move(a);
        if(b.size() > ans.size()) ans = std::move(b);
        if(c.size() > ans.size()) ans = std::move(c);
        if(d.size() > ans.size()) ans = std::move(d);
    }

    cout << ans.size() << '\n';
    for(auto &p: ans) cout << p.first << ' ' << p.second << '\n';
    return 0;
}
