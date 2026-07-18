#include <bits/stdc++.h>
using namespace std;

static vector<int> primes_up_to(int n){
    vector<int> ps;
    for(int x=2;x<=n;x++){
        bool ok=true;
        for(int d=2;d*d<=x;d++) if(x%d==0){ok=false;break;}
        if(ok) ps.push_back(x);
    }
    return ps;
}

static void add_line_inc(int lid, int q, int C, int row, vector<pair<int,int>>& out){
    int q2=q*q;
    auto addp = [&](int pid){ if(pid>=0 && pid<C) out.push_back({row+1, pid+1}); };
    if(lid < q2){
        int a = lid / q, b = lid % q;
        for(int x=0;x<q;x++){
            int y = ( (long long)a*x + b ) % q;
            addp(x*q + y);
        }
        addp(q2 + a);
    }else if(lid < q2 + q){
        int c = lid - q2;
        for(int y=0;y<q;y++) addp(c*q + y);
        addp(q2 + q);
    }else{
        for(int a=0;a<q;a++) addp(q2 + a);
        addp(q2 + q);
    }
}

static vector<pair<int,int>> projective_candidate(int R, int C, int q){
    int v = q*q + q + 1;
    int rr = min(R, v);
    vector<pair<int,int>> out;
    out.reserve((size_t)rr * (q+1));
    for(int lid=0; lid<rr; ++lid) add_line_inc(lid, q, C, lid, out);
    return out;
}

static vector<pair<int,int>> star_candidate(int n, int m){
    vector<pair<int,int>> out;
    out.reserve((size_t)n + m);
    if(n==0 || m==0) return out;
    for(int c=1;c<=m;c++) out.push_back({1,c});
    for(int r=2;r<=n;r++) out.push_back({r,1});
    return out;
}

// Use columns as distinct edges of K_n (plus singleton filler columns).  Then any
// two rows occur together in at most one column, exactly the rectangle-free
// condition.  This dominates the star construction on very skinny instances
// where a projective-plane block is too short in the wide direction.
static vector<pair<int,int>> edge_pack_candidate(int n, int m){
    vector<pair<int,int>> out;
    long long maxEdges = 1LL*n*(n-1)/2;
    int eLim = (int)min<long long>(m, maxEdges);
    out.reserve((size_t)m + eLim);
    int col = 1;
    for(int i=1; i<=n && col<=eLim; ++i){
        for(int j=i+1; j<=n && col<=eLim; ++j){
            out.push_back({i, col});
            out.push_back({j, col});
            ++col;
        }
    }
    for(; col<=m; ++col) out.push_back({1, col});
    return out;
}

static uint64_t row_pair_key(int a, int b){
    if(a > b) swap(a,b);
    return (uint64_t)(uint32_t)a << 32 | (uint32_t)b;
}

// Deterministically make the selected algebraic/edge construction maximal by
// adding any safe missing cells.  This is cheap because n*m <= 100000, and it
// often recovers points lost by rectangular clipping of a projective plane.
static void greedy_augment(int n, int m, vector<pair<int,int>>& pts){
    if(n <= 1 || m <= 1){
        vector<char> seen((size_t)n*m, 0);
        vector<pair<int,int>> clean;
        clean.reserve((size_t)n*m);
        for(auto p: pts){
            size_t id = (size_t)(p.first-1)*m + (p.second-1);
            if(!seen[id]){ seen[id]=1; clean.push_back(p); }
        }
        pts.swap(clean);
        for(int r=1;r<=n;r++) for(int c=1;c<=m;c++){
            size_t id=(size_t)(r-1)*m+(c-1);
            if(!seen[id]){ seen[id]=1; pts.push_back({r,c}); }
        }
        return;
    }

    vector<char> occ((size_t)n*m, 0);
    vector<vector<int>> colRows(m+1);
    vector<pair<int,int>> clean;
    clean.reserve(pts.size() + 1024);
    for(auto p: pts){
        if(p.first < 1 || p.first > n || p.second < 1 || p.second > m) continue;
        size_t id = (size_t)(p.first-1)*m + (p.second-1);
        if(!occ[id]){
            occ[id] = 1;
            clean.push_back(p);
            colRows[p.second].push_back(p.first);
        }
    }
    pts.swap(clean);

    long long pairEstimate = 0;
    for(int c=1;c<=m;c++){
        long long d = (long long)colRows[c].size();
        pairEstimate += d*(d-1)/2;
        if(pairEstimate > 3000000LL) return; // avoid pathological skinny stars
    }

    unordered_set<uint64_t> used;
    used.reserve((size_t)pairEstimate*2 + 1024);
    for(int c=1;c<=m;c++){
        auto &v = colRows[c];
        for(int i=0;i<(int)v.size();i++) for(int j=i+1;j<(int)v.size();j++)
            used.insert(row_pair_key(v[i], v[j]));
    }

    auto try_add = [&](int r, int c)->bool{
        size_t id = (size_t)(r-1)*m + (c-1);
        if(occ[id]) return false;
        for(int rr: colRows[c]) if(used.find(row_pair_key(r, rr)) != used.end()) return false;
        occ[id] = 1;
        for(int rr: colRows[c]) used.insert(row_pair_key(r, rr));
        colRows[c].push_back(r);
        pts.push_back({r,c});
        return true;
    };

    // Two complementary sweeps; the second can exploit pairs opened by the first
    // while remaining deterministic and validator-safe.
    for(int r=1;r<=n;r++) for(int c=1;c<=m;c++) try_add(r,c);
    for(int c=1;c<=m;c++) for(int r=n;r>=1;r--) try_add(r,c);
}

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n,m;
    if(!(cin>>n>>m)) return 0;

    vector<pair<int,int>> best = star_candidate(n,m);
    auto edge = edge_pack_candidate(n,m);
    if(edge.size() > best.size()) best.swap(edge);
    auto edgeT = edge_pack_candidate(m,n);
    for(auto &p: edgeT) swap(p.first, p.second);
    if(edgeT.size() > best.size()) best.swap(edgeT);

    int mn = min(n,m), mx = max(n,m);
    if(mn >= 20){
        int lim = (int)(sqrt(mx) + 80);
        lim = max(lim, 20);
        lim = min(lim, 700);
        vector<int> ps = primes_up_to(lim);
        for(int q: ps){
            // Incidence graph of the projective plane PG(2,q), clipped to the requested rectangle.
            auto cand = projective_candidate(n,m,q);
            if(cand.size() > best.size()) best.swap(cand);
            auto tc = projective_candidate(m,n,q);
            for(auto &p: tc) swap(p.first, p.second);
            if(tc.size() > best.size()) best.swap(tc);
        }
    }

    greedy_augment(n, m, best);

    cout << best.size() << '\n';
    for(auto &p: best) cout << p.first << ' ' << p.second << '\n';
    return 0;
}
