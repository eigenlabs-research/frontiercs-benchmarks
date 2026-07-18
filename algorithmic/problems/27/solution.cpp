#include <bits/stdc++.h>
using namespace std;

using pii = pair<int,int>;

static bool isPrime(int x){
    if(x < 2) return false;
    for(int d=2; d*d<=x; ++d) if(x%d==0) return false;
    return true;
}

static int nextPrimeWithProjectiveCapacity(long long need){
    for(int p=2;;++p){
        if(isPrime(p) && 1LL*p*p + p + 1 >= need) return p;
    }
}

static int nextBinaryPowerWithProjectiveCapacity(long long need){
    for(int q=2; q<=512; q<<=1){
        if(1LL*q*q + q + 1 >= need) return q;
    }
    return 1024;
}

static int gf2Poly(int q){
    switch(q){
        case 2: return 0x7;      // x^2+x+1 (reduction for one-bit elements)
        case 4: return 0x7;      // x^2+x+1
        case 8: return 0xb;      // x^3+x+1
        case 16: return 0x13;    // x^4+x+1
        case 32: return 0x25;    // x^5+x^2+1
        case 64: return 0x43;    // x^6+x+1
        case 128: return 0x83;   // x^7+x+1
        case 256: return 0x11d;  // AES/Rabin irreducible polynomial
        case 512: return 0x211;  // x^9+x^4+1
    }
    return 0;
}

static int gf2Mul(int a, int b, int q){
    int poly = gf2Poly(q);
    int res = 0;
    while(b){
        if(b & 1) res ^= a;
        b >>= 1;
        a <<= 1;
        if(a & q) a ^= poly;
    }
    return res & (q - 1);
}

// Columns are all 2-subsets of rows, then singleton columns.  This is optimal
// for the very-wide regime by the standard row-pair counting bound.
static vector<pii> pairColumns(int n, int m){
    vector<pii> ans;
    ans.reserve((size_t)min<long long>(1LL*n*m, 2LL*m));
    if(n == 1){
        for(int c=1;c<=m;c++) ans.push_back({1,c});
        return ans;
    }
    int c = 1;
    for(int i=1; i<=n && c<=m; ++i){
        for(int j=i+1; j<=n && c<=m; ++j){
            ans.push_back({i,c});
            ans.push_back({j,c});
            ++c;
        }
    }
    for(; c<=m; ++c) ans.push_back({1,c});
    return ans;
}

static vector<pii> pairRows(int n, int m){
    vector<pii> t = pairColumns(m, n);
    vector<pii> ans;
    ans.reserve(t.size());
    for(auto &e: t) ans.push_back({e.second, e.first});
    return ans;
}

// Incidence graph of a projective plane of prime order p, restricted to the
// requested prefix of points and lines.  Any two lines meet in at most one
// retained point, so no rectangle can occur.  Points are ordered as affine
// points (x,y), then points at infinity for each slope, then vertical infinity.
static vector<pii> projective(int n, int m){
    int p = nextPrimeWithProjectiveCapacity(max(n,m));
    int affine = p*p;
    auto pointIdAffine = [p](int x, int y){ return x*p + y; };          // 0-based
    auto pointIdInfSlope = [affine](int s){ return affine + s; };
    int pointIdInfVert = affine + p;

    vector<pii> ans;
    ans.reserve((size_t)min<long long>(1LL*n*m, 1LL*(p+1)*m));

    auto addPoint = [&](int pid, int col){
        if(pid < n) ans.push_back({pid + 1, col + 1});
    };

    for(int col=0; col<m; ++col){
        if(col < p*p){
            // non-vertical line y = s*x + t
            int s = col / p, t = col % p;
            for(int x=0; x<p; ++x){
                int y = (int)((1LL*s*x + t) % p);
                addPoint(pointIdAffine(x,y), col);
            }
            addPoint(pointIdInfSlope(s), col);
        }else if(col < p*p + p){
            // vertical line x = a
            int a = col - p*p;
            for(int y=0; y<p; ++y) addPoint(pointIdAffine(a,y), col);
            addPoint(pointIdInfVert, col);
        }else{
            // line at infinity
            for(int s=0; s<p; ++s) addPoint(pointIdInfSlope(s), col);
            addPoint(pointIdInfVert, col);
        }
    }
    return ans;
}

static vector<pii> projectiveTransposed(int n, int m){
    vector<pii> t = projective(m, n);
    vector<pii> ans;
    ans.reserve(t.size());
    for(auto &e: t) ans.push_back({e.second, e.first});
    return ans;
}

// Same projective-plane incidence construction over GF(2^a).  This avoids the
// sometimes-large jump to the next prime order (for example q=32 instead of
// q=37 around one thousand rows/columns), improving density while preserving
// the "two lines meet once" property.
static vector<pii> projectiveBinary(int n, int m){
    int q = nextBinaryPowerWithProjectiveCapacity(max(n,m));
    if(q > 512) return {};
    int affine = q*q;
    auto pointIdAffine = [q](int x, int y){ return x*q + y; };
    auto pointIdInfSlope = [affine](int s){ return affine + s; };
    int pointIdInfVert = affine + q;

    vector<pii> ans;
    ans.reserve((size_t)min<long long>(1LL*n*m, 1LL*(q+1)*m));
    auto addPoint = [&](int pid, int col){
        if(pid < n) ans.push_back({pid + 1, col + 1});
    };

    for(int col=0; col<m; ++col){
        if(col < q*q){
            int s = col / q, t = col % q;
            for(int x=0; x<q; ++x){
                int y = gf2Mul(s, x, q) ^ t;
                addPoint(pointIdAffine(x,y), col);
            }
            addPoint(pointIdInfSlope(s), col);
        }else if(col < q*q + q){
            int a = col - q*q;
            for(int y=0; y<q; ++y) addPoint(pointIdAffine(a,y), col);
            addPoint(pointIdInfVert, col);
        }else{
            for(int s=0; s<q; ++s) addPoint(pointIdInfSlope(s), col);
            addPoint(pointIdInfVert, col);
        }
    }
    return ans;
}

static vector<pii> projectiveBinaryTransposed(int n, int m){
    vector<pii> t = projectiveBinary(m, n);
    vector<pii> ans;
    ans.reserve(t.size());
    for(auto &e: t) ans.push_back({e.second, e.first});
    return ans;
}

static uint64_t rowPairKey(int a, int b){
    if(a > b) swap(a,b);
    return (uint64_t)(uint32_t)a << 32 | (uint32_t)b;
}

// Safe monotone repair/augmentation pass.  Starting from any valid C4-free set,
// try to insert currently empty cells; a cell (r,c) is legal exactly when row r
// has not already shared a different column with any row that is present in c.
static vector<pii> greedyAugment(int n, int m, const vector<pii>& base, int mode){
    if((long long)n * m == (long long)base.size()) return base;

    vector<unsigned char> occ((size_t)n * m, 0);
    vector<vector<int>> col(m + 1);
    unordered_set<uint64_t> shared;
    shared.reserve(base.size() * 8 + 1024);

    for(auto &e: base){
        int r = e.first, c = e.second;
        size_t id = (size_t)(r - 1) * m + (c - 1);
        if(!occ[id]){
            for(int rr: col[c]) shared.insert(rowPairKey(r, rr));
            occ[id] = 1;
            col[c].push_back(r);
        }
    }

    vector<pii> ans = base;
    ans.reserve((size_t)n * m);

    auto tryAdd = [&](int r, int c){
        size_t id = (size_t)(r - 1) * m + (c - 1);
        if(occ[id]) return;
        for(int rr: col[c]){
            if(shared.find(rowPairKey(r, rr)) != shared.end()) return;
        }
        for(int rr: col[c]) shared.insert(rowPairKey(r, rr));
        occ[id] = 1;
        col[c].push_back(r);
        ans.push_back({r,c});
    };

    if(mode == 0){
        for(int r=1; r<=n; ++r) for(int c=1; c<=m; ++c) tryAdd(r,c);
    }else if(mode == 1){
        for(int c=1; c<=m; ++c) for(int r=1; r<=n; ++r) tryAdd(r,c);
    }else if(mode == 2){
        int step = max(1, m / max(1, n));
        for(int r=1; r<=n; ++r){
            int start = (int)((1LL * (r - 1) * step) % m) + 1;
            for(int t=0; t<m; ++t){
                int c = ((start - 1 + t) % m) + 1;
                tryAdd(r,c);
            }
        }
    }else{
        int step = max(1, n / max(1, m));
        for(int c=1; c<=m; ++c){
            int start = (int)((1LL * (c - 1) * step) % n) + 1;
            for(int t=0; t<n; ++t){
                int r = ((start - 1 + t) % n) + 1;
                tryAdd(r,c);
            }
        }
    }
    return ans;
}

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n, m;
    if(!(cin >> n >> m)) return 0;

    vector<vector<pii>> cand;
    cand.push_back(pairColumns(n,m));
    cand.push_back(pairRows(n,m));
    cand.push_back(projective(n,m));
    cand.push_back(projectiveTransposed(n,m));
    cand.push_back(projectiveBinary(n,m));
    cand.push_back(projectiveBinaryTransposed(n,m));

    int initial = (int)cand.size();
    for(int i=0; i<initial; ++i){
        cand.push_back(greedyAugment(n, m, cand[i], 0));
        cand.push_back(greedyAugment(n, m, cand[i], 1));
        cand.push_back(greedyAugment(n, m, cand[i], 2));
        cand.push_back(greedyAugment(n, m, cand[i], 3));
    }

    int best = 0;
    for(int i=1;i<(int)cand.size();++i){
        if(cand[i].size() > cand[best].size()) best = i;
    }

    const auto &ans = cand[best];
    cout << ans.size() << '\n';
    for(auto &e: ans) cout << e.first << ' ' << e.second << '\n';
    return 0;
}
