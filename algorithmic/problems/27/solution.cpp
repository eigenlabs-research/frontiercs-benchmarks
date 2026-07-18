#include <bits/stdc++.h>
using namespace std;

struct Cand {
    vector<pair<int,int>> pts;
};

static bool isPrime(int x){
    if(x<2) return false;
    for(int d=2; d*d<=x; ++d) if(x%d==0) return false;
    return true;
}

static Cand greedyPacking(int n, int m, uint64_t seed){
    // Work with the smaller side as "rows" of a packing of column subsets.
    bool tr = false;
    int a=n, b=m;
    if(a>b){ tr=true; swap(a,b); }
    Cand res;
    vector<vector<unsigned char>> used(a, vector<unsigned char>(a, 0));
    long long remPairs = 1LL*a*(a-1)/2;
    vector<int> order(a);
    iota(order.begin(), order.end(), 0);
    mt19937_64 rng(seed);

    for(int col=0; col<b; ++col){
        int remCols = b-col;
        int base = 1;
        while(base+1<=a && 1LL*(base+1)*base/2 <= remPairs / max(1,remCols)) base++;
        int target = base;
        long long cbase = 1LL*base*(base-1)/2;
        long long cnext = 1LL*(base+1)*base/2;
        if(base+1<=a && cnext + 1LL*(remCols-1)*cbase <= remPairs) target = base+1;
        target = max(1, target);

        shuffle(order.begin(), order.end(), rng);
        vector<int> sub;
        sub.reserve(target);
        // Several randomized passes; adding a row is legal iff it shares no already-used pair
        // with the rows already chosen for this column.
        for(int pass=0; pass<3 && (int)sub.size()<target; ++pass){
            if(pass) shuffle(order.begin(), order.end(), rng);
            for(int r: order){
                if((int)sub.size()>=target) break;
                bool already=false, ok=true;
                for(int x: sub){
                    if(x==r){ already=true; break; }
                    int u=min(x,r), v=max(x,r);
                    if(used[u][v]){ ok=false; break; }
                }
                if(!already && ok) sub.push_back(r);
            }
        }
        if(sub.empty()) sub.push_back(col%a);
        for(int i=0;i<(int)sub.size();++i) for(int j=i+1;j<(int)sub.size();++j){
            int u=min(sub[i],sub[j]), v=max(sub[i],sub[j]);
            if(!used[u][v]){ used[u][v]=1; remPairs--; }
        }
        for(int r: sub){
            if(!tr) res.pts.push_back({r+1, col+1});
            else    res.pts.push_back({col+1, r+1});
        }
    }
    return res;
}

static Cand affineConstruction(int n, int m, int q, bool tr){
    // Incidence between points (x,y) and non-vertical lines y=s*x+t over F_q.
    // q is prime.  Truncating points/lines preserves K_{2,2}-freeness.
    int R = tr ? m : n;
    int C = tr ? n : m;
    Cand res;
    for(int i=0;i<R;i++){
        int x=i/q, y=i%q;
        for(int j=0;j<C;j++){
            int s=j/q, t=j%q;
            if(( (long long)s*x + t ) % q == y){
                if(!tr) res.pts.push_back({i+1,j+1});
                else    res.pts.push_back({j+1,i+1});
            }
        }
    }
    return res;
}

static Cand bestAffineLines(int n, int m, int q, bool tr){
    // Same affine plane, but allow the q vertical lines too and choose the C
    // line-columns with largest intersection with the truncated point set.
    // This is a targeted improvement over taking the first C non-vertical lines:
    // for arbitrary n,m the boundary of the q by q point square is uneven, so
    // vertical/high-degree lines can add points while preserving the linear-space
    // property (two points determine at most one selected line).
    int R = tr ? m : n;
    int C = tr ? n : m;
    if(R > q*q || C > q*q + q) return Cand{};

    vector<vector<int>> lines;
    lines.reserve(q*q + q);
    for(int s=0; s<q; ++s){
        for(int t=0; t<q; ++t){
            vector<int> v;
            v.reserve(q);
            for(int x=0; x<q; ++x){
                int y = ((long long)s*x + t) % q;
                int id = x*q + y;
                if(id < R) v.push_back(id);
            }
            lines.push_back(move(v));
        }
    }
    for(int x=0; x<q; ++x){
        vector<int> v;
        v.reserve(q);
        for(int y=0; y<q; ++y){
            int id = x*q + y;
            if(id < R) v.push_back(id);
        }
        lines.push_back(move(v));
    }

    vector<int> ord(lines.size());
    iota(ord.begin(), ord.end(), 0);
    sort(ord.begin(), ord.end(), [&](int a, int b){
        if(lines[a].size() != lines[b].size()) return lines[a].size() > lines[b].size();
        return a < b;
    });

    Cand res;
    for(int newc=0; newc<C; ++newc){
        for(int id: lines[ord[newc]]){
            if(!tr) res.pts.push_back({id+1, newc+1});
            else    res.pts.push_back({newc+1, id+1});
        }
    }
    return res;
}

static bool validQuick(const Cand& c, int n, int m){
    vector<vector<int>> cols(m+1);
    vector<unsigned long long> enc; enc.reserve(c.pts.size());
    for(auto [r,cc]: c.pts){
        if(r<1||r>n||cc<1||cc>m) return false;
        enc.push_back((unsigned long long)(r-1)*(unsigned long long)m + (cc-1));
        cols[cc].push_back(r);
    }
    sort(enc.begin(), enc.end());
    for(size_t i=1;i<enc.size();++i) if(enc[i]==enc[i-1]) return false;
    unordered_set<unsigned long long> seen;
    seen.reserve(c.pts.size()*4+10);
    for(auto &v: cols){
        sort(v.begin(), v.end());
        for(int i=0;i<(int)v.size();++i) for(int j=i+1;j<(int)v.size();++j){
            unsigned long long key=((unsigned long long)v[i]<<32) | (unsigned int)v[j];
            if(!seen.insert(key).second) return false;
        }
    }
    return true;
}

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n,m;
    if(!(cin>>n>>m)) return 0;

    Cand best;
    auto consider = [&](Cand &&c){
        if(c.pts.size() > best.pts.size() && validQuick(c,n,m)) best = move(c);
    };

    // Main general mechanism: randomized partial linear-space packing on the smaller side.
    for(int s=0; s<24; ++s) consider(greedyPacking(n,m,0x9e3779b97f4a7c15ULL + 1000003ULL*s + (uint64_t)n*9176 + m));

    // Algebraic linear-space candidates; often stronger for square-ish prime-order regions.
    int need = (int)ceil(sqrt((double)max(n,m)));
    int tried=0;
    for(int q=max(2,need); q<=450 && tried<18; ++q){
        if(!isPrime(q)) continue;
        if(1LL*q*q < max(n,m)) continue;
        tried++;
        consider(affineConstruction(n,m,q,false));
        consider(affineConstruction(n,m,q,true));
        // The degree-sorted line selection is O(q^3), so run it only for the
        // nearest field size; larger q rarely offsets the lower density.
        if(tried == 1){
            consider(bestAffineLines(n,m,q,false));
            consider(bestAffineLines(n,m,q,true));
        }
    }

    cout << best.pts.size() << '\n';
    for(auto [r,c]: best.pts) cout << r << ' ' << c << '\n';
    return 0;
}
