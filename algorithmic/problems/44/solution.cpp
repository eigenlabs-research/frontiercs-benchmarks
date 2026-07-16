#include <bits/stdc++.h>
using namespace std;

struct Pt { long long x, y; };

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n;
    if (!(cin >> n)) return 0;
    vector<Pt> a(n);
    long long xmin = LLONG_MAX, xmax = LLONG_MIN, ymin = LLONG_MAX, ymax = LLONG_MIN;
    for (auto &p : a) {
        cin >> p.x >> p.y;
        xmin = min(xmin, p.x); xmax = max(xmax, p.x);
        ymin = min(ymin, p.y); ymax = max(ymax, p.y);
    }
    vector<char> prime(max(2, n), true);
    prime[0] = false; prime[1] = false;
    for (int i = 2; 1LL * i * i < n; ++i) if (prime[i])
        for (int j = i * i; j < n; j += i) prime[j] = false;

    auto distance = [&](int u, int v) {
        return hypot((double)a[u].x - a[v].x, (double)a[u].y - a[v].y);
    };
    auto edgeCost = [&](const vector<int>& p, int e) {
        double z = distance(p[e], p[e + 1]);
        if ((e + 1) % 10 == 0 && !prime[p[e]]) z *= 1.1;
        return z;
    };
    auto score = [&](const vector<int>& p) {
        double z = 0;
        for (int e = 0; e < n; ++e) z += edgeCost(p, e);
        return z;
    };
    vector<int> best;
    double bestScore = numeric_limits<double>::infinity();
    auto consider = [&](vector<int> ord) {
        for (int rev = 0; rev < 2; ++rev) {
            if (rev) reverse(ord.begin(), ord.end());
            vector<int> p; p.reserve(n + 1); p.push_back(0);
            for (int v : ord) p.push_back(v);
            p.push_back(0);
            double z = score(p);
            if (z < bestScore) { bestScore = z; best.swap(p); }
            if (rev) reverse(ord.begin(), ord.end());
        }
    };

    // Horizontal strip sweeps: particularly effective when nearby cities form rows.
    vector<int> byY;
    for (int i = 1; i < n; ++i) byY.push_back(i);
    sort(byY.begin(), byY.end(), [&](int u, int v) {
        if (a[u].y != a[v].y) return a[u].y < a[v].y;
        return a[u].x < a[v].x;
    });
    vector<int> bands = {1, 2, 4, 8, 16, 32, 64, 128};
    for (int b : bands) {
        if (b > (int)byY.size() || b == 0) continue;
        vector<int> ord;
        for (int k = 0; k < b; ++k) {
            int l = (long long)k * byY.size() / b;
            int r = (long long)(k + 1) * byY.size() / b;
            vector<int> part(byY.begin() + l, byY.begin() + r);
            sort(part.begin(), part.end()); // IDs are increasing x
            if (k & 1) reverse(part.begin(), part.end());
            ord.insert(ord.end(), part.begin(), part.end());
        }
        consider(move(ord));
    }

    // Vertical strip sweeps complement horizontal ones on column-like instances.
    vector<int> ids;
    for (int i = 1; i < n; ++i) ids.push_back(i);
    for (int b : bands) {
        if (b > (int)ids.size() || b == 0) continue;
        vector<int> ord;
        for (int k = 0; k < b; ++k) {
            int l = (long long)k * ids.size() / b;
            int r = (long long)(k + 1) * ids.size() / b;
            vector<int> part(ids.begin() + l, ids.begin() + r);
            sort(part.begin(), part.end(), [&](int u, int v) {
                if (a[u].y != a[v].y) return a[u].y < a[v].y;
                return u < v;
            });
            if (k & 1) reverse(part.begin(), part.end());
            ord.insert(ord.end(), part.begin(), part.end());
        }
        consider(move(ord));
    }

    // A Morton traversal supplies a scale-independent fallback for irregular clouds.
    auto morton = [&](int id) {
        auto norm = [](long long v, long long lo, long long hi) -> unsigned long long {
            if (hi == lo) return 0;
            return (unsigned long long)((__int128)(v - lo) * ((1ULL << 21) - 1) / (hi - lo));
        };
        unsigned long long x = norm(a[id].x, xmin, xmax), y = norm(a[id].y, ymin, ymax), z = 0;
        for (int k = 0; k < 21; ++k) z |= ((x >> k) & 1ULL) << (2 * k), z |= ((y >> k) & 1ULL) << (2 * k + 1);
        return z;
    };
    vector<unsigned long long> mkey(n);
    for (int i = 1; i < n; ++i) mkey[i] = morton(i);
    vector<int> mo = ids;
    sort(mo.begin(), mo.end(), [&](int u, int v) { return mkey[u] < mkey[v]; });
    consider(move(mo));

    // Exact greedy nearest-neighbor construction.  The kd tree stores an alive
    // count per subtree, so deleted cities are skipped without rebuilding it.
    // This is deliberately a separate seed rather than an approximation to a
    // space-filling curve: it reacts directly to local cluster geometry.
    vector<int> kdOrder = ids, kdLeft(n, -1), kdRight(n, -1), kdParent(n, -1), kdAt(n, -1), where(n, -1), kdAlive(n, 0);
    vector<long long> lx(n), rx(n), ly(n), ry(n);
    int kdNodes = 0;
    function<int(int,int,int,int)> build = [&](int l, int r, int dep, int par) -> int {
        if (l >= r) return -1;
        int mid = (l + r) / 2, axis = dep & 1;
        nth_element(kdOrder.begin() + l, kdOrder.begin() + mid, kdOrder.begin() + r, [&](int u, int v) {
            long long au = axis ? a[u].y : a[u].x, av = axis ? a[v].y : a[v].x;
            return au == av ? u < v : au < av;
        });
        int z = kdNodes++, v = kdOrder[mid];
        kdAt[z] = v; where[v] = z; kdParent[z] = par;
        kdLeft[z] = build(l, mid, dep + 1, z);
        kdRight[z] = build(mid + 1, r, dep + 1, z);
        lx[z] = rx[z] = a[v].x; ly[z] = ry[z] = a[v].y; kdAlive[z] = 1;
        for (int c : {kdLeft[z], kdRight[z]}) if (c != -1) {
            lx[z] = min(lx[z], lx[c]); rx[z] = max(rx[z], rx[c]);
            ly[z] = min(ly[z], ly[c]); ry[z] = max(ry[z], ry[c]);
            kdAlive[z] += kdAlive[c];
        }
        return z;
    };
    int kdRoot = build(0, (int)kdOrder.size(), 0, -1);
    auto boxDist = [&](int z, int v) -> double {
        double dx = 0, dy = 0;
        if (a[v].x < lx[z]) dx = (double)lx[z] - a[v].x;
        else if (a[v].x > rx[z]) dx = (double)a[v].x - rx[z];
        if (a[v].y < ly[z]) dy = (double)ly[z] - a[v].y;
        else if (a[v].y > ry[z]) dy = (double)a[v].y - ry[z];
        return dx * dx + dy * dy;
    };
    vector<char> used(n, false); used[0] = true;
    vector<int> nn; nn.reserve(n - 1);
    int cur = 0;
    for (int step = 1; step < n; ++step) {
        int near = -1;
        double bestD = numeric_limits<double>::infinity();
        function<void(int)> query = [&](int z) {
            if (z == -1 || kdAlive[z] == 0 || boxDist(z, cur) >= bestD) return;
            int v = kdAt[z];
            if (!used[v]) {
                double dx = (double)a[v].x - a[cur].x, dy = (double)a[v].y - a[cur].y;
                double d = dx * dx + dy * dy;
                if (d < bestD || (d == bestD && v < near)) bestD = d, near = v;
            }
            int c1 = kdLeft[z], c2 = kdRight[z];
            if (c1 != -1 && c2 != -1 && boxDist(c2, cur) < boxDist(c1, cur)) swap(c1, c2);
            query(c1); query(c2);
        };
        query(kdRoot);
        // A nonempty root guarantees a result; this guard also keeps the route
        // valid should an implementation platform ever encounter an anomaly.
        if (near == -1) for (int v = 1; v < n; ++v) if (!used[v]) { near = v; break; }
        nn.push_back(near); used[near] = true; cur = near;
        for (int z = where[near]; z != -1; z = kdParent[z]) --kdAlive[z];
    }
    consider(move(nn));

    // At a tenth-step source, try nearby prime cities.  The exact affected edges
    // are evaluated, so a swap is accepted only if it really lowers carrot cost.
    for (int pos = 9; pos < n; pos += 10) {
        if (prime[best[pos]]) continue;
        int chosen = -1;
        double gain = 0;
        for (int q = max(1, pos - 16); q <= min(n - 1, pos + 16); ++q) {
            if (!prime[best[q]]) continue;
            vector<int> es = {pos - 1, pos, q - 1, q};
            sort(es.begin(), es.end()); es.erase(unique(es.begin(), es.end()), es.end());
            double before = 0, after = 0;
            for (int e : es) if (e >= 0 && e < n) before += edgeCost(best, e);
            swap(best[pos], best[q]);
            for (int e : es) if (e >= 0 && e < n) after += edgeCost(best, e);
            swap(best[pos], best[q]);
            if (before - after > gain) gain = before - after, chosen = q;
        }
        if (chosen != -1) swap(best[pos], best[chosen]);
    }

    // Ablation: exact short 2-opt reversals.  Unlike ordinary Euclidean 2-opt,
    // reversing changes which city is the source at tenth steps, so evaluate
    // every edge in the (small) reversed interval with its actual multiplier.
    // The bounded neighborhood keeps this linear-time in practice on 200k nodes.
    const int WINDOW = 10;
    for (int i = 0; i + 2 < n; ++i) {
        for (int len = 2; len <= WINDOW && i + len < n; ++len) {
            int j = i + len;
            double before = 0.0, after = 0.0;
            for (int e = i; e <= j; ++e) before += edgeCost(best, e);
            reverse(best.begin() + i + 1, best.begin() + j + 1);
            for (int e = i; e <= j; ++e) after += edgeCost(best, e);
            if (after + 1e-9 < before) {
                // Keep this strictly improving move; subsequent candidates see it.
            } else {
                reverse(best.begin() + i + 1, best.begin() + j + 1);
            }
        }
    }

    cout << n + 1 << '\n';
    for (int v : best) cout << v << '\n';
}
