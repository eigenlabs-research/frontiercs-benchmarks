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

    cout << n + 1 << '\n';
    for (int v : best) cout << v << '\n';
}
