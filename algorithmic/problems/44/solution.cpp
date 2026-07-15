#include <bits/stdc++.h>
using namespace std;

struct Point { long long x, y; };

// Hilbert index of a point in a 2^30 by 2^30 square.  Unlike sorting on one
// coordinate, this keeps nearby points together in both coordinate directions.
static unsigned long long hilbert(long long x, long long y) {
    unsigned long long d = 0;
    for (long long s = 1LL << 29; s; s >>= 1) {
        long long rx = (x & s) != 0, ry = (y & s) != 0;
        d += (unsigned long long)(s * s) * (unsigned long long)((3 * rx) ^ ry);
        if (!ry) {
            if (rx) { x = s - 1 - x; y = s - 1 - y; }
            swap(x, y);
        }
    }
    return d;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n;
    if (!(cin >> n)) return 0;
    vector<Point> a(n);
    long long minx = LLONG_MAX, maxx = LLONG_MIN, miny = LLONG_MAX, maxy = LLONG_MIN;
    for (auto &p : a) {
        cin >> p.x >> p.y;
        minx = min(minx, p.x); maxx = max(maxx, p.x);
        miny = min(miny, p.y); maxy = max(maxy, p.y);
    }
    const long long S = (1LL << 30) - 1;
    long long rangex = max(1LL, maxx - minx), rangey = max(1LL, maxy - miny);
    // City 0 fixes where the cyclic traversal is cut.  A Hilbert curve has
    // orientation-dependent seams, so retain several dihedral orientations
    // rather than assuming the conventional one is suitable for that cut.
    vector<vector<pair<unsigned long long,int>>> orders(4);
    for (auto &ord : orders) ord.reserve(n);
    for (int i = 0; i < n; ++i) {
        // Products fit in signed 64 bits for the specified coordinate range.
        long long xx = (a[i].x - minx) * S / rangex;
        long long yy = (a[i].y - miny) * S / rangey;
        long long tx[4] = {xx, yy, S - xx, xx};
        long long ty[4] = {yy, xx, yy, S - yy};
        for (int q = 0; q < 4; ++q) orders[q].push_back({hilbert(tx[q], ty[q]), i});
    }
    for (auto &ord : orders) sort(ord.begin(), ord.end());

    vector<char> prime(n, true);
    if (n > 0) prime[0] = false;
    if (n > 1) prime[1] = false;
    for (int i = 2; 1LL * i * i < n; ++i)
        if (prime[i]) for (int j = i * i; j < n; j += i) prime[j] = false;

    auto edge = [&](int u, int v, int step) -> long double {
        long double dx = (long double)a[u].x - a[v].x;
        long double dy = (long double)a[u].y - a[v].y;
        long double z = sqrtl(dx * dx + dy * dy);
        return (step % 10 == 0 && !prime[u]) ? z * 1.1L : z;
    };
    auto routeCost = [&](const vector<int>& p) {
        long double ans = 0;
        for (int t = 1; t <= n; ++t) ans += edge(p[t-1], p[t], t);
        return ans;
    };

    // A Hilbert traversal is cyclic here; rotate every orientation at the
    // required endpoint and score both directions with the exact carrot rule.
    vector<int> path;
    long double selectedCost = numeric_limits<long double>::infinity();
    for (const auto &ord : orders) {
        int at0 = 0;
        while (ord[at0].second != 0) ++at0;
        vector<int> candidate;
        candidate.reserve(n + 1);
        candidate.push_back(0);
        for (int k = 1; k < n; ++k) candidate.push_back(ord[(at0 + k) % n].second);
        candidate.push_back(0);
        long double c = routeCost(candidate);
        if (c < selectedCost) {
            selectedCost = c;
            path = candidate;
        }
        reverse(candidate.begin() + 1, candidate.end() - 1);
        c = routeCost(candidate);
        if (c < selectedCost) {
            selectedCost = c;
            path = candidate;
        }
    }

    // A raster tour is a deliberately different family from a space-filling
    // curve.  It is especially useful for data arranged in horizontal layers:
    // consecutive bands are crossed only once, while each band is swept in the
    // opposite x direction.  Keep several resolutions and let exact scoring
    // decide whether this decomposition beats Hilbert on this particular map.
    for (int bands : {2, 4, 8, 16, 32}) {
        vector<vector<int>> bucket(bands);
        for (int i = 0; i < n; ++i) {
            long long yy = (a[i].y - miny) * S / rangey;
            int b = min(bands - 1, (int)(yy * bands / (S + 1)));
            bucket[b].push_back(i); // IDs already have increasing x order.
        }
        vector<int> cyclic;
        cyclic.reserve(n);
        for (int b = 0; b < bands; ++b) {
            if (b & 1) {
                for (auto it = bucket[b].rbegin(); it != bucket[b].rend(); ++it) cyclic.push_back(*it);
            } else {
                cyclic.insert(cyclic.end(), bucket[b].begin(), bucket[b].end());
            }
        }
        int at0 = (int)(find(cyclic.begin(), cyclic.end(), 0) - cyclic.begin());
        vector<int> candidate;
        candidate.reserve(n + 1);
        candidate.push_back(0);
        for (int k = 1; k < n; ++k) candidate.push_back(cyclic[(at0 + k) % n]);
        candidate.push_back(0);
        long double c = routeCost(candidate);
        if (c < selectedCost) { selectedCost = c; path = candidate; }
        reverse(candidate.begin() + 1, candidate.end() - 1);
        c = routeCost(candidate);
        if (c < selectedCost) { selectedCost = c; path = candidate; }
    }

    // For a reversal, all non-penalized internal edges have exactly the same
    // length before and after.  Thus only its two boundary edges and internal
    // steps divisible by ten need recomputing.  This supports a substantially
    // wider local repair without turning it into an O(N*width^2) scan.
    auto improve = [&](int longest) {
        for (int l = 1; l < n; ++l) {
            int bestR = -1;
            long double bestDelta = 0;
            for (int len = 2; len <= longest && l + len - 1 < n; ++len) {
                int r = l + len - 1;
                long double oldc = edge(path[l-1], path[l], l) + edge(path[r], path[r+1], r+1);
                long double newc = edge(path[l-1], path[r], l) + edge(path[l], path[r+1], r+1);
                for (int t = ((l + 9) / 10) * 10; t <= r; t += 10) {
                    if (t <= l) continue;
                    int j = l + r - t + 1;
                    oldc += edge(path[t-1], path[t], t);
                    newc += edge(path[j], path[j-1], t);
                }
                if (newc - oldc < bestDelta - 1e-7L) {
                    bestDelta = newc - oldc;
                    bestR = r;
                }
            }
            if (bestR != -1) reverse(path.begin() + l, path.begin() + bestR + 1);
        }
    };
    improve(16);
    improve(5);

    // Reversals retain the order of neither side of a cut.  A different
    // short-range repair is to remove one city and put it a few positions
    // away.  This fixes a misplaced point without undoing an otherwise good
    // local chain.  Score every affected step explicitly: this is important
    // because the shifted (geometrically identical) edges can move onto or
    // off a carrot step.
    auto relocate = [&](int radius) {
        for (int i = 1; i < n; ++i) {
            int chosen = -1;
            long double chosenDelta = 0;
            for (int k = max(1, i - radius); k <= min(n - 1, i + radius); ++k) {
                if (k == i) continue;
                int lo = min(i, k), hi = max(i, k);
                auto after = [&](int pos) -> int {
                    if (i < k) {
                        if (pos < i || pos > k) return path[pos];
                        if (pos == k) return path[i];
                        return path[pos + 1];
                    }
                    if (pos < k || pos > i) return path[pos];
                    if (pos == k) return path[i];
                    return path[pos - 1];
                };
                long double oldc = 0, newc = 0;
                for (int t = lo; t <= hi + 1; ++t) {
                    oldc += edge(path[t - 1], path[t], t);
                    newc += edge(after(t - 1), after(t), t);
                }
                if (newc - oldc < chosenDelta - 1e-7L) {
                    chosenDelta = newc - oldc;
                    chosen = k;
                }
            }
            if (chosen != -1) {
                // erase/insert would shift the entire suffix on every accepted
                // move, making this bounded local search quadratic.  Rotate
                // only the affected interval (whose length is <= radius + 1).
                if (i < chosen)
                    rotate(path.begin() + i, path.begin() + i + 1, path.begin() + chosen + 1);
                else
                    rotate(path.begin() + chosen, path.begin() + i, path.begin() + i + 1);
            }
        }
    };
    relocate(4);
    // A relocation can expose a newly profitable short reversal.
    improve(5);

    cout << n + 1 << '\n';
    for (int v : path) cout << v << '\n';
}
