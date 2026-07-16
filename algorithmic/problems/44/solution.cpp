#include <bits/stdc++.h>
using namespace std;

struct Point {
    long long x, y;
};

// Hilbert index on a 2^31 by 2^31 square.  Consecutive indices are spatially local.
static uint64_t hilbert(uint32_t x, uint32_t y) {
    uint64_t d = 0;
    for (uint32_t s = 1u << 30; s; s >>= 1) {
        uint32_t rx = (x & s) != 0;
        uint32_t ry = (y & s) != 0;
        d += uint64_t(s) * uint64_t(s) * ((3u * rx) ^ ry);
        if (!ry) {
            if (rx) {
                uint32_t hi = (s << 1) - 1;
                x = hi - x;
                y = hi - y;
            }
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
    vector<Point> p(n);
    long long minx = LLONG_MAX, maxx = LLONG_MIN, miny = LLONG_MAX, maxy = LLONG_MIN;
    for (auto &q : p) {
        cin >> q.x >> q.y;
        minx = min(minx, q.x); maxx = max(maxx, q.x);
        miny = min(miny, q.y); maxy = max(maxy, q.y);
    }

    vector<char> prime(n, true);
    prime[0] = prime[1] = false;
    for (int i = 2; 1LL * i * i < n; ++i)
        if (prime[i])
            for (int j = i * i; j < n; j += i) prime[j] = false;

    auto edge = [&](int a, int b) {
        double dx = (double)p[a].x - p[b].x;
        double dy = (double)p[a].y - p[b].y;
        return sqrt(dx * dx + dy * dy);
    };
    auto score = [&](const vector<int>& route) {
        double ans = 0;
        for (int t = 1; t <= n; ++t) {
            double mul = (t % 10 == 0 && !prime[route[t - 1]]) ? 1.1 : 1.0;
            ans += mul * edge(route[t - 1], route[t]);
        }
        return ans;
    };

    // Keep the supplied monotone order as a safe fallback.
    vector<int> best;
    best.reserve(n + 1);
    for (int i = 0; i < n; ++i) best.push_back(i);
    best.push_back(0);
    double bestCost = score(best);

    const uint32_t LIM = (1u << 31) - 1;
    unsigned long long span = max<unsigned long long>(static_cast<unsigned long long>(maxx - minx),
                                                        static_cast<unsigned long long>(maxy - miny));
    vector<uint32_t> ux(n), uy(n);
    for (int i = 0; i < n; ++i) {
        // One common scale preserves the geometry's aspect ratio.
        ux[i] = uint32_t(static_cast<unsigned long long>(p[i].x - minx) * LIM / span);
        uy[i] = uint32_t(static_cast<unsigned long long>(p[i].y - miny) * LIM / span);
    }

    // The eight square symmetries change where Hilbert's coarse cells are cut.
    for (int tr = 0; tr < 8; ++tr) {
        vector<pair<uint64_t, int>> order;
        order.reserve(n);
        for (int i = 0; i < n; ++i) {
            uint32_t a = ux[i], b = uy[i];
            uint32_t x, y;
            switch (tr) {
                case 0: x = a;       y = b;       break;
                case 1: x = LIM - a; y = b;       break;
                case 2: x = a;       y = LIM - b; break;
                case 3: x = LIM - a; y = LIM - b; break;
                case 4: x = b;       y = a;       break;
                case 5: x = LIM - b; y = a;       break;
                case 6: x = b;       y = LIM - a; break;
                default:x = LIM - b; y = LIM - a; break;
            }
            order.push_back({hilbert(x, y), i});
        }
        sort(order.begin(), order.end());
        int at = 0;
        while (order[at].second != 0) ++at;

        for (int backwards = 0; backwards < 2; ++backwards) {
            vector<int> route;
            route.reserve(n + 1);
            route.push_back(0);
            for (int step = 1; step < n; ++step) {
                int pos = backwards ? (at - step + n) % n : (at + step) % n;
                route.push_back(order[pos].second);
            }
            route.push_back(0);
            double cur = score(route);
            if (cur < bestCost) {
                bestCost = cur;
                best.swap(route);
            }
        }
    }

    // Discriminating local search: a normal 2-opt delta is not enough here,
    // since reversing a segment changes the city which starts every tenth edge.
    // For short segments, account for each such edge exactly.  Consequently an
    // accepted move is a strict improvement in the checker objective, not just
    // in unpenalized Euclidean length.
    auto refine2opt = [&](vector<int>& route) {
        vector<double> d(n + 1);
        auto rebuild = [&](int lo, int hi) {
            lo = max(lo, 1);
            hi = min(hi, n);
            for (int t = lo; t <= hi; ++t) d[t] = edge(route[t - 1], route[t]);
        };
        rebuild(1, n);
        constexpr int WINDOW = 24;
        for (int pass = 0; pass < 2; ++pass) {
            bool changed = false;
            for (int l = 1; l < n - 1; ++l) {
                int last = min(n - 1, l + WINDOW);
                for (int r = l + 1; r <= last; ++r) {
                    // Ignoring multipliers, all interior undirected edges are
                    // retained by a reversal; only its two boundary edges vary.
                    double ndL = edge(route[l - 1], route[r]);
                    double ndR = edge(route[l], route[r + 1]);
                    double delta = ndL + ndR - d[l] - d[r + 1];

                    // Add the complete difference for penalized positions.
                    // At an interior step t, the new edge is old edge q,
                    // q = l+r-t+1, traversed in the opposite direction.
                    int firstPenalty = ((l + 9) / 10) * 10;
                    for (int t = firstPenalty; t <= r + 1; t += 10) {
                        double oldExtra = prime[route[t - 1]] ? 0.0 : 0.1 * d[t];
                        int newSource;
                        double newDist;
                        if (t == l) {
                            newSource = route[l - 1];
                            newDist = ndL;
                        } else if (t == r + 1) {
                            newSource = route[l];
                            newDist = ndR;
                        } else {
                            int q = l + r - t + 1;
                            newSource = route[q];
                            newDist = d[q];
                        }
                        double newExtra = prime[newSource] ? 0.0 : 0.1 * newDist;
                        delta += newExtra - oldExtra;
                    }
                    if (delta < -1e-7) {
                        reverse(route.begin() + l, route.begin() + r + 1);
                        rebuild(l, r + 1);
                        changed = true;
                    }
                }
            }
            if (!changed) break;
        }
    };
    refine2opt(best);

    cout << n + 1 << '\n';
    for (int v : best) cout << v << '\n';
    return 0;
}
