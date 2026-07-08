#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <utility>
#include <vector>
using namespace std;

static const int HBITS = 21;
static const uint32_t HN = 1u << HBITS;
static const uint32_t HMASK = HN - 1u;

struct Solver {
    int N = 0;
    vector<long long> x, y;
    vector<char> prime;

    static uint64_t hilbert_rec(uint32_t x, uint32_t y, int pow, int rot) {
        if (pow == 0) return 0;
        uint32_t h = 1u << (pow - 1);
        int seg = (x < h) ? ((y < h) ? 0 : 3) : ((y < h) ? 1 : 2);
        seg = (seg + rot) & 3;
        static const int rotate_delta[4] = {3, 0, 0, 1};
        uint32_t nx = x & (h - 1);
        uint32_t ny = y & (h - 1);
        int nrot = (rot + rotate_delta[seg]) & 3;
        uint64_t sub = uint64_t(1) << (2 * pow - 2);
        uint64_t add = hilbert_rec(nx, ny, pow - 1, nrot);
        return uint64_t(seg) * sub + ((seg == 1 || seg == 2) ? add : (sub - add - 1));
    }

    static uint64_t hilbert_index(uint32_t x, uint32_t y) {
        return hilbert_rec(x, y, HBITS, 0);
    }

    static uint64_t morton_index(uint32_t x, uint32_t y) {
        uint64_t d = 0;
        for (int b = 0; b < HBITS; ++b) {
            d |= uint64_t((x >> b) & 1u) << (2 * b);
            d |= uint64_t((y >> b) & 1u) << (2 * b + 1);
        }
        return d;
    }

    double dist_id(int a, int b) const {
        return hypot((double)x[a] - (double)x[b], (double)y[a] - (double)y[b]);
    }

    int city_at(const vector<int> &route, int pos) const {
        if (pos == 0 || pos == N) return 0;
        return route[pos - 1];
    }

    double step_cost(const vector<int> &route, int t) const {
        if (t < 1 || t > N) return 0.0;
        int a = city_at(route, t - 1);
        int b = city_at(route, t);
        double m = (t % 10 == 0 && !prime[a]) ? 1.1 : 1.0;
        return m * dist_id(a, b);
    }

    double route_cost(const vector<int> &route) const {
        double total = 0.0;
        int prev = 0;
        for (int i = 0; i < (int)route.size(); ++i) {
            int t = i + 1;
            double m = (t % 10 == 0 && !prime[prev]) ? 1.1 : 1.0;
            total += m * dist_id(prev, route[i]);
            prev = route[i];
        }
        double m = (N % 10 == 0 && !prime[prev]) ? 1.1 : 1.0;
        total += m * dist_id(prev, 0);
        return total;
    }

    double swap_delta(vector<int> &route, int i, int j) const {
        if (i == j) return 0.0;
        if (i > j) swap(i, j);
        int steps[8] = {i + 1, i + 2, j + 1, j + 2, i, i + 3, j, j + 3};
        sort(begin(steps), end(steps));
        double before = 0.0;
        int last = -1;
        for (int t : steps) {
            if (t != last && 1 <= t && t <= N) before += step_cost(route, t);
            last = t;
        }
        swap(route[i], route[j]);
        double after = 0.0;
        last = -1;
        for (int t : steps) {
            if (t != last && 1 <= t && t <= N) after += step_cost(route, t);
            last = t;
        }
        swap(route[i], route[j]);
        return after - before;
    }

    double reverse_delta(const vector<int> &route, int l, int r) const {
        if (l >= r) return 0.0;
        auto after_city = [&](int pos) -> int {
            if (pos == 0 || pos == N) return 0;
            int idx = pos - 1;
            if (idx < l || idx > r) return route[idx];
            return route[l + r - idx];
        };
        auto after_step = [&](int t) -> double {
            int a = after_city(t - 1);
            int b = after_city(t);
            double m = (t % 10 == 0 && !prime[a]) ? 1.1 : 1.0;
            return m * dist_id(a, b);
        };
        auto extra_before = [&](int t) -> double {
            int a = city_at(route, t - 1);
            if (prime[a]) return 0.0;
            int b = city_at(route, t);
            return 0.1 * dist_id(a, b);
        };
        auto extra_after = [&](int t) -> double {
            int a = after_city(t - 1);
            if (prime[a]) return 0.0;
            int b = after_city(t);
            return 0.1 * dist_id(a, b);
        };

        double before = 0.0, after = 0.0;
        int t1 = l + 1;
        int t2 = r + 2;
        if (1 <= t1 && t1 <= N) {
            before += step_cost(route, t1);
            after += after_step(t1);
        }
        if (t2 != t1 && 1 <= t2 && t2 <= N) {
            before += step_cost(route, t2);
            after += after_step(t2);
        }
        int first = ((l + 2 + 9) / 10) * 10;
        for (int t = first; t <= r + 1; t += 10) {
            before += extra_before(t);
            after += extra_after(t);
        }
        return after - before;
    }

    void improve_two_opt(vector<int> &route) const {
        int M = (int)route.size();
        if (M < 3) return;
        int passes, span;
        if (M <= 350) {
            passes = 8;
            span = M;
        } else if (M <= 1200) {
            passes = 4;
            span = 140;
        } else if (M <= 6000) {
            passes = 2;
            span = 60;
        } else {
            passes = 1;
            span = 14;
        }
        for (int pass = 0; pass < passes; ++pass) {
            bool changed = false;
            for (int l = 0; l < M - 1; ++l) {
                int hi = min(M - 1, l + span);
                for (int r = l + 1; r <= hi; ++r) {
                    double d = reverse_delta(route, l, r);
                    if (d < -1e-7) {
                        reverse(route.begin() + l, route.begin() + r + 1);
                        changed = true;
                    }
                }
            }
            if (!changed) break;
        }
    }

    void improve_prime_slots(vector<int> &route) const {
        int M = (int)route.size();
        if (M <= 1) return;
        int window = (N <= 5000) ? 120 : 55;
        for (int source_pos = 9; source_pos <= N - 1; source_pos += 10) {
            int i = source_pos - 1;
            if (i < 0 || i >= M || prime[route[i]]) continue;
            double best = 0.0;
            int best_j = -1;
            int lo = max(0, i - window);
            int hi = min(M - 1, i + window);
            for (int j = lo; j <= hi; ++j) {
                if (j == i || !prime[route[j]]) continue;
                double d = swap_delta(route, i, j);
                if (d < best) {
                    best = d;
                    best_j = j;
                }
            }
            if (best_j >= 0) swap(route[i], route[best_j]);
        }
    }

    vector<vector<int>> cycle_breaks(const vector<int> &order, int limit) const {
        int M = (int)order.size();
        if (M <= 1) return {order};
        vector<pair<double, int>> top;
        top.reserve(limit);
        for (int k = 0; k < M; ++k) {
            int a = order[(k + M - 1) % M];
            int b = order[k];
            double v = dist_id(0, b) + dist_id(a, 0) - dist_id(a, b);
            if ((int)top.size() < limit || v < top.back().first) {
                top.push_back({v, k});
                sort(top.begin(), top.end());
                if ((int)top.size() > limit) top.pop_back();
            }
        }
        vector<vector<int>> routes;
        routes.reserve(top.size());
        for (auto item : top) {
            int best_k = item.second;
            vector<int> route;
            route.reserve(M);
            for (int c = 0; c < M; ++c) route.push_back(order[(best_k + c) % M]);
            routes.push_back(std::move(route));
        }
        return routes;
    }

    vector<uint32_t> scaled_coords(bool swap_axes, bool flip_x, bool flip_y) const {
        long long minx = x[0], maxx = x[0], miny = y[0], maxy = y[0];
        for (int i = 1; i < N; ++i) {
            minx = min(minx, x[i]);
            maxx = max(maxx, x[i]);
            miny = min(miny, y[i]);
            maxy = max(maxy, y[i]);
        }
        auto scale = [](long long v, long long lo, long long hi) -> uint32_t {
            if (hi == lo) return HMASK / 2;
            __int128 num = (__int128)(v - lo) * HMASK;
            return (uint32_t)(num / (hi - lo));
        };
        vector<uint32_t> out(2 * N);
        for (int i = 0; i < N; ++i) {
            uint32_t sx = scale(x[i], minx, maxx);
            uint32_t sy = scale(y[i], miny, maxy);
            if (flip_x) sx = HMASK - sx;
            if (flip_y) sy = HMASK - sy;
            if (swap_axes) swap(sx, sy);
            out[2 * i] = sx;
            out[2 * i + 1] = sy;
        }
        return out;
    }

    vector<int> key_order(const vector<uint64_t> &key) const {
        vector<int> ids;
        ids.reserve(max(0, N - 1));
        for (int i = 1; i < N; ++i) ids.push_back(i);
        stable_sort(ids.begin(), ids.end(), [&](int a, int b) {
            if (key[a] != key[b]) return key[a] < key[b];
            return a < b;
        });
        return ids;
    }

    vector<int> x_strip_snake(int block) const {
        vector<int> route;
        route.reserve(max(0, N - 1));
        vector<int> part;
        for (int l = 1, band = 0; l < N; l += block, ++band) {
            int r = min(N, l + block);
            part.clear();
            for (int i = l; i < r; ++i) part.push_back(i);
            stable_sort(part.begin(), part.end(), [&](int a, int b) {
                if (y[a] != y[b]) return band % 2 == 0 ? y[a] < y[b] : y[a] > y[b];
                return band % 2 == 0 ? a < b : a > b;
            });
            route.insert(route.end(), part.begin(), part.end());
        }
        return route;
    }

    vector<int> y_band_snake(int block) const {
        vector<int> ids;
        ids.reserve(max(0, N - 1));
        for (int i = 1; i < N; ++i) ids.push_back(i);
        stable_sort(ids.begin(), ids.end(), [&](int a, int b) {
            if (y[a] != y[b]) return y[a] < y[b];
            return a < b;
        });
        vector<int> route;
        route.reserve(ids.size());
        for (int l = 0, band = 0; l < (int)ids.size(); l += block, ++band) {
            int r = min((int)ids.size(), l + block);
            if (band % 2 == 0) {
                stable_sort(ids.begin() + l, ids.begin() + r);
            } else {
                stable_sort(ids.begin() + l, ids.begin() + r, greater<int>());
            }
            route.insert(route.end(), ids.begin() + l, ids.begin() + r);
        }
        return route;
    }

    vector<int> angle_order(bool far_first) const {
        vector<int> ids;
        ids.reserve(max(0, N - 1));
        vector<double> ang(N), rad(N);
        for (int i = 1; i < N; ++i) {
            ids.push_back(i);
            double dx = (double)x[i] - (double)x[0];
            double dy = (double)y[i] - (double)y[0];
            ang[i] = atan2(dy, dx);
            rad[i] = dx * dx + dy * dy;
        }
        stable_sort(ids.begin(), ids.end(), [&](int a, int b) {
            if (ang[a] != ang[b]) return ang[a] < ang[b];
            if (rad[a] != rad[b]) return far_first ? rad[a] > rad[b] : rad[a] < rad[b];
            return a < b;
        });
        return ids;
    }

    vector<int> angle_radius_snake(int block) const {
        vector<int> ids;
        ids.reserve(max(0, N - 1));
        vector<double> ang(N), rad(N);
        for (int i = 1; i < N; ++i) {
            ids.push_back(i);
            double dx = (double)x[i] - (double)x[0];
            double dy = (double)y[i] - (double)y[0];
            ang[i] = atan2(dy, dx);
            rad[i] = dx * dx + dy * dy;
        }
        stable_sort(ids.begin(), ids.end(), [&](int a, int b) {
            if (ang[a] != ang[b]) return ang[a] < ang[b];
            return rad[a] < rad[b];
        });
        vector<int> route;
        route.reserve(ids.size());
        for (int l = 0, band = 0; l < (int)ids.size(); l += block, ++band) {
            int r = min((int)ids.size(), l + block);
            stable_sort(ids.begin() + l, ids.begin() + r, [&](int a, int b) {
                if (rad[a] != rad[b]) return band % 2 == 0 ? rad[a] < rad[b] : rad[a] > rad[b];
                return band % 2 == 0 ? a < b : a > b;
            });
            route.insert(route.end(), ids.begin() + l, ids.begin() + r);
        }
        return route;
    }

    vector<int> modulo_lane_snake(int lanes) const {
        vector<int> route;
        route.reserve(max(0, N - 1));
        vector<int> part;
        for (int lane = 0; lane < lanes; ++lane) {
            part.clear();
            for (int i = 1; i < N; ++i) {
                if (i % lanes == lane) part.push_back(i);
            }
            if (lane % 2 == 1) reverse(part.begin(), part.end());
            route.insert(route.end(), part.begin(), part.end());
        }
        return route;
    }

    vector<int> nearest_neighbor_small() const {
        int M = N - 1;
        vector<int> route;
        route.reserve(M);
        vector<char> used(N, 0);
        used[0] = 1;
        int cur = 0;
        for (int step = 0; step < M; ++step) {
            int best = -1;
            long double bd = 0;
            for (int i = 1; i < N; ++i) {
                if (used[i]) continue;
                long double dx = (long double)x[cur] - x[i];
                long double dy = (long double)y[cur] - y[i];
                long double d = dx * dx + dy * dy;
                if (best < 0 || d < bd) {
                    bd = d;
                    best = i;
                }
            }
            used[best] = 1;
            route.push_back(best);
            cur = best;
        }
        return route;
    }

    vector<int> exact_small() const {
        int M = N - 1;
        if (M <= 0) return {};
        if (M > 17) return {};
        int S = 1 << M;
        const double INF = numeric_limits<double>::infinity();
        vector<double> dp((size_t)S * M, INF);
        vector<int16_t> par((size_t)S * M, -1);
        auto at = [&](int mask, int last) -> size_t { return (size_t)mask * M + last; };

        for (int i = 0; i < M; ++i) dp[at(1 << i, i)] = dist_id(0, i + 1);
        for (int mask = 1; mask < S; ++mask) {
            int k = __builtin_popcount((unsigned)mask);
            int t = k + 1;
            for (int last = 0; last < M; ++last) {
                double cur = dp[at(mask, last)];
                if (!isfinite(cur)) continue;
                int a = last + 1;
                double mult = (t % 10 == 0 && !prime[a]) ? 1.1 : 1.0;
                int rem = (S - 1) ^ mask;
                while (rem) {
                    int bit = rem & -rem;
                    int nxt = __builtin_ctz((unsigned)bit);
                    int nmask = mask | bit;
                    double cand = cur + mult * dist_id(a, nxt + 1);
                    size_t idx = at(nmask, nxt);
                    if (cand < dp[idx]) {
                        dp[idx] = cand;
                        par[idx] = (int16_t)last;
                    }
                    rem -= bit;
                }
            }
        }

        int full = S - 1;
        double best = INF;
        int last_best = 0;
        for (int last = 0; last < M; ++last) {
            int a = last + 1;
            double mult = (N % 10 == 0 && !prime[a]) ? 1.1 : 1.0;
            double cand = dp[at(full, last)] + mult * dist_id(a, 0);
            if (cand < best) {
                best = cand;
                last_best = last;
            }
        }

        vector<int> route(M);
        int mask = full;
        int last = last_best;
        for (int pos = M - 1; pos >= 0; --pos) {
            route[pos] = last + 1;
            int prev = par[at(mask, last)];
            mask ^= 1 << last;
            last = prev;
        }
        return route;
    }

    vector<int> solve() {
        prime.assign(max(2, N), true);
        prime[0] = false;
        prime[1] = false;
        for (int p = 2; 1LL * p * p < N; ++p) {
            if (prime[p]) {
                for (long long q = 1LL * p * p; q < N; q += p) prime[(int)q] = false;
            }
        }

        vector<int> best;
        double best_cost = numeric_limits<double>::infinity();
        auto consider = [&](vector<int> route) {
            if ((int)route.size() != N - 1) return;
            double c = route_cost(route);
            if (c < best_cost) {
                best_cost = c;
                best = std::move(route);
            }
        };
        auto consider_path_variants = [&](vector<int> route, int breaks) {
            consider(route);
            vector<int> rev = route;
            reverse(rev.begin(), rev.end());
            consider(rev);
            if (breaks > 0) {
                auto routes = cycle_breaks(route, breaks);
                for (auto &r : routes) consider(std::move(r));
                reverse(route.begin(), route.end());
                routes = cycle_breaks(route, breaks);
                for (auto &r : routes) consider(std::move(r));
            }
        };

        vector<int> identity;
        identity.reserve(max(0, N - 1));
        for (int i = 1; i < N; ++i) identity.push_back(i);
        consider(identity);
        reverse(identity.begin(), identity.end());
        consider(identity);

        int root = max(8, (int)(sqrt((double)max(1, N - 1)) + 0.5));
        vector<int> blocks = {root / 2, root, root * 2, 128, 384, 1024};
        sort(blocks.begin(), blocks.end());
        blocks.erase(unique(blocks.begin(), blocks.end()), blocks.end());
        for (int b : blocks) {
            if (b <= 1) continue;
            consider_path_variants(x_strip_snake(b), 2);
            consider_path_variants(y_band_snake(b), 2);
        }
        consider_path_variants(angle_order(false), 2);
        consider_path_variants(angle_order(true), 2);
        for (int b : blocks) {
            if (b > 1) consider_path_variants(angle_radius_snake(b), 2);
        }
        for (int lanes = 2; lanes <= 10; ++lanes) {
            consider_path_variants(modulo_lane_snake(lanes), 2);
        }

        const bool configs[8][3] = {
            {false, false, false},
            {false, true, false},
            {false, false, true},
            {false, true, true},
            {true, false, false},
            {true, true, false},
            {true, false, true},
            {true, true, true}
        };
        for (int ci = 0; ci < 8; ++ci) {
            const bool *cfg = configs[ci];
            vector<uint32_t> sc = scaled_coords(cfg[0], cfg[1], cfg[2]);
            vector<uint64_t> key(N);
            for (int i = 0; i < N; ++i) key[i] = hilbert_index(sc[2 * i], sc[2 * i + 1]);
            vector<int> order = key_order(key);
            auto routes = cycle_breaks(order, 3);
            for (auto &r : routes) consider(std::move(r));
            reverse(order.begin(), order.end());
            routes = cycle_breaks(order, 3);
            for (auto &r : routes) consider(std::move(r));
            if (ci < 4) {
                for (int i = 0; i < N; ++i) key[i] = morton_index(sc[2 * i], sc[2 * i + 1]);
                order = key_order(key);
                routes = cycle_breaks(order, 2);
                for (auto &r : routes) consider(std::move(r));
                reverse(order.begin(), order.end());
                routes = cycle_breaks(order, 2);
                for (auto &r : routes) consider(std::move(r));
            }
        }

        if (N <= 5000) consider(nearest_neighbor_small());
        if (N <= 18) consider(exact_small());

        improve_two_opt(best);
        improve_prime_slots(best);
        if (N <= 350) improve_two_opt(best);
        return best;
    }
};

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    Solver s;
    if (!(cin >> s.N)) return 0;
    s.x.resize(s.N);
    s.y.resize(s.N);
    for (int i = 0; i < s.N; ++i) cin >> s.x[i] >> s.y[i];

    vector<int> route = s.solve();
    cout << s.N + 1 << '\n';
    cout << 0 << '\n';
    for (int v : route) cout << v << '\n';
    cout << 0 << '\n';
    return 0;
}
