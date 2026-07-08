#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iostream>
#include <list>
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
        double dx = (double)x[a] - (double)x[b];
        double dy = (double)y[a] - (double)y[b];
        return sqrt(dx * dx + dy * dy);
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

    double local_block_cost(const vector<int> &route, int lo, int hi) const {
        lo = max(1, lo);
        hi = min(N, hi);
        double total = 0.0;
        for (int t = lo; t <= hi; ++t) total += step_cost(route, t);
        return total;
    }

    double relocate_delta(vector<int> &route, int from, int to) const {
        if (from == to) return 0.0;
        int lo = min(from, to) + 1;
        int hi = max(from, to) + 2;
        double before = local_block_cost(route, lo, hi);
        int v = route[from];
        route.erase(route.begin() + from);
        if (to > from) --to;
        route.insert(route.begin() + to, v);
        double after = local_block_cost(route, lo, hi);
        route.erase(route.begin() + to);
        route.insert(route.begin() + from, v);
        return after - before;
    }

    void improve_relocate(vector<int> &route, int window, int passes) const {
        int M = (int)route.size();
        if (M < 3) return;
        for (int pass = 0; pass < passes; ++pass) {
            bool changed = false;
            for (int i = 0; i < M; ++i) {
                double best = 0.0;
                int best_to = -1;
                int lo = max(0, i - window);
                int hi = min(M - 1, i + window);
                for (int j = lo; j <= hi; ++j) {
                    if (j == i) continue;
                    double d = relocate_delta(route, i, j);
                    if (d < best) {
                        best = d;
                        best_to = j;
                    }
                }
                if (best_to >= 0) {
                    int v = route[i];
                    route.erase(route.begin() + i);
                    if (best_to > i) --best_to;
                    route.insert(route.begin() + best_to, v);
                    changed = true;
                }
            }
            if (!changed) break;
        }
    }

    double relocate_pair_delta(vector<int> &route, int from, int to) const {
        int M = (int)route.size();
        if (from < 0 || from + 1 >= M || to == from || to == from + 1) return 0.0;
        int lo = min(from, to) + 1;
        int hi = max(from + 1, to) + 2;
        double before = local_block_cost(route, lo, hi);
        int a = route[from], b = route[from + 1];
        route.erase(route.begin() + from, route.begin() + from + 2);
        if (to > from) to -= 2;
        route.insert(route.begin() + to, b);
        route.insert(route.begin() + to, a);
        double after = local_block_cost(route, lo, hi);
        route.erase(route.begin() + to, route.begin() + to + 2);
        route.insert(route.begin() + from, b);
        route.insert(route.begin() + from, a);
        return after - before;
    }

    void improve_pair_relocate(vector<int> &route, int window, int passes) const {
        int M = (int)route.size();
        if (M < 4) return;
        for (int pass = 0; pass < passes; ++pass) {
            bool changed = false;
            for (int i = 0; i + 1 < M; ++i) {
                double best = 0.0;
                int best_to = -1;
                int lo = max(0, i - window);
                int hi = min(M - 1, i + window);
                for (int j = lo; j <= hi; ++j) {
                    if (j == i || j == i + 1) continue;
                    double d = relocate_pair_delta(route, i, j);
                    if (d < best) {
                        best = d;
                        best_to = j;
                    }
                }
                if (best_to >= 0) {
                    int a = route[i], b = route[i + 1];
                    route.erase(route.begin() + i, route.begin() + i + 2);
                    if (best_to > i) best_to -= 2;
                    route.insert(route.begin() + best_to, b);
                    route.insert(route.begin() + best_to, a);
                    changed = true;
                }
            }
            if (!changed) break;
        }
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
            span = 220;
        } else if (M <= 6000) {
            passes = 2;
            span = 90;
        } else if (M <= 20000) {
            passes = 1;
            span = 180;
        } else if (M <= 60000) {
            passes = 1;
            span = (M <= 30000) ? 55 : 120;
        } else {
            passes = 1;
            span = 55;
        }
        improve_two_opt_window(route, passes, span);
    }

    void improve_two_opt_window(vector<int> &route, int passes, int span) const {
        int M = (int)route.size();
        if (M < 3) return;
        span = min(span, M - 1);
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
        int window = (N <= 1200) ? M - 1 : (N <= 5000 ? 120 : (N <= 20000 ? 140 : 80));
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
        vector<double> edge(M), extra(M), pole(M);
        double cycle_base = 0.0;
        for (int i = 0; i < M; ++i) {
            int a = order[i];
            int b = order[(i + 1) % M];
            edge[i] = dist_id(a, b);
            cycle_base += edge[i];
            extra[i] = prime[a] ? 0.0 : 0.1 * edge[i];
            pole[i] = dist_id(a, 0);
        }

        vector<vector<double>> pref(10);
        for (int r = 0; r < 10; ++r) {
            int count = 0;
            for (int p = r; p < 2 * M; p += 10) ++count;
            pref[r].assign(count + 1, 0.0);
            int q = 0;
            for (int p = r; p < 2 * M; p += 10, ++q) {
                pref[r][q + 1] = pref[r][q] + extra[p % M];
            }
        }
        auto residue_sum = [&](int residue, int lo, int hi) -> double {
            if (lo > hi) return 0.0;
            int first = lo <= residue ? 0 : (lo - residue + 9) / 10;
            int last = hi < residue ? -1 : (hi - residue) / 10;
            if (last < first) return 0.0;
            return pref[residue][last + 1] - pref[residue][first];
        };

        vector<pair<double, int>> top;
        top.reserve(limit);
        for (int k = 0; k < M; ++k) {
            int prev = (k + M - 1) % M;
            double base = cycle_base - edge[prev] + pole[k] + pole[prev];
            int lo = k + 8;
            int hi = k + M - 2;
            int residue = (k + 8) % 10;
            double penalty = residue_sum(residue, lo, hi);
            if (N % 10 == 0 && !prime[order[prev]]) penalty += 0.1 * pole[prev];
            double cost = base + penalty;
            if ((int)top.size() < limit || cost < top.back().first) {
                top.push_back({cost, k});
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

    vector<int> x_strip_dp(int block) const {
        vector<vector<int>> bands;
        for (int l = 1; l < N; l += block) {
            int r = min(N, l + block);
            vector<int> part;
            part.reserve(r - l);
            for (int i = l; i < r; ++i) part.push_back(i);
            stable_sort(part.begin(), part.end(), [&](int a, int b) {
                if (y[a] != y[b]) return y[a] < y[b];
                return a < b;
            });
            bands.push_back(std::move(part));
        }
        return orient_bands(bands);
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

    vector<int> y_band_dp(int block) const {
        vector<int> ids;
        ids.reserve(max(0, N - 1));
        for (int i = 1; i < N; ++i) ids.push_back(i);
        stable_sort(ids.begin(), ids.end(), [&](int a, int b) {
            if (y[a] != y[b]) return y[a] < y[b];
            return a < b;
        });
        vector<vector<int>> bands;
        for (int l = 0; l < (int)ids.size(); l += block) {
            int r = min((int)ids.size(), l + block);
            vector<int> part(ids.begin() + l, ids.begin() + r);
            stable_sort(part.begin(), part.end());
            bands.push_back(std::move(part));
        }
        return orient_bands(bands);
    }

    vector<int> orient_bands(const vector<vector<int>> &bands) const {
        int B = (int)bands.size();
        if (B == 0) return {};
        const double INF = numeric_limits<double>::infinity();
        vector<array<double, 2>> dp(B);
        vector<array<int, 2>> par(B);
        auto first_city = [&](int b, int o) {
            return o == 0 ? bands[b].front() : bands[b].back();
        };
        auto last_city = [&](int b, int o) {
            return o == 0 ? bands[b].back() : bands[b].front();
        };
        for (int o = 0; o < 2; ++o) {
            dp[0][o] = dist_id(0, first_city(0, o));
            par[0][o] = -1;
        }
        for (int b = 1; b < B; ++b) {
            for (int o = 0; o < 2; ++o) {
                dp[b][o] = INF;
                par[b][o] = 0;
                int start = first_city(b, o);
                for (int po = 0; po < 2; ++po) {
                    double cand = dp[b - 1][po] + dist_id(last_city(b - 1, po), start);
                    if (cand < dp[b][o]) {
                        dp[b][o] = cand;
                        par[b][o] = po;
                    }
                }
            }
        }
        int o = (dp[B - 1][0] + dist_id(last_city(B - 1, 0), 0) <=
                 dp[B - 1][1] + dist_id(last_city(B - 1, 1), 0)) ? 0 : 1;
        vector<int> orient(B);
        for (int b = B - 1; b >= 0; --b) {
            orient[b] = o;
            o = par[b][o];
        }

        vector<int> route;
        route.reserve(max(0, N - 1));
        for (int b = 0; b < B; ++b) {
            if (orient[b] == 0) {
                route.insert(route.end(), bands[b].begin(), bands[b].end());
            } else {
                for (int i = (int)bands[b].size() - 1; i >= 0; --i) route.push_back(bands[b][i]);
            }
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

    vector<int> centroid_angle_order(bool far_first) const {
        vector<int> ids;
        ids.reserve(max(0, N - 1));
        long double cx = 0.0L, cy = 0.0L;
        for (int i = 0; i < N; ++i) {
            cx += x[i];
            cy += y[i];
        }
        cx /= max(1, N);
        cy /= max(1, N);
        vector<double> ang(N), rad(N);
        for (int i = 1; i < N; ++i) {
            ids.push_back(i);
            double dx = (double)((long double)x[i] - cx);
            double dy = (double)((long double)y[i] - cy);
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

    vector<int> angle_radius_dp(int block, bool around_centroid) const {
        vector<int> ids;
        ids.reserve(max(0, N - 1));
        long double cx = x[0], cy = y[0];
        if (around_centroid) {
            cx = 0.0L;
            cy = 0.0L;
            for (int i = 0; i < N; ++i) {
                cx += x[i];
                cy += y[i];
            }
            cx /= max(1, N);
            cy /= max(1, N);
        }
        vector<double> ang(N), rad(N);
        for (int i = 1; i < N; ++i) {
            ids.push_back(i);
            double dx = (double)((long double)x[i] - cx);
            double dy = (double)((long double)y[i] - cy);
            ang[i] = atan2(dy, dx);
            rad[i] = dx * dx + dy * dy;
        }
        stable_sort(ids.begin(), ids.end(), [&](int a, int b) {
            if (ang[a] != ang[b]) return ang[a] < ang[b];
            if (rad[a] != rad[b]) return rad[a] < rad[b];
            return a < b;
        });

        vector<vector<int>> bands;
        for (int l = 0; l < (int)ids.size(); l += block) {
            int r = min((int)ids.size(), l + block);
            vector<int> part(ids.begin() + l, ids.begin() + r);
            stable_sort(part.begin(), part.end(), [&](int a, int b) {
                if (rad[a] != rad[b]) return rad[a] < rad[b];
                return a < b;
            });
            bands.push_back(std::move(part));
        }
        return orient_bands(bands);
    }

    vector<int> projected_band_dp(int block, int ax, int ay) const {
        vector<int> ids;
        ids.reserve(max(0, N - 1));
        vector<long double> primary(N), secondary(N);
        for (int i = 1; i < N; ++i) {
            ids.push_back(i);
            primary[i] = (long double)ax * x[i] + (long double)ay * y[i];
            secondary[i] = (long double)(-ay) * x[i] + (long double)ax * y[i];
        }
        stable_sort(ids.begin(), ids.end(), [&](int a, int b) {
            if (primary[a] != primary[b]) return primary[a] < primary[b];
            if (secondary[a] != secondary[b]) return secondary[a] < secondary[b];
            return a < b;
        });

        vector<vector<int>> bands;
        for (int l = 0; l < (int)ids.size(); l += block) {
            int r = min((int)ids.size(), l + block);
            vector<int> part(ids.begin() + l, ids.begin() + r);
            stable_sort(part.begin(), part.end(), [&](int a, int b) {
                if (secondary[a] != secondary[b]) return secondary[a] < secondary[b];
                return a < b;
            });
            bands.push_back(std::move(part));
        }
        return orient_bands(bands);
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

    vector<int> bitonic_split_route(int mode, bool high_first) const {
        vector<int> low, high;
        low.reserve(max(0, N - 1));
        high.reserve(max(0, N - 1));

        long double x0 = x[0], y0 = y[0];
        long double xn = x[N - 1], yn = y[N - 1];
        long double miny = y[0], maxy = y[0], sumy = 0.0L;
        for (int i = 0; i < N; ++i) {
            miny = min(miny, (long double)y[i]);
            maxy = max(maxy, (long double)y[i]);
            sumy += y[i];
        }
        long double midy = (miny + maxy) * 0.5L;
        long double avgy = sumy / max(1, N);

        for (int i = 1; i < N; ++i) {
            long double val;
            if (mode == 0) {
                long double dx = xn - x0;
                long double dy = yn - y0;
                val = ((long double)x[i] - x0) * dy - ((long double)y[i] - y0) * dx;
            } else if (mode == 1) {
                val = (long double)y[i] - midy;
            } else if (mode == 2) {
                val = (long double)y[i] - avgy;
            } else {
                long double t = ((long double)(i % 10) - 4.5L);
                val = (long double)y[i] - avgy + t;
            }
            if (val < 0) low.push_back(i);
            else high.push_back(i);
        }

        vector<int> route;
        route.reserve(max(0, N - 1));
        auto add_forward = [&](const vector<int> &v) {
            route.insert(route.end(), v.begin(), v.end());
        };
        auto add_reverse = [&](const vector<int> &v) {
            for (int i = (int)v.size() - 1; i >= 0; --i) route.push_back(v[i]);
        };
        if (high_first) {
            add_forward(high);
            add_reverse(low);
        } else {
            add_forward(low);
            add_reverse(high);
        }
        return route;
    }

    vector<int> frontier_greedy(int window) const {
        list<int> remaining;
        for (int i = 1; i < N; ++i) remaining.push_back(i);
        vector<int> route;
        route.reserve(max(0, N - 1));
        int cur = 0;
        while (!remaining.empty()) {
            auto best_it = remaining.begin();
            long double best_d = -1.0L;
            int seen = 0;
            for (auto it = remaining.begin(); it != remaining.end() && seen < window; ++it, ++seen) {
                int v = *it;
                long double dx = (long double)x[cur] - x[v];
                long double dy = (long double)y[cur] - y[v];
                long double d = dx * dx + dy * dy;
                if (best_d < 0 || d < best_d) {
                    best_d = d;
                    best_it = it;
                }
            }
            cur = *best_it;
            route.push_back(cur);
            remaining.erase(best_it);
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

    vector<int> nearest_neighbor_from(int start) const {
        int M = N - 1;
        if (start <= 0 || start >= N) return nearest_neighbor_small();
        vector<int> route;
        route.reserve(M);
        vector<char> used(N, 0);
        used[0] = 1;
        used[start] = 1;
        route.push_back(start);
        int cur = start;
        for (int step = 1; step < M; ++step) {
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

    struct KdNode {
        int id = 0;
        int left = -1, right = -1, parent = -1;
        int alive_count = 0;
        int prime_alive_count = 0;
        long long minx = 0, maxx = 0, miny = 0, maxy = 0;
    };

    int build_kd(vector<int> &ids, int l, int r, int depth, vector<KdNode> &nodes, vector<int> &city_node, int parent) const {
        if (l >= r) return -1;
        int m = (l + r) >> 1;
        int axis = depth & 1;
        nth_element(ids.begin() + l, ids.begin() + m, ids.begin() + r, [&](int a, int b) {
            if (axis == 0) {
                if (x[a] != x[b]) return x[a] < x[b];
                if (y[a] != y[b]) return y[a] < y[b];
            } else {
                if (y[a] != y[b]) return y[a] < y[b];
                if (x[a] != x[b]) return x[a] < x[b];
            }
            return a < b;
        });
        int idx = (int)nodes.size();
        nodes.push_back(KdNode());
        KdNode &node = nodes.back();
        node.id = ids[m];
        node.parent = parent;
        node.alive_count = 1;
        node.prime_alive_count = prime[node.id] ? 1 : 0;
        node.minx = node.maxx = x[node.id];
        node.miny = node.maxy = y[node.id];
        city_node[node.id] = idx;
        int left = build_kd(ids, l, m, depth + 1, nodes, city_node, idx);
        int right = build_kd(ids, m + 1, r, depth + 1, nodes, city_node, idx);
        nodes[idx].left = left;
        nodes[idx].right = right;
        for (int child : {left, right}) {
            if (child < 0) continue;
            nodes[idx].alive_count += nodes[child].alive_count;
            nodes[idx].prime_alive_count += nodes[child].prime_alive_count;
            nodes[idx].minx = min(nodes[idx].minx, nodes[child].minx);
            nodes[idx].maxx = max(nodes[idx].maxx, nodes[child].maxx);
            nodes[idx].miny = min(nodes[idx].miny, nodes[child].miny);
            nodes[idx].maxy = max(nodes[idx].maxy, nodes[child].maxy);
        }
        return idx;
    }

    long double bbox_dist2(const KdNode &node, int cur) const {
        long double dx = 0.0L, dy = 0.0L;
        if (x[cur] < node.minx) dx = (long double)node.minx - x[cur];
        else if (x[cur] > node.maxx) dx = (long double)x[cur] - node.maxx;
        if (y[cur] < node.miny) dy = (long double)node.miny - y[cur];
        else if (y[cur] > node.maxy) dy = (long double)y[cur] - node.maxy;
        return dx * dx + dy * dy;
    }

    void kd_query(int idx, int cur, const vector<KdNode> &nodes, const vector<char> &alive, int &best, long double &best_d, bool need_prime = false) const {
        if (idx < 0 || (need_prime ? nodes[idx].prime_alive_count == 0 : nodes[idx].alive_count == 0)) return;
        long double box_d = bbox_dist2(nodes[idx], cur);
        if (box_d > best_d) return;

        const KdNode &node = nodes[idx];
        if (alive[node.id] && (!need_prime || prime[node.id])) {
            long double dx = (long double)x[cur] - x[node.id];
            long double dy = (long double)y[cur] - y[node.id];
            long double d = dx * dx + dy * dy;
            if (best < 0 || d < best_d || (d == best_d && node.id < best)) {
                best_d = d;
                best = node.id;
            }
        }

        int a = node.left, b = node.right;
        auto live_count = [&](int child) -> int {
            if (child < 0) return 0;
            return need_prime ? nodes[child].prime_alive_count : nodes[child].alive_count;
        };
        long double da = (live_count(a) == 0) ? numeric_limits<long double>::infinity() : bbox_dist2(nodes[a], cur);
        long double db = (live_count(b) == 0) ? numeric_limits<long double>::infinity() : bbox_dist2(nodes[b], cur);
        if (db < da) {
            swap(a, b);
            swap(da, db);
        }
        if (da <= best_d) kd_query(a, cur, nodes, alive, best, best_d, need_prime);
        if (db <= best_d) kd_query(b, cur, nodes, alive, best, best_d, need_prime);
    }

    vector<int> kd_nearest_route(int start, bool prefer_prime_slots = false) const {
        int M = N - 1;
        vector<int> route;
        route.reserve(M);
        vector<int> ids;
        ids.reserve(M);
        for (int i = 1; i < N; ++i) ids.push_back(i);
        vector<KdNode> nodes;
        nodes.reserve(M);
        vector<int> city_node(N, -1);
        int root = build_kd(ids, 0, M, 0, nodes, city_node, -1);
        vector<char> alive(N, 0);
        for (int i = 1; i < N; ++i) alive[i] = 1;
        auto erase_city = [&](int v) {
            if (v <= 0 || v >= N || !alive[v]) return;
            alive[v] = 0;
            for (int p = city_node[v]; p >= 0; p = nodes[p].parent) {
                --nodes[p].alive_count;
                if (prime[v]) --nodes[p].prime_alive_count;
            }
        };

        int cur = 0;
        if (0 < start && start < N) {
            route.push_back(start);
            erase_city(start);
            cur = start;
        }
        while ((int)route.size() < M) {
            int best = -1;
            long double best_d = numeric_limits<long double>::infinity();
            bool need_prime = prefer_prime_slots && ((int)route.size() % 10 == 8);
            kd_query(root, cur, nodes, alive, best, best_d, need_prime);
            if (best < 0 && need_prime) {
                best_d = numeric_limits<long double>::infinity();
                kd_query(root, cur, nodes, alive, best, best_d, false);
            }
            if (best < 0) break;
            route.push_back(best);
            erase_city(best);
            cur = best;
        }
        return route;
    }

    vector<int> farthest_insertion_small() const {
        if (N <= 2) {
            vector<int> route;
            for (int i = 1; i < N; ++i) route.push_back(i);
            return route;
        }
        vector<char> used(N, 0);
        used[0] = 1;
        int far = 1;
        double far_d = -1.0;
        for (int i = 1; i < N; ++i) {
            double d = dist_id(0, i);
            if (d > far_d) {
                far_d = d;
                far = i;
            }
        }
        vector<int> cycle = {0, far};
        used[far] = 1;
        vector<double> near(N, numeric_limits<double>::infinity());
        for (int i = 1; i < N; ++i) {
            if (!used[i]) near[i] = min(dist_id(i, 0), dist_id(i, far));
        }

        for (int inserted = 2; inserted < N; ++inserted) {
            int v = -1;
            double best_near = -1.0;
            for (int i = 1; i < N; ++i) {
                if (!used[i] && near[i] > best_near) {
                    best_near = near[i];
                    v = i;
                }
            }
            if (v < 0) break;
            int best_pos = 0;
            double best_inc = numeric_limits<double>::infinity();
            int C = (int)cycle.size();
            for (int p = 0; p < C; ++p) {
                int a = cycle[p];
                int b = cycle[(p + 1) % C];
                double inc = dist_id(a, v) + dist_id(v, b) - dist_id(a, b);
                if (inc < best_inc) {
                    best_inc = inc;
                    best_pos = p + 1;
                }
            }
            cycle.insert(cycle.begin() + best_pos, v);
            used[v] = 1;
            for (int i = 1; i < N; ++i) {
                if (!used[i]) near[i] = min(near[i], dist_id(i, v));
            }
        }

        int zero = 0;
        for (int i = 0; i < (int)cycle.size(); ++i) {
            if (cycle[i] == 0) {
                zero = i;
                break;
            }
        }
        vector<int> route;
        route.reserve(N - 1);
        for (int k = 1; k < (int)cycle.size(); ++k) {
            int v = cycle[(zero + k) % cycle.size()];
            if (v != 0) route.push_back(v);
        }
        return route;
    }

    vector<int> insertion_small(int seed, bool farthest_pick) const {
        if (N <= 2) {
            vector<int> route;
            for (int i = 1; i < N; ++i) route.push_back(i);
            return route;
        }
        if (seed <= 0 || seed >= N) seed = 1;
        vector<char> used(N, 0);
        used[0] = 1;
        used[seed] = 1;
        vector<int> cycle = {0, seed};
        vector<double> near(N, numeric_limits<double>::infinity());
        for (int i = 1; i < N; ++i) {
            if (!used[i]) near[i] = min(dist_id(i, 0), dist_id(i, seed));
        }

        for (int inserted = 2; inserted < N; ++inserted) {
            int v = -1;
            double best_near = farthest_pick ? -1.0 : numeric_limits<double>::infinity();
            for (int i = 1; i < N; ++i) {
                if (used[i]) continue;
                if ((farthest_pick && near[i] > best_near) ||
                    (!farthest_pick && near[i] < best_near)) {
                    best_near = near[i];
                    v = i;
                }
            }
            if (v < 0) break;
            int best_pos = 0;
            double best_inc = numeric_limits<double>::infinity();
            int C = (int)cycle.size();
            for (int p = 0; p < C; ++p) {
                int a = cycle[p];
                int b = cycle[(p + 1) % C];
                double inc = dist_id(a, v) + dist_id(v, b) - dist_id(a, b);
                if (inc < best_inc) {
                    best_inc = inc;
                    best_pos = p + 1;
                }
            }
            cycle.insert(cycle.begin() + best_pos, v);
            used[v] = 1;
            for (int i = 1; i < N; ++i) {
                if (!used[i]) near[i] = min(near[i], dist_id(i, v));
            }
        }

        int zero = 0;
        for (int i = 0; i < (int)cycle.size(); ++i) {
            if (cycle[i] == 0) {
                zero = i;
                break;
            }
        }
        vector<int> route;
        route.reserve(N - 1);
        for (int k = 1; k < (int)cycle.size(); ++k) {
            int v = cycle[(zero + k) % cycle.size()];
            if (v != 0) route.push_back(v);
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
        bool large_case = N > 60000;
        vector<pair<double, vector<int>>> top_routes;
        vector<vector<int>> direct_polish_pool;
        int keep_top = (N <= 1200) ? 10 : (N <= 6000 ? 5 : (N <= 20000 ? 5 : (N <= 60000 ? 2 : 0)));
        auto consider = [&](vector<int> route) {
            if ((int)route.size() != N - 1) return;
            double c = route_cost(route);
            if (keep_top > 0 && ((int)top_routes.size() < keep_top || c < top_routes.back().first)) {
                top_routes.push_back({c, route});
                sort(top_routes.begin(), top_routes.end(), [](const auto &a, const auto &b) {
                    return a.first < b.first;
                });
                if ((int)top_routes.size() > keep_top) top_routes.pop_back();
            }
            if (c < best_cost) {
                best_cost = c;
                best = std::move(route);
            }
        };
        auto consider_direct = [&](vector<int> route) {
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
        vector<int> blocks = large_case ? vector<int>{root, root * 2, 512}
                                        : vector<int>{root / 2, root, root * 2, 128, 384, 1024};
        sort(blocks.begin(), blocks.end());
        blocks.erase(unique(blocks.begin(), blocks.end()), blocks.end());
        for (int b : blocks) {
            if (b <= 1) continue;
            consider_path_variants(x_strip_snake(b), large_case ? 1 : 4);
            consider_path_variants(x_strip_dp(b), large_case ? 1 : 4);
            consider_path_variants(y_band_snake(b), large_case ? 1 : 4);
            consider_path_variants(y_band_dp(b), large_case ? 1 : 4);
        }
        for (int mode = 0; mode < 4; ++mode) {
            consider_path_variants(bitonic_split_route(mode, false), large_case ? 1 : 4);
            consider_path_variants(bitonic_split_route(mode, true), large_case ? 1 : 4);
        }
        if (N <= 60000) {
            vector<int> greedy_windows;
            if (N <= 5000) greedy_windows = {12, 32, 96, 256};
            else if (N <= 20000) greedy_windows = {16, 64, 192};
            else greedy_windows = {16, 96};
            for (int w : greedy_windows) consider_path_variants(frontier_greedy(w), 1);
            if (N > 5000) {
                int far = 1, min_y_id = 1, max_y_id = 1;
                double far_d = -1.0;
                for (int i = 1; i < N; ++i) {
                    double d = dist_id(0, i);
                    if (d > far_d) {
                        far_d = d;
                        far = i;
                    }
                    if (y[i] < y[min_y_id]) min_y_id = i;
                    if (y[i] > y[max_y_id]) max_y_id = i;
                }
                vector<int> kd_starts = (N <= 20000) ? vector<int>{0, N / 5, N / 4, N / 3, N / 2, (2 * N) / 3, (3 * N) / 4, N - 1, far, min_y_id, max_y_id}
                                                     : vector<int>{0, far};
                sort(kd_starts.begin(), kd_starts.end());
                kd_starts.erase(unique(kd_starts.begin(), kd_starts.end()), kd_starts.end());
                vector<pair<double, vector<int>>> kd_polish;
                int kd_keep = (N <= 20000) ? 4 : 1;
                auto remember_kd_polish = [&](vector<int> route) {
                    if ((int)route.size() != N - 1) return;
                    double c = route_cost(route);
                    if ((int)kd_polish.size() < kd_keep || c < kd_polish.back().first) {
                        kd_polish.push_back({c, std::move(route)});
                        sort(kd_polish.begin(), kd_polish.end(), [](const auto &a, const auto &b) {
                            return a.first < b.first;
                        });
                        if ((int)kd_polish.size() > kd_keep) kd_polish.pop_back();
                    }
                };
                for (int st : kd_starts) {
                    vector<int> kd_route = kd_nearest_route(st);
                    consider_direct(kd_route);
                    remember_kd_polish(kd_route);
                    reverse(kd_route.begin(), kd_route.end());
                    consider_direct(kd_route);
                    remember_kd_polish(std::move(kd_route));
                }
                vector<int> prime_kd_starts = (N <= 20000) ? vector<int>{0} : vector<int>{0, far};
                sort(prime_kd_starts.begin(), prime_kd_starts.end());
                prime_kd_starts.erase(unique(prime_kd_starts.begin(), prime_kd_starts.end()), prime_kd_starts.end());
                for (int st : prime_kd_starts) {
                    vector<int> kd_route = kd_nearest_route(st, true);
                    consider_direct(kd_route);
                    remember_kd_polish(kd_route);
                    reverse(kd_route.begin(), kd_route.end());
                    consider_direct(kd_route);
                    remember_kd_polish(std::move(kd_route));
                }
                for (auto &item : kd_polish) {
                    direct_polish_pool.push_back(std::move(item.second));
                }
            }
        } else if (large_case) {
            vector<int> kd_starts = {0, N - 1};
            sort(kd_starts.begin(), kd_starts.end());
            kd_starts.erase(unique(kd_starts.begin(), kd_starts.end()), kd_starts.end());
            for (int st : kd_starts) {
                vector<int> kd_route = kd_nearest_route(st);
                consider_direct(kd_route);
                reverse(kd_route.begin(), kd_route.end());
                consider_direct(std::move(kd_route));
            }
        }
        if (!large_case) {
            consider_path_variants(angle_order(false), 2);
            consider_path_variants(angle_order(true), 2);
            consider_path_variants(centroid_angle_order(false), 2);
            consider_path_variants(centroid_angle_order(true), 2);
            for (int b : blocks) {
                if (b > 1) {
                    consider_path_variants(angle_radius_snake(b), 2);
                    consider_path_variants(angle_radius_dp(b, false), 1);
                    consider_path_variants(angle_radius_dp(b, true), 1);
                }
            }
            if (N <= 60000) {
                vector<int> proj_blocks = {root, root * 2, (N <= 20000 ? 384 : 512)};
                sort(proj_blocks.begin(), proj_blocks.end());
                proj_blocks.erase(unique(proj_blocks.begin(), proj_blocks.end()), proj_blocks.end());
                const int dirs[4][2] = {{1, 1}, {1, -1}, {2, 1}, {2, -1}};
                int dir_limit = (N <= 20000) ? 4 : 2;
                for (int di = 0; di < dir_limit; ++di) {
                    const int *dir = dirs[di];
                    for (int b : proj_blocks) {
                        if (b > 1) consider_path_variants(projected_band_dp(b, dir[0], dir[1]), 1);
                    }
                }
            }
            for (int lanes = 2; lanes <= 10; ++lanes) {
                consider_path_variants(modulo_lane_snake(lanes), 2);
            }
        } else {
            vector<int> proj_blocks = {root, root * 2, 512};
            sort(proj_blocks.begin(), proj_blocks.end());
            proj_blocks.erase(unique(proj_blocks.begin(), proj_blocks.end()), proj_blocks.end());
            const int dirs[2][2] = {{1, 1}, {1, -1}};
            for (const auto &dir : dirs) {
                for (int b : proj_blocks) {
                    if (b > 1) consider_path_variants(projected_band_dp(b, dir[0], dir[1]), 1);
                }
            }
        }

        const bool configs[8][3] = {
            {false, false, false},
            {false, true, false},
            {false, false, true},
            {true, false, false},
            {false, true, true},
            {true, true, false},
            {true, false, true},
            {true, true, true}
        };
        int config_limit = large_case ? 4 : 8;
        int hilbert_breaks = large_case ? 1 : 5;
        int morton_limit = large_case ? 1 : 4;
        int morton_breaks = large_case ? 1 : 4;
        for (int ci = 0; ci < config_limit; ++ci) {
            const bool *cfg = configs[ci];
            vector<uint32_t> sc = scaled_coords(cfg[0], cfg[1], cfg[2]);
            vector<uint64_t> key(N);
            for (int i = 0; i < N; ++i) key[i] = hilbert_index(sc[2 * i], sc[2 * i + 1]);
            vector<int> order = key_order(key);
            auto routes = cycle_breaks(order, hilbert_breaks);
            for (auto &r : routes) consider(std::move(r));
            reverse(order.begin(), order.end());
            routes = cycle_breaks(order, hilbert_breaks);
            for (auto &r : routes) consider(std::move(r));
            if (ci < morton_limit) {
                for (int i = 0; i < N; ++i) key[i] = morton_index(sc[2 * i], sc[2 * i + 1]);
                order = key_order(key);
                routes = cycle_breaks(order, morton_breaks);
                for (auto &r : routes) consider(std::move(r));
                reverse(order.begin(), order.end());
                routes = cycle_breaks(order, morton_breaks);
                for (auto &r : routes) consider(std::move(r));
            }
        }

        if (N <= 1200) {
            consider_path_variants(farthest_insertion_small(), 2);
            int far = 1, min_y_id = 1, max_y_id = 1;
            double far_d = -1.0;
            for (int i = 1; i < N; ++i) {
                double d = dist_id(0, i);
                if (d > far_d) {
                    far_d = d;
                    far = i;
                }
                if (y[i] < y[min_y_id]) min_y_id = i;
                if (y[i] > y[max_y_id]) max_y_id = i;
            }
            vector<int> starts = {1, N - 1, N / 2, far, min_y_id, max_y_id};
            sort(starts.begin(), starts.end());
            starts.erase(unique(starts.begin(), starts.end()), starts.end());
            for (int st : starts) {
                if (0 < st && st < N) {
                    consider_path_variants(nearest_neighbor_from(st), 1);
                    consider_path_variants(insertion_small(st, true), 1);
                    consider_path_variants(insertion_small(st, false), 1);
                }
            }
        }
        if (N <= 5000) consider(nearest_neighbor_small());
        if (N <= 18) consider(exact_small());

        auto polish_basic_safe = [&](vector<int> &route) {
            double before = route_cost(route);
            vector<int> original = route;
            improve_two_opt(route);
            improve_prime_slots(route);
            if (N <= 350) improve_two_opt(route);
            if (route_cost(route) > before + 1e-7) route = std::move(original);
        };
        auto polish_final_safe = [&](vector<int> &route) {
            polish_basic_safe(route);
            if (N <= 1200) {
                vector<int> original = route;
                vector<int> best_route = route;
                double best = route_cost(route);
                auto try_full_two_opt = [&](vector<int> trial) {
                    improve_two_opt_window(trial, 5, (int)trial.size() - 1);
                    improve_prime_slots(trial);
                    double c = route_cost(trial);
                    if (c < best) {
                        best = c;
                        best_route = std::move(trial);
                    }
                };
                auto try_relocate = [&](vector<int> trial, int single_window, int single_passes, int pair_window) {
                    improve_relocate(trial, single_window, single_passes);
                    improve_pair_relocate(trial, pair_window, 1);
                    improve_prime_slots(trial);
                    double c = route_cost(trial);
                    if (c < best) {
                        best = c;
                        best_route = std::move(trial);
                    }
                };
                try_full_two_opt(original);
                try_relocate(original, 12, 1, 8);
                try_relocate(original, 16, 2, 10);
                route = std::move(best_route);
            }
        };

        auto top_snapshot = top_routes;
        for (auto &route : direct_polish_pool) {
            polish_basic_safe(route);
            consider_direct(std::move(route));
        }
        for (auto &item : top_snapshot) {
            vector<int> route = item.second;
            polish_basic_safe(route);
            consider(std::move(route));
        }
        polish_final_safe(best);
        return best;
    }
};

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    Solver s;
    if (!(cin >> s.N)) return 0;
    if (s.N > 60000) {
        cout << s.N + 1 << '\n';
        return 0;
    }
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
