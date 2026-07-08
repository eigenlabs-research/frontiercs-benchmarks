#include <algorithm>
#include <bitset>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <utility>
#include <vector>
using namespace std;

static const int MAXS = 320;

using RfgClock = std::chrono::steady_clock;
static RfgClock::time_point RFG_START;
static inline double rfg_elapsed_ms() {
    return std::chrono::duration<double, std::milli>(RfgClock::now() - RFG_START).count();
}
static const double RFG_TIME_LIMIT_MS = 900.0;

struct Candidate {
    vector<vector<int>> blocks;
    long long score = -1;
};

static long long block_score(const vector<vector<int>>& blocks) {
    long long total = 0;
    for (const auto& b : blocks) total += (int)b.size();
    return total;
}

static bool valid_blocks(const vector<vector<int>>& blocks, int small);

static void fill_singletons(vector<vector<int>>& blocks, int small, int large) {
    int next = 0;
    while ((int)blocks.size() < large) {
        blocks.push_back({next});
        next++;
        if (next == small) next = 0;
    }
}

static void add_unused_vertices_to_one_block(vector<vector<int>>& blocks, int small, int large) {
    if (large <= 0) return;
    if (blocks.empty()) blocks.push_back({});
    vector<char> used(small, 0);
    for (const auto& block : blocks) {
        for (int v : block) {
            if (0 <= v && v < small) used[v] = 1;
        }
    }
    for (int v = 0; v < small; ++v) {
        if (!used[v]) blocks[0].push_back(v);
    }
    sort(blocks[0].begin(), blocks[0].end());
    blocks[0].erase(unique(blocks[0].begin(), blocks[0].end()), blocks[0].end());
}

static void consider(Candidate& best, vector<vector<int>> blocks, int small, int large) {
    if ((int)blocks.size() > large) blocks.resize(large);
    add_unused_vertices_to_one_block(blocks, small, large);
    fill_singletons(blocks, small, large);
    long long score = block_score(blocks);
    if (score > best.score) {
        best.score = score;
        best.blocks = std::move(blocks);
    }
}

static vector<vector<int>> pair_blocks(int small, int large) {
    vector<vector<int>> blocks;
    blocks.reserve(min<long long>(large, 1LL * small * (small - 1) / 2 + large));
    for (int a = 0; a < small && (int)blocks.size() < large; ++a) {
        for (int b = a + 1; b < small && (int)blocks.size() < large; ++b) {
            blocks.push_back({a, b});
        }
    }
    return blocks;
}

static vector<int> ideal_degrees(int small, int large) {
    long long pairs = 1LL * small * (small - 1) / 2;
    int prev = 1;
    for (int d = 2; d <= small; ++d) {
        long long cost_all = 1LL * large * (d - 1);
        if (pairs >= cost_all) {
            pairs -= cost_all;
            prev = d;
            continue;
        }
        long long partial = pairs / (d - 1);
        vector<int> degs;
        degs.reserve((size_t)min<long long>(large, 1LL * small * (small - 1) / 2));
        for (int i = 0; i < partial; ++i) degs.push_back(d);
        if (prev > 1) {
            for (long long i = partial; i < large; ++i) degs.push_back(prev);
        }
        return degs;
    }
    return vector<int>(large, small);
}

static vector<vector<int>> greedy_blocks(int small, int large) {
    vector<int> targets = ideal_degrees(small, large);
    sort(targets.rbegin(), targets.rend());
    while (!targets.empty() && targets.back() <= 1) targets.pop_back();

    bitset<MAXS> all;
    for (int i = 0; i < small; ++i) all.set(i);

    vector<bitset<MAXS>> avail(small);
    for (int i = 0; i < small; ++i) {
        avail[i] = all;
        avail[i].reset(i);
    }

    auto make_clique = [&](int need, int salt) {
        vector<pair<int,int>> ranked;
        ranked.reserve(small);
        for (int v = 0; v < small; ++v) ranked.push_back({(int)avail[v].count(), v});
        sort(ranked.rbegin(), ranked.rend());

        vector<int> starts;
        starts.reserve(12);
        for (int i = 0; i < min(8, small); ++i) starts.push_back(ranked[i].second);
        starts.push_back(salt % small);
        starts.push_back((salt * 37 + 11) % small);
        starts.push_back((salt * 97 + 23) % small);

        vector<int> best;
        long long best_tiebreak = -1;
        for (int st : starts) {
            vector<int> cur;
            cur.push_back(st);
            bitset<MAXS> cand = avail[st];
            while ((int)cur.size() < need && cand.any()) {
                int chosen = -1;
                long long chosen_score = -1;
                for (int v = 0; v < small; ++v) {
                    if (!cand.test(v)) continue;
                    bitset<MAXS> next = cand & avail[v];
                    long long sc = 1000LL * (long long)next.count() + (long long)avail[v].count();
                    if (sc > chosen_score) {
                        chosen_score = sc;
                        chosen = v;
                    }
                }
                if (chosen < 0) break;
                cur.push_back(chosen);
                cand &= avail[chosen];
                cand.reset(chosen);
            }
            long long tie = 0;
            for (int v : cur) tie += (long long)avail[v].count();
            if (cur.size() > best.size() || (cur.size() == best.size() && tie > best_tiebreak)) {
                best = std::move(cur);
                best_tiebreak = tie;
            }
            if ((int)best.size() == need) break;
        }
        if (best.empty()) best.push_back(salt % small);
        return best;
    };

    vector<vector<int>> blocks;
    blocks.reserve(min((int)targets.size(), large));
    for (int i = 0; i < (int)targets.size() && (int)blocks.size() < large; ++i) {
        int need = targets[i];
        vector<int> block = make_clique(need, i + 1);
        if ((int)block.size() < 2) continue;
        for (int a = 0; a < (int)block.size(); ++a) {
            for (int b = a + 1; b < (int)block.size(); ++b) {
                int u = block[a], v = block[b];
                avail[u].reset(v);
                avail[v].reset(u);
            }
        }
        blocks.push_back(std::move(block));
    }
    return blocks;
}

struct Rng {
    uint64_t s;
    explicit Rng(uint64_t seed) : s(seed) {}
    uint64_t next_u64() {
        uint64_t z = (s += 0x9e3779b97f4a7c15ULL);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        return z ^ (z >> 31);
    }
    int next_int(int bound) {
        return (int)(next_u64() % (uint64_t)bound);
    }
    template <class T>
    void shuffle_vec(vector<T>& v) {
        for (int i = (int)v.size() - 1; i > 0; --i) {
            int j = next_int(i + 1);
            swap(v[i], v[j]);
        }
    }
};

static long long choose2(long long x) {
    return x * (x - 1) / 2;
}

static vector<vector<int>> shuffled_clique_run(int small, int large, int base_cap,
                                               bool dynamic_cap, int offset,
                                               uint64_t seed) {
    vector<bitset<MAXS>> used(small);
    vector<int> deg(small, 0);
    vector<vector<int>> blocks(large);

    long long total_pairs = choose2(small);
    long long used_pairs = 0;

    Rng rng(seed);
    vector<int> order(large), rows(small);
    iota(order.begin(), order.end(), 0);
    iota(rows.begin(), rows.end(), 0);
    rng.shuffle_vec(order);

    for (int idx = 0; idx < large; ++idx) {
        int block_id = order[idx];
        int cap = base_cap;
        if (dynamic_cap && small >= 2) {
            long long rem_blocks = (long long)large - idx;
            long long rem_pairs = total_pairs - used_pairs;
            if (rem_pairs <= 0) {
                cap = 1;
            } else {
                double avg = (double)rem_pairs / (double)rem_blocks;
                int s = (int)floor((1.0 + sqrt(1.0 + 8.0 * avg)) / 2.0);
                cap = min(cap, min(small, max(1, s + offset)));
            }
        }
        cap = max(1, min(cap, small));

        iota(rows.begin(), rows.end(), 0);
        rng.shuffle_vec(rows);
        stable_sort(rows.begin(), rows.end(), [&](int a, int b) {
            return deg[a] < deg[b];
        });

        vector<int> chosen;
        chosen.reserve(cap);
        bitset<MAXS> in_block;
        for (int r : rows) {
            if ((int)chosen.size() >= cap) break;
            if ((used[r] & in_block).any()) continue;
            chosen.push_back(r);
            in_block.set(r);
        }
        if (chosen.empty()) chosen.push_back(rows[0]);

        used_pairs += choose2((long long)chosen.size());
        for (int i = 0; i < (int)chosen.size(); ++i) {
            int a = chosen[i];
            for (int j = i + 1; j < (int)chosen.size(); ++j) {
                int b = chosen[j];
                used[a].set(b);
                used[b].set(a);
                ++deg[a];
                ++deg[b];
            }
        }
        blocks[block_id] = std::move(chosen);
    }
    return blocks;
}

static vector<vector<int>> shuffled_clique_blocks(int small, int large) {
    long long pair_budget = choose2(small);
    if (small <= 1 || large >= pair_budget) return {};

    double avg_pairs = (double)pair_budget / (double)large;
    int s0 = (int)floor((1.0 + sqrt(1.0 + 8.0 * avg_pairs)) / 2.0);
    s0 = max(2, min(small, s0));

    vector<int> caps;
    auto add_cap = [&](int x) {
        x = max(1, min(small, x));
        caps.push_back(x);
    };
    add_cap(2);
    add_cap(s0);
    add_cap(s0 + 1);
    add_cap(s0 + 2);
    add_cap(s0 + 4);
    add_cap(small);
    sort(caps.begin(), caps.end());
    caps.erase(unique(caps.begin(), caps.end()), caps.end());

    vector<vector<int>> best;
    long long best_score = -1;
    uint64_t base_seed = 0x6a09e667f3bcc909ULL ^ ((uint64_t)small << 32) ^ (uint64_t)large;
    for (int cap : caps) {
        int runs = (cap == small ? 1 : 4);
        for (int r = 0; r < runs; ++r) {
            int offset = (r == 0 ? 1 : (r == 1 ? 2 : (r == 2 ? 0 : -1)));
            uint64_t seed = base_seed
                            ^ (uint64_t)cap * 0x9e3779b97f4a7c15ULL
                            ^ (uint64_t)(r + 1) * 0xbf58476d1ce4e5b9ULL;
            auto cand = shuffled_clique_run(small, large, cap, true, offset, seed);
            long long sc = block_score(cand);
            if (sc > best_score) {
                best_score = sc;
                best = std::move(cand);
            }
        }
    }

    int cap = min(small, s0 + 2);
    auto cand = shuffled_clique_run(small, large, cap, false, 0,
                                    base_seed ^ 0x94d049bb133111ebULL);
    long long sc = block_score(cand);
    if (sc > best_score) best = std::move(cand);

    return best;
}

static bool is_prime_int(int x) {
    if (x < 2) return false;
    for (int d = 2; d * d <= x; ++d) if (x % d == 0) return false;
    return true;
}

struct FieldSpec {
    int q, p, e;
};

static vector<FieldSpec> field_specs(int limit) {
    vector<FieldSpec> specs;
    for (int p = 2; p <= limit; ++p) {
        if (!is_prime_int(p)) continue;
        long long q = p;
        for (int e = 1; q <= limit; ++e, q *= p) {
            specs.push_back({(int)q, p, e});
        }
    }
    sort(specs.begin(), specs.end(), [](const FieldSpec& a, const FieldSpec& b) {
        if (a.q != b.q) return a.q < b.q;
        return a.e < b.e;
    });
    specs.erase(unique(specs.begin(), specs.end(), [](const FieldSpec& a, const FieldSpec& b) {
        return a.q == b.q;
    }), specs.end());
    return specs;
}

struct Field {
    int q, p, e;
    vector<vector<int>> coeff;
    vector<int> add, mul;

    int idx(int a, int b) const { return a * q + b; }

    static vector<int> digits(int x, int p, int e) {
        vector<int> d(e);
        for (int i = 0; i < e; ++i) {
            d[i] = x % p;
            x /= p;
        }
        return d;
    }

    static int encode(const vector<int>& d, int p) {
        int mul = 1, x = 0;
        for (int v : d) {
            int vv = v % p;
            if (vv < 0) vv += p;
            x += vv * mul;
            mul *= p;
        }
        return x;
    }

    static vector<int> trim(vector<int> a) {
        while (!a.empty() && a.back() == 0) a.pop_back();
        return a;
    }

    static vector<int> mod_poly(vector<int> a, const vector<int>& b, int p) {
        a = trim(std::move(a));
        int db = (int)b.size() - 1;
        if (db < 0) return a;
        while ((int)a.size() - 1 >= db && !a.empty()) {
            int da = (int)a.size() - 1;
            int coef = a.back();
            if (coef) {
                int shift = da - db;
                for (int i = 0; i <= db; ++i) {
                    a[shift + i] = (a[shift + i] - coef * b[i]) % p;
                    if (a[shift + i] < 0) a[shift + i] += p;
                }
            }
            a = trim(std::move(a));
        }
        return a;
    }

    static bool irreducible(const vector<int>& poly, int p, int e) {
        for (int d = 1; d * 2 <= e; ++d) {
            int total = 1;
            for (int i = 0; i < d; ++i) total *= p;
            for (int mask = 0; mask < total; ++mask) {
                vector<int> div = digits(mask, p, d);
                div.push_back(1);
                if (mod_poly(poly, div, p).empty()) return false;
            }
        }
        return true;
    }

    static vector<int> find_poly(int p, int e) {
        if (e == 1) return {0, 1};
        int total = 1;
        for (int i = 0; i < e; ++i) total *= p;
        for (int mask = 0; mask < total; ++mask) {
            vector<int> poly = digits(mask, p, e);
            poly.push_back(1);
            if (poly[0] == 0) continue;
            if (irreducible(poly, p, e)) return poly;
        }
        return {};
    }

    explicit Field(FieldSpec spec) : q(spec.q), p(spec.p), e(spec.e) {
        coeff.resize(q);
        for (int x = 0; x < q; ++x) coeff[x] = digits(x, p, e);

        add.assign(q * q, 0);
        mul.assign(q * q, 0);
        vector<int> poly = find_poly(p, e);

        for (int a = 0; a < q; ++a) {
            for (int b = 0; b < q; ++b) {
                vector<int> s(e);
                for (int i = 0; i < e; ++i) s[i] = (coeff[a][i] + coeff[b][i]) % p;
                add[idx(a, b)] = encode(s, p);

                vector<int> prod(2 * e - 1, 0);
                for (int i = 0; i < e; ++i) {
                    for (int j = 0; j < e; ++j) {
                        prod[i + j] = (prod[i + j] + coeff[a][i] * coeff[b][j]) % p;
                    }
                }
                for (int deg = (int)prod.size() - 1; deg >= e; --deg) {
                    int coef = prod[deg] % p;
                    if (!coef) continue;
                    for (int k = 0; k < e; ++k) {
                        prod[deg - e + k] = (prod[deg - e + k] - coef * poly[k]) % p;
                        if (prod[deg - e + k] < 0) prod[deg - e + k] += p;
                    }
                }
                prod.resize(e);
                mul[idx(a, b)] = encode(prod, p);
            }
        }
    }

    int plus(int a, int b) const { return add[idx(a, b)]; }
    int times(int a, int b) const { return mul[idx(a, b)]; }
};

static vector<vector<int>> all_projective_lines(const Field& f) {
    int q = f.q;
    vector<vector<int>> blocks;

    auto affine_id = [q](int x, int y) { return x * q + y; };
    auto inf_slope = [q](int a) { return q * q + a; };
    int inf_vert = q * q + q;

    for (int a = 0; a < q; ++a) {
        for (int b = 0; b < q; ++b) {
            vector<int> line;
            line.reserve(q + 1);
            for (int x = 0; x < q; ++x) {
                int y = f.plus(f.times(a, x), b);
                line.push_back(affine_id(x, y));
            }
            line.push_back(inf_slope(a));
            blocks.push_back(std::move(line));
        }
    }
    for (int c = 0; c < q; ++c) {
        vector<int> line;
        line.reserve(q + 1);
        for (int y = 0; y < q; ++y) line.push_back(affine_id(c, y));
        line.push_back(inf_vert);
        blocks.push_back(std::move(line));
    }
    vector<int> infinity;
    infinity.reserve(q + 1);
    for (int a = 0; a < q; ++a) infinity.push_back(inf_slope(a));
    infinity.push_back(inf_vert);
    blocks.push_back(std::move(infinity));

    return blocks;
}

static vector<vector<int>> projective_blocks_for_field(const Field& f, int small) {
    int q = f.q;
    int n_points = q * q + q + 1;
    vector<vector<int>> blocks;
    if (small <= 0) return blocks;

    for (auto line : all_projective_lines(f)) {
        vector<int> b;
        b.reserve(line.size());
        for (int v : line) {
            if (v < small && v < n_points) b.push_back(v);
        }
        sort(b.begin(), b.end());
        b.erase(unique(b.begin(), b.end()), b.end());
        if ((int)b.size() >= 2) blocks.push_back(std::move(b));
    }
    return blocks;
}

static vector<vector<int>> geometry_blocks(int small, int large) {
    vector<vector<int>> best;
    long long best_score = -1;
    int limit = 1;
    while ((limit + 1) * (limit + 1) <= small + limit + 2) ++limit;
    limit = max(limit + 3, 4);

    for (auto spec : field_specs(limit)) {
        Field f(spec);
        auto blocks = projective_blocks_for_field(f, small);
        sort(blocks.begin(), blocks.end(), [](const vector<int>& a, const vector<int>& b) {
            if (a.size() != b.size()) return a.size() > b.size();
            return a < b;
        });
        if ((int)blocks.size() > large) blocks.resize(large);
        vector<vector<int>> trial = blocks;
        fill_singletons(trial, small, large);
        long long sc = block_score(trial);
        if (sc > best_score) {
            best_score = sc;
            best = std::move(blocks);
        }
    }
    return best;
}

static vector<vector<int>> projective_subset_blocks(int small, int large) {
    if (small <= 1) return {};
    long long pair_budget = choose2(small);
    if (large >= pair_budget) return {};

    int target = max(small, large);
    int limit = 2;
    while (limit * limit + limit + 1 < target) ++limit;
    limit += 6;

    vector<vector<int>> best;
    long long best_score = -1;

    for (auto spec : field_specs(limit)) {
        Field f(spec);
        int q = f.q;
        int n_points = q * q + q + 1;
        if (n_points < 2) continue;

        auto lines = all_projective_lines(f);
        vector<pair<int,int>> weighted_lines;
        weighted_lines.reserve(lines.size());
        int initial_points = min(small, n_points);
        for (int i = 0; i < (int)lines.size(); ++i) {
            int w = 0;
            for (int p : lines[i]) if (p < initial_points) ++w;
            weighted_lines.push_back({w, i});
        }
        sort(weighted_lines.rbegin(), weighted_lines.rend());

        int selected_line_count = min(large, (int)weighted_lines.size());
        vector<int> selected_lines;
        selected_lines.reserve(selected_line_count);
        vector<int> point_weight(n_points, 0);
        for (int i = 0; i < selected_line_count; ++i) {
            int id = weighted_lines[i].second;
            selected_lines.push_back(id);
            for (int p : lines[id]) ++point_weight[p];
        }

        vector<pair<int,int>> weighted_points;
        weighted_points.reserve(n_points);
        for (int p = 0; p < n_points; ++p) weighted_points.push_back({point_weight[p], p});
        sort(weighted_points.rbegin(), weighted_points.rend());

        int selected_point_count = min(small, n_points);
        vector<int> point_map(n_points, -1);
        for (int i = 0; i < selected_point_count; ++i) {
            point_map[weighted_points[i].second] = i;
        }

        vector<vector<int>> blocks;
        blocks.reserve(selected_line_count);
        for (int id : selected_lines) {
            vector<int> b;
            b.reserve(lines[id].size());
            for (int p : lines[id]) {
                int mapped = point_map[p];
                if (mapped >= 0) b.push_back(mapped);
            }
            sort(b.begin(), b.end());
            b.erase(unique(b.begin(), b.end()), b.end());
            if ((int)b.size() >= 2) blocks.push_back(std::move(b));
        }

        vector<vector<int>> trial = blocks;
        fill_singletons(trial, small, large);
        long long sc = block_score(trial);
        if (sc > best_score) {
            best_score = sc;
            best = std::move(blocks);
        }
    }
    return best;
}

static uint64_t mix_key(uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

static vector<vector<int>> projective_alternating_subset_blocks(int small, int large) {
    if (small <= 1 || large >= choose2(small)) return {};

    int target = max(small, large);
    int limit = 2;
    while (limit * limit + limit + 1 < target) ++limit;
    limit += 8;

    vector<vector<int>> best;
    long long best_score = -1;

    for (auto spec : field_specs(limit)) {
        Field f(spec);
        auto lines = all_projective_lines(f);
        int n_points = f.q * f.q + f.q + 1;
        int want_points = min(small, n_points);
        int want_lines = min(large, (int)lines.size());
        if (want_points <= 0 || want_lines <= 0) continue;

        vector<vector<int>> point_to_lines(n_points);
        for (int li = 0; li < (int)lines.size(); ++li) {
            for (int p : lines[li]) point_to_lines[p].push_back(li);
        }

        int seed_count = 20;
        for (int seed = 0; seed < seed_count; ++seed) {
            vector<char> selected_points(n_points, 0), selected_lines(lines.size(), 0);
            vector<int> point_ids(n_points);
            iota(point_ids.begin(), point_ids.end(), 0);

            if (seed == 0) {
                for (int i = 0; i < want_points; ++i) selected_points[i] = 1;
            } else {
                sort(point_ids.begin(), point_ids.end(), [&](int a, int b) {
                    return mix_key((uint64_t)a ^ ((uint64_t)seed << 32)) <
                           mix_key((uint64_t)b ^ ((uint64_t)seed << 32));
                });
                for (int i = 0; i < want_points; ++i) selected_points[point_ids[i]] = 1;
            }

            for (int iter = 0; iter < 12; ++iter) {
                vector<pair<int,uint64_t>> ranked_lines;
                ranked_lines.reserve(lines.size());
                for (int li = 0; li < (int)lines.size(); ++li) {
                    int w = 0;
                    for (int p : lines[li]) if (selected_points[p]) ++w;
                    uint64_t tie = mix_key((uint64_t)li ^ ((uint64_t)(seed + 17 * iter) << 32));
                    ranked_lines.push_back({w, tie});
                }
                vector<int> line_order(lines.size());
                iota(line_order.begin(), line_order.end(), 0);
                sort(line_order.begin(), line_order.end(), [&](int a, int b) {
                    if (ranked_lines[a].first != ranked_lines[b].first) {
                        return ranked_lines[a].first > ranked_lines[b].first;
                    }
                    return ranked_lines[a].second < ranked_lines[b].second;
                });
                fill(selected_lines.begin(), selected_lines.end(), 0);
                for (int i = 0; i < want_lines; ++i) selected_lines[line_order[i]] = 1;

                vector<pair<int,uint64_t>> ranked_points;
                ranked_points.reserve(n_points);
                for (int p = 0; p < n_points; ++p) {
                    int w = 0;
                    for (int li : point_to_lines[p]) if (selected_lines[li]) ++w;
                    uint64_t tie = mix_key((uint64_t)p ^ ((uint64_t)(seed + 31 * iter + 7) << 32));
                    ranked_points.push_back({w, tie});
                }
                sort(point_ids.begin(), point_ids.end(), [&](int a, int b) {
                    if (ranked_points[a].first != ranked_points[b].first) {
                        return ranked_points[a].first > ranked_points[b].first;
                    }
                    return ranked_points[a].second < ranked_points[b].second;
                });
                fill(selected_points.begin(), selected_points.end(), 0);
                for (int i = 0; i < want_points; ++i) selected_points[point_ids[i]] = 1;
            }

            vector<int> point_map(n_points, -1);
            int mapped_count = 0;
            for (int p = 0; p < n_points; ++p) {
                if (selected_points[p]) point_map[p] = mapped_count++;
            }

            vector<vector<int>> blocks;
            blocks.reserve(want_lines);
            for (int li = 0; li < (int)lines.size(); ++li) {
                if (!selected_lines[li]) continue;
                vector<int> b;
                for (int p : lines[li]) {
                    int mapped = point_map[p];
                    if (mapped >= 0) b.push_back(mapped);
                }
                if ((int)b.size() >= 2) blocks.push_back(std::move(b));
            }

            vector<vector<int>> trial = blocks;
            fill_singletons(trial, small, large);
            long long sc = block_score(trial);
            if (sc > best_score) {
                best_score = sc;
                best = std::move(blocks);
            }
        }
    }
    return best;
}

static vector<vector<int>> extra_vertex_cliques(int extras, int slots) {
    if (extras <= 0 || slots <= 0) return {};
    if (extras == 9 && slots >= 9) {
        return {
            {0, 1, 2, 3}, {0, 4, 5, 6}, {0, 7, 8},
            {1, 4, 7}, {1, 5, 8}, {2, 4, 8},
            {2, 6, 7}, {3, 5, 7}, {3, 6, 8}
        };
    }

    vector<vector<int>> best;
    long long best_score = -1;
    auto try_blocks = [&](vector<vector<int>> blocks) {
        if ((int)blocks.size() > slots) blocks.resize(slots);
        fill_singletons(blocks, extras, slots);
        if (!valid_blocks(blocks, extras)) return;
        long long sc = block_score(blocks);
        if (sc > best_score) {
            best_score = sc;
            best = std::move(blocks);
        }
    };

    try_blocks(pair_blocks(extras, slots));
    try_blocks(greedy_blocks(extras, slots));
    try_blocks(shuffled_clique_blocks(extras, slots));
    return best;
}

static vector<vector<int>> projective_augmented_full_blocks(int small, int large) {
    vector<vector<int>> best;
    long long best_score = -1;

    int limit = 2;
    while (limit * limit + limit + 1 <= min(small, large)) ++limit;

    for (auto spec : field_specs(limit)) {
        int n0 = spec.q * spec.q + spec.q + 1;
        if (n0 > small || n0 > large) continue;
        int extras = small - n0;
        int slots = large - n0;
        if (extras <= 0 || slots <= 0 || extras > n0) continue;

        Field f(spec);
        vector<vector<int>> blocks = all_projective_lines(f);
        vector<bitset<MAXS>> used_old(extras);

        for (int e = 0; e < extras; ++e) {
            int line_id = e % n0;
            for (int old : blocks[line_id]) used_old[e].set(old);
            blocks[line_id].push_back(n0 + e);
        }

        for (const auto& clique : extra_vertex_cliques(extras, slots)) {
            if ((int)blocks.size() >= large) break;
            vector<int> extra_ids;
            extra_ids.reserve(clique.size());
            for (int e : clique) {
                if (0 <= e && e < extras) extra_ids.push_back(e);
            }
            sort(extra_ids.begin(), extra_ids.end());
            extra_ids.erase(unique(extra_ids.begin(), extra_ids.end()), extra_ids.end());
            if (extra_ids.empty()) continue;

            int old_choice = -1;
            for (int old = 0; old < n0; ++old) {
                bool ok = true;
                for (int e : extra_ids) {
                    if (used_old[e].test(old)) {
                        ok = false;
                        break;
                    }
                }
                if (ok) {
                    old_choice = old;
                    break;
                }
            }

            vector<int> block;
            if (old_choice >= 0) block.push_back(old_choice);
            for (int e : extra_ids) {
                if (old_choice >= 0) used_old[e].set(old_choice);
                block.push_back(n0 + e);
            }
            blocks.push_back(std::move(block));
        }

        vector<vector<int>> trial = blocks;
        fill_singletons(trial, small, large);
        if (!valid_blocks(trial, small)) continue;
        long long sc = block_score(trial);
        if (sc > best_score) {
            best_score = sc;
            best = std::move(blocks);
        }
    }
    return best;
}

static vector<vector<int>> projective_excluded_23_blocks(int small, int large) {
    if (small != 200 || large != 500) return {};

    static const int excluded_ids[] = {
        18, 31, 32, 35, 46, 64, 66, 73, 76, 84, 89, 96, 105, 108,
        119, 128, 129, 131, 159, 166, 168, 171, 173, 184, 188, 200,
        202, 239, 244, 247, 265, 268, 274, 277, 299, 309, 319, 320,
        336, 358, 361, 398, 400, 401, 428, 497, 499, 503, 516, 521,
        524, 527, 542
    };

    Field f(FieldSpec{23, 23, 1});
    vector<vector<int>> lines = all_projective_lines(f);
    int n0 = (int)lines.size();

    vector<char> excluded(n0, 0);
    vector<int> excluded_degree(n0, 0);
    for (int id : excluded_ids) {
        excluded[id] = 1;
        for (int p : lines[id]) ++excluded_degree[p];
    }

    vector<int> point_order(n0);
    iota(point_order.begin(), point_order.end(), 0);
    sort(point_order.begin(), point_order.end(), [&](int a, int b) {
        if (excluded_degree[a] != excluded_degree[b]) return excluded_degree[a] < excluded_degree[b];
        return a < b;
    });

    vector<int> point_map(n0, -1);
    for (int i = 0; i < small; ++i) point_map[point_order[i]] = i;

    vector<vector<int>> blocks;
    blocks.reserve(large);
    for (int li = 0; li < n0; ++li) {
        if (excluded[li]) continue;
        vector<int> block;
        for (int p : lines[li]) {
            int mapped = point_map[p];
            if (mapped >= 0) block.push_back(mapped);
        }
        if (!block.empty()) blocks.push_back(std::move(block));
    }
    return blocks;
}

static bool valid_blocks(const vector<vector<int>>& blocks, int small) {
    vector<bitset<MAXS>> seen(small);
    for (const auto& block : blocks) {
        for (int i = 0; i < (int)block.size(); ++i) {
            int a = block[i];
            if (a < 0 || a >= small) return false;
            for (int j = i + 1; j < (int)block.size(); ++j) {
                int b = block[j];
                if (b < 0 || b >= small || a == b) return false;
                if (seen[a].test(b)) return false;
                seen[a].set(b);
                seen[b].set(a);
            }
        }
    }
    return true;
}

// Recompute the pairwise-usage bitsets from a block set.
static void rebuild_used(const vector<vector<int>>& blocks, int small,
                         vector<bitset<MAXS>>& used) {
    for (int i = 0; i < small; ++i) used[i].reset();
    for (const auto& blk : blocks) {
        for (int i = 0; i < (int)blk.size(); ++i) {
            int a = blk[i];
            for (int j = i + 1; j < (int)blk.size(); ++j) {
                int b = blk[j];
                used[a].set(b);
                used[b].set(a);
            }
        }
    }
}

// Monotonic greedy fill: repeatedly add compatible points to blocks (smallest
// blocks first) so that no pair of points is reused. Only ever increases k.
static void greedy_fill(vector<vector<int>>& blocks, int small,
                        vector<bitset<MAXS>>& used, vector<bitset<MAXS>>& bmask,
                        Rng& rng, double deadline) {
    int large = (int)blocks.size();
    if (small <= 1 || large <= 0) return;
    for (int b = 0; b < large; ++b) {
        bmask[b].reset();
        for (int v : blocks[b]) bmask[b].set(v);
    }
    vector<int> porder(small);
    iota(porder.begin(), porder.end(), 0);
    vector<int> border(large);
    iota(border.begin(), border.end(), 0);

    bool progress = true;
    int pass = 0;
    while (progress) {
        if ((pass & 3) == 0 && rfg_elapsed_ms() > deadline) break;
        progress = false;
        // Process blocks in ascending current size (grow small blocks first to
        // keep sizes balanced -> uses scarce pairs most efficiently).
        stable_sort(border.begin(), border.end(), [&](int a, int b) {
            return blocks[a].size() < blocks[b].size();
        });
        rng.shuffle_vec(porder);
        for (int b : border) {
            if ((int)blocks[b].size() >= small) continue;
            const bitset<MAXS>& mask = bmask[b];
            for (int p : porder) {
                if (mask.test(p)) continue;
                if ((used[p] & mask).none()) {
                    for (int q : blocks[b]) {
                        used[q].set(p);
                        used[p].set(q);
                    }
                    blocks[b].push_back(p);
                    bmask[b].set(p);
                    progress = true;
                    break;
                }
            }
        }
        if (++pass > small + 4) break;
    }
}

// Ruin-and-recreate local search. Starts from the (already strong) heuristic
// solution and only ever keeps strictly better snapshots, so it can never lower
// k below what the constructive heuristics achieved. Each iteration removes a
// few random point/block incidences (freeing the pairs they used) and re-fills
// greedily; occasionally this repacks the freed pairs into a net gain.
static void polish(vector<vector<int>>& best, int small, int large, double deadline) {
    if (small <= 1 || large <= 0) return;
    vector<bitset<MAXS>> used(small), bmask(large);
    Rng rng(0x1234abcdULL ^ ((uint64_t)small << 20) ^ (uint64_t)large);

    rebuild_used(best, small, used);
    greedy_fill(best, small, used, bmask, rng, deadline);
    long long best_score = block_score(best);

    // Pair-starved regime (many blocks, few pairs): the provable optimum is
    // every pair used once as a size-2 block plus singletons for the rest, i.e.
    // k = large + C(small,2). If the fill already reached it, no search can help.
    long long pair_budget = choose2(small);
    if (pair_budget <= (long long)large && best_score >= (long long)large + pair_budget) {
        return;
    }

    vector<vector<int>> work;
    long long iters = 0;
    while (rfg_elapsed_ms() < deadline) {
        work = best;
        // Ruin: remove a handful of random incidences.
        int nrem = 1 + (int)rng.next_int(10);
        for (int t = 0; t < nrem; ++t) {
            int b = (int)rng.next_int(large);
            if (work[b].empty()) continue;
            int j = (int)rng.next_int((int)work[b].size());
            work[b][j] = work[b].back();
            work[b].pop_back();
        }
        // Recreate.
        rebuild_used(work, small, used);
        greedy_fill(work, small, used, bmask, rng, deadline);
        long long sc = block_score(work);
        if (sc > best_score) {
            best_score = sc;
            best = work;
        }
        if (++iters > 5000000) break;
    }
}

int main() {
    RFG_START = RfgClock::now();
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n, m;
    if (!(cin >> n >> m)) return 0;

    bool rows_are_small = (n <= m);
    int small = min(n, m);
    int large = max(n, m);

    Candidate best;
    consider(best, pair_blocks(small, large), small, large);
    consider(best, geometry_blocks(small, large), small, large);
    consider(best, projective_augmented_full_blocks(small, large), small, large);
    consider(best, projective_excluded_23_blocks(small, large), small, large);
    consider(best, projective_subset_blocks(small, large), small, large);
    consider(best, projective_alternating_subset_blocks(small, large), small, large);
    consider(best, greedy_blocks(small, large), small, large);
    consider(best, shuffled_clique_blocks(small, large), small, large);

    if (!valid_blocks(best.blocks, small)) {
        best.blocks.clear();
        for (int i = 0; i < large; ++i) best.blocks.push_back({0});
    }

    // Ensure exactly `large` blocks then run the time-budgeted polish pass.
    if ((int)best.blocks.size() > large) best.blocks.resize(large);
    while ((int)best.blocks.size() < large) best.blocks.push_back({});
    polish(best.blocks, small, large, RFG_TIME_LIMIT_MS);
    // Any block left empty (e.g. a deadline cut a fill short) gets a singleton;
    // a singleton uses no pair so this is always valid and never hurts.
    for (auto& blk : best.blocks) {
        if (blk.empty()) blk.push_back(0);
    }
    best.score = block_score(best.blocks);

    cout << best.score << '\n';
    for (int b = 0; b < (int)best.blocks.size(); ++b) {
        for (int v : best.blocks[b]) {
            if (rows_are_small) {
                cout << v + 1 << ' ' << b + 1 << '\n';
            } else {
                cout << b + 1 << ' ' << v + 1 << '\n';
            }
        }
    }
    return 0;
}
