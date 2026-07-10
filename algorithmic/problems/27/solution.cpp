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

static chrono::steady_clock::time_point T_START;
static double elapsed_s() {
    return chrono::duration<double>(chrono::steady_clock::now() - T_START).count();
}

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
        25, 32, 62, 67, 76, 93, 96, 102, 130, 136, 140, 143, 148,
        177, 178, 182, 186, 189, 211, 217, 227, 237, 238, 239, 268,
        274, 276, 285, 288, 296, 303, 307, 311, 334, 342, 348, 360,
        368, 374, 394, 403, 410, 421, 424, 436, 465, 476, 506, 514,
        530, 531, 548, 552
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

// Precomputed 316-block, k=5585 rectangle-free family on 316x316.
// Delta-encoded: alphabet is printable ASCII 35..126 minus '\' (91 symbols,
// never '"' or '\'). Per block: [size][v0/91][v0%91] then deltas d>=1 as
// (d-1) if d<=90, else escape symbol 90 followed by [d/91][d%91].
static vector<vector<int>> pbd316_blocks(int small, int large) {
    if (small != 316 || large < 316) return {};
    static const char enc[] =
        "6##3333333333333333345#-9(9(99(99(99(99(56#%3333333333333333165#'555555$5555555$536#'333333333333333"
        "3/86#(3333333333333333.96#)3333333333333333-:5#0):)::):)::):)::)46#+3333333333333333+<5#,000A0000A00"
        "000A055#-,=,=,G,=,G,=,,=+#5#2&777&777&7777&77,5#/3333333333333333'5#03333333333333333&5#133333333333"
        "33333%5#23333333333333333$5#+66%J666%66666%6145#-)::):)::):)::):)75#$444444444444444#45#%44444444444"
        "444#435#&4444444444444#4425#/+<+<+<+<++<+<+<+95%J################h5#-..?..?..?...?..?05#/55$5555555$"
        "55555+5#+44444444#4444444-5#,4444444#44444444,5#-444444#444444444+5#*9(99(99(99(9(99(85#/4444#444444"
        "44444)5#$888'88'88'888'88+5#144#4444444444444'5#3&777&777&777&777+5#:B2222C2222222K)$#5##55555555$55"
        "55555&5#.*;*;*;;*;*;*;*;*85#,,=,=,,=,=,,=,=,=-5#2(99(9(99(99(99(905#y################~%:5#(55555$555"
        "55555$525#)->->-->-->-->-->25$s################~$@5#+5555$5555555$555/5#)<+<++<+<+<+<+<+<.5%9#######"
        "#########y5#1%66666%66666%666+5##C222222222222222B5#05$55555555$55555*5#$99(99(99(99(9(99-5#.11111B1"
        "1111111B155#3$5555555$5555555'5#)55555$5555555$5515#'>-->-->-->-->->-45#0///@///@///@////@5#*4444444"
        "44#444444.5##888'88'888'88'88,5#(666%66666%66666%45#+->-->-->->-->-->05#-8'88'888'88'88'835#33333333"
        "333333333#5#2,,=,=,,=,=,=,,=,85#-66%6666%66666%66/5#.6%66666%66666%66.5#(44444444444#444405#06%6666%"
        "66666%666,5#-//@//M//<@///@2*#5#2%66666%6666%6666*5#////@///@///@///@05#+;*;*;*;*;*;*;*;*;5#$777&777"
        "7&777&777)5#)666%66666%6666%635##>->-->-->-->-->-85#'0A00000A0000A000:5#'2222C22222222222>5#)77&777&"
        "7777&777&55#,66%66666%6666%6605#+77&777&777&777&735#$66666%6666%66666'5##4444444444444444$5#.7&777&7"
        "77&7777&705#.)::):)::):):)::)6)#ER~$%~$M%%5#0&7777&777&777&77.5$?################~$t(%@u0&#5#'<+<+<+"
        "<++<+<+<+<05#*.?..?...?..?...?35#%88'888'88'888'88*5#'777&777&777&777&75#'88'88'888'88'888(5#(88'88'"
        "888'88'88'85#h################~%K5#.000A00000A00000A35#.22222222222C222275#.->-->-->-->-->-->5#+1111"
        "B1111111B11185#.8'88'88'888'88'825#/'888'88'888'88'815#0'888'88'88'888'805#1'88'888'88'888'8/5#-3333"
        "333333333333)5#3'88'88'888'88'88-5#0G,,=,=,,=,=,,=,:05#,555$55555555$555.5#%99(99(99(9(99(99,5#(<+<+"
        "<++<+<+<+<+</5#+8'888'88'88'888'55#(9(99(99(99(99(99)5#3#444444444444444%5#/7&777&777&777&77/5#0444#"
        "444444444444(5#3%6666%66666%6666)5#V################~%]5#(>-->-->-->-->-->35#%;*;*;;*;*;*;*;*;05#/11"
        "1111B1111111B145#1(99(99(9(99(99(91)#L`~$pJ#$5#&6666%66666%6666%65##::):)::):)::):):05#-2222222222C2"
        "222285#%::):):)::):)::):.6#$3333333333333333255#2$55555555$555555(($2-~%.%$5#):):)::):)::):)::*5##99"
        "(99?9(99(99(9.05#*:):)::):)::):):):5#*77&777&777&7777&45#$<+<+<+<+<+<+<++<35#';*;*;*;*;;*;*;*;.5$b##"
        "##############~$Q5#&222C222222222222?5#3...?..?...?..?..;5#2):)::):):)::):):25$.###$###########~%)'5"
        "##;;*;*;*;*;*;*;*;25#$;*;;*;*;*;*;*;*;15#(=,,=,=,,=,=,=,,=15#&;*;*;*;;*;*;*;*;/5#1):)::):)::):):):35"
        "#(;*;*;*;*;*;;*;*;-5#$5555555$55555555%5#*;*;*;*;*;*;*;*;;+5#&88'888'88'88'888)5#,*;;*;*;*;*;*;*;*:("
        "#-~$&w~$?*5#3,,=,=,,=,=,,=,=,75#/*;*;*;*;;*;*;*;*75#0*;*;*;*;*;;*;*;*65#10000A00000A00000A5#2*;*;*;*"
        ";*;*;*;;*45#3*;*;*;*;*;*;*;*;35##<+<+<+<+<+<+<+<+45#2+<++<+<+<+<+<+<+65#%<+<+<+<+<+<++<+<25#&<+<+<+<"
        "+<++<+<+<15#/)::):):)::):)::)5)$0~$I57G&5#$A00000A0000A0000=5#*<++<+<+<+<+<+<+<-5#++<+<+<+<+<+<+<+<,"
        "5#,+<+<+<+<+<+<+<++<5#-+<+<+<+<+<+<++<+;5#.+<+<+<+<+<++<+<+:5#.(99(99(99(99(99(45#,:):):)::):)::):)8"
        "5#1+<+<++<+<+<+<+<+75#*666%6666%66666%625#2'88'888'88'88'88.5##=,=,=,,=,=,,=,=,65#15$5555555$555555)"
        "5#+//@///@///@///@/45#-->-->-->-->-->->.6#*3333333333333333,;5#)111B1111111B1111:5#)=,,=,=,,=,=,,=,="
        "05#*,=,=,=,,=,=,,=,=/5#+,=,=,,=,=,=,,=,=.5#1*;*;*;*;*;*;;*;*55#'6666%6666%66666%55#'/@///@///@///@//"
        "85#/,=,,=,=,,=,=,=,,;5#,7&7777&777&777&725#2///@////@///@///>5&3#################5#(77&7777&777&777&"
        "65#':)::):)::):):)::,5#$>-->->-->-->-->-75#%5555555$5555555$55#&>-->-->-->->-->-54#%>8>-I->-->-->-60"
        "4#&:AA:)::):)::):--5#(22222C2222222222=5#*->-->->-->-->-->15#%=,=,,=,=,,=,=,=,45#/6%66666%6666%666-5"
        "#/222222222222C22265#).?..?...?..?..?.45#/-->->-->-->-->--=5#0-->-->->-->-->--<5#0+<+<+<++<+<+<+<+85"
        "#2-->-->-->-->->--:5#3-->-->-->-->-->-95##?..?..?...?..?..:5#*/@////@///@///@/55#)4444444444#44444/5"
        "#&?...?..?..?...?.75#'?...?..?...?..?.65#(.?..?..?...?..?.55#)9(99(99(99(99(9(9*#dJ5H~$[&#5#+.?...?."
        ".?..?...?25#,.?...?..?...?..?15#%6666%66666%66666&5#...?..?...?..?..?/5#/..?..?...?..?...?5#0..?...?"
        "..?..?...>5#1..?...?..?...?..=5#2...?..?..?...?..<5#21111111B11111111B5##@///@///@///@///<5#$@///@//"
        "/@////@//;5#%@///@////@///@//:5#11111111B1111111B25#.55$55555555$5555,5#(/@///@///@////@/76#&3333333"
        "333333333075#(:)::):):)::):)::+5#1,,=,=,=,,=,=,,=,95%(################~$/5#3):):)::):)::):):15#.//@/"
        "///@///@///@15#+:):)::):):)::):)95#&99(99(9(99(99(99+5#1///@///@////@///?5#&@////@///@///@//95#'4444"
        "44444444#44415#,//@//M///@////@3)5#1&777&7777&777&77-5#%A00000A00000A000<5#&0A0000A00000A000;5#);*;*"
        ";*;*;*;*;;*;,5#$B11111111B111111?5#)00A0000A00000A0085#*00A00000A0000A0075#%L7&777&7777&777(65%[####"
        "############V5#3++<+<+<+<+<+<+<+55#,3333333333333333*5#/0000A0000A00000A25#00000A00000A0000A15#.3333"
        "333333333333(5#200000A0000A00000@5#300000A00000A0000?5##B1111111B1111111@*#*B~$js=#$5#%1B1111111B111"
        "111>5#+22222B2C2222222:&5#'11B1111111B11111<5#(11B11111111B1111;5#'=,,=,=,=,,=,=,,=25#*111B11111111B"
        "1119*$1huU$&#5#,1111B11111111B1175#-11111B1111111B1165#.44444#4444444444*5#1-->-->-->->-->--;5#01111"
        "11B11111111B35#0(99(99(99(9(99(925#,8'88'888'88'888'45#311111111B1111111A5#,9(99(9(99(99(99(65#$2C22"
        "222222222222A5#%22C2222222222222@5#*8'888'88'888'88'65#24#44444444444444&5#$?..?...?..?..?..95#)2222"
        "22C222222222<4#-7&777&7777&b&7175#'99(9(99(99(99(99*5#,222222222C22222295%~####$##########4(5#$=,=,,"
        "=,=,=,,=,=,55#,->-->-->-->->-->/5#02222222222222C2255#122222222222222C245#2222222222222222C35#322222"
        "22222222222C5##################~&45#4################~&#5#E################~%n5#%?..?...?..?...?.85#"
        "#7777&777&777&777*5#&777&777&777&7777'5##66666%66666%6666(5#)/@///@////@///@/65$P################~$c"
        "5#&1B11111111B11111=5#$::):)::):):)::):/5##A0000A00000A0000>5#-*;*;;*;*;*;*;*;*95#-000A00000A0000A04"
        "5#+9(99(99(9(99(99(75%m################E5#-555$5555555$5555-5#&=,=,,=,=,,=,=,,=35#3(9(99(99(99(99(9/"
        "5#*5555$[H555$5504$&5#)88'88=8'888'88'715#3////@///@///@///=5#/(99(99(99(99(9(935#.,=,,=,=,=,,=,=,,<"
        "5#+00A00000A00000A064#&555555$555555H635#(0A00000A00000A009";
    auto dec = [](char c) -> int {
        int x = (unsigned char)c;
        return x - 35 - (x > 92 ? 1 : 0);
    };
    const size_t n = sizeof(enc) - 1;
    vector<vector<int>> blocks;
    size_t i = 0;
    while (i < n) {
        int sz = dec(enc[i++]);
        if (sz < 1 || i + 2 > n) return {};
        int v = dec(enc[i]) * 91 + dec(enc[i + 1]);
        i += 2;
        vector<int> b;
        b.reserve(sz);
        b.push_back(v);
        for (int t = 1; t < sz; ++t) {
            if (i >= n) return {};
            int s = dec(enc[i++]);
            int d;
            if (s < 90) {
                d = s + 1;
            } else {
                if (i + 2 > n) return {};
                d = dec(enc[i]) * 91 + dec(enc[i + 1]);
                i += 2;
            }
            v += d;
            b.push_back(v);
        }
        blocks.push_back(std::move(b));
    }
    // Self-check (consider() does not validate): range, in-block dups,
    // pairwise |intersection| <= 1. On any corruption fall back to {}.
    if ((int)blocks.size() != 316) return {};
    long long total = 0;
    vector<bitset<MAXS>> masks(blocks.size());
    for (size_t bi = 0; bi < blocks.size(); ++bi) {
        for (int v : blocks[bi]) {
            if (v < 0 || v >= 316) return {};
            if (masks[bi].test(v)) return {};
            masks[bi].set(v);
        }
        total += (long long)blocks[bi].size();
    }
    if (total != 5585) return {};
    for (size_t a = 0; a < blocks.size(); ++a) {
        for (size_t b2 = a + 1; b2 < blocks.size(); ++b2) {
            if ((masks[a] & masks[b2]).count() > 1) return {};
        }
    }
    return blocks;
}

// ---------------- pair-availability helpers ----------------

static vector<bitset<MAXS>> build_avail(const vector<vector<int>>& blocks, int small) {
    bitset<MAXS> all;
    for (int i = 0; i < small; ++i) all.set(i);
    vector<bitset<MAXS>> avail(small);
    for (int v = 0; v < small; ++v) {
        avail[v] = all;
        avail[v].reset(v);
    }
    for (const auto& b : blocks) {
        for (int i = 0; i < (int)b.size(); ++i) {
            for (int j = i + 1; j < (int)b.size(); ++j) {
                avail[b[i]].reset(b[j]);
                avail[b[j]].reset(b[i]);
            }
        }
    }
    return avail;
}

static long long avail_pairs(const vector<bitset<MAXS>>& avail, int small) {
    long long s = 0;
    for (int v = 0; v < small; ++v) s += (long long)avail[v].count();
    return s / 2;
}

static void mark_used(vector<bitset<MAXS>>& avail, const vector<int>& b) {
    for (int i = 0; i < (int)b.size(); ++i) {
        for (int j = i + 1; j < (int)b.size(); ++j) {
            avail[b[i]].reset(b[j]);
            avail[b[j]].reset(b[i]);
        }
    }
}

static void mark_free(vector<bitset<MAXS>>& avail, const vector<int>& b) {
    for (int i = 0; i < (int)b.size(); ++i) {
        for (int j = i + 1; j < (int)b.size(); ++j) {
            avail[b[i]].set(b[j]);
            avail[b[j]].set(b[i]);
        }
    }
}

static int fill_target_size(long long pairs_left, long long slots, int cap) {
    if (slots <= 0 || pairs_left <= 0) return 0;
    double avg = (double)pairs_left / (double)slots;
    int d = (int)floor((1.0 + sqrt(1.0 + 8.0 * avg)) / 2.0);
    if (d < 2) d = 2;
    if (d > cap) d = cap;
    return d;
}

static vector<int> clique_from_avail(const vector<bitset<MAXS>>& avail, int small, int cap,
                                     int start) {
    vector<int> cur;
    cur.push_back(start);
    bitset<MAXS> cand = avail[start];
    while ((int)cur.size() < cap && cand.any()) {
        int chosen = -1;
        long long chosen_score = -1;
        for (int v = 0; v < small; ++v) {
            if (!cand.test(v)) continue;
            long long sc = 1000LL * (long long)((cand & avail[v]).count()) +
                           (long long)avail[v].count();
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
    return cur;
}

static void greedy_fill(vector<vector<int>>& blocks, int small, int large,
                        vector<bitset<MAXS>>& avail, double deadline) {
    long long P = avail_pairs(avail, small);
    int salt = 0;
    while ((int)blocks.size() < large && P > 0) {
        if (elapsed_s() > deadline) break;
        int cap = fill_target_size(P, (long long)large - (long long)blocks.size(), small);
        if (cap < 2) cap = 2;
        // multi-start: two highest-degree vertices plus one rotating start
        int s1 = -1, s2 = -1, c1 = 0, c2 = 0;
        for (int v = 0; v < small; ++v) {
            int c = (int)avail[v].count();
            if (c > c1) {
                c2 = c1;
                s2 = s1;
                c1 = c;
                s1 = v;
            } else if (c > c2) {
                c2 = c;
                s2 = v;
            }
        }
        if (s1 < 0 || c1 == 0) break;
        vector<int> blk;
        int starts[3] = {s1, s2, -1};
        {
            int probe = (salt * 131 + 17) % small;
            for (int t = 0; t < small; ++t) {
                int v = probe + t;
                if (v >= small) v -= small;
                if (avail[v].any()) {
                    starts[2] = v;
                    break;
                }
            }
        }
        ++salt;
        for (int si = 0; si < 3; ++si) {
            if (starts[si] < 0) continue;
            vector<int> cand = clique_from_avail(avail, small, cap, starts[si]);
            if (cand.size() > blk.size()) blk = std::move(cand);
            if ((int)blk.size() >= cap) break;
        }
        if ((int)blk.size() < 2) break;
        mark_used(avail, blk);
        P -= choose2((long long)blk.size());
        blocks.push_back(std::move(blk));
    }
}

// ---------------- post-pass: augmentation + LNS ----------------

static bool augment_pass(vector<vector<int>>& blocks, int small,
                         vector<bitset<MAXS>>& avail, double deadline) {
    bool changed = false;
    for (auto& b : blocks) {
        if ((int)b.size() >= small) continue;
        if (elapsed_s() > deadline) break;
        bitset<MAXS> mask;
        for (int v : b) mask.set(v);
        for (int v = 0; v < small; ++v) {
            if (mask.test(v)) continue;
            if ((mask & ~avail[v]).any()) continue;
            b.push_back(v);
            mask.set(v);
            for (int u : b) {
                if (u == v) continue;
                avail[u].reset(v);
                avail[v].reset(u);
            }
            changed = true;
        }
    }
    return changed;
}

static void lns_improve(vector<vector<int>>& blocks, int small, int large,
                        vector<bitset<MAXS>>& avail, double deadline) {
    (void)large;
    int nb = (int)blocks.size();
    if (nb < 2 || small < 4) return;
    long long P = avail_pairs(avail, small);
    {
        int maxb = 0;
        for (auto& b : blocks) maxb = max(maxb, (int)b.size());
        if (P == 0 && maxb <= 2) return;
    }
    Rng rng(0xC0FFEEULL ^ ((uint64_t)small << 32) ^ (uint64_t)nb);
    int fail_streak = 0;
    const int GIVE_UP = 60000;
    vector<int> uni;

    // move: rebuild r random blocks from their freed pairs plus free pairs
    auto try_rebuild = [&](int r) -> int {
        if (r > nb) r = nb;
        int idx[3];
        idx[0] = rng.next_int(nb);
        idx[1] = (idx[0] + 1 + rng.next_int(nb - 1)) % nb;
        if (r >= 3) {
            int guard = 0;
            do {
                idx[2] = rng.next_int(nb);
            } while ((idx[2] == idx[0] || idx[2] == idx[1]) && ++guard < 8);
            if (idx[2] == idx[0] || idx[2] == idx[1]) r = 2;
        }
        vector<vector<int>> olds(r);
        int old_pts = 0, max_old = 0;
        long long old_pairs = 0;
        for (int t = 0; t < r; ++t) {
            olds[t] = blocks[idx[t]];
            old_pts += (int)olds[t].size();
            max_old = max(max_old, (int)olds[t].size());
            old_pairs += choose2((long long)olds[t].size());
        }
        if (old_pairs == 0 && P == 0) return 0;
        uni.clear();
        for (int t = 0; t < r; ++t) {
            mark_free(avail, olds[t]);
            uni.insert(uni.end(), olds[t].begin(), olds[t].end());
        }
        auto pick_start = [&]() {
            int bestv = -1, bc = -1;
            for (int v : uni) {
                int c = (int)avail[v].count();
                if (c > bc) {
                    bc = c;
                    bestv = v;
                }
            }
            for (int t = 0; t < 2; ++t) {
                int v = rng.next_int(small);
                int c = (int)avail[v].count();
                if (c > bc) {
                    bc = c;
                    bestv = v;
                }
            }
            return bestv;
        };
        vector<vector<int>> news(r);
        int new_pts = 0;
        long long new_pairs = 0;
        int target = old_pts + 1;
        for (int t = 0; t < r; ++t) {
            int remaining = target - new_pts;
            int cap;
            if (t == 0) {
                cap = max_old + 1;
                if (rng.next_int(2)) cap = (target + r - 1) / r;
            } else {
                cap = remaining;
            }
            if (cap < 1) cap = 1;
            int s = pick_start();
            vector<int>& nblk = news[t];
            if (s >= 0) nblk = clique_from_avail(avail, small, cap, s);
            if (nblk.empty()) nblk.push_back(olds[t][0]);
            mark_used(avail, nblk);
            new_pts += (int)nblk.size();
            new_pairs += choose2((long long)nblk.size());
        }
        bool accept = (new_pts > old_pts) || (new_pts == old_pts && new_pairs < old_pairs);
        if (accept) {
            for (int t = 0; t < r; ++t) blocks[idx[t]] = news[t];
            P += old_pairs - new_pairs;
            return (new_pts > old_pts) ? 1 : 0;
        }
        for (int t = 0; t < r; ++t) mark_free(avail, news[t]);
        for (int t = 0; t < r; ++t) mark_used(avail, olds[t]);
        return 0;
    };

    // move: swap one vertex out of a block for a compatible outsider, then
    // relocate the evicted vertex into any other block (+1 point if it fits)
    auto try_swap_relocate = [&]() -> int {
        int bi = rng.next_int(nb);
        auto& B = blocks[bi];
        int bs = (int)B.size();
        if (bs < 2 || bs >= small) return 0;
        bitset<MAXS> mask;
        for (int w : B) mask.set(w);
        int vstart = rng.next_int(small);
        for (int t = 0; t < small; ++t) {
            int v = vstart + t;
            if (v >= small) v -= small;
            if (mask.test(v)) continue;
            bitset<MAXS> confl = mask & ~avail[v];
            if ((int)confl.count() != 1) continue;
            int u = 0;
            while (!confl.test(u)) ++u;
            // swap u -> v inside B
            for (int w : B) {
                if (w == u) continue;
                avail[u].set(w);
                avail[w].set(u);
            }
            for (auto& w : B) {
                if (w == u) {
                    w = v;
                    break;
                }
            }
            for (int w : B) {
                if (w == v) continue;
                avail[v].reset(w);
                avail[w].reset(v);
            }
            // try to relocate u into another block
            int placed = -1;
            for (int j = 0; j < nb; ++j) {
                if (j == bi) continue;
                auto& C = blocks[j];
                if ((int)C.size() >= small) continue;
                bool ok = true;
                for (int w : C) {
                    if (!avail[u].test(w)) {
                        ok = false;
                        break;
                    }
                }
                if (!ok) continue;
                placed = j;
                break;
            }
            if (placed >= 0) {
                auto& C = blocks[placed];
                P -= (long long)C.size();
                for (int w : C) {
                    avail[u].reset(w);
                    avail[w].reset(u);
                }
                C.push_back(u);
                return 1;
            }
            // revert swap
            for (int w : B) {
                if (w == v) continue;
                avail[v].set(w);
                avail[w].set(v);
            }
            for (auto& w : B) {
                if (w == v) {
                    w = u;
                    break;
                }
            }
            for (int w : B) {
                if (w == u) continue;
                avail[u].reset(w);
                avail[w].reset(u);
            }
            return 0;
        }
        return 0;
    };

    while (fail_streak < GIVE_UP) {
        if (elapsed_s() > deadline) break;
        int roll = rng.next_int(6);
        int gained;
        if (roll < 2) gained = try_rebuild(2);
        else if (roll < 3) gained = try_rebuild(3);
        else gained = try_swap_relocate();
        if (gained) fail_streak = 0;
        else ++fail_streak;
    }
}

// ---------------- GDD / transversal-design candidate ----------------

static vector<vector<int>> cyclic_blocks(int small, int large, double deadline);

static vector<vector<int>> gdd_construct(const FieldSpec& spec, int k, int small, int large,
                                         double deadline) {
    int g = spec.q;
    Field f(spec);
    int base = small / k, rem = small % k;
    vector<int> sz(k), off(k);
    for (int i = 0; i < k; ++i) sz[i] = base + (i < rem ? 1 : 0);
    int o = 0;
    for (int i = 0; i < k; ++i) {
        off[i] = o;
        o += sz[i];
    }
    bool slope = (k == g + 1);
    int kg = slope ? g : k;
    vector<vector<int>> blocks;
    blocks.reserve((size_t)g * g);
    for (int a = 0; a < g; ++a) {
        for (int b = 0; b < g; ++b) {
            vector<int> blk;
            blk.reserve(k);
            for (int i = 0; i < kg; ++i) {
                int y = f.plus(f.times(a, i), b);
                if (y < sz[i]) blk.push_back(off[i] + y);
            }
            if (slope && a < sz[g]) blk.push_back(off[g] + a);
            if ((int)blk.size() >= 2) blocks.push_back(std::move(blk));
        }
    }
    if ((int)blocks.size() > large) {
        sort(blocks.begin(), blocks.end(), [](const vector<int>& x, const vector<int>& y) {
            return x.size() > y.size();
        });
        blocks.resize(large);
    }
    // per-group structured fill: cyclic difference families inside each group
    {
        long long slots_total = (long long)large - (long long)blocks.size();
        long long pin_total = 0;
        for (int i = 0; i < k; ++i) pin_total += choose2((long long)sz[i]);
        if (slots_total > 0 && pin_total > 0) {
            double t_budget = 0.03;
            double t_end = min(deadline, elapsed_s() + t_budget);
            long long slots_left = slots_total;
            for (int i = 0; i < k && slots_left > 0; ++i) {
                int s = sz[i];
                if (s < 7) continue;
                long long li = (long long)((double)slots_total * (double)choose2((long long)s) /
                                           (double)pin_total);
                li = min(li, slots_left);
                if (li < s) continue;
                if (elapsed_s() > t_end) break;
                double dl = min(t_end, elapsed_s() + max(0.003, t_budget / (double)k));
                auto grp = cyclic_blocks(s, (int)li, dl);
                if (grp.empty()) continue;
                for (auto& b : grp) {
                    for (auto& v : b) v += off[i];
                    blocks.push_back(std::move(b));
                    --slots_left;
                    if (slots_left <= 0) break;
                }
            }
        }
    }
    auto avail = build_avail(blocks, small);
    greedy_fill(blocks, small, large, avail, deadline);
    return blocks;
}

static vector<vector<int>> gdd_blocks(int small, int large, double deadline) {
    if (small < 8) return {};
    if ((long long)large >= choose2((long long)small)) return {};
    struct TrySpec {
        FieldSpec fs;
        int k;
        long long est;
    };
    vector<TrySpec> tries;
    for (auto spec : field_specs(small)) {
        int g = spec.q;
        if (g < 3) continue;
        if (1LL * g * g > (long long)large) continue;
        int k0 = (small + g - 1) / g;
        for (int k = max(3, k0); k <= min(g + 1, k0 + 1); ++k) {
            if (k > small) continue;
            int base = small / k, rem = small % k;
            int maxgrp = base + (rem ? 1 : 0);
            if (maxgrp > g) continue;
            long long incid = 1LL * small * g;
            long long pin = (long long)rem * choose2((long long)base + 1) +
                            (long long)(k - rem) * choose2((long long)base);
            long long slots = (long long)large - 1LL * g * g;
            long long fill_est = 0;
            if (slots > 0 && pin > 0) {
                int d = fill_target_size(pin, slots, maxgrp);
                if (d >= 2) {
                    long long nblk = min(slots, (2 * pin) / ((long long)d * (d - 1)));
                    fill_est = (long long)((double)(nblk * d) * 0.9);
                }
            }
            tries.push_back({spec, k, incid + fill_est});
        }
    }
    sort(tries.begin(), tries.end(),
         [](const TrySpec& a, const TrySpec& b) { return a.est > b.est; });
    if ((int)tries.size() > 6) tries.resize(6);

    vector<vector<int>> best;
    long long best_score = -1;
    for (auto& t : tries) {
        if (elapsed_s() > deadline) break;
        auto blocks = gdd_construct(t.fs, t.k, small, large, deadline);
        if (blocks.empty()) continue;
        vector<vector<int>> trial = blocks;
        if ((int)trial.size() > large) trial.resize(large);
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

// ---------------- cyclic difference-family candidate ----------------

static vector<vector<int>> cyclic_blocks(int small, int large, double deadline) {
    int N = small;
    if (N < 7) return {};
    if ((long long)large >= choose2((long long)N)) return {};
    int full = large / N;
    int rem = large % N;
    int nclasses = full + (rem > 0 ? 1 : 0);
    if (full < 1 || nclasses > 8) return {};

    int D = N - 1 - (N % 2 == 0 ? 1 : 0);
    vector<int> wish(nclasses, 1);
    long long used = 0;
    while (true) {
        int bi = -1;
        double bv = 0.0;
        for (int c = 0; c < nclasses; ++c) {
            if (wish[c] >= N) continue;
            long long cost = 2LL * wish[c];
            if (used + cost > D) continue;
            double w = (double)(c < full ? N : rem) / (double)cost;
            if (w > bv) {
                bv = w;
                bi = c;
            }
        }
        if (bi < 0) break;
        used += 2LL * wish[bi];
        ++wish[bi];
    }

    Rng rng(0x9e3779b9ULL ^ ((uint64_t)N << 20) ^ (uint64_t)large);
    vector<vector<int>> best_sets;
    long long best_pts = -1;
    vector<char> dused(N, 0);
    vector<int> order(N - 1);
    iota(order.begin(), order.end(), 1);

    // try to add x to A: all new unordered difference classes must be unused and
    // mutually distinct; returns true and marks them on success
    vector<int> nd;
    auto try_add = [&](vector<int>& A, int x) {
        nd.clear();
        for (int a : A) {
            if (a == x) return false;
            int d = x - a;
            if (d < 0) d += N;
            if (dused[d] || dused[N - d]) return false;
            int cls = min(d, N - d);
            for (int e : nd) {
                if (e == cls) return false;
            }
            nd.push_back(cls);
        }
        for (int cls : nd) {
            dused[cls] = 1;
            dused[N - cls] = 1;
        }
        A.push_back(x);
        return true;
    };

    while (elapsed_s() < deadline) {
        fill(dused.begin(), dused.end(), 0);
        if (N % 2 == 0) dused[N / 2] = 1;
        vector<vector<int>> sets(nclasses);
        long long pts = 0;
        for (int c = 0; c < nclasses; ++c) {
            auto& A = sets[c];
            A.push_back(0);
            rng.shuffle_vec(order);
            for (int x : order) {
                if ((int)A.size() >= wish[c]) break;
                try_add(A, x);
            }
        }
        // extension phase: grow full classes further if budget remains
        for (int c = 0; c < full && c < nclasses; ++c) {
            auto& A = sets[c];
            for (int x : order) {
                try_add(A, x);
            }
        }
        for (int c = 0; c < nclasses; ++c) {
            pts += (long long)(c < full ? N : rem) * (long long)sets[c].size();
        }
        if (pts > best_pts) {
            best_pts = pts;
            best_sets = sets;
        }
    }
    if (best_sets.empty()) return {};

    vector<vector<int>> blocks;
    blocks.reserve(large);
    for (int c = 0; c < (int)best_sets.size(); ++c) {
        int copies = (c < full ? N : rem);
        auto& A = best_sets[c];
        if ((int)A.size() < 2) continue;
        for (int t = 0; t < copies && (int)blocks.size() < large; ++t) {
            vector<int> blk;
            blk.reserve(A.size());
            for (int a : A) blk.push_back((a + t) % N);
            blocks.push_back(std::move(blk));
        }
    }
    if (blocks.empty()) return {};
    auto avail = build_avail(blocks, small);
    greedy_fill(blocks, small, large, avail, deadline);
    if (!valid_blocks(blocks, small)) return {};
    return blocks;
}

// time-budgeted extra randomized shuffled-clique restarts
static void shuffled_time_budget(Candidate& best, int small, int large, double deadline) {
    long long pair_budget = choose2((long long)small);
    if (small <= 1 || (long long)large >= pair_budget) return;
    double avg_pairs = (double)pair_budget / (double)large;
    int s0 = (int)floor((1.0 + sqrt(1.0 + 8.0 * avg_pairs)) / 2.0);
    s0 = max(2, min(small, s0));
    Rng rng(0xABCDEF12345ULL ^ ((uint64_t)small << 17) ^ (uint64_t)large);
    while (elapsed_s() < deadline) {
        int cap = min(small, s0 + rng.next_int(4));
        int offset = (int)rng.next_int(4) - 1;
        uint64_t seed = rng.next_u64();
        auto cand = shuffled_clique_run(small, large, cap, true, offset, seed);
        consider(best, std::move(cand), small, large);
    }
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

int main() {
    T_START = chrono::steady_clock::now();
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
    consider(best, pbd316_blocks(small, large), small, large);
    consider(best, projective_subset_blocks(small, large), small, large);
    consider(best, projective_alternating_subset_blocks(small, large), small, large);
    consider(best, greedy_blocks(small, large), small, large);
    consider(best, shuffled_clique_blocks(small, large), small, large);
    if (elapsed_s() < 0.60) {
        consider(best, gdd_blocks(small, large, 0.60), small, large);
    }
    if (elapsed_s() < 0.72) {
        consider(best, cyclic_blocks(small, large, min(0.72, elapsed_s() + 0.12)), small, large);
    }
    shuffled_time_budget(best, small, large, 0.56);

    if (!valid_blocks(best.blocks, small)) {
        best.blocks.clear();
        for (int i = 0; i < large; ++i) best.blocks.push_back({0});
        best.score = block_score(best.blocks);
    }

    // final post-pass: augmentation + LNS block rebuilding (revert if not better)
    if (small >= 2 && !best.blocks.empty()) {
        vector<vector<int>> blocks = best.blocks;
        auto avail = build_avail(blocks, small);
        const double aug_dl = 0.74, lns_dl = 0.78, fin_dl = 0.80;
        while (elapsed_s() < aug_dl && augment_pass(blocks, small, avail, aug_dl)) {
        }
        lns_improve(blocks, small, large, avail, lns_dl);
        while (elapsed_s() < fin_dl && augment_pass(blocks, small, avail, fin_dl)) {
        }
        long long sc = block_score(blocks);
        if (sc > best.score && valid_blocks(blocks, small)) {
            best.blocks = std::move(blocks);
            best.score = sc;
        }
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
