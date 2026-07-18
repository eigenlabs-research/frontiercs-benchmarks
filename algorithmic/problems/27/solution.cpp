#include <bits/stdc++.h>
using namespace std;

/*
  We view the smaller grid dimension as a set of "rows" and the larger
  dimension as columns.  A valid output is then a family of subsets B_j of the
  small side such that no unordered pair of small-side indices occurs in two
  different subsets.  This is exactly enough to forbid two grid rows sharing two
  selected columns.  The construction below greedily builds a large packing of
  such subsets; in very unbalanced cases it uses the optimal cheap pattern of
  one point in every column plus one distinct extra row-pair when possible.
*/

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n, m;
    if (!(cin >> n >> m)) return 0;

    bool transposed = false;
    int R = n, C = m;
    if (R > C) { swap(R, C); transposed = true; }

    vector<vector<int>> block(C + 1); // columns on the larger side, rows in [1..R]

    if (R == 1) {
        for (int c = 1; c <= C; ++c) block[c].push_back(1);
    } else {
        long long P = 1LL * R * (R - 1) / 2;

        // When there are at least as many columns as row-pairs, degree 2 is the
        // most pair-efficient way to add points after putting one point per col.
        if (C >= P) {
            int col = 1;
            for (int a = 1; a <= R; ++a) {
                for (int b = a + 1; b <= R; ++b) {
                    block[col].push_back(a);
                    block[col].push_back(b);
                    ++col;
                }
            }
            for (; col <= C; ++col) block[col].push_back(1);
        } else {
            // Target degree suggested by the pair-count upper bound:
            // C * d(d-1)/2 <= R(R-1)/2.
            int target = 1;
            while (target < R && 1LL * C * (target + 1) * target / 2 <= P) ++target;
            target = max(2, target);

            auto attempt = [&](uint32_t seed) {
                mt19937 rng(seed);
                vector<vector<int>> cur(C + 1);
                vector<vector<unsigned char>> used(R + 1, vector<unsigned char>(R + 1, 0));
                vector<int> deg(R + 1, 0), cols(C);
                iota(cols.begin(), cols.end(), 1);
                shuffle(cols.begin(), cols.end(), rng);

                auto can_add = [&](const vector<int>& v, int x) {
                    for (int y : v) {
                        int a = min(x, y), b = max(x, y);
                        if (used[a][b]) return false;
                    }
                    return true;
                };
                auto do_add = [&](vector<int>& v, int x) {
                    for (int y : v) {
                        int a = min(x, y), b = max(x, y);
                        used[a][b] = 1;
                    }
                    v.push_back(x);
                    deg[x]++;
                };

                vector<int> rows(R);
                iota(rows.begin(), rows.end(), 1);

                // First pass: try to give every column about 'target' rows.
                for (int c : cols) {
                    shuffle(rows.begin(), rows.end(), rng);
                    while ((int)cur[c].size() < target) {
                        int best = -1, bestScore = INT_MAX;
                        for (int x : rows) {
                            bool already = false;
                            for (int y : cur[c]) if (y == x) { already = true; break; }
                            if (already || !can_add(cur[c], x)) continue;
                            int score = deg[x] * 1000 + int(rng() % 1000);
                            if (score < bestScore) bestScore = score, best = x;
                        }
                        if (best == -1) break;
                        do_add(cur[c], best);
                    }
                    if (cur[c].empty()) do_add(cur[c], int(rng() % R) + 1);
                }

                // Improvement pass: greedily spend remaining unused pairs where possible.
                bool changed = true;
                int rounds = 0;
                while (changed && rounds++ < 3) {
                    changed = false;
                    shuffle(cols.begin(), cols.end(), rng);
                    for (int c : cols) {
                        int best = -1, bestScore = INT_MAX;
                        for (int x = 1; x <= R; ++x) {
                            bool already = false;
                            for (int y : cur[c]) if (y == x) { already = true; break; }
                            if (already || !can_add(cur[c], x)) continue;
                            int score = deg[x] * 1000 + int(rng() % 1000);
                            if (score < bestScore) bestScore = score, best = x;
                        }
                        if (best != -1) {
                            do_add(cur[c], best);
                            changed = true;
                        }
                    }
                }
                return cur;
            };

            long long bestK = -1;
            vector<vector<int>> best;
            int tries = (R * C <= 30000 ? 13 : 9);
            for (int t = 0; t < tries; ++t) {
                auto cand = attempt(1234567u + 1009u * t + 17u * R + 31u * C);
                long long k = 0;
                for (int c = 1; c <= C; ++c) k += cand[c].size();
                if (k > bestK) { bestK = k; best.swap(cand); }
            }
            block.swap(best);
        }
    }

    // A second, structurally different candidate: finite-field lines.  Extend
    // the previous prime-field affine construction to the small prime powers
    // that fit under n*m<=100000, and use projective planes when the extra
    // q+1 rows are available.  Any two such lines meet in at most one point.
    auto count_points = [&](const vector<vector<int>>& v) {
        long long s = 0;
        for (int c = 1; c <= C; ++c) s += (int)v[c].size();
        return s;
    };
    auto field_tables = [&](int q) {
        vector<vector<int>> add(q, vector<int>(q)), mul(q, vector<int>(q));
        auto is_prime = [](int x) {
            if (x < 2) return false;
            for (int d = 2; d * d <= x; ++d) if (x % d == 0) return false;
            return true;
        };
        if (is_prime(q)) {
            for (int a = 0; a < q; ++a) for (int b = 0; b < q; ++b) {
                add[a][b] = (a + b) % q;
                mul[a][b] = (a * b) % q;
            }
        } else if (q == 9) {
            for (int a = 0; a < q; ++a) for (int b = 0; b < q; ++b) {
                int a0 = a % 3, a1 = a / 3, b0 = b % 3, b1 = b / 3;
                int c0 = (a0 * b0 + 2 * a1 * b1) % 3; // x^2 = -1
                int c1 = (a0 * b1 + a1 * b0) % 3;
                add[a][b] = ((a0 + b0) % 3) + 3 * ((a1 + b1) % 3);
                mul[a][b] = c0 + 3 * c1;
            }
        } else { // binary fields: q=4,8,16 with x^k+x+1
            int poly = (q == 4 ? 0b111 : (q == 8 ? 0b1011 : 0b10011));
            int k = (q == 4 ? 2 : (q == 8 ? 3 : 4));
            for (int a = 0; a < q; ++a) for (int b = 0; b < q; ++b) {
                add[a][b] = a ^ b;
                int p = 0;
                for (int i = 0; i < k; ++i) if (b & (1 << i)) p ^= a << i;
                for (int i = 2 * k - 2; i >= k; --i) if (p & (1 << i)) p ^= poly << (i - k);
                mul[a][b] = p & (q - 1);
            }
        }
        return pair<vector<vector<int>>, vector<vector<int>>>(add, mul);
    };

    vector<int> orders = {2,3,4,5,7,8,9,11,13,16,17};
    vector<vector<int>> algebraic(C + 1);
    long long algebraicK = -1;
    for (int q : orders) {
        if (q * q > R) continue;
        auto tables = field_tables(q);
        auto &add = tables.first;
        auto &mul = tables.second;

        vector<vector<int>> cand(C + 1);
        int col = 1;
        bool projective = (q * q + q + 1 <= R);
        for (int a = 0; a < q && col <= C; ++a) {
            for (int b = 0; b < q && col <= C; ++b) {
                for (int x = 0; x < q; ++x) {
                    int y = add[mul[a][x]][b];
                    cand[col].push_back(x * q + y + 1);
                }
                if (projective) cand[col].push_back(q * q + a + 1);
                ++col;
            }
        }
        for (int x0 = 0; x0 < q && col <= C; ++x0) {
            for (int y = 0; y < q; ++y) cand[col].push_back(x0 * q + y + 1);
            if (projective) cand[col].push_back(q * q + q + 1);
            ++col;
        }
        if (projective && col <= C) {
            for (int a = 0; a < q; ++a) cand[col].push_back(q * q + a + 1);
            cand[col].push_back(q * q + q + 1);
            ++col;
        }
        for (; col <= C; ++col) cand[col].push_back(1);
        long long k = count_points(cand);
        if (k > algebraicK) { algebraicK = k; algebraic.swap(cand); }
    }
    if (algebraicK > count_points(block)) block.swap(algebraic);

    long long K = count_points(block);
    cout << K << '\n';
    for (int c = 1; c <= C; ++c) {
        for (int r : block[c]) {
            if (!transposed) cout << r << ' ' << c << '\n';
            else cout << c << ' ' << r << '\n';
        }
    }
    return 0;
}
