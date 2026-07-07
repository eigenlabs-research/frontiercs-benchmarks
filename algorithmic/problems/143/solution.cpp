#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>
using namespace std;

struct Card {
    int s = -1;
    int v = -1;
};

static mt19937_64 rng(0x9e3779b97f4a7c15ULL);
static int match_hands = 0;

static int used_val[14];
static bool used_card[14][4];

struct Score {
    int type = 0;
    long long key = 0;
};

static Score make_score(int type, const vector<int>& vals) {
    Score s;
    s.type = type;
    for (int v : vals) s.key = s.key * 20 + v;
    return s;
}

static bool less_score(const Score& a, const Score& b) {
    if (a.type != b.type) return a.type < b.type;
    return a.key < b.key;
}

static Score eval7(const array<Card, 7>& cards) {
    memset(used_val, 0, sizeof(used_val));
    memset(used_card, 0, sizeof(used_card));
    for (const auto& c : cards) {
        used_val[c.v]++;
        used_card[c.v][c.s] = true;
    }
    used_val[0] = used_val[13];
    for (int s = 0; s < 4; ++s) used_card[0][s] = used_card[13][s];

    for (int suit = 0; suit < 4; ++suit) {
        for (int lo = 9; lo >= 0; --lo) {
            bool ok = true;
            for (int j = 0; j < 5; ++j) ok &= used_card[lo + j][suit];
            if (ok) {
                vector<int> a;
                for (int j = 4; j >= 0; --j) a.push_back(lo + j);
                return make_score(8, a);
            }
        }
    }
    for (int v = 13; v >= 1; --v) {
        if (used_val[v] == 4) {
            vector<int> a(4, v);
            for (int k = 13; k >= 1; --k) {
                if (k != v && used_val[k]) {
                    a.push_back(k);
                    return make_score(7, a);
                }
            }
        }
    }
    for (int t = 13; t >= 1; --t) {
        if (used_val[t] >= 3) {
            for (int p = 13; p >= 1; --p) {
                if (p != t && used_val[p] >= 2) {
                    vector<int> a = {t, t, t, p, p};
                    return make_score(6, a);
                }
            }
        }
    }
    for (int suit = 0; suit < 4; ++suit) {
        vector<int> a;
        for (int v = 13; v >= 1; --v) {
            if (used_card[v][suit]) a.push_back(v);
            if ((int)a.size() == 5) return make_score(5, a);
        }
    }
    for (int lo = 9; lo >= 0; --lo) {
        bool ok = true;
        for (int j = 0; j < 5; ++j) ok &= used_val[lo + j] > 0;
        if (ok) {
            vector<int> a;
            for (int j = 4; j >= 0; --j) a.push_back(lo + j);
            return make_score(4, a);
        }
    }
    for (int t = 13; t >= 1; --t) {
        if (used_val[t] >= 3) {
            vector<int> a = {t, t, t};
            for (int v = 13; v >= 1; --v) {
                if (v != t && used_val[v]) a.push_back(v);
                if ((int)a.size() == 5) return make_score(3, a);
            }
        }
    }
    for (int hi = 13; hi >= 1; --hi) {
        if (used_val[hi] >= 2) {
            for (int lo = hi - 1; lo >= 1; --lo) {
                if (used_val[lo] >= 2) {
                    vector<int> a = {hi, hi, lo, lo};
                    for (int k = 13; k >= 1; --k) {
                        if (k != hi && k != lo && used_val[k]) {
                            a.push_back(k);
                            return make_score(2, a);
                        }
                    }
                }
            }
        }
    }
    for (int p = 13; p >= 1; --p) {
        if (used_val[p] >= 2) {
            vector<int> a = {p, p};
            for (int v = 13; v >= 1; --v) {
                if (v != p && used_val[v]) a.push_back(v);
                if ((int)a.size() == 5) return make_score(1, a);
            }
        }
    }
    vector<int> a;
    for (int v = 13; v >= 1; --v) {
        if (used_val[v]) a.push_back(v);
        if ((int)a.size() == 5) return make_score(0, a);
    }
    return {};
}

static int compare_showdown(const Card alice[2], const Card bob[2], const array<Card, 5>& board) {
    array<Card, 7> a = {alice[0], alice[1], board[0], board[1], board[2], board[3], board[4]};
    array<Card, 7> b = {bob[0], bob[1], board[0], board[1], board[2], board[3], board[4]};
    Score sa = eval7(a), sb = eval7(b);
    if (less_score(sa, sb)) return -1;
    if (less_score(sb, sa)) return 1;
    return 0;
}

static int code_card(const Card& c) {
    return (c.v - 1) * 4 + c.s;
}

static Card decode_card(int x) {
    return Card{x % 4, x / 4 + 1};
}

static vector<int> deck_without(const vector<Card>& known) {
    bool used[52] = {};
    for (auto c : known) {
        if (c.s >= 0) used[code_card(c)] = true;
    }
    vector<int> d;
    d.reserve(52 - known.size());
    for (int i = 0; i < 52; ++i) {
        if (!used[i]) d.push_back(i);
    }
    return d;
}

static void draw_prefix(vector<int>& d, int need) {
    for (int i = 0; i < need; ++i) {
        uniform_int_distribution<int> dist(i, (int)d.size() - 1);
        int j = dist(rng);
        swap(d[i], d[j]);
    }
}

static uint64_t q_key(int round, const vector<Card>& board, Card b0, Card b1) {
    int c0 = code_card(b0), c1 = code_card(b1);
    if (c1 < c0) swap(c0, c1);
    vector<int> bs;
    for (auto c : board) bs.push_back(code_card(c));
    sort(bs.begin(), bs.end());
    uint64_t key = (uint64_t)round;
    key = key * 64 + c0;
    key = key * 64 + c1;
    key = key * 8 + bs.size();
    for (int x : bs) key = key * 64 + x;
    return key;
}

static unordered_map<uint64_t, double> bob_equity_cache;

static double estimate_random_alice_equity_for_bob(int round, const vector<Card>& board, Card b0, Card b1) {
    uint64_t key = q_key(round, board, b0, b1);
    auto it = bob_equity_cache.find(key);
    if (it != bob_equity_cache.end()) return it->second;

    vector<Card> known = {b0, b1};
    known.insert(known.end(), board.begin(), board.end());
    vector<int> base = deck_without(known);
    int missing_board = 5 - (int)board.size();

    int trials;
    bool large_match = match_hands >= 8000;
    if (round == 1) trials = 140;
    else if (round == 2) trials = 40;
    else if (round == 3) trials = large_match ? 62 : 72;
    else trials = large_match ? 92 : 110;

    double sum = 0.0;
    int total = 0;

    auto score_one = [&](Card a0, Card a1, array<Card, 5> full_board) {
        Card aa[2] = {a0, a1};
        Card bb[2] = {b0, b1};
        int cmp = compare_showdown(aa, bb, full_board);
        if (cmp > 0) sum += 1.0;
        else if (cmp == 0) sum += 0.5;
        total++;
    };

    for (int t = 0; t < trials; ++t) {
        vector<int> d = base;
        draw_prefix(d, 2 + missing_board);
        Card a0 = decode_card(d[0]);
        Card a1 = decode_card(d[1]);
        array<Card, 5> full_board;
        for (int i = 0; i < (int)board.size(); ++i) full_board[i] = board[i];
        for (int i = 0; i < missing_board; ++i) full_board[board.size() + i] = decode_card(d[2 + i]);
        score_one(a0, a1, full_board);
    }

    double q = total ? sum / total : 0.5;
    if (bob_equity_cache.size() < 250000) bob_equity_cache.emplace(key, q);
    return q;
}

struct World {
    Card bob[2];
    array<Card, 5> board;
    int outcome;
    double bob_q;
};

static double smooth_fold_from_q(double q, int pot, int x, double sigma = 0.030) {
    double threshold = (double)(pot + x) / (double)(pot + 2 * x);
    double fold_prob = 0.5 * erfc((threshold - q) / (sqrt(2.0) * sigma));
    return min(1.0, max(0.0, fold_prob));
}

static vector<World> sample_worlds(const Card alice[2], const vector<Card>& board, int round, int samples) {
    vector<Card> known = {alice[0], alice[1]};
    known.insert(known.end(), board.begin(), board.end());
    vector<int> base = deck_without(known);
    int missing_board = 5 - (int)board.size();
    vector<World> worlds;
    worlds.reserve(samples);
    for (int t = 0; t < samples; ++t) {
        vector<int> d = base;
        draw_prefix(d, 2 + missing_board);
        World w;
        w.bob[0] = decode_card(d[0]);
        w.bob[1] = decode_card(d[1]);
        for (int i = 0; i < (int)board.size(); ++i) w.board[i] = board[i];
        for (int i = 0; i < missing_board; ++i) w.board[board.size() + i] = decode_card(d[2 + i]);
        w.outcome = compare_showdown(alice, w.bob, w.board);
        w.bob_q = estimate_random_alice_equity_for_bob(round, board, w.bob[0], w.bob[1]);
        worlds.push_back(w);
    }
    return worlds;
}

static double showdown_delta_from_outcome(int outcome, int a, int pot) {
    if (outcome > 0) return (a + pot) - 100;
    if (outcome == 0) return (a + pot / 2.0) - 100;
    return a - 100;
}

static pair<int, long long> score_pair(const Score& s) {
    return {s.type, s.key};
}

static int compare_scores(const Score& a, const Score& b) {
    if (less_score(a, b)) return -1;
    if (less_score(b, a)) return 1;
    return 0;
}

static int choose_river_action(const Card alice[2], const vector<Card>& board_vec, int a, int pot) {
    array<Card, 5> board;
    for (int i = 0; i < 5; ++i) board[i] = board_vec[i];

    vector<Card> board_known = board_vec;
    vector<int> no_board = deck_without(board_known);
    vector<pair<int, long long>> all_scores;
    all_scores.reserve(no_board.size() * no_board.size() / 2);

    for (int i = 0; i < (int)no_board.size(); ++i) {
        for (int j = i + 1; j < (int)no_board.size(); ++j) {
            Card h[2] = {decode_card(no_board[i]), decode_card(no_board[j])};
            array<Card, 7> cards = {h[0], h[1], board[0], board[1], board[2], board[3], board[4]};
            all_scores.push_back(score_pair(eval7(cards)));
        }
    }
    sort(all_scores.begin(), all_scores.end());

    array<Card, 7> alice_cards = {alice[0], alice[1], board[0], board[1], board[2], board[3], board[4]};
    Score alice_score = eval7(alice_cards);

    vector<Card> known = {alice[0], alice[1]};
    known.insert(known.end(), board_vec.begin(), board_vec.end());
    vector<int> bob_deck = deck_without(known);

    struct RiverWorld {
        int outcome;
        double alice_win;
        double draw;
    };
    vector<RiverWorld> worlds;
    worlds.reserve(bob_deck.size() * bob_deck.size() / 2);
    double check_ev = 0.0;
    double equity = 0.0;

    for (int i = 0; i < (int)bob_deck.size(); ++i) {
        for (int j = i + 1; j < (int)bob_deck.size(); ++j) {
            Card b0 = decode_card(bob_deck[i]);
            Card b1 = decode_card(bob_deck[j]);
            array<Card, 7> bob_cards = {b0, b1, board[0], board[1], board[2], board[3], board[4]};
            Score bob_score = eval7(bob_cards);
            int outcome = compare_scores(alice_score, bob_score);
            auto sp = score_pair(bob_score);
            auto lo = lower_bound(all_scores.begin(), all_scores.end(), sp);
            auto hi = upper_bound(all_scores.begin(), all_scores.end(), sp);
            double greater = all_scores.end() - hi;
            double equal = hi - lo;
            double n = (double)all_scores.size();
            worlds.push_back({outcome, greater / n, equal / n});
            check_ev += showdown_delta_from_outcome(outcome, a, pot);
            if (outcome > 0) equity += 1.0;
            else if (outcome == 0) equity += 0.5;
        }
    }

    check_ev /= worlds.size();
    equity /= worlds.size();

    vector<int> cand = {1, 2, 3, 5, 8, 10, 13, 16, 20, 25, 32, 40, 50, 64, 80, 100, a};
    cand.push_back(max(1, pot / 2));
    cand.push_back(max(1, pot));
    cand.push_back(max(1, min(a, pot * 2)));
    cand.push_back(max(1, min(a, pot * 4)));
    sort(cand.begin(), cand.end());
    cand.erase(unique(cand.begin(), cand.end()), cand.end());

    int best_x = 0;
    double best_ev = check_ev;
    for (int x : cand) {
        if (x < 1 || x > a) continue;
        double threshold = (double)(pot + x) / (double)(pot + 2 * x);
        double ev = 0.0;
        int folds = 0;
        for (const auto& w : worlds) {
            double lose = w.alice_win;
            double draw = w.draw;
            double win = max(0.0, 1.0 - lose - draw);
            double v_lose = -x;
            double v_draw = pot / 2.0;
            double v_win = pot + x;
            double mean = lose * v_lose + draw * v_draw + win * v_win;
            double second = lose * v_lose * v_lose + draw * v_draw * v_draw + win * v_win * v_win;
            double var = max(0.0, second - mean * mean);
            double fold_prob;
            if (var < 1e-9) {
                fold_prob = mean <= 0.0 ? 1.0 : 0.0;
            } else {
                double z = (0.0 - mean) / sqrt(var / 100.0);
                fold_prob = 0.5 * erfc(-z / sqrt(2.0));
            }
            fold_prob = min(1.0, max(0.0, fold_prob));
            double fold_delta = a + pot - 100;
            double call_delta = showdown_delta_from_outcome(w.outcome, a - x, pot + 2 * x);
            ev += fold_prob * fold_delta + (1.0 - fold_prob) * call_delta;
            if (fold_prob >= 0.5) folds++;
        }
        ev /= worlds.size();

        double fold_rate = (double)folds / worlds.size();
        double margin = 0.45;
        if (equity > 0.62) margin -= 0.25;
        if (equity < 0.38 && fold_rate < 0.45) margin += 0.45;
        if (x >= a && equity < 0.44 && fold_rate < 0.62) margin += 0.80;
        if (ev > best_ev + margin) {
            best_ev = ev;
            best_x = x;
        }
    }
    return best_x;
}

static int choose_action(const Card alice[2], const vector<Card>& board, int round, int a, int pot) {
    if (round == 1) {
        if (alice[0].v == alice[1].v && alice[0].v >= 12) return min(a, 18);
        return 0;
    }
    if (round == 4) return choose_river_action(alice, board, a, pot);

    int samples;
    bool large_match = match_hands >= 8000;
    if (round == 1) samples = 28;
    else if (round == 2) samples = 34;
    else if (round == 3) samples = large_match ? 52 : 64;
    else samples = large_match ? 76 : 96;

    vector<World> worlds = sample_worlds(alice, board, round, samples);
    double check_ev = 0.0;
    double equity = 0.0;
    for (const auto& w : worlds) {
        check_ev += showdown_delta_from_outcome(w.outcome, a, pot);
        if (w.outcome > 0) equity += 1.0;
        else if (w.outcome == 0) equity += 0.5;
    }
    check_ev /= worlds.size();
    equity /= worlds.size();

    vector<int> cand = {1, 2, 3, 5, 8, 10, 13, 16, 20, 25, 32, 40, 50, 64, 80, 100, a};
    cand.push_back(max(1, pot / 2));
    cand.push_back(max(1, pot));
    cand.push_back(max(1, min(a, pot * 2)));
    cand.push_back(max(1, min(a, pot * 4)));
    sort(cand.begin(), cand.end());
    cand.erase(unique(cand.begin(), cand.end()), cand.end());

    int best_x = 0;
    double best_ev = check_ev;

    for (int x : cand) {
        if (x < 1 || x > a) continue;
        double threshold = (double)(pot + x) / (double)(pot + 2 * x);
        double ev = 0.0;
        int folds = 0;
        for (const auto& w : worlds) {
            double fold_prob = smooth_fold_from_q(w.bob_q, pot, x);
            double fold_delta = a + pot - 100;
            double call_delta = showdown_delta_from_outcome(w.outcome, a - x, pot + 2 * x);
            ev += fold_prob * fold_delta + (1.0 - fold_prob) * call_delta;
            if (fold_prob >= 0.5) folds++;
        }
        ev /= worlds.size();

        double fold_rate = (double)folds / worlds.size();
        double margin = 0.45;
        if (round == 1) margin = 1000000.00;
        else if (round == 2) margin = 20.00;
        else if (round == 3) margin = 1.50;

        if (equity > 0.62) margin -= 0.25;
        if (equity < 0.38 && fold_rate < 0.45) margin += 0.45;
        if (x >= a && equity < 0.44 && fold_rate < 0.62) margin += 0.80;

        if (ev > best_ev + margin) {
            best_ev = ev;
            best_x = x;
        }
    }

    return best_x;
}

static bool read_token(string& tok) {
    if (!(cin >> tok)) return false;
    if (tok == "-1") exit(0);
    return true;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int G;
    if (!(cin >> G)) return 0;
    match_hands = G;

    for (int hand = 1; hand <= G; ++hand) {
        bool ended = false;
        while (!ended) {
            string tok;
            if (!read_token(tok)) return 0;
            if (tok == "SCORE") return 0;
            if (tok != "STATE") return 0;

            int h, r, a, b, pot, k;
            cin >> h >> r >> a >> b >> pot >> k;

            string tag;
            Card alice[2];
            cin >> tag >> alice[0].s >> alice[0].v >> alice[1].s >> alice[1].v;

            cin >> tag;
            vector<Card> board(k);
            for (int i = 0; i < k; ++i) cin >> board[i].s >> board[i].v;

            uint64_t seed = 1469598103934665603ULL;
            auto mix = [&](uint64_t x) {
                seed ^= x + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
            };
            mix(h); mix(r); mix(a); mix(pot);
            mix(code_card(alice[0])); mix(code_card(alice[1]));
            for (auto c : board) mix(code_card(c));
            rng.seed(seed);

            int raise = choose_action(alice, board, r, a, pot);
            if (raise <= 0) {
                cout << "ACTION CHECK" << endl;
            } else {
                cout << "ACTION RAISE " << raise << endl;
            }

            string opp;
            if (!read_token(opp)) return 0;
            if (opp != "OPP") return 0;
            string oa;
            cin >> oa;
            if (oa == "CALL") {
                int x;
                cin >> x;
            }

            if (oa == "FOLD" || r == 4) {
                string res;
                if (!read_token(res)) return 0;
                if (res != "RESULT") return 0;
                int delta;
                cin >> delta;
                ended = true;
            }
        }
    }

    string tok;
    if (cin >> tok && tok == "SCORE") {
        double w;
        cin >> w;
    }
    return 0;
}
