#include <chrono>
#include <climits>
#include <cstdio>
#include <random>
#include <utility>
#include <vector>
using namespace std;

// v3 with occasional perturbation restart from best known solution.
int main() {
    int n, m;
    if (scanf("%d %d", &n, &m) != 2) return 0;

    vector<int> deg(n, 0);
    vector<pair<int,int>> edges(m);
    for (int i = 0; i < m; i++) {
        int u, v; scanf("%d %d", &u, &v);
        --u; --v;
        edges[i] = {u, v};
        deg[u]++; deg[v]++;
    }
    vector<int> adjStart(n + 1, 0);
    for (int v = 0; v < n; v++) adjStart[v + 1] = adjStart[v] + deg[v];
    vector<int> adj(2 * m);
    {
        vector<int> pos(adjStart.begin(), adjStart.end());
        for (auto &e : edges) {
            adj[pos[e.first]++] = e.second;
            adj[pos[e.second]++] = e.first;
        }
    }

    auto t0 = chrono::steady_clock::now();
    auto elapsed = [&] { return chrono::duration<double>(chrono::steady_clock::now() - t0).count(); };
    const double BUDGET = 0.90;

    mt19937 rng(0x9E3779B9u);

    vector<int> side(n, 0), best(n, 0);
    vector<int> crossCnt(n, 0), gain(n, 0);
    vector<int> tabuUntil(n, 0);
    long long bestCut = -1;

    auto rebuild = [&]() -> long long {
        long long cut = 0;
        for (int v = 0; v < n; v++) {
            int c = 0;
            for (int i = adjStart[v]; i < adjStart[v + 1]; i++)
                c += side[adj[i]] != side[v];
            crossCnt[v] = c;
            gain[v] = deg[v] - 2 * c;
            cut += c;
        }
        return cut / 2;
    };

    auto greedyInit = [&]() {
        vector<int> order(n);
        for (int v = 0; v < n; v++) order[v] = v;
        for (int v = n - 1; v > 0; v--) swap(order[v], order[rng() % (v + 1)]);
        vector<char> placed(n, 0);
        vector<int> net(n, 0);
        for (int v = 0; v < n; v++) side[v] = 0;
        for (int idx = 0; idx < n; idx++) {
            int v = order[idx];
            int s = (net[v] > 0) ? 0 : (net[v] < 0 ? 1 : (int)(rng() & 1));
            side[v] = s;
            placed[v] = 1;
            for (int i = adjStart[v]; i < adjStart[v + 1]; i++) {
                int w = adj[i];
                if (!placed[w]) net[w] += (s == 1 ? 1 : -1);
            }
        }
    };

    auto randomInit = [&]() {
        for (int v = 0; v < n; v++) side[v] = (int)(rng() & 1);
    };

    auto perturbBestInit = [&]() {
        side = best;
        int flips = max(1, n / 20);
        for (int i = 0; i < flips; i++) {
            int v = rng() % n;
            side[v] ^= 1;
        }
    };

    auto flip = [&](int v, long long &cut) {
        cut += gain[v];
        int sv = side[v] ^= 1;
        gain[v] = -gain[v];
        for (int i = adjStart[v]; i < adjStart[v + 1]; i++) {
            int w = adj[i];
            if (side[w] != sv) {
                crossCnt[w]++;
                gain[w] -= 2;
            } else {
                crossCnt[w]--;
                gain[w] += 2;
            }
        }
    };

    vector<int> ties;
    ties.reserve(64);
    const int CHECK_MASK = 1023;
    const int TABU_MIN = 20;
    const int TABU_MAX = 60;

    int restart = 0;
    while (elapsed() < BUDGET) {
        if (bestCut < 0) greedyInit();
        else {
            int r = restart % 6;
            if (r == 0 || r == 1) greedyInit();
            else if (r == 2 || r == 3) randomInit();
            else perturbBestInit();
        }
        restart++;

        long long cut = rebuild();
        for (int v = 0; v < n; v++) tabuUntil[v] = 0;
        if (cut > bestCut) { bestCut = cut; best = side; }

        long long localBest = cut;
        int iter = 0;
        int sinceImprove = 0;
        int stagnationLimit = 2 * n + 500;

        while (true) {
            if ((iter & CHECK_MASK) == 0 && elapsed() >= BUDGET) break;
            iter++;

            int bestGain = INT_MIN;
            ties.clear();
            for (int v = 0; v < n; v++) {
                bool allowed = tabuUntil[v] <= iter || (cut + gain[v] > bestCut);
                if (!allowed) continue;
                if (gain[v] > bestGain) {
                    bestGain = gain[v];
                    ties.clear();
                    ties.push_back(v);
                } else if (gain[v] == bestGain) {
                    ties.push_back(v);
                }
            }
            if (ties.empty()) break;

            int v = ties[rng() % ties.size()];
            flip(v, cut);
            int tenure = TABU_MIN + (int)(rng() % (TABU_MAX - TABU_MIN + 1));
            tabuUntil[v] = iter + tenure;

            if (cut > bestCut) { bestCut = cut; best = side; }
            if (cut > localBest) { localBest = cut; sinceImprove = 0; }
            else if (++sinceImprove >= stagnationLimit) break;
        }
    }

    if (bestCut < 0) best.assign(n, 0);
    for (int v = 0; v < n; v++) printf("%d%c", best[v], v == n - 1 ? '\n' : ' ');
    return 0;
}
