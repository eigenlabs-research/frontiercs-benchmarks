#include <bits/stdc++.h>
using namespace std;

struct Point { long long x, y; };
struct Node {
    int id, left = -1, right = -1, parent = -1, alive = 0;
    long long lx, rx, ly, ry;
};

int n;
vector<Point> p;
vector<int> ord;
vector<Node> tree;
vector<char> used, primeId;

int build(int lo, int hi, int depth, int parent) {
    if (lo >= hi) return -1;
    int mid = (lo + hi) >> 1;
    int axis = depth & 1;
    nth_element(ord.begin() + lo, ord.begin() + mid, ord.begin() + hi,
        [axis](int a, int b) {
            if (axis == 0) return p[a].x != p[b].x ? p[a].x < p[b].x : p[a].y < p[b].y;
            return p[a].y != p[b].y ? p[a].y < p[b].y : p[a].x < p[b].x;
        });
    int at = (int)tree.size();
    tree.push_back(Node());
    tree[at].id = ord[mid];
    tree[at].parent = parent;
    // Recursive calls may grow tree and reallocate it, so keep only indices
    // across them rather than references into the vector.
    int left = build(lo, mid, depth + 1, at);
    int right = build(mid + 1, hi, depth + 1, at);
    tree[at].left = left;
    tree[at].right = right;
    Node &q = tree[at];
    q.alive = hi - lo;
    q.lx = q.rx = p[q.id].x;
    q.ly = q.ry = p[q.id].y;
    for (int ch : {q.left, q.right}) if (ch != -1) {
        q.lx = min(q.lx, tree[ch].lx); q.rx = max(q.rx, tree[ch].rx);
        q.ly = min(q.ly, tree[ch].ly); q.ry = max(q.ry, tree[ch].ry);
    }
    return at;
}

inline double pointDist(int a, int b) {
    double dx = (double)p[a].x - p[b].x, dy = (double)p[a].y - p[b].y;
    return hypot(dx, dy);
}
inline double boxDist(int v, int id) {
    const Node &q = tree[v];
    double dx = 0, dy = 0;
    if (p[id].x < q.lx) dx = (double)q.lx - p[id].x;
    else if (p[id].x > q.rx) dx = (double)p[id].x - q.rx;
    if (p[id].y < q.ly) dy = (double)q.ly - p[id].y;
    else if (p[id].y > q.ry) dy = (double)p[id].y - q.ry;
    return dx * dx + dy * dy;
}

void nearest(int v, int from, int &answer, double &best2) {
    if (v == -1 || tree[v].alive == 0 || boxDist(v, from) >= best2) return;
    int id = tree[v].id;
    if (!used[id]) {
        double dx = (double)p[from].x - p[id].x, dy = (double)p[from].y - p[id].y;
        double d = dx * dx + dy * dy;
        if (d < best2 || (d == best2 && id < answer)) best2 = d, answer = id;
    }
    int a = tree[v].left, b = tree[v].right;
    double da = a == -1 ? numeric_limits<double>::infinity() : boxDist(a, from);
    double db = b == -1 ? numeric_limits<double>::infinity() : boxDist(b, from);
    if (da > db) swap(a, b);
    nearest(a, from, answer, best2);
    nearest(b, from, answer, best2);
}

void erasePoint(int id, int nodeOfId[]) {
    used[id] = 1;
    for (int v = nodeOfId[id]; v != -1; v = tree[v].parent) --tree[v].alive;
}

// Exact static k-nearest query.  Unlike the construction query this does not
// use `alive`: it supplies geometric candidate edges for later 2-opt moves.
vector<int> spatialNeighbors(int from, int want) {
    using Entry = pair<double, int>;
    priority_queue<Entry, vector<Entry>, greater<Entry>> pending;
    priority_queue<Entry> best;
    pending.push({boxDist(0, from), 0});
    while (!pending.empty()) {
        auto [bound, v] = pending.top(); pending.pop();
        if ((int)best.size() == want && bound > best.top().first) break;
        int id = tree[v].id;
        if (id != from) {
            double dx = (double)p[from].x - p[id].x;
            double dy = (double)p[from].y - p[id].y;
            double d = dx * dx + dy * dy;
            if ((int)best.size() < want) best.push({d, id});
            else if (d < best.top().first || (d == best.top().first && id < best.top().second)) {
                best.pop(); best.push({d, id});
            }
        }
        for (int ch : {tree[v].left, tree[v].right}) if (ch != -1) {
            double b = boxDist(ch, from);
            if ((int)best.size() < want || b <= best.top().first) pending.push({b, ch});
        }
    }
    vector<int> ans;
    while (!best.empty()) ans.push_back(best.top().second), best.pop();
    return ans;
}

double routeCost(const vector<int> &r) {
    double z = 0;
    for (int t = 1; t <= n; ++t) {
        double mult = (t % 10 == 0 && !primeId[r[t - 1]]) ? 1.1 : 1.0;
        z += mult * pointDist(r[t - 1], r[t]);
    }
    return z;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    if (!(cin >> n)) return 0;
    p.resize(n);
    for (auto &a : p) cin >> a.x >> a.y;
    primeId.assign(n, true);
    if (n > 0) primeId[0] = false;
    if (n > 1) primeId[1] = false;
    for (int i = 2; 1LL * i * i < n; ++i) if (primeId[i])
        for (int j = i * i; j < n; j += i) primeId[j] = false;

    ord.resize(n); iota(ord.begin(), ord.end(), 0);
    tree.reserve(n);
    int root = build(0, n, 0, -1);
    vector<int> nodeOfId(n);
    for (int i = 0; i < n; ++i) nodeOfId[tree[i].id] = i;
    // Spatial, rather than input-rank, neighbors are the candidate edges for
    // local search.  The input order is x-sorted, but it is not generally a
    // good proxy for proximity when y varies or clusters are present.
    const int SPATIAL_NEIGHBORS = 10;
    vector<vector<int>> close(n);
    for (int id = 0; id < n; ++id) {
        close[id] = spatialNeighbors(id, SPATIAL_NEIGHBORS);
        // Retain the incumbent's x-rank candidates, so this broader search
        // cannot lose an already available improving 2-opt move.
        for (int d = 1; d <= 4 && id + d < n; ++d) close[id].push_back(id + d);
    }

    used.assign(n, 0);
    vector<int> route; route.reserve(n + 1); route.push_back(0);
    erasePoint(0, nodeOfId.data());
    int cur = 0;
    for (int take = 1; take < n; ++take) {
        int nxt = -1;
        double best2 = numeric_limits<double>::infinity();
        nearest(root, cur, nxt, best2);
        // The tree contains every unused city, so nxt is always found.
        route.push_back(nxt);
        erasePoint(nxt, nodeOfId.data());
        cur = nxt;
    }
    route.push_back(0);

    // Traversing the same geometric cycle in the opposite direction can move
    // non-prime cities away from the tenth-step sources at no distance cost.
    vector<int> rev(n + 1); rev[0] = rev[n] = 0;
    for (int i = 1; i < n; ++i) rev[i] = route[n - i];
    if (routeCost(rev) < routeCost(route)) route.swap(rev);

    // A full 2-opt reversal keeps all its internal Euclidean edge lengths,
    // but not necessarily their carrot multipliers.  The two prefix systems
    // below evaluate that change exactly: oldPenalty indexes old step sources,
    // while reversePenalty indexes the source after an edge is traversed back.
    // Candidate new edges are from an exact spatial k-nearest graph.  This
    // catches crossings and incompatible cluster exits that x-rank pairs miss.
    // The previous fixed cap could terminate while this exact candidate descent
    // still had improving moves.  Spend more moves on small instances (where
    // each full candidate scan is cheap), but retain a conservative cap at the
    // largest scale for the 2.5 second limit.
    const int ROUNDS = n <= 30000 ? 80 : (n <= 80000 ? 40 : 20);
    for (int round = 0; round < ROUNDS; ++round) {
        vector<int> pos(n);
        // Keep city 0 at its rooted position, rather than overwriting it
        // with the duplicate closing occurrence at route[n].
        for (int k = 0; k < n; ++k) pos[route[k]] = k;
        vector<double> edge(n + 1);
        for (int k = 1; k <= n; ++k) edge[k] = pointDist(route[k - 1], route[k]);
        vector<double> oldPenalty(n + 1, 0.0);
        vector<vector<double>> reversePenalty(10, vector<double>(n + 1, 0.0));
        for (int k = 1; k <= n; ++k) {
            oldPenalty[k] = oldPenalty[k - 1];
            if (k % 10 == 0 && !primeId[route[k - 1]]) oldPenalty[k] += 0.1 * edge[k];
            for (int r = 0; r < 10; ++r) reversePenalty[r][k] = reversePenalty[r][k - 1];
            if (!primeId[route[k]]) reversePenalty[k % 10][k] += 0.1 * edge[k];
        }
        auto stepCost = [&](int t, int a, int b) {
            return (1.0 + ((t % 10 == 0 && !primeId[a]) ? 0.1 : 0.0)) * pointDist(a, b);
        };
        double bestDelta = 0.0;
        int bestI = -1, bestJ = -1;
        for (int u = 0; u < n; ++u) for (int v : close[u]) {
            if (v == 0) continue;
            int l = pos[u], r = pos[v];
            if (l > r) swap(l, r);
            // Make the candidate cities themselves adjacent: reversing
            // route[l+1..r] creates the edge route[l]-route[r].
            if (r - l < 2) continue;
            int i = l + 1, j = r;
            double before = stepCost(i, route[i - 1], route[i]) +
                            stepCost(j + 1, route[j], route[j + 1]);
            double after = stepCost(i, route[i - 1], route[j]) +
                           stepCost(j + 1, route[i], route[j + 1]);
            // For a new tenth step t, its source was formerly at index
            // i+j-t+1, hence all relevant old edge indices share this residue.
            int residue = (i + j + 1) % 10;
            double oldInside = oldPenalty[j] - oldPenalty[i];
            double newInside = reversePenalty[residue][j] - reversePenalty[residue][i];
            double delta = after - before + newInside - oldInside;
            if (delta < bestDelta) bestDelta = delta, bestI = i, bestJ = j;
        }
        if (bestI == -1) break;
        reverse(route.begin() + bestI, route.begin() + bestJ + 1);
    }

    cout << n + 1 << '\n';
    for (int id : route) cout << id << '\n';
}
