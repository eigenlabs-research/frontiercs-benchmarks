#include <bits/stdc++.h>
using namespace std;

struct Candidate {
    double r = -1;
    vector<array<double, 3>> p;
};

static Candidate make_candidate(int n, int type) {
    // type 0: Z^3, type 1: FCC (coordinate sum even), type 2: BCC
    int d = 1;
    auto capacity = [&](int q) -> long long {
        long long Q = q;
        if (type == 0) return Q * Q * Q;
        if (type == 1) return (Q * Q * Q + 1) / 2;
        long long even = (q + 1) / 2, odd = q / 2;
        return even * even * even + odd * odd * odd;
    };
    while (capacity(d) < n) ++d;

    vector<array<int, 3>> a;
    a.reserve(n);
    for (int x = 0; x < d && (int)a.size() < n; ++x)
        for (int y = 0; y < d && (int)a.size() < n; ++y)
            for (int z = 0; z < d && (int)a.size() < n; ++z) {
                bool take;
                if (type == 0) take = true;
                else if (type == 1) take = ((x + y + z) & 1) == 0;
                else take = ((x & 1) == (y & 1) && (y & 1) == (z & 1));
                if (take) a.push_back({x, y, z});
            }

    int lo[3] = {INT_MAX, INT_MAX, INT_MAX};
    int hi[3] = {INT_MIN, INT_MIN, INT_MIN};
    for (auto q : a) for (int k = 0; k < 3; ++k) {
        lo[k] = min(lo[k], q[k]); hi[k] = max(hi[k], q[k]);
    }
    int span = max({hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]});
    long long mind2 = LLONG_MAX;
    for (int i = 0; i < n; ++i) for (int j = 0; j < i; ++j) {
        long long dx = a[i][0] - a[j][0], dy = a[i][1] - a[j][1], dz = a[i][2] - a[j][2];
        mind2 = min(mind2, dx * dx + dy * dy + dz * dz);
    }
    double md = sqrt((double)mind2);
    // With lattice scale s, r is min(s*md/2, (1-s*span)/2).
    // Center the occupied lattice bounding box and make the two bounds equal.
    double scale = 1.0 / (span + md);
    Candidate out;
    out.r = md * scale * 0.5;
    out.p.reserve(n);
    for (auto q : a)
        out.p.push_back({out.r + scale * (q[0] - lo[0]),
                         out.r + scale * (q[1] - lo[1]),
                         out.r + scale * (q[2] - lo[2])});
    return out;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n;
    if (!(cin >> n)) return 0;
    Candidate best;
    for (int type = 0; type < 3; ++type) {
        Candidate cur = make_candidate(n, type);
        if (cur.r > best.r) best = move(cur);
    }
    cout << setprecision(17);
    for (const auto &q : best.p) cout << q[0] << ' ' << q[1] << ' ' << q[2] << '\n';
}
