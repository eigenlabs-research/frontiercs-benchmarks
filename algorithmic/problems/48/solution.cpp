#include <bits/stdc++.h>
using namespace std;

// Integer points with i+j+k even form an FCC lattice.  Its nearest-neighbour
// distance is sqrt(2), versus 1 for the ordinary cubic lattice at equal pitch.
int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n;
    if (!(cin >> n)) return 0;

    // A finite FCC crop is not optimal for two spheres: its closest pair is
    // only a face diagonal.  Put the pair on the cube body diagonal instead.
    // Here the face clearance r equals half their distance.
    if (n == 2) {
        const double q = sqrt(3.0);
        const double r = q / (2.0 * (1.0 + q));
        cout << setprecision(17)
             << r << ' ' << r << ' ' << r << '\n'
             << 1.0-r << ' ' << 1.0-r << ' ' << 1.0-r << '\n';
        return 0;
    }

    auto fccCount = [](long long side) -> long long {
        // side is the number of integer coordinates on an edge (0..side-1).
        long long even = (side + 1) / 2, odd = side / 2;
        return even * even * even + 3 * even * odd * odd;
    };

    long long side = 1;
    while (fccCount(side) < n) ++side;

    vector<array<int, 3>> p;
    p.reserve(n);
    for (int x = 0; x < side && (int)p.size() < n; ++x)
        for (int y = 0; y < side && (int)p.size() < n; ++y)
            for (int z = 0; z < side && (int)p.size() < n; ++z)
                if (((x + y + z) & 1) == 0)
                    p.push_back({x, y, z});

    int lo[3] = {INT_MAX, INT_MAX, INT_MAX};
    int hi[3] = {INT_MIN, INT_MIN, INT_MIN};
    for (const auto &a : p)
        for (int d = 0; d < 3; ++d) {
            lo[d] = min(lo[d], a[d]);
            hi[d] = max(hi[d], a[d]);
        }
    int span = max({hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]});

    // Nearest FCC sites are sqrt(2) pitches apart.  Leaving sqrt(2)/2
    // pitches at either face makes the face and pair constraints equal.
    const double root2 = sqrt(2.0);
    const double pitch = 1.0 / (span + root2);
    const double margin = 0.5 * root2 * pitch;

    cout << setprecision(17);
    for (const auto &a : p) {
        cout << margin + pitch * (a[0] - lo[0]) << ' '
             << margin + pitch * (a[1] - lo[1]) << ' '
             << margin + pitch * (a[2] - lo[2]) << '\n';
    }
    return 0;
}
