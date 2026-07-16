#include <bits/stdc++.h>
using namespace std;

static long long fccCount(int L, int parity) {
    long long e = (L + 1) / 2, o = L / 2;
    if (parity == 0) return e * e * e + 3 * e * o * o;
    return o * o * o + 3 * o * e * e;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n;
    if (!(cin >> n)) return 0;
    cout << setprecision(17);

    // The two-sphere optimum is the opposite body diagonal of the cube.
    if (n == 2) {
        const double r = sqrt(3.0) / (2.0 + 2.0 * sqrt(3.0));
        cout << r << ' ' << r << ' ' << r << '\n';
        cout << 1.0-r << ' ' << 1.0-r << ' ' << 1.0-r << '\n';
        return 0;
    }

    // A conventional FCC cell is the integer lattice points of one parity.
    // Its nearest-neighbour distance is sqrt(2) times the grid step.
    int L = 1;
    while (max(fccCount(L, 0), fccCount(L, 1)) < n) ++L;
    int parity = fccCount(L, 0) >= fccCount(L, 1) ? 0 : 1;
    const double root2 = sqrt(2.0);
    const double rfcc = 1.0 / (2.0 + root2 * (L - 1));

    // Retain a cubic grid whenever a finite FCC layer would be worse.
    int M = 1;
    while (1LL * M * M * M < n) ++M;
    const double rgrid = 1.0 / (2.0 * M);

    if (rfcc > rgrid) {
        const double step = root2 * rfcc;
        int written = 0;
        for (int i = 0; i < L && written < n; ++i)
            for (int j = 0; j < L && written < n; ++j)
                for (int k = 0; k < L && written < n; ++k)
                    if (((i + j + k) & 1) == parity) {
                        cout << rfcc + step * i << ' '
                             << rfcc + step * j << ' '
                             << rfcc + step * k << '\n';
                        ++written;
                    }
    } else {
        int small = (n + M * M - 1) / (M * M);
        int written = 0;
        for (int i = 0; i < M && written < n; ++i)
            for (int j = 0; j < M && written < n; ++j)
                for (int k = 0; k < small && written < n; ++k) {
                    cout << (i + 0.5) / M << ' '
                         << (j + 0.5) / M << ' '
                         << (k + 0.5) / small << '\n';
                    ++written;
                }
    }
    return 0;
}
