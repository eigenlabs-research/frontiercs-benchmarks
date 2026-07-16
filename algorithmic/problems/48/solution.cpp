#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n;
    if (!(cin >> n)) return 0;
    cout << setprecision(17);

    // Two opposite corners (pulled equally away from every face) are better
    // than either lattice construction for this exceptional case.
    if (n == 2) {
        long double r = sqrtl(3.0L) / (2.0L * (1.0L + sqrtl(3.0L)));
        double a = (double)r, b = (double)(1.0L - r);
        cout << a << ' ' << a << ' ' << a << '\n';
        cout << b << ' ' << b << ' ' << b << '\n';
        return 0;
    }

    int gridM = 1;
    while (1LL * gridM * gridM * gridM < n) ++gridM;
    long double gridR = 1.0L / (2.0L * gridM);

    // Integer triples of even coordinate sum form an FCC lattice.  Its
    // nearest points differ by sqrt(2)*s; choosing margin s/sqrt(2) makes
    // both the face and pair constraints tight.
    int fccM = 1;
    auto fccCount = [](long long m) { return (m * m * m + 1) / 2; };
    while (fccCount(fccM) < n) ++fccM;
    long double fccR = 1.0L / (sqrtl(2.0L) * (fccM - 1) + 2.0L);

    if (fccR > gridR) {
        long double s = 1.0L / (fccM - 1 + sqrtl(2.0L));
        long double margin = s / sqrtl(2.0L);
        int made = 0;
        for (int x = 0; x < fccM && made < n; ++x)
            for (int y = 0; y < fccM && made < n; ++y)
                for (int z = 0; z < fccM && made < n; ++z)
                    if (((x + y + z) & 1) == 0) {
                        cout << (double)(margin + x * s) << ' '
                             << (double)(margin + y * s) << ' '
                             << (double)(margin + z * s) << '\n';
                        ++made;
                    }
    } else {
        int made = 0;
        for (int x = 0; x < gridM && made < n; ++x)
            for (int y = 0; y < gridM && made < n; ++y)
                for (int z = 0; z < gridM && made < n; ++z) {
                    cout << (double)((x + 0.5L) / gridM) << ' '
                         << (double)((y + 0.5L) / gridM) << ' '
                         << (double)((z + 0.5L) / gridM) << '\n';
                    ++made;
                }
    }
    return 0;
}
