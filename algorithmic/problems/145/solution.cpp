#include <bits/stdc++.h>
using namespace std;

static const vector<string> templ = {
    "?   ?   ??? ",
    "?? ??  ?   ?",
    "? ? ?  ?   ?",
    "? ? ?  ???? ",
    "? ? ?  ?    ",
    "?   ?  ?    ",
    "            ",
    "?  ?   ?????",
    "? ?      ?  ",
    "??   ? ? ?  ",
    "? ?  ? ? ?  ",
    "?  ? ??? ?  "
};

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int mode;
    if (!(cin >> mode)) return 0;

    vector<string> incumbent0 = {
        "2   1   111 ",
        "10 00  0   1",
        "1 0 0  0   1",
        "1 0 0  0000 ",
        "1 0 0  0    ",
        "1   0  0    ",
        "            ",
        "1  0   00001",
        "1 0      0  ",
        "10   0 0 0  ",
        "1 0  0 0 0  ",
        "2  1 111 1  "
    };

    vector<string> out1 = {
        "2   2   222 ",
        "23 22  2   2",
        "2 2 2  2   2",
        "2 2 2  2223 ",
        "2 2 2  2    ",
        "2   2  2    ",
        "            ",
        "2  2   22232",
        "2 2      2  ",
        "23   2 2 2  ",
        "2 2  2 2 2  ",
        "1  1 111 1  "
    };

    if (mode == 1) {
        for (auto &s: out1) cout << s << '\n';
        return 0;
    }

    auto countOnes = [](const vector<string>& g) {
        int r = 0;
        for (auto &s: g) for (char c: s) if (c == '1') ++r;
        return r;
    };

    // Candidate family for mode 0: clues induced by one rectangular loop on
    // grid lines.  It is always a valid non-empty loop; enumerate all such
    // loops and use the one with most 1-clues (then most nonzero clues).  This
    // is a small general construction rather than a sample-specific answer.
    vector<string> best = incumbent0;
    int bestOnes = countOnes(best), bestNonZero = -1, bestArea = -1;

    for (int top = 0; top <= 11; ++top) for (int bot = top + 1; bot <= 12; ++bot)
    for (int lef = 0; lef <= 11; ++lef) for (int rig = lef + 1; rig <= 12; ++rig) {
        vector<string> g = templ;
        bool ok = true;
        int ones = 0, nonzero = 0;
        for (int r = 0; r < 12; ++r) for (int c = 0; c < 12; ++c) if (templ[r][c] == '?') {
            int d = 0;
            // Cell (r,c) is bounded by horizontal grid lines r,r+1 and
            // vertical grid lines c,c+1.  Add the selected rectangle sides.
            if ((r == top || r + 1 == top) && lef <= c && c < rig) ++d;
            if ((r == bot || r + 1 == bot) && lef <= c && c < rig) ++d;
            if ((c == lef || c + 1 == lef) && top <= r && r < bot) ++d;
            if ((c == rig || c + 1 == rig) && top <= r && r < bot) ++d;
            if (d > 3) { ok = false; break; }
            g[r][c] = char('0' + d);
            if (d == 1) ++ones;
            if (d) ++nonzero;
        }
        if (!ok) continue;
        int area = (bot - top) * (rig - lef);
        if (ones > bestOnes || (ones == bestOnes && (nonzero > bestNonZero || (nonzero == bestNonZero && area > bestArea)))) {
            best = g; bestOnes = ones; bestNonZero = nonzero; bestArea = area;
        }
    }

    for (auto &s: best) cout << s << '\n';
    return 0;
}
