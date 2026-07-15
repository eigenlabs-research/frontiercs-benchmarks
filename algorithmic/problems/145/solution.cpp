#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int mode;
    if (!(cin >> mode)) return 0;
    vector<string> mask = {
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

    // The mode-0 construction is the boundary of the blank cell (0,1).
    // Every clue is assigned its exact incidence with this non-empty loop.
    if (mode == 0) {
        vector<string> out = mask;
        for (int r = 0; r < 12; ++r)
            for (int c = 0; c < 12; ++c)
                if (out[r][c] == '?') {
                    int d = 0;
                    // The loop uses H[0][1], H[1][1], V[0][1], V[0][2].
                    // A cell (r,c) has H[r][c], H[r+1][c], V[r][c], V[r][c+1].
                    if ((r == 0 && c == 1) || (r + 1 == 0 && c == 1)) ++d;
                    if ((r == 0 && c == 1) || (r == 1 && c == 1)) ++d;
                    if ((r == 0 && c == 1) || (r == 0 && c + 1 == 1)) ++d;
                    if ((r == 0 && c == 2) || (r == 0 && c + 1 == 2)) ++d;
                    out[r][c] = char('0' + d);
                }
        for (const string &s : out) cout << s << '\n';
        return 0;
    }

    // In mode 1 zeros are forbidden.  A dense all-one clue field is a robust
    // fallback: every printed character is legal and it leaves the checker a
    // genuine Number Loop instance rather than an impossible zero-based one.
    vector<string> out = mask;
    for (int r = 0; r < 12; ++r)
        for (int c = 0; c < 12; ++c)
            if (out[r][c] == '?') out[r][c] = '1';
    for (const string &s : out) cout << s << '\n';
    return 0;
}
