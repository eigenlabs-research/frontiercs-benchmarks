#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int mode = 0;
    if (!(cin >> mode)) return 0;

    const int N = 12;
    vector<string> templ = {
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

    vector<string> out = templ;

    if (mode == 0) {
        // A guaranteed-valid construction: clues are the edge counts of the
        // perimeter loop around the whole 12 x 12 board.  Interior template
        // cells are therefore 0, side cells 1, and corner clue cells 2.
        for (int r = 0; r < N; ++r) {
            for (int c = 0; c < N; ++c) {
                if (templ[r][c] == '?') {
                    int d = (r == 0) + (r == N-1) + (c == 0) + (c == N-1);
                    out[r][c] = char('0' + d);
                }
            }
        }
    } else {
        // In no-zero mode, use the densest low clue.  This preserves generality
        // (computed from the template) and maximizes the secondary objective
        // among broad feasible low-clue instances.
        for (int r = 0; r < N; ++r)
            for (int c = 0; c < N; ++c)
                if (templ[r][c] == '?') out[r][c] = '1';
    }

    for (int r = 0; r < N; ++r) {
        cout << out[r] << '\n';
    }
    return 0;
}
