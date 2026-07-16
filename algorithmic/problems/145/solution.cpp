#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int mode = 0;
    if (!(cin >> mode)) return 0;

    const int N = 12;
    vector<string> tmpl = {
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

    /*
      Construct clues from the boundary of one fixed simply-connected polyomino.
      The polyomino is a comb: four vertical one-cell-wide bars, joined through
      the clue-free middle row.  Every clue cell lies either on this polyomino or
      immediately next to it, so all induced clue counts are in {1,2,3}; this is
      valid in both input modes and guarantees at least the polyomino boundary as
      a legal non-empty single loop.
    */
    bool in[N][N] = {};
    for (int r = 0; r < N; ++r) {
        in[r][1] = true;
        in[r][4] = true;
        in[r][7] = true;
        in[r][10] = true;
    }
    for (int c = 1; c <= 10; ++c) in[6][c] = true;

    vector<string> out(N, string(N, ' '));
    const int dr[4] = {-1, 1, 0, 0};
    const int dc[4] = {0, 0, -1, 1};
    for (int r = 0; r < N; ++r) {
        for (int c = 0; c < N; ++c) {
            if (tmpl[r][c] != '?') continue;
            int cnt = 0;
            for (int k = 0; k < 4; ++k) {
                int nr = r + dr[k], nc = c + dc[k];
                bool here = in[r][c];
                bool there = (0 <= nr && nr < N && 0 <= nc && nc < N) ? in[nr][nc] : false;
                if (here != there) ++cnt;
            }
            // The chosen comb dominates all clue cells, hence cnt is never zero.
            // It also has no isolated clue cell with four boundary sides.
            if (mode == 1 && cnt == 0) cnt = 1; // defensive; not reached
            if (cnt > 3) cnt = 3;              // defensive; not reached
            out[r][c] = char('0' + cnt);
        }
    }

    for (int r = 0; r < N; ++r) {
        cout << out[r] << '\n';
    }
    return 0;
}
