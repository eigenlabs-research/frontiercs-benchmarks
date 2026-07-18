#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int mode = 0;
    if (!(cin >> mode)) return 0;

    const int N = 12;
    const vector<string> templ = {
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

    auto is_clue = [&](int r, int c) { return templ[r][c] == '?'; };

    vector<string> ans = templ;
    for (int r = 0; r < N; ++r)
        for (int c = 0; c < N; ++c)
            if (ans[r][c] == '?') ans[r][c] = '1';

    if (mode == 0) {
        // Clues induced by the boundary of a connected polyomino: all cells to
        // the left of a two-plateau vertical cut.  Compared with the plain
        // outside border, this moves a long forced boundary through the dense
        // clue columns (7/8/9) while retaining useful zero clues elsewhere.
        // The boundary of this simply connected set is a concrete valid loop.
        auto inS = [&](int r, int c) -> bool {
            if (r < 0 || r >= N || c < 0 || c >= N) return false;
            int cut = (r <= 6 ? 7 : 8); // cells 0..cut are inside
            return c <= cut;
        };
        const int dr[4] = {-1, 1, 0, 0};
        const int dc[4] = {0, 0, -1, 1};
        for (int r = 0; r < N; ++r) {
            for (int c = 0; c < N; ++c) if (is_clue(r, c)) {
                bool here = inS(r, c);
                int d = 0;
                for (int k = 0; k < 4; ++k)
                    if (here != inS(r + dr[k], c + dc[k])) ++d;
                ans[r][c] = char('0' + d);
            }
        }
    } else {
        // Use the boundary of a one-cell-wide horizontal snake polyomino.
        // Every row of clues is adjacent to the snake, so no clue is zero;
        // all produced values are in {1,2,3} and the boundary itself is a
        // concrete valid loop for the printed puzzle.
        auto inS = [&](int r, int c) -> bool {
            if (r < 0 || r >= N || c < 0 || c >= N) return false;
            if (r % 2 == 0) return true;
            if (r == 1 && c == 11) return true;
            if (r == 3 && c == 0) return true;
            if (r == 5 && c == 11) return true;
            if (r == 7 && c == 0) return true;
            if (r == 9 && c == 11) return true;
            return false;
        };
        const int dr[4] = {-1, 1, 0, 0};
        const int dc[4] = {0, 0, -1, 1};
        for (int r = 0; r < N; ++r) {
            for (int c = 0; c < N; ++c) if (is_clue(r, c)) {
                bool here = inS(r, c);
                int d = 0;
                for (int k = 0; k < 4; ++k) {
                    int nr = r + dr[k], nc = c + dc[k];
                    if (here != inS(nr, nc)) ++d;
                }
                if (d < 1) d = 1;
                if (d > 3) d = 3;
                ans[r][c] = char('0' + d);
            }
        }
    }

    for (const string &row : ans) cout << row << '\n';
    return 0;
}
