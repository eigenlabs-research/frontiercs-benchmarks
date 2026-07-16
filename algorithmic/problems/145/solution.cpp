#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int mode;
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

    /*
       Use one explicit loop: the boundary of a connected, hole-free polyomino.
       The polyomino is a horizontal snake consisting of every even row, joined
       at alternating ends.  Its boundary is therefore a single closed loop.
       Every cell of the board is either in the snake or adjacent to it, and no
       cell is completely surrounded by the same side, so all clue counts are in
       {1,2,3}; this is valid for both input modes.
    */
    auto inside = [&](int r, int c) -> bool {
        if (r % 2 == 0 && r <= 10) return true;
        if (r == 1 || r == 5 || r == 9) return c == 11;
        if (r == 3 || r == 7 || r == 11) return c == 0;
        return false;
    };

    vector<string> out(N, string(N, ' '));
    for (int r = 0; r < N; ++r) {
        for (int c = 0; c < N; ++c) {
            if (templ[r][c] != '?') continue;
            bool s = inside(r, c);
            int cnt = 0;
            const int dr[4] = {-1, 1, 0, 0};
            const int dc[4] = {0, 0, -1, 1};
            for (int k = 0; k < 4; ++k) {
                int nr = r + dr[k], nc = c + dc[k];
                bool t = (0 <= nr && nr < N && 0 <= nc && nc < N) ? inside(nr, nc) : false;
                if (t != s) ++cnt;
            }
            // The construction above yields 1..3 on all clue cells.
            if (mode == 1 && cnt == 0) cnt = 2; // defensive; should not occur
            out[r][c] = char('0' + cnt);
        }
    }

    for (const string &row : out) cout << row << '\n';
    return 0;
}
