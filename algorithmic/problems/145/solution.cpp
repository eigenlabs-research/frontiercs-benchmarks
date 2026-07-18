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

    vector<string> out = templ;

    if (mode == 0) {
        // Conservative repair for the 0-allowed case: encode the perimeter loop
        // of the whole 12x12 board.  This gives an explicit known valid loop:
        // each clue is just the number of outside-border edges around that cell.
        for (int r = 0; r < N; ++r) {
            for (int c = 0; c < N; ++c) {
                if (templ[r][c] != '?') continue;
                int cnt = 0;
                if (r == 0) cnt++;
                if (r == N - 1) cnt++;
                if (c == 0) cnt++;
                if (c == N - 1) cnt++;
                out[r][c] = char('0' + cnt);
            }
        }
    } else {
        // For the 1-only mode, use a one-cell-wide serpentine polyomino and
        // output the boundary counts of its single outer boundary loop.  The
        // snake runs on every other row and is connected alternately at the two
        // sides, so all template clues touch the boundary and receive 1..3.
        bool in[N][N] = {};
        for (int r = 0; r < N; ++r) {
            for (int c = 0; c < N; ++c) {
                if (r % 2 == 0 && r <= 10) in[r][c] = true;
                if (r == 1 || r == 5 || r == 9) in[r][c] = (c == 11);
                if (r == 3 || r == 7) in[r][c] = (c == 0);
            }
        }
        auto inside = [&](int r, int c) -> bool {
            if (r < 0 || r >= N || c < 0 || c >= N) return false;
            return in[r][c];
        };
        for (int r = 0; r < N; ++r) {
            for (int c = 0; c < N; ++c) {
                if (templ[r][c] != '?') continue;
                int cnt = 0;
                cnt += (inside(r, c) != inside(r - 1, c));
                cnt += (inside(r, c) != inside(r + 1, c));
                cnt += (inside(r, c) != inside(r, c - 1));
                cnt += (inside(r, c) != inside(r, c + 1));
                if (cnt < 1) cnt = 1;
                if (cnt > 3) cnt = 3;
                out[r][c] = char('0' + cnt);
            }
        }
    }

    for (const string &s : out) cout << s << '\n';
    return 0;
}
