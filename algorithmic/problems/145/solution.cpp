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

    vector<string> out(N, string(N, ' '));

    if (mode == 0) {
        // Clues induced by the outer-border loop.  This uses many zeroes to
        // forbid interior activity in the zero-allowed version while retaining
        // a guaranteed valid single-loop witness.
        for (int r = 0; r < N; ++r) {
            for (int c = 0; c < N; ++c) if (templ[r][c] == '?') {
                int d = 0;
                if (r == 0) ++d;
                if (r == N - 1) ++d;
                if (c == 0) ++d;
                if (c == N - 1) ++d;
                out[r][c] = char('0' + d);
            }
        }
    } else {
        // A connected, hole-free thin polyomino: full odd rows joined at the
        // right edge.  Its boundary is one simple loop, and every template clue
        // cell is adjacent to this boundary, so all emitted clues are in 1..3.
        bool in[N][N] = {};
        for (int r : {1, 3, 5, 7, 9, 11})
            for (int c = 0; c < N; ++c) in[r][c] = true;
        for (int r : {2, 4, 6, 8, 10}) in[r][N - 1] = true;

        for (int r = 0; r < N; ++r) {
            for (int c = 0; c < N; ++c) if (templ[r][c] == '?') {
                int d = 0;
                const int dr[4] = {-1, 1, 0, 0};
                const int dc[4] = {0, 0, -1, 1};
                for (int k = 0; k < 4; ++k) {
                    int nr = r + dr[k], nc = c + dc[k];
                    bool nb = (0 <= nr && nr < N && 0 <= nc && nc < N) ? in[nr][nc] : false;
                    if (in[r][c] != nb) ++d;
                }
                if (d == 0) d = 2; // should not occur for this template; safety only
                out[r][c] = char('0' + d);
            }
        }
    }

    for (int r = 0; r < N; ++r) {
        for (int c = 0; c < N; ++c) if (templ[r][c] == ' ') out[r][c] = ' ';
        cout << out[r] << '\n';
    }
    return 0;
}
