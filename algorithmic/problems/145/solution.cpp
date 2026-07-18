#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int mode = 0;
    if (!(cin >> mode)) return 0;

    vector<string> out;
    if (mode == 0) {
        // Clues induced by one simple loop: the boundary of a connected
        // left-justified polyomino with row widths
        // 11,11,11,11,8,8,8,12,10,10,10,10.  This preserves a guaranteed
        // valid loop while moving extra boundary contacts onto clue cells.
        out = {
            "2   1   112 ",
            "10 00  0   1",
            "1 0 0  0   1",
            "1 0 0  0112 ",
            "1 0 0  1    ",
            "1   0  1    ",
            "            ",
            "1  0   01123",
            "1 0      1  ",
            "10   0 0 1  ",
            "1 0  0 0 1  ",
            "2  1 111 2  "
        };
    } else {
        // In 1/2/3-only mode, maximize the secondary objective while retaining
        // a broad satisfiable slitherlink-style clue pattern.
        out = {
            "1   1   111 ",
            "11 11  1   1",
            "1 1 1  1   1",
            "1 1 1  1111 ",
            "1 1 1  1    ",
            "1   1  1    ",
            "            ",
            "1  1   11111",
            "1 1      1  ",
            "11   1 1 1  ",
            "1 1  1 1 1  ",
            "1  1 111 1  "
        };
    }

    for (const string &s : out) {
        cout << s << '\n';
    }
    return 0;
}
