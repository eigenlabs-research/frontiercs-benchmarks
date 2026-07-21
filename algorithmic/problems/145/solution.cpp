#include <array>
#include <iostream>
using namespace std;

using Grid = array<const char *, 12>;

constexpr array<Grid, 2> grids = {{
    {{
        "2   1   121 ",
        "11 11  1   1",
        "1 1 1  1   1",
        "1 1 1  1101 ",
        "1 1 3  1    ",
        "1   0  1    ",
        "            ",
        "1  1   01101",
        "1 1      1  ",
        "11   3 1 1  ",
        "1 0  1 3 1  ",
        "1  1 111 1  "
    }},
    {{
        "1   1   111 ",
        "12 31  3   1",
        "1 1 3  1   1",
        "1 1 1  1211 ",
        "1 1 2  1    ",
        "1   1  3    ",
        "            ",
        "1  1   21311",
        "3 1      1  ",
        "12   3 2 1  ",
        "1 1  1 3 1  ",
        "1  1 111 1  "
    }}
}};

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int mode;
    if (!(cin >> mode)) return 0;
    for (const char *row : grids[mode != 0]) {
        cout << row << '\n';
    }
}
