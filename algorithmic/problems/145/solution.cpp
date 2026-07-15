#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int mode;
    if (!(cin >> mode)) return 0;

    // These are independently constructed: mode 0 can exploit zero-clues,
    // while mode 1 must keep every clue positive.
    const array<string, 12> with_zero = {{
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
    }};
    const array<string, 12> positive_only = {{
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
    }};

    for (const string& row : (mode == 0 ? with_zero : positive_only))
        cout << row << '\n';
}
