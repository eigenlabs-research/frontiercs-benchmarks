#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int mode = 0;
    if (!(cin >> mode)) return 0;

    // Candidate A001: maximize the secondary objective while using the same
    // admissible construction in both modes.  Every clue position is filled
    // with '1', which is allowed for mode 0 and mode 1; blanks are preserved
    // exactly from the fixed template.
    const vector<string> out = {
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
    for (const string &s : out) cout << s << '\n';
    return 0;
}
