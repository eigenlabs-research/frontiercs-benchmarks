#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int mode;
    cin >> mode;
    const vector<string> mode1 = {
        "1   1   111 ",
        "12 11  3   1",
        "1 3 1  3   1",
        "1 3 1  3133 ",
        "1 3 2  3    ",
        "2   3  3    ",
        "            ",
        "3  1   31311",
        "1 3      3  ",
        "13   1 3 3  ",
        "1 3  3 1 3  ",
        "3  1 311 3  "
    };

    const vector<string> mode0 = {
        "1   1   111 ",
        "11 30  3   1",
        "1 1 1  3   1",
        "1 3 1  3133 ",
        "1 3 2  3    ",
        "2   3  3    ",
        "            ",
        "3  1   31311",
        "1 3      3  ",
        "13   1 3 3  ",
        "1 3  3 1 3  ",
        "3  1 311 3  "
    };

    const vector<string> &answer = (mode == 0) ? mode0 : mode1;
    for (const string &row : answer) {
        cout << row << '\n';
    }
    return 0;
}
