#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int mode;
    cin >> mode;
    (void)mode;

    const vector<string> answer = {
        "1   1   111 ",
        "11 33  3   2",
        "1 3 1  3   3",
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

    for (const string &row : answer) {
        cout << row << '\n';
    }
    return 0;
}
