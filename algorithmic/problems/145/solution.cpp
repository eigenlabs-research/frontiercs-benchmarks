// Number Loop Construction - problem 145
// Mode 0: 46 ones, 1 solution -> score 0.911
// Mode 1: 42 ones, 1 solution -> score 0.875
#include <iostream>
#include <string>
#include <vector>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int mode;
    cin >> mode;
    const vector<string> mode0 = {
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
    };

    const vector<string> mode1 = {
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
    };

    const vector<string> &answer = (mode == 0 ? mode0 : mode1);
    for (const string &row : answer) {
        cout << row << '\n';
    }
    return 0;
}
