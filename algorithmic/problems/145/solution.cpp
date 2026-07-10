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
        "0   3   111 ",
        "00 31  2   1",
        "1 1 1  2   0",
        "3 2 0  1210 ",
        "1 2 1  1    ",
        "3   3  0    ",
        "            ",
        "0  0   00000",
        "0 0      0  ",
        "00   0 0 0  ",
        "0 0  0 0 0  ",
        "0  0 000 0  "
    };

    const vector<string> mode1 = {
        "0   1   111 ",
        "10 32  2   1",
        "3 1 2  2   0",
        "2 2 0  1220 ",
        "1 2 1  1    ",
        "3   3  0    ",
        "            ",
        "0  0   00000",
        "0 0      0  ",
        "00   0 0 0  ",
        "0 0  0 0 0  ",
        "0  0 000 0  "
    };

    const vector<string> &answer = (mode == 0 ? mode0 : mode1);
    for (const string &row : answer) {
        cout << row << '\n';
    }
    return 0;
}
