#include <iostream>
#include <string>
using namespace std;
int main() {
    int allowed_zero;
    if (!(cin >> allowed_zero)) return 0;
    string grid[12] = {
        "0   0   000 ",
        "00 00  0   0",
        "0 0 0  0   0",
        "0 0 0  0000 ",
        "0 0 0  0    ",
        "0   0  0    ",
        "            ",
        "0  0   00000",
        "0 0      0  ",
        "00   0 0 0  ",
        "0 0  0 0 0  ",
        "0  0 000 0  "
    };
    if(allowed_zero == 1) {
        for(int r=0; r<12; ++r) {
            for(int c=0; c<12; ++c) {
                if(grid[r][c] == '0') grid[r][c] = '1';
            }
        }
    }
    for(int r=0; r<12; ++r) cout << grid[r] << "\n";
    return 0;
}
