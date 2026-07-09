#include <iostream>
using namespace std;
int main() {
    int J, M;
    if (!(cin >> J >> M)) return 0;
    for(int i=0; i<M; ++i) {
        for(int j=0; j<J; ++j) {
            cout << j << (j==J-1 ? "" : " ");
        }
        cout << "\n";
    }
    return 0;
}
