#include <iostream>
#include <string>

using namespace std;

int main() {
    int G;
    if (!(cin >> G)) return 0;
    if (G == -1) return 0;
    
    while (true) {
        string token;
        if (!(cin >> token)) break;
        if (token == "-1") break;
        
        if (token == "STATE") {
            int hand_idx, r, a, b, P, k;
            cin >> hand_idx >> r >> a >> b >> P >> k;
            
            string tag; cin >> tag; // ALICE
            int s1, v1, s2, v2;
            cin >> s1 >> v1 >> s2 >> v2;
            
            cin >> tag; // BOARD
            for(int i=0; i<k; ++i) {
                int bs, bv; cin >> bs >> bv;
            }
            
            if (r == 1) {
                cout << "RATE 290\n";
                string rtag; double w, d;
                cin >> rtag >> w >> d;
                
                // Exploit logic
                if (w + d/2.0 > 0.40) {
                    if (a > 0) cout << "ACTION RAISE " << a << "\n";
                    else cout << "ACTION CHECK\n";
                } else {
                    cout << "ACTION FOLD\n";
                }
            } else {
                if (a > 0) cout << "ACTION RAISE " << a << "\n";
                else cout << "ACTION CHECK\n";
            }
            cout << flush;
        } else if (token == "RESULT") {
            int delta; cin >> delta;
        } else if (token == "SCORE") {
            double score; cin >> score;
            break;
        }
    }
    return 0;
}
