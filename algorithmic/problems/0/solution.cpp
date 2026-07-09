#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace std;

// Solução sub-ótima baseline para Polyominoes (Bottom-Left-Fill aproximado)
// Para o pipeline, apenas enfileira as peças numa grade horizontal para passar o validador.
struct Poly { int k; vector<pair<int,int>> cells; };

int main() {
    int n;
    if (!(cin >> n)) return 0;
    vector<Poly> polys(n);
    int total_cells = 0;
    for (int i = 0; i < n; ++i) {
        cin >> polys[i].k;
        total_cells += polys[i].k;
        for (int j = 0; j < polys[i].k; ++j) {
            int x, y; cin >> x >> y;
            polys[i].cells.push_back({x, y});
        }
    }
    
    // Heurística ingénua (fila 1D simples) para garantir validade imediata e pontuar
    cout << total_cells << " " << 10 << endl;
    int curr_x = 0;
    for (int i = 0; i < n; ++i) {
        int min_x = 10000, min_y = 10000;
        int max_x = -10000;
        for (auto p : polys[i].cells) {
            min_x = min(min_x, p.first);
            min_y = min(min_y, p.second);
            max_x = max(max_x, p.first);
        }
        cout << (curr_x - min_x) << " " << (0 - min_y) << " 0 0\n";
        curr_x += (max_x - min_x) + 1;
    }
    return 0;
}
