#include <bits/stdc++.h>
using namespace std;

struct Polyomino {
    int k;
    vector<pair<int, int>> cells;
};

struct Placement {
    int x, y, r, f;
};

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(nullptr);

    int n;
    cin >> n;
    vector<Polyomino> polyominoes(n);
    for (int i = 0; i < n; ++i) {
        cin >> polyominoes[i].k;
        polyominoes[i].cells.resize(polyominoes[i].k);
        for (int j = 0; j < polyominoes[i].k; ++j) {
            cin >> polyominoes[i].cells[j].first >> polyominoes[i].cells[j].second;
        }
    }

    // Place all polyominoes in a simple grid layout
    int W = 100, H = 100;
    vector<Placement> placements(n);
    int x = 0, y = 0;
    for (int i = 0; i < n; ++i) {
        placements[i] = {x, y, 0, 0};
        x += 2;
        if (x >= W) {
            x = 0;
            y += 2;
        }
    }

    // Output the result
    cout << W << " " << H << "\n";
    for (int i = 0; i < n; ++i) {
        cout << placements[i].x << " " << placements[i].y << " " << placements[i].r << " " << placements[i].f << "\n";
    }

    return 0;
}