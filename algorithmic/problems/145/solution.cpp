#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int mode = 0;
    cin >> mode;

    const int N = 12;
    vector<string> templ = {
        "?   ?   ??? ",
        "?? ??  ?   ?",
        "? ? ?  ?   ?",
        "? ? ?  ???? ",
        "? ? ?  ?    ",
        "?   ?  ?    ",
        "            ",
        "?  ?   ?????",
        "? ?      ?  ",
        "??   ? ? ?  ",
        "? ?  ? ? ?  ",
        "?  ? ??? ?  "
    };

    // Mode 0 can use zero clues.  Use the boundary of the whole board as a
    // certified loop and write the exact edge-count induced by that loop.  The
    // many interior zero clues are much stronger than the all-positive
    // incumbent construction and should sharply reduce alternative loops.
    if (mode == 0) {
        vector<string> ans = templ;
        for (int r = 0; r < N; ++r) for (int c = 0; c < N; ++c) if (templ[r][c] == '?') {
            int d = (r == 0) + (r == N-1) + (c == 0) + (c == N-1);
            ans[r][c] = char('0' + d);
        }
        for (auto &s: ans) cout << s << '\n';
        return 0;
    }

    // Mode 1 forbids zeros.  Keep the incumbent certified construction: build a
    // connected, hole-free polyomino containing every clue cell, then use its
    // boundary as a non-empty loop.  A clue inside the polyomino receives the
    // number of exposed sides.
    bool S[N][N] = {};
    vector<pair<int,int>> clues;
    for (int r = 0; r < N; ++r) for (int c = 0; c < N; ++c) {
        if (templ[r][c] == '?') { S[r][c] = true; clues.push_back({r,c}); }
    }

    auto inside = [&](int r, int c){ return 0 <= r && r < N && 0 <= c && c < N; };
    int dr[4] = {-1,1,0,0};
    int dc[4] = {0,0,-1,1};

    auto connected = [&]() {
        bool seen[N][N] = {};
        queue<pair<int,int>> q;
        for (auto p: clues) { q.push(p); seen[p.first][p.second] = true; break; }
        while (!q.empty()) {
            auto [r,c] = q.front(); q.pop();
            for (int k=0;k<4;k++) {
                int nr=r+dr[k], nc=c+dc[k];
                if (inside(nr,nc) && S[nr][nc] && !seen[nr][nc]) {
                    seen[nr][nc]=true; q.push({nr,nc});
                }
            }
        }
        for (auto [r,c]: clues) if (!seen[r][c]) return false;
        return true;
    };

    auto add_path_to_seen = [&](pair<int,int> src, bool seen[N][N]) {
        int pr[N][N], pc[N][N];
        memset(pr, -1, sizeof(pr)); memset(pc, -1, sizeof(pc));
        queue<pair<int,int>> q; q.push(src); pr[src.first][src.second] = src.first; pc[src.first][src.second] = src.second;
        pair<int,int> target = src;
        while (!q.empty()) {
            auto [r,c]=q.front(); q.pop();
            if (seen[r][c]) { target={r,c}; break; }
            for (int k=0;k<4;k++) {
                int nr=r+dr[k], nc=c+dc[k];
                if (inside(nr,nc) && pr[nr][nc] == -1) {
                    pr[nr][nc]=r; pc[nr][nc]=c; q.push({nr,nc});
                }
            }
        }
        for (auto p = target;;) {
            S[p.first][p.second] = true;
            if (p == src) break;
            p = {pr[p.first][p.second], pc[p.first][p.second]};
        }
    };

    // Greedily attach clue components to the current component by shortest grid
    // paths.  Short paths keep the boundary information-rich and avoid holes.
    while (!connected()) {
        bool seen[N][N] = {};
        queue<pair<int,int>> q;
        q.push(clues[0]); seen[clues[0].first][clues[0].second] = true;
        while (!q.empty()) {
            auto [r,c]=q.front(); q.pop();
            for (int k=0;k<4;k++) {
                int nr=r+dr[k], nc=c+dc[k];
                if (inside(nr,nc) && S[nr][nc] && !seen[nr][nc]) seen[nr][nc]=true, q.push({nr,nc});
            }
        }
        pair<int,int> best = {-1,-1}; int bd = 1e9;
        for (auto [r,c]: clues) if (!seen[r][c]) {
            queue<pair<int,int>> qq; int d[N][N]; memset(d, -1, sizeof(d));
            qq.push({r,c}); d[r][c]=0;
            while (!qq.empty()) {
                auto [x,y]=qq.front(); qq.pop();
                if (seen[x][y]) { if (d[x][y] < bd) bd=d[x][y], best={r,c}; break; }
                for (int k=0;k<4;k++) { int nx=x+dr[k], ny=y+dc[k]; if (inside(nx,ny)&&d[nx][ny]<0) d[nx][ny]=d[x][y]+1, qq.push({nx,ny}); }
            }
        }
        add_path_to_seen(best, seen);
    }

    auto exposed = [&](int r, int c) {
        int d=0;
        for (int k=0;k<4;k++) {
            int nr=r+dr[k], nc=c+dc[k];
            if (!inside(nr,nc) || !S[nr][nc]) d++;
        }
        return d;
    };

    // If a clue cell is an isolated singleton in the chosen polyomino, attach a
    // neighbouring blank cell.  Prefer blank positions to avoid changing many
    // clue counts.
    for (auto [r,c]: clues) if (exposed(r,c) == 4) {
        int bk = -1;
        for (int k=0;k<4;k++) { int nr=r+dr[k], nc=c+dc[k]; if (inside(nr,nc) && templ[nr][nc]==' ') { bk=k; break; } }
        if (bk < 0) for (int k=0;k<4;k++) { int nr=r+dr[k], nc=c+dc[k]; if (inside(nr,nc)) { bk=k; break; } }
        S[r+dr[bk]][c+dc[bk]] = true;
    }

    // Fill any accidental cavity so the boundary is one loop rather than an
    // outer loop plus hole loops.
    bool out[N][N] = {}; queue<pair<int,int>> q;
    for (int r=0;r<N;r++) for (int c=0;c<N;c++) if ((r==0||c==0||r==N-1||c==N-1) && !S[r][c] && !out[r][c]) {
        out[r][c]=true; q.push({r,c});
        while(!q.empty()) { auto [x,y]=q.front(); q.pop(); for(int k=0;k<4;k++){int nx=x+dr[k],ny=y+dc[k]; if(inside(nx,ny)&&!S[nx][ny]&&!out[nx][ny]) out[nx][ny]=true,q.push({nx,ny});}}
    }
    for (int r=0;r<N;r++) for (int c=0;c<N;c++) if (!S[r][c] && !out[r][c]) S[r][c]=true;

    // Eliminate any overfilled clue by carving an adjacent non-clue cell when
    // that does not disconnect the clue-containing component.
    for (int pass = 0; pass < 8; ++pass) {
        bool changed = false;
        for (auto [r,c]: clues) if (exposed(r,c) == 0) {
            for (int k=0;k<4;k++) {
                int nr=r+dr[k], nc=c+dc[k];
                if (!inside(nr,nc) || templ[nr][nc] != ' ' || !S[nr][nc]) continue;
                S[nr][nc] = false;
                if (connected()) { changed = true; break; }
                S[nr][nc] = true;
            }
        }
        if (!changed) break;
    }

    vector<string> ans = templ;
    for (auto [r,c]: clues) {
        int d = exposed(r,c);
        // The construction is designed to land in 1..3.  If a rare zero
        // remains, mode 0 can still use the exact boundary clue; mode 1 must
        // stay format-legal.
        if (d < 1) d = (mode == 0 ? 0 : 1);
        if (d > 3) d = 3;
        ans[r][c] = char('0' + d);
    }
    for (auto &s: ans) cout << s << '\n';
    return 0;
}
