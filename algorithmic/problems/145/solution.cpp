#include <bits/stdc++.h>
using namespace std;

static const vector<string> mask = {
    "?   ?   ??? ", "?? ??  ?   ?", "? ? ?  ?   ?", "? ? ?  ???? ",
    "? ? ?  ?    ", "?   ?  ?    ", "            ", "?  ?   ?????",
    "? ?      ?  ", "??   ? ? ?  ", "? ?  ? ? ?  ", "?  ? ??? ?  "
};

int val(const vector<vector<int>>& a, int r, int c) {
    int z=0;
    const int dr[4]={-1,1,0,0}, dc[4]={0,0,-1,1};
    for(int k=0;k<4;k++) {
        int x=r+dr[k], y=c+dc[k];
        if(x<0||x>=12||y<0||y>=12 || a[x][y]!=a[r][c]) ++z;
    }
    return z;
}
int ones(const vector<vector<int>>& a) {
    int s=0;
    for(int r=0;r<12;r++) for(int c=0;c<12;c++)
        if(mask[r][c]=='?' && val(a,r,c)==1) ++s;
    return s;
}
bool usable(const vector<vector<int>>& a) {
    for(int r=0;r<12;r++) for(int c=0;c<12;c++) if(mask[r][c]=='?') {
        int x=val(a,r,c); if(x<1 || x>3) return false;
    }
    return true;
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int mode; if(!(cin>>mode)) return 0;
    vector<vector<int>> best(12,vector<int>(12,0));
    if(mode==0) {
        // The boundary of the whole board is one loop.  All interior clues are
        // zeros, which pins almost every non-boundary edge.
    } else {
        int bs=-1;
        // A collection of alternating complete rows (or columns), joined at
        // alternate ends, is an induced path of cells.  Its boundary is one
        // simple loop and it touches every cell of the board.
        for(int orient=0;orient<2;orient++) for(int parity=0;parity<2;parity++)
            for(int first=0;first<2;first++) {
                vector<vector<int>> a(12,vector<int>(12));
                if(!orient) {
                    for(int r=parity;r<12;r+=2) for(int c=0;c<12;c++) a[r][c]=1;
                    int t=0;
                    for(int r=parity+1;r<11;r+=2) { a[r][(t++%2)?(11-first):first]=1; }
                } else {
                    for(int c=parity;c<12;c+=2) for(int r=0;r<12;r++) a[r][c]=1;
                    int t=0;
                    for(int c=parity+1;c<11;c+=2) { a[(t++%2)?(11-first):first][c]=1; }
                }
                if(usable(a) && ones(a)>bs) bs=ones(a),best=a;
            }

        // Search further among induced trees of cells.  The boundary of any
        // nonempty tree polyomino is a single loop.  Growing only at a cell
        // with one existing neighbour preserves that property; this gives
        // substantially more branch (and hence '1') clues than a plain snake.
        mt19937 rng(712367);
        for(int it=0;it<1800;it++) {
            vector<vector<int>> a(12,vector<int>(12));
            a[rng()%12][rng()%12]=1;
            for(int step=0;step<80;step++) {
                vector<pair<int,int>> q;
                int high=-1000000;
                for(int r=0;r<12;r++) for(int c=0;c<12;c++) if(!a[r][c]) {
                    int nr=0, pr=-1,pc=-1;
                    const int dr[4]={-1,1,0,0},dc[4]={0,0,-1,1};
                    for(int k=0;k<4;k++){int x=r+dr[k],y=c+dc[k]; if(x>=0&&x<12&&y>=0&&y<12&&a[x][y]) nr++,pr=x,pc=y;}
                    if(nr!=1) continue;
                    // Do not turn a numbered cell into a four-way interior cell.
                    if(mask[pr][pc]=='?' && val(a,pr,pc)==1) continue;
                    int w=(int)(rng()%17);
                    if(mask[r][c]=='?' && val(a,r,c)==0) w+=35;
                    for(int k=0;k<4;k++){int x=r+dr[k],y=c+dc[k]; if(x>=0&&x<12&&y>=0&&y<12&&mask[x][y]=='?'&&val(a,x,y)==0) w+=12;}
                    if(w>high) high=w,q.clear();
                    if(w==high) q.push_back({r,c});
                }
                if(q.empty()) break;
                auto p=q[rng()%q.size()]; a[p.first][p.second]=1;
                if(step>15 && usable(a)) { int z=ones(a); if(z>bs) bs=z,best=a; }
            }
        }
    }
    for(int r=0;r<12;r++) {
        string out=mask[r];
        for(int c=0;c<12;c++) if(out[c]=='?') out[c]=char('0'+val(best,r,c));
        cout<<out<<'\n';
    }
}
