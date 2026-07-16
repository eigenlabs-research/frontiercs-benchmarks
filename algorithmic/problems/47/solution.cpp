#include <bits/stdc++.h>
using namespace std;
struct Item { string id; int w,h; long long v, lim; };
struct Put { int t,x,y,rot,w,h; };
struct Ans { long long val=0; vector<Put> p; };

static bool getnum(const string& s,const string& key,long long& z){
    regex r("\\\""+key+"\\\"\\s*:\\s*(-?[0-9]+)"); smatch m;
    if(!regex_search(s,m,r)) return false; z=stoll(m[1]); return true;
}
int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    string in((istreambuf_iterator<char>(cin)),{});
    long long qW=0,qH=0; getnum(in,"W",qW); getnum(in,"H",qH);
    bool rotate = false;
    { regex r("\\\"allow_rotate\\\"\\s*:\\s*(true|false)"); smatch m; if(regex_search(in,m,r)) rotate=m[1]=="true"; }
    int W=(int)qW,H=(int)qH;
    vector<Item>a;
    size_t pos=0;
    while((pos=in.find("\"type\"",pos))!=string::npos){
        size_t l=in.rfind('{',pos), r=in.find('}',pos); if(l==string::npos||r==string::npos) break;
        string o=in.substr(l,r-l+1); smatch m;
        regex rt("\\\"type\\\"\\s*:\\s*\\\"([^\\\"]*)\\\"");
        if(!regex_search(o,m,rt)){pos=r+1;continue;}
        Item z; z.id=m[1]; long long x;
        if(!getnum(o,"w",x)){pos=r+1;continue;} z.w=x;
        if(!getnum(o,"h",x)){pos=r+1;continue;} z.h=x;
        if(!getnum(o,"v",x)){pos=r+1;continue;} z.v=x;
        if(!getnum(o,"limit",x)){pos=r+1;continue;} z.lim=x;
        a.push_back(z); pos=r+1;
    }
    Ans best;
    vector<int> heights;
    for(auto &z:a){ if(z.h<=H) heights.push_back(z.h); if(rotate&&z.w<=H) heights.push_back(z.w); }
    sort(heights.begin(),heights.end()); heights.erase(unique(heights.begin(),heights.end()),heights.end());
    auto start=chrono::steady_clock::now();
    // A candidate is a stack of shelves.  Within a shelf, left-to-right greedy is a bounded 1-D packing.
    for(int mode=0;mode<3;mode++){
        vector<long long> used(a.size()); Ans cur; int y=0;
        while(y<H && cur.p.size()<100000){
            if(chrono::duration<double>(chrono::steady_clock::now()-start).count()>0.82) break;
            long double bestRate=-1; long long bestRowVal=-1; int chooseH=0; vector<pair<int,int>> choose;
            for(int sh:heights){
                if(sh>H-y) break;
                vector<long long> take(a.size()); vector<pair<int,int>> row; int rem=W; long long rv=0;
                while(rem>0){
                    int bt=-1, br=0, bw=0; long double bk=-1;
                    for(int i=0;i<(int)a.size();i++) if(used[i]+take[i]<a[i].lim){
                        for(int rr=0;rr<=(rotate?1:0);rr++){
                            int iw=rr?a[i].h:a[i].w, ih=rr?a[i].w:a[i].h;
                            if(iw>rem||ih>sh) continue;
                            long double k;
                            if(mode==0) k=(long double)a[i].v/(iw*1.0L*ih); // area density
                            else if(mode==1) k=(long double)a[i].v/iw;       // row profit density
                            else k=(long double)a[i].v;                      // high-value first
                            // Small width tie break improves the chance of closing a row.
                            if(k>bk+1e-18L || (fabsl(k-bk)<1e-18L && iw<bw)){bk=k;bt=i;br=rr;bw=iw;}
                        }
                    }
                    if(bt<0) break;
                    row.push_back({bt,br}); take[bt]++; rem-=bw; rv+=a[bt].v;
                    if(row.size()>100000) break;
                }
                long double rate=(long double)rv/sh;
                if(!row.empty() && (rate>bestRate+1e-15L || (fabsl(rate-bestRate)<1e-15L && rv>bestRowVal))){
                    bestRate=rate; bestRowVal=rv; chooseH=sh; choose.swap(row);
                }
            }
            if(chooseH==0) break;
            int x=0;
            for(auto [i,rr]:choose){
                int iw=rr?a[i].h:a[i].w, ih=rr?a[i].w:a[i].h;
                cur.p.push_back({i,x,y,rr,iw,ih}); x+=iw; used[i]++; cur.val+=a[i].v;
            }
            y+=chooseH;
        }
        if(cur.val>best.val) best=move(cur);
    }
    cout << "{\"placements\":[";
    for(size_t i=0;i<best.p.size();i++){
        if(i) cout<<','; auto &z=best.p[i];
        cout << "{\"type\":\"" << a[z.t].id << "\",\"x\":" << z.x << ",\"y\":" << z.y << ",\"rot\":" << z.rot << '}';
    }
    cout << "]}";
}
