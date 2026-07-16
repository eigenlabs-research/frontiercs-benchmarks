#include <bits/stdc++.h>
using namespace std;

struct Item { string id; int w,h; long long v; int lim; };
struct P { int t,x,y,r,w,h; };

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    string s((istreambuf_iterator<char>(cin)), {});
    smatch m;
    regex binre("\\\"W\\\"\\s*:\\s*([0-9]+)\\s*,\\s*\\\"H\\\"\\s*:\\s*([0-9]+)\\s*,\\s*\\\"allow_rotate\\\"\\s*:\\s*(true|false)");
    if(!regex_search(s,m,binre)) { cout << "{\"placements\":[]}"; return 0; }
    int W=stoi(m[1]), H=stoi(m[2]); bool rot=(m[3]=="true");
    regex ire("\\{\\s*\\\"type\\\"\\s*:\\s*\\\"([^\\\"]*)\\\"\\s*,\\s*\\\"w\\\"\\s*:\\s*([0-9]+)\\s*,\\s*\\\"h\\\"\\s*:\\s*([0-9]+)\\s*,\\s*\\\"v\\\"\\s*:\\s*([0-9]+)\\s*,\\s*\\\"limit\\\"\\s*:\\s*([0-9]+)\\s*\\}");
    vector<Item> a;
    for(sregex_iterator it(s.begin(),s.end(),ire), e; it!=e; ++it) {
        auto q=*it; a.push_back({q[1],stoi(q[2]),stoi(q[3]),stoll(q[4]),stoi(q[5])});
    }
    if(W<=0 || H<=0 || a.empty()) { cout << "{\"placements\":[]}"; return 0; }

    vector<P> best; long long bestv=-1;
    // A number of differently perturbed density orders makes the skyline construction less
    // dependent on one arbitrary ordering of similarly useful rectangles.
    mt19937 rng(712367);
    int trials=18;
    for(int z=0; z<trials; ++z){
        vector<int> ord(a.size()); iota(ord.begin(),ord.end(),0);
        vector<double> key(a.size());
        for(int i=0;i<(int)a.size();++i){
            double den=(double)a[i].v/((double)a[i].w*a[i].h);
            if(z==1) den=(double)a[i].v;
            else if(z==2) den=den*sqrt((double)min(a[i].w,a[i].h));
            else if(z>=3) den*=0.78 + (rng()%10001)/10000.0*0.44;
            key[i]=den;
        }
        sort(ord.begin(),ord.end(),[&](int i,int j){
            if(key[i]!=key[j]) return key[i]>key[j];
            return a[i].v>a[j].v;
        });
        // One integer per x-coordinate is a compact skyline.  The generated coordinates are
        // integral, and assigning the maximum height over an item's interval guarantees no overlap.
        vector<int> sky(W,0); vector<P> out; long long val=0;
        for(int ti:ord){
            const Item &it=a[ti];
            for(int copy=0;copy<it.lim;copy++){
                int bx=-1,by=INT_MAX,br=0,btop=INT_MAX,bw=0,bh=0;
                for(int r=0;r<=(rot?1:0);r++){
                    int w=r?it.h:it.w, h=r?it.w:it.h;
                    if(w>W || h>H) continue;
                    // Sliding maximum computes every possible bottom-left position in O(W).
                    deque<int> dq;
                    for(int x=0;x<W;x++){
                        while(!dq.empty() && sky[dq.back()]<=sky[x]) dq.pop_back();
                        dq.push_back(x);
                        while(!dq.empty() && dq.front()<=x-w) dq.pop_front();
                        if(x>=w-1){
                            int left=x-w+1, y=sky[dq.front()];
                            if(y+h>H) continue;
                            int top=y+h;
                            // Bottom-left preference; alternate tie direction between trials.
                            bool better=(y<by || (y==by && top<btop) ||
                              (y==by && top==btop && ((z&1) ? left>bx : (bx<0||left<bx))));
                            if(better) bx=left,by=y,br=r,btop=top,bw=w,bh=h;
                        }
                    }
                }
                if(bx<0) break;
                for(int x=bx;x<bx+bw;x++) sky[x]=by+bh;
                out.push_back({ti,bx,by,br,bw,bh}); val+=it.v;
            }
        }
        if(val>bestv){ bestv=val; best.swap(out); }
    }
    cout << "{\"placements\":[";
    for(size_t i=0;i<best.size();++i){
        if(i) cout << ',';
        const P&p=best[i];
        cout << "{\"type\":\"" << a[p.t].id << "\",\"x\":" << p.x
             << ",\"y\":" << p.y << ",\"rot\":" << p.r << '}';
    }
    cout << "]}";
    return 0;
}
