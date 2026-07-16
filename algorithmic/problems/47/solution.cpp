#include <bits/stdc++.h>
using namespace std;
struct Item { string id; int w,h; long long v; int lim; };
struct P { int t,x,y,r; };
struct Ans { long long val=0; vector<P> p; };
int W,H; bool rotok; vector<Item> a;

// The input objects have fixed keys; this deliberately small parser accepts JSON whitespace.
string allin(){ return string((istreambuf_iterator<char>(cin)),{}); }
long long keynum(const string&s,const string& k){
    regex q("\\\""+k+"\\\"\\s*:\\s*([0-9]+)"); smatch m;
    return regex_search(s,m,q)?stoll(m[1]):-1;
}
Ans shelves(int mode, bool trans){
    int BW=trans?H:W, BH=trans?W:H;
    vector<int> used(a.size()); Ans out;
    int yy=0, step=0;
    while(yy<BH){
        vector<int> heights;
        for(int i=0;i<(int)a.size();++i) if(used[i]<a[i].lim){
            if(a[i].h<=BH-yy) heights.push_back(a[i].h);
            if(rotok && a[i].w<=BH-yy) heights.push_back(a[i].w);
        }
        sort(heights.begin(),heights.end()); heights.erase(unique(heights.begin(),heights.end()),heights.end());
        long double bestfit=-1; vector<P> bestrow; vector<int> bestadd; int besth=0; long long bestv=0;
        for(int sh:heights){
            vector<int> cnt(a.size()); vector<P> row; int xx=0; long long vv=0;
            while(true){
                int bi=-1, br=0; long double bs=-1;
                for(int i=0;i<(int)a.size();++i) if(used[i]+cnt[i]<a[i].lim){
                    for(int r=0;r<=(rotok?1:0);++r){
                        int iw=r?a[i].h:a[i].w, ih=r?a[i].w:a[i].h;
                        // Dimensions in the transposed construction are swapped.
                        if(trans) swap(iw,ih);
                        if(ih>sh || iw>BW-xx) continue;
                        long double s;
                        int kind=mode%5;
                        if(kind==0) s=(long double)a[i].v/iw;
                        else if(kind==1) s=(long double)a[i].v/(iw*ih);
                        else if(kind==2) s=(long double)a[i].v;
                        else if(kind==3) s=(long double)a[i].v/(iw+0.18L*ih);
                        else s=(long double)a[i].v/(ih+0.18L*iw);
                        // deterministic mild tie diversification, not input-specific tuning
                        s*=1.0L + ((i*37+step*13+mode*19+r*7)%23-11)*0.003L;
                        if(s>bs){bs=s;bi=i;br=r;}
                    }
                }
                if(bi<0) break;
                int iw=br?a[bi].h:a[bi].w, ih=br?a[bi].w:a[bi].h;
                if(trans) swap(iw,ih);
                row.push_back({bi,xx,yy,br}); xx+=iw; cnt[bi]++; vv+=a[bi].v;
            }
            if(row.empty()) continue;
            long double fit;
            int sel=mode/5;
            if(sel==0) fit=(long double)vv/sh;
            else if(sel==1) fit=(long double)vv/(sh*(0.65L+0.35L*(long double)xx/BW));
            else fit=(long double)vv/(sh+0.12L*(BH-yy));
            if(fit>bestfit){bestfit=fit;bestrow=row;bestadd=cnt;besth=sh;bestv=vv;}
        }
        if(besth==0) break;
        for(int i=0;i<(int)a.size();++i) used[i]+=bestadd[i];
        out.val+=bestv;
        for(auto q:bestrow){
            if(trans) { // (base x,base y), base dimensions are original dimensions swapped
                out.p.push_back({q.t,q.y,q.x,q.r});
            } else out.p.push_back(q);
        }
        yy+=besth; ++step;
    }
    return out;
}
int main(){
    ios::sync_with_stdio(false);cin.tie(nullptr);
    string s=allin();
    smatch m;
    regex bin("\\\"bin\\\"\\s*:\\s*\\{\\s*\\\"W\\\"\\s*:\\s*([0-9]+)\\s*,\\s*\\\"H\\\"\\s*:\\s*([0-9]+)\\s*,\\s*\\\"allow_rotate\\\"\\s*:\\s*(true|false)");
    if(!regex_search(s,m,bin)) { cout<<"{\"placements\":[]}"; return 0; }
    W=stoi(m[1]);H=stoi(m[2]);rotok=m[3]=="true";
    regex it("\\{\\s*\\\"type\\\"\\s*:\\s*\\\"([^\\\"]*)\\\"\\s*,\\s*\\\"w\\\"\\s*:\\s*([0-9]+)\\s*,\\s*\\\"h\\\"\\s*:\\s*([0-9]+)\\s*,\\s*\\\"v\\\"\\s*:\\s*([0-9]+)\\s*,\\s*\\\"limit\\\"\\s*:\\s*([0-9]+)\\s*\\}");
    for(sregex_iterator z(s.begin(),s.end(),it),e;z!=e;++z){ auto x=*z; a.push_back({x[1],stoi(x[2]),stoi(x[3]),stoll(x[4]),stoi(x[5])}); }
    Ans best;
    // Three strip-selection policies and five within-strip priorities, in each axis.
    for(int mode=0;mode<15;mode++) for(int tr=0;tr<2;tr++){
        Ans q=shelves(mode,tr); if(q.val>best.val) best=move(q);
    }
    cout<<"{\"placements\":[";
    for(int j=0;j<(int)best.p.size();j++){
        auto q=best.p[j]; if(j) cout<<',';
        cout<<"{\"type\":\""<<a[q.t].id<<"\",\"x\":"<<q.x<<",\"y\":"<<q.y<<",\"rot\":"<<q.r<<'}';
    }
    cout<<"]}";
}
