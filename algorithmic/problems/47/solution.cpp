#include <bits/stdc++.h>
using namespace std;

struct Item { string id; int w,h; long long v,lim; };
struct Seg { int l,r,y; };
struct P { int t,x,y,rot; };

// The input grammar is deliberately small, but this reader accepts arbitrary JSON whitespace
// and field order rather than depending on the presentation of a particular instance.
struct Json {
    string s; size_t p=0;
    Json(){ s.assign(istreambuf_iterator<char>(cin),{}); }
    void ws(){ while(p<s.size() && isspace((unsigned char)s[p])) ++p; }
    char get(){ ws(); return p<s.size()?s[p++]:0; }
    string str(){ ws(); string r; if(p<s.size()&&s[p]=='\"') ++p; while(p<s.size()&&s[p]!='\"'){ if(s[p]=='\\' && p+1<s.size()) ++p; r+=s[p++]; } if(p<s.size())++p; return r; }
    long long num(){ ws(); int sg=1; if(s[p]=='-')sg=-1,++p; long long x=0; while(p<s.size()&&isdigit((unsigned char)s[p])) x=x*10+s[p++]-'0'; return sg*x; }
    bool boolean(){ ws(); bool b=s.compare(p,4,"true")==0; p+=b?4:5; return b; }
};

static int W,H; static bool canrot; static vector<Item> a;

static void mergeSeg(vector<Seg>& z){
    sort(z.begin(),z.end(),[](const Seg&A,const Seg&B){return A.l<B.l;});
    vector<Seg> q;
    for(auto e:z) if(e.l<e.r) {
        if(!q.empty() && q.back().r==e.l && q.back().y==e.y) q.back().r=e.r;
        else q.push_back(e);
    }
    z.swap(q);
}
static int topAt(const vector<Seg>& s,int x,int width){
    int y=0, e=x+width;
    for(const auto& q:s) if(q.r>x && q.l<e) y=max(y,q.y);
    return y;
}
static void raise(vector<Seg>& s,int x,int width,int ny){
    int e=x+width; vector<Seg> z; z.reserve(s.size()+3);
    for(auto q:s){
        if(q.r<=x || q.l>=e) z.push_back(q);
        else { if(q.l<x) z.push_back({q.l,x,q.y}); if(q.r>e) z.push_back({e,q.r,q.y}); }
    }
    z.push_back({x,e,ny}); mergeSeg(z); s.swap(z);
}

struct Pos { bool ok=false; int x=0,y=0,rot=0; };
// Bottom-left skyline placement.  The final tie breaker changes between restarts,
// giving the same simple mechanism several useful ways to resolve near ties.
static Pos bestPos(const vector<Seg>& s, const Item& it, int orientMode){
    Pos ans; int bestTop=INT_MAX;
    for(int rot=0;rot<=(canrot?1:0);rot++){
        int iw=rot?it.h:it.w, ih=rot?it.w:it.h;
        if(iw>W||ih>H) continue;
        // Sweep candidate left edges.  The multiset is the maximum skyline height in
        // [x,x+iw), avoiding a quadratic re-scan of the skyline for every candidate.
        multiset<int> active; int add=0;
        for(int i=0;i<(int)s.size();i++){
            int x=s[i].l;
            while(add<(int)s.size() && s[add].l<x+iw) active.insert(s[add++].y);
            if(i) { auto q=active.find(s[i-1].y); if(q!=active.end()) active.erase(q); }
            if(x+iw>W || active.empty()) continue;
            int y=*active.rbegin(); if(y+ih>H) continue;
            int nt=y+ih;
            bool take=!ans.ok || y<ans.y || (y==ans.y && nt<bestTop);
            if(y==ans.y && nt==bestTop){
                if(orientMode&1) take=(rot>ans.rot) || (rot==ans.rot && x>ans.x);
                else take=(rot<ans.rot) || (rot==ans.rot && x<ans.x);
            }
            if(take) ans={true,x,y,rot},bestTop=nt;
        }
    }
    return ans;
}

static vector<P> pack(vector<int> order,int mode){
    vector<Seg> sky={{0,W,0}}; vector<P> out;
    for(int t:order){
        long long usable=min(a[t].lim, (long long)W*H/(1LL*a[t].w*a[t].h)+2);
        for(long long k=0;k<usable;k++){
            Pos z=bestPos(sky,a[t],mode);
            if(!z.ok) break;
            int ih=z.rot?a[t].w:a[t].h;
            out.push_back({t,z.x,z.y,z.rot});
            raise(sky,z.x,z.rot?a[t].h:a[t].w,z.y+ih);
        }
    }
    return out;
}
static long long score(const vector<P>& p){ long long r=0; for(auto x:p)r+=a[x.t].v; return r; }

int main(){
    Json j; j.get(); vector<Item> items;
    for(int root=0;root<2;root++){
        string key=j.str(); j.get();
        if(key=="bin"){
            j.get(); for(int k=0;k<3;k++){ string z=j.str();j.get(); if(z=="W")W=j.num(); else if(z=="H")H=j.num(); else canrot=j.boolean(); if(k<2)j.get(); } j.get();
        } else {
            j.get(); j.ws();
            while(j.p<j.s.size()&&j.s[j.p]!=']'){
                j.get(); Item it;
                for(int k=0;k<5;k++){ string z=j.str();j.get(); if(z=="type")it.id=j.str(); else if(z=="w")it.w=j.num(); else if(z=="h")it.h=j.num(); else if(z=="v")it.v=j.num(); else it.lim=j.num(); if(k<4)j.get(); }
                j.get(); items.push_back(it); j.ws(); if(j.s[j.p]==',')j.p++,j.ws();
            } j.get();
        }
        j.ws(); if(j.p<j.s.size()&&j.s[j.p]==',')j.p++;
    }
    a=items; int n=a.size(); vector<P> answer; long long best=-1;
    // A small deterministic portfolio of priorities; all candidates use the identical skyline
    // feasibility construction, so selecting the most valuable output is safe.
    for(int mode=0;mode<10;mode++){
        vector<int> ord(n); iota(ord.begin(),ord.end(),0);
        sort(ord.begin(),ord.end(),[&](int i,int k){
            long double di=(long double)a[i].v/(a[i].w*a[i].h), dk=(long double)a[k].v/(a[k].w*a[k].h);
            long double xi, xk;
            if(mode==0) xi=di,xk=dk;
            else if(mode==1) xi=a[i].v,xk=a[k].v;
            else if(mode==2) xi=di*(1.0L+0.16L*(a[i].lim<80)),xk=dk*(1.0L+0.16L*(a[k].lim<80));
            else if(mode==3) xi=di/(min(a[i].w,a[i].h)),xk=dk/(min(a[k].w,a[k].h));
            else if(mode==4) xi=di*min(a[i].w,a[i].h),xk=dk*min(a[k].w,a[k].h);
            else { unsigned long long hi=(unsigned)(i*1103515245u+mode*12345u), hk=(unsigned)(k*1103515245u+mode*12345u); xi=di*(0.82L+(hi%1000)/2500.0L); xk=dk*(0.82L+(hk%1000)/2500.0L); }
            if(xi!=xk)return xi>xk; return a[i].v>a[k].v;
        });
        auto got=pack(ord,mode); long long val=score(got); if(val>best) best=val,answer.swap(got);
    }
    cout<<"{\"placements\":[";
    for(size_t i=0;i<answer.size();i++){ if(i)cout<<','; auto q=answer[i]; cout<<"{\"type\":\""<<a[q.t].id<<"\",\"x\":"<<q.x<<",\"y\":"<<q.y<<",\"rot\":"<<q.rot<<'}'; }
    cout<<"]}\n";
}
