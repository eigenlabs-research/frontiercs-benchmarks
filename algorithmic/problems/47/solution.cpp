#include <bits/stdc++.h>
using namespace std;

struct Item { string id; int w,h,lim; long long v; };
struct Pl { int t,x,y,r; };
struct Ans { long long val=0; vector<Pl> p; };

// Small JSON reader: the task input only contains objects, arrays, strings, integers and booleans.
struct Json {
    string s; size_t q=0;
    Json(){ s.assign(istreambuf_iterator<char>(cin),{}); }
    void ws(){ while(q<s.size() && isspace((unsigned char)s[q])) ++q; }
    char ch(){ ws(); return q<s.size()?s[q]:0; }
    void need(char c){ ws(); if(q<s.size() && s[q]==c) ++q; }
    string str(){ need('"'); string r; while(q<s.size() && s[q]!='"'){ if(s[q]=='\\' && q+1<s.size()) ++q; r+=s[q++]; } need('"'); return r; }
    long long num(){ ws(); int sg=1; if(s[q]=='-') sg=-1,++q; long long x=0; while(q<s.size()&&isdigit((unsigned char)s[q])) x=x*10+s[q++]-'0'; return sg*x; }
    bool boolean(){ ws(); bool z=s.compare(q,4,"true")==0; q+=z?4:5; return z; }
};

static pair<int,int> bestSpot(const vector<int>& sky,int w,int h,int H){
    int W=(int)sky.size(); if(w>W) return {INT_MAX,-1};
    deque<int> d; int best=INT_MAX,bx=-1;
    for(int i=0;i<W;i++){
        while(!d.empty() && sky[d.back()]<=sky[i]) d.pop_back();
        d.push_back(i);
        while(!d.empty() && d.front()<=i-w) d.pop_front();
        if(i>=w-1){ int y=sky[d.front()]; if(y+h<=H && y+h<best) best=y+h,bx=i-w+1; }
    }
    return {best,bx};
}

static Ans pack(const vector<Item>& a,int W,int H,bool rot, vector<int> ord, bool flipTie){
    int n=a.size(); vector<int> used(n), dead(n); vector<int> sky(W,0); Ans z;
    while(true){
        int t=-1;
        for(int k:ord) if(!dead[k] && used[k]<a[k].lim){ t=k; break; }
        if(t<0) break;
        int bw=0,bh=0,br=0,bx=-1, by=INT_MAX;
        for(int r=0;r<= (rot && a[t].w!=a[t].h);r++){
            int w=r?a[t].h:a[t].w, h=r?a[t].w:a[t].h;
            auto e=bestSpot(sky,w,h,H);
            if(e.second<0) continue;
            int y=e.first-h;
            // Lowest resulting skyline is a robust tie breaker for either orientation.
            if(e.first<by+bh || (e.first==by+bh && ((flipTie&&r>br)||(!flipTie&&r<br))))
                bw=w,bh=h,br=r,bx=e.second,by=y;
        }
        if(bx<0){ dead[t]=1; continue; } // skyline heights only rise, so it can never fit later
        for(int x=bx;x<bx+bw;x++) sky[x]=by+bh;
        used[t]++; z.val+=a[t].v; z.p.push_back({t,bx,by,br});
    }
    return z;
}

static string esc(const string& s){ string r; for(char c:s){ if(c=='"'||c=='\\') r+='\\'; r+=c; } return r; }
int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    Json j; vector<Item>a; int W=0,H=0; bool allow=false;
    j.need('{');
    while(j.ch()!='}' && j.ch()){
        string key=j.str(); j.need(':');
        if(key=="bin"){
            j.need('{'); while(j.ch()!='}') { string k=j.str(); j.need(':'); if(k=="W")W=j.num(); else if(k=="H")H=j.num(); else if(k=="allow_rotate")allow=j.boolean(); if(j.ch()==',')j.need(','); } j.need('}');
        } else if(key=="items"){
            j.need('['); while(j.ch()!=']'){
                j.need('{'); Item it; while(j.ch()!='}') { string k=j.str(); j.need(':'); if(k=="type")it.id=j.str(); else if(k=="w")it.w=j.num(); else if(k=="h")it.h=j.num(); else if(k=="v")it.v=j.num(); else if(k=="limit")it.lim=j.num(); if(j.ch()==',')j.need(','); } j.need('}'); a.push_back(it); if(j.ch()==',')j.need(',');
            } j.need(']');
        }
        if(j.ch()==',')j.need(',');
    }
    int n=a.size(); vector<int> base(n); iota(base.begin(),base.end(),0);
    vector<pair<long double,int>> den(n);
    for(int i=0;i<n;i++) den[i]={ (long double)a[i].v/(a[i].w*a[i].h),i};
    Ans best; uint32_t seed=1234567u+W*239u+H;
    // Independent density perturbations test whether strict fractional-density priority causes bad skyline fragmentation.
    for(int run=0;run<14;run++){
        vector<int> ord=base;
        mt19937 rng(seed+run*7919u);
        vector<long double> key(n);
        for(int i=0;i<n;i++){
            long double d=den[i].first;
            if(run==1) key[i]=d*(1.0L+0.12L*(long double)min(a[i].w,a[i].h)/max(a[i].w,a[i].h));
            else if(run==2) key[i]=d*(1.0L+0.08L*(long double)(a[i].w*a[i].h)/(W*H));
            // Density-only perturbations cannot test whether large awkward pieces must establish the skyline first.
            else if(run==3) key[i]=(long double)a[i].w*a[i].h;
            else if(run==4) key[i]=max(a[i].w,a[i].h);
            else if(run==5) key[i]=min(a[i].w,a[i].h);
            else if(run==6) key[i]=a[i].v;
            else { long double noise=(long double)((int)(rng()%2001)-1000)/1000.0L*(run<10?0.06L:0.18L); key[i]=d*(1+noise); }
        }
        sort(ord.begin(),ord.end(),[&](int x,int y){ if(key[x]!=key[y]) return key[x]>key[y]; return a[x].v>a[y].v; });
        Ans cur=pack(a,W,H,allow,ord,run&1);
        if(cur.val>best.val) best=move(cur);
    }
    cout << "{\"placements\":[";
    for(size_t i=0;i<best.p.size();i++){ auto p=best.p[i]; if(i) cout<<','; cout<<"{\"type\":\""<<esc(a[p.t].id)<<"\",\"x\":"<<p.x<<",\"y\":"<<p.y<<",\"rot\":"<<p.r<<'}'; }
    cout << "]}\n";
}
