#include <bits/stdc++.h>
using namespace std;

struct Item { string id; int w,h; long long v; int lim; };
struct R { int x,y,w,h; };
struct P { int t,x,y,rot,w,h; };

struct Parser {
    string s; size_t p=0;
    Parser(string z):s(move(z)){}
    void ws(){ while(p<s.size() && isspace((unsigned char)s[p])) ++p; }
    char ch(){ ws(); return p<s.size()?s[p]:0; }
    void take(char c){ ws(); if(p<s.size() && s[p]==c) ++p; }
    string str(){ ws(); ++p; string r; while(p<s.size() && s[p]!='"'){ if(s[p]=='\\' && p+1<s.size()) ++p; r+=s[p++]; } if(p<s.size()) ++p; return r; }
    long long num(){ ws(); int sg=1; if(s[p]=='-'){sg=-1;++p;} long long x=0; while(p<s.size()&&isdigit((unsigned char)s[p])) x=x*10+s[p++]-'0'; return sg*x; }
    bool boolean(){ ws(); bool r=s.compare(p,4,"true")==0; p+=r?4:5; return r; }
};

static bool intersects(const R&a,int x,int y,int w,int h){
    return x<a.x+a.w && a.x<x+w && y<a.y+a.h && a.y<y+h;
}
static void prune(vector<R>& a){
    vector<char> bad(a.size());
    for(int i=0;i<(int)a.size();++i) if(!bad[i]) for(int j=0;j<(int)a.size();++j) if(i!=j && !bad[j]) {
        if(a[i].x>=a[j].x && a[i].y>=a[j].y && a[i].x+a[i].w<=a[j].x+a[j].w && a[i].y+a[i].h<=a[j].y+a[j].h) { bad[i]=1; break; }
    }
    vector<R> b; b.reserve(a.size());
    for(int i=0;i<(int)a.size();++i) if(!bad[i] && a[i].w>0 && a[i].h>0) b.push_back(a[i]);
    a.swap(b);
}

static vector<P> pack(const vector<Item>& it,int W,int H,bool allow,const vector<int>& order,int mode){
    vector<R> freeR(1,{0,0,W,H}); vector<P> out;
    for(int t:order){
        for(int cnt=0;cnt<it[t].lim;cnt++){
            int bi=-1, br=0; long long bk=LLONG_MAX, bk2=LLONG_MAX;
            for(int i=0;i<(int)freeR.size();++i) for(int rot=0;rot<=allow;rot++){
                int w=rot?it[t].h:it[t].w, h=rot?it[t].w:it[t].h;
                R q=freeR[i]; if(w>q.w||h>q.h) continue;
                long long a,b;
                if(mode==0){ a=min(q.w-w,q.h-h); b=max(q.w-w,q.h-h); }
                else if(mode==1){ a=1LL*(q.w-w)*(q.h-h); b=min(q.w-w,q.h-h); }
                else { a=max(q.w-w,q.h-h); b=1LL*(q.w-w)*(q.h-h); }
                if(a<bk || (a==bk && (b<bk2 || (b==bk2 && (q.y<freeR[bi<0?i:bi].y || (q.y==freeR[bi<0?i:bi].y && q.x<freeR[bi<0?i:bi].x)))))) bi=i,br=rot,bk=a,bk2=b;
            }
            if(bi<0) break;
            int x=freeR[bi].x,y=freeR[bi].y,w=br?it[t].h:it[t].w,h=br?it[t].w:it[t].h;
            out.push_back({t,x,y,br,w,h});
            vector<R> nf; nf.reserve(freeR.size()+4);
            for(R q:freeR){
                if(!intersects(q,x,y,w,h)){ nf.push_back(q); continue; }
                if(x>q.x) nf.push_back({q.x,q.y,x-q.x,q.h});
                if(x+w<q.x+q.w) nf.push_back({x+w,q.y,q.x+q.w-x-w,q.h});
                if(y>q.y) nf.push_back({q.x,q.y,q.w,y-q.y});
                if(y+h<q.y+q.h) nf.push_back({q.x,y+h,q.w,q.y+q.h-y-h});
            }
            freeR.swap(nf); prune(freeR);
            // A very fragmented residual has negligible utility; retaining its largest regions
            // keeps the construction bounded on adversarial tiny-rectangle inputs.
            if(freeR.size()>1800){
                nth_element(freeR.begin(),freeR.begin()+1200,freeR.end(),[](const R&a,const R&b){return 1LL*a.w*a.h>1LL*b.w*b.h;});
                freeR.resize(1200);
            }
        }
    }
    return out;
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    string input((istreambuf_iterator<char>(cin)),{}); Parser q(input);
    int W=0,H=0; bool allow=false; vector<Item> it;
    q.take('{');
    while(q.ch()!='}' && q.ch()){
        string key=q.str(); q.take(':');
        if(key=="bin"){
            q.take('{'); while(q.ch()!='}') { string k=q.str(); q.take(':'); if(k=="W") W=q.num(); else if(k=="H") H=q.num(); else allow=q.boolean(); if(q.ch()==',')q.take(','); } q.take('}');
        } else {
            q.take('['); while(q.ch()!=']'){
                q.take('{'); Item z; while(q.ch()!='}') { string k=q.str(); q.take(':'); if(k=="type") z.id=q.str(); else if(k=="w")z.w=q.num(); else if(k=="h")z.h=q.num(); else if(k=="v")z.v=q.num(); else z.lim=q.num(); if(q.ch()==',')q.take(','); } q.take('}'); it.push_back(z); if(q.ch()==',')q.take(',');
            } q.take(']');
        }
        if(q.ch()==',')q.take(',');
    }
    vector<vector<int>> orders;
    vector<int> a(it.size()); iota(a.begin(),a.end(),0);
    auto makeorder=[&](int kind){ vector<int> z=a; sort(z.begin(),z.end(),[&](int i,int j){
        long long ai=1LL*it[i].w*it[i].h, aj=1LL*it[j].w*it[j].h;
        if(kind==0) { __int128 l=(__int128)it[i].v*aj, r=(__int128)it[j].v*ai; if(l!=r)return l>r; }
        if(kind==1 && it[i].v!=it[j].v)return it[i].v>it[j].v;
        if(kind==2) { long long li=it[i].v*min(it[i].w,it[i].h), lj=it[j].v*min(it[j].w,it[j].h); if(li!=lj)return li>lj; }
        if(kind==3 && ai!=aj)return ai<aj;
        if(ai!=aj)return ai>aj; return i<j;
    }); return z; };
    for(int k=0;k<4;k++) orders.push_back(makeorder(k));
    vector<P> best; long long bv=-1;
    for(int k=0;k<(int)orders.size();k++) for(int mode=0;mode<3;mode++){
        // The density order is tried with all placement tie-breaks; alternate orders receive
        // the strongest short-side rule, limiting runtime while retaining diverse layouts.
        if(k && mode) continue;
        vector<P> cur=pack(it,W,H,allow,orders[k],mode); long long val=0; for(auto&p:cur) val+=it[p.t].v;
        if(val>bv) bv=val,best.swap(cur);
    }
    cout << "{\"placements\":[";
    for(size_t i=0;i<best.size();++i){ if(i) cout<<','; auto&p=best[i]; cout<<"{\"type\":\""<<it[p.t].id<<"\",\"x\":"<<p.x<<",\"y\":"<<p.y<<",\"rot\":"<<p.rot<<'}'; }
    cout << "]}\n";
}
