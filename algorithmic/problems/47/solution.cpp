#include <bits/stdc++.h>
using namespace std;

struct Item { string id; int w,h,lim; long long v; };
struct Put { int t,x,y,r; };
struct Answer { vector<Put> p; long long value=0; };

static size_t afterKey(const string& s, const string& key, size_t from=0) {
    size_t q=s.find("\""+key+"\"",from);
    if(q==string::npos) return q;
    q=s.find(':',q+key.size()+2);
    return q==string::npos?q:q+1;
}
static void ws(const string&s,size_t& p){ while(p<s.size() && isspace((unsigned char)s[p])) ++p; }
static long long numberAt(const string&s, size_t p) { ws(s,p); bool neg=false; if(s[p]=='-') neg=true,++p; long long z=0; while(p<s.size() && isdigit((unsigned char)s[p])) z=z*10+s[p++]-'0'; return neg?-z:z; }
static string stringAt(const string&s,size_t p) {
    ws(s,p); string r; if(p>=s.size()||s[p++]!='\"') return r;
    while(p<s.size() && s[p]!='\"') { if(s[p]=='\\' && p+1<s.size()) { ++p; char c=s[p++]; if(c=='n')r+='\n'; else if(c=='t')r+='\t'; else r+=c; } else r+=s[p++]; } return r;
}
static string esc(const string& s) { string r; for(char c:s) { if(c=='\"'||c=='\\') r+='\\'; r+=c; } return r; }

// Returns the lowest skyline position for a rectangle.  The deque makes each
// placement search linear in the (at most 2000) bin width, rather than w*W.
static pair<int,int> bestPlace(const vector<int>& sky, int w, int policy) {
    int W=(int)sky.size(), best=INT_MAX, bx=-1;
    deque<int> q;
    for(int i=0;i<W;i++) {
        while(!q.empty() && sky[q.back()]<=sky[i]) q.pop_back();
        q.push_back(i);
        while(!q.empty() && q.front()<=i-w) q.pop_front();
        if(i>=w-1) {
            int x=i-w+1, y=sky[q.front()];
            bool take=y<best;
            if(y==best) {
                if(policy==1) take=x>bx;
                else if(policy==2) take=abs(2*x+w-W)<abs(2*bx+w-W);
            }
            if(take) best=y,bx=x;
        }
    }
    return {bx,best};
}

static Answer makeAnswer(const vector<Item>& a,int W,int H,bool rot,int mode) {
    int style=mode%6;
    int n=a.size(); vector<int> ord(n); iota(ord.begin(),ord.end(),0);
    auto density=[&](int i){ return (long double)a[i].v/((long double)a[i].w*a[i].h); };
    sort(ord.begin(),ord.end(),[&](int i,int j) {
        long double di=density(i), dj=density(j);
        if(style==3) { if(a[i].v!=a[j].v) return a[i].v>a[j].v; }
        else if(style==4) { // value divided by sqrt(area): favors valuable compact pieces
            long double si=(long double)a[i].v/sqrt((long double)a[i].w*a[i].h);
            long double sj=(long double)a[j].v/sqrt((long double)a[j].w*a[j].h);
            if(fabsl(si-sj)>1e-12L) return si>sj;
        } else if(fabsl(di-dj)>1e-15L) return di>dj;
        long long ai=1LL*a[i].w*a[i].h, aj=1LL*a[j].w*a[j].h;
        if(style==1 && ai!=aj) return ai<aj;
        if(style==2 && ai!=aj) return ai>aj;
        if(style==5 && max(a[i].w,a[i].h)!=max(a[j].w,a[j].h)) return max(a[i].w,a[i].h)<max(a[j].w,a[j].h);
        return a[i].v>a[j].v;
    });
    vector<int> sky(W,0); Answer z;
    int policy=(mode/6)%3;
    for(int id:ord) {
        for(int copy=0;copy<a[id].lim;copy++) {
            struct Choice { int x=-1,y=INT_MAX,w=0,h=0,r=0; } c;
            for(int r=0;r<=(rot?1:0);r++) {
                int w=r?a[id].h:a[id].w, h=r?a[id].w:a[id].h;
                if(w>W || h>H) continue;
                auto [x,y]=bestPlace(sky,w,policy);
                if(y+h>H) continue;
                // Lower bottoms preserve low shelves; then prefer the lower top.
                if(y<c.y || (y==c.y && y+h<c.y+c.h) ||
                   (y==c.y && y+h==c.y+c.h && w>c.w)) c={x,y,w,h,r};
            }
            // Skyline heights only rise, so an item which cannot be placed now
            // cannot become placeable after more rectangles are added.
            if(c.x<0) break;
            for(int x=c.x;x<c.x+c.w;x++) sky[x]=c.y+c.h;
            z.p.push_back({id,c.x,c.y,c.r}); z.value+=a[id].v;
        }
    }
    return z;
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    string s((istreambuf_iterator<char>(cin)),{});
    size_t p=afterKey(s,"W"); if(p==string::npos) return 0; int W=(int)numberAt(s,p);
    p=afterKey(s,"H"); int H=(int)numberAt(s,p);
    p=afterKey(s,"allow_rotate"); ws(s,p); bool rot=s.compare(p,4,"true")==0;
    vector<Item> a; size_t at=0;
    while(true) {
        at=s.find("\"type\"",at); if(at==string::npos) break;
        size_t q=afterKey(s,"type",at); Item e; e.id=stringAt(s,q);
        q=afterKey(s,"w",at); e.w=(int)numberAt(s,q);
        q=afterKey(s,"h",at); e.h=(int)numberAt(s,q);
        q=afterKey(s,"v",at); e.v=numberAt(s,q);
        q=afterKey(s,"limit",at); e.lim=(int)numberAt(s,q);
        a.push_back(e); at=q+1;
    }
    Answer best;
    // Different sort and tie-break choices are cheap diversified skyline packings.
    for(int mode=0;mode<18;mode++) {
        Answer cur=makeAnswer(a,W,H,rot,mode);
        if(cur.value>best.value) best=move(cur);
    }
    cout << "{\"placements\":[";
    for(size_t i=0;i<best.p.size();i++) { const Put& e=best.p[i]; if(i) cout<<',';
        cout << "{\"type\":\""<<esc(a[e.t].id)<<"\",\"x\":"<<e.x<<",\"y\":"<<e.y<<",\"rot\":"<<e.r<<'}';
    }
    cout << "]}\n";
}
