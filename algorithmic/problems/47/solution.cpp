#include <bits/stdc++.h>
using namespace std;

struct Item { string id; long long w,h,v,lim; };
struct Pl { int t; long long x,y; int rot; };
struct Seg { long long l,r,y; };

static string getString(const string& s, const string& key) {
    size_t p=s.find("\""+key+"\"");
    if(p==string::npos) return "";
    p=s.find(':',p)+1; while(p<s.size() && isspace((unsigned char)s[p])) ++p;
    if(p>=s.size() || s[p]!='\"') return "";
    ++p; string z;
    while(p<s.size() && s[p]!='\"') { if(s[p]=='\\' && p+1<s.size()) ++p; z+=s[p++]; }
    return z;
}
static long long getNum(const string& s, const string& key) {
    size_t p=s.find("\""+key+"\"");
    if(p==string::npos) return 0;
    p=s.find(':',p)+1; while(p<s.size() && isspace((unsigned char)s[p])) ++p;
    long long sign=1, x=0; if(s[p]=='-') sign=-1,++p;
    while(p<s.size() && isdigit((unsigned char)s[p])) x=x*10+(s[p++]-'0');
    return sign*x;
}
static bool getBool(const string& s, const string& key) {
    size_t p=s.find("\""+key+"\""); if(p==string::npos) return false;
    p=s.find(':',p)+1; while(p<s.size() && isspace((unsigned char)s[p])) ++p;
    return s.compare(p,4,"true")==0;
}

// A skyline is a compact representation of the already occupied region.  Each new
// rectangle is put at its lowest possible left edge, so all produced layouts are valid.
struct Skyline {
    long long W,H; vector<Seg> a;
    Skyline(long long w,long long h):W(w),H(h) { a.push_back({0,w,0}); }
    bool place(long long w,long long h,long long &bx,long long &by) const {
        bool ok=false; bx=by=0;
        for(int i=0;i<(int)a.size();++i) {
            long long x=a[i].l; if(x+w>W) continue;
            long long top=0;
            for(int j=i;j<(int)a.size() && a[j].l<x+w;++j) top=max(top,a[j].y);
            if(top+h>H) continue;
            // The second tie break tends to retain broad, low shelves.
            if(!ok || top<by || (top==by && x<bx)) ok=true,bx=x,by=top;
        }
        return ok;
    }
    void add(long long x,long long y,long long w,long long h) {
        long long q=x+w; vector<Seg> b;
        for(auto s:a) {
            if(s.r<=x || s.l>=q) b.push_back(s);
            else {
                if(s.l<x) b.push_back({s.l,x,s.y});
                if(s.r>q) b.push_back({q,s.r,s.y});
            }
        }
        b.push_back({x,q,y+h});
        sort(b.begin(),b.end(),[](const Seg&u,const Seg&v){return u.l<v.l;});
        a.clear();
        for(auto s:b) {
            if(s.l>=s.r) continue;
            if(!a.empty() && a.back().r==s.l && a.back().y==s.y) a.back().r=s.r;
            else a.push_back(s);
        }
    }
};

static string esc(const string& s) {
    string r; for(char c:s) { if(c=='\"'||c=='\\') r+='\\'; r+=c; } return r;
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    string in((istreambuf_iterator<char>(cin)), {});
    size_t bp=in.find("\"bin\"");
    size_t b0=in.find('{',bp), b1=in.find('}',b0);
    string bin=in.substr(b0,b1-b0+1);
    long long W=getNum(bin,"W"), H=getNum(bin,"H"); bool rotate=getBool(bin,"allow_rotate");
    vector<Item> it;
    size_t ip=in.find("\"items\""); if(ip!=string::npos) {
        size_t p=in.find('[',ip)+1;
        while(p<in.size()) {
            while(p<in.size() && (isspace((unsigned char)in[p])||in[p]==',')) ++p;
            if(p>=in.size()||in[p]==']') break;
            if(in[p]!='{') {++p; continue;}
            size_t st=p; int dep=0;
            do { if(in[p]=='{') ++dep; else if(in[p]=='}') --dep; ++p; } while(p<in.size()&&dep);
            string o=in.substr(st,p-st);
            it.push_back({getString(o,"type"),getNum(o,"w"),getNum(o,"h"),getNum(o,"v"),getNum(o,"limit")});
        }
    }
    vector<Pl> answer; long long best=-1;
    // Different density exponents are still the same skyline-greedy mechanism; they
    // make it robust when value is correlated with either area or number of pieces.
    vector<double> powers={0.0,.20,.40,.60,.80,1.0,1.25,1.5,1.8,2.1};
    for(int pass=0;pass<(int)powers.size();++pass) {
        double pw=powers[pass]; vector<int> ord(it.size()); iota(ord.begin(),ord.end(),0);
        sort(ord.begin(),ord.end(),[&](int u,int v){
            double au=(double)it[u].w*it[u].h, av=(double)it[v].w*it[v].h;
            double ku=it[u].v/pow(max(1.0,au),pw), kv=it[v].v/pow(max(1.0,av),pw);
            if(fabs(ku-kv)>1e-9) return ku>kv;
            // Alternate harmless deterministic ties, giving distinct layouts on equal densities.
            if(pass&1) return max(it[u].w,it[u].h)<max(it[v].w,it[v].h);
            return it[u].v>it[v].v;
        });
        Skyline sky(W,H); vector<Pl> cur; long long val=0;
        for(int t:ord) for(long long k=0;k<it[t].lim;++k) {
            long long x0,y0,x1,y1; bool f0=sky.place(it[t].w,it[t].h,x0,y0), f1=false;
            if(rotate && it[t].w!=it[t].h) f1=sky.place(it[t].h,it[t].w,x1,y1);
            if(!f0 && !f1) break; // a skyline only rises, so this type cannot fit later either
            int r=0; long long x=x0,y=y0,w=it[t].w,h=it[t].h;
            if(f1 && (!f0 || y1<y0 || (y1==y0 && x1<x0))) r=1,x=x1,y=y1,w=it[t].h,h=it[t].w;
            sky.add(x,y,w,h); cur.push_back({t,x,y,r}); val+=it[t].v;
        }
        if(val>best) best=val,answer.swap(cur);
    }
    cout << "{\"placements\":[";
    for(size_t i=0;i<answer.size();++i) { if(i) cout<<','; auto p=answer[i];
        cout<<"{\"type\":\""<<esc(it[p.t].id)<<"\",\"x\":"<<p.x<<",\"y\":"<<p.y<<",\"rot\":"<<p.rot<<'}';
    }
    cout << "]}\n";
}
