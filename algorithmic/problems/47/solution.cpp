#include <bits/stdc++.h>
using namespace std;

struct JS {
    string s; size_t p=0;
    JS(){ s.assign(istreambuf_iterator<char>(cin), {}); }
    void ws(){ while(p<s.size() && isspace((unsigned char)s[p])) ++p; }
    void ch(char c){ ws(); if(p<s.size() && s[p]==c) ++p; }
    string str(){ ws(); string r; if(p<s.size() && s[p]=='"') ++p; while(p<s.size() && s[p]!='"'){ if(s[p]=='\\' && p+1<s.size()) { ++p; r+=s[p++]; } else r+=s[p++]; } if(p<s.size()) ++p; return r; }
    long long num(){ ws(); long long z=1,r=0; if(p<s.size()&&s[p]=='-'){z=-1;++p;} while(p<s.size()&&isdigit((unsigned char)s[p])) r=r*10+s[p++]-'0'; return z*r; }
    bool boolean(){ ws(); if(s.compare(p,4,"true")==0){p+=4;return true;} p+=5;return false; }
};
struct Item { string name; int w,h,lim; long long v; };
struct Seg { int x,w,h; };
struct Put { int t,x,y,rot; };
struct Ans { long long val=-1; vector<Put> a; };

// The profile is a disjoint horizontal partition of the bin.  Raising an interval
// makes every recorded placement lie below the profile, hence overlap is impossible.
static bool position(const vector<Seg>& q, int W, int iw, int ih, int &bx, int &by) {
    bool ok=false; int bestY=INT_MAX, bestX=INT_MAX;
    for(int i=0;i<(int)q.size();++i){
        int x=q[i].x;
        if(x+iw>W) continue;
        int need=x+iw, top=0;
        for(int j=i;j<(int)q.size() && q[j].x<need;++j) top=max(top,q[j].h);
        if(top<bestY || (top==bestY && x<bestX)) bestY=top,bestX=x,ok=true;
    }
    if(ok) bx=bestX,by=bestY;
    return ok;
}
static void raiseProfile(vector<Seg>& q,int x,int w,int nh){
    int r=x+w; vector<Seg> z;
    for(auto e:q){
        int a=e.x,b=e.x+e.w;
        if(a<x) z.push_back({a,min(b,x)-a,e.h});
        int l=max(a,x), u=min(b,r);
        if(l<u) z.push_back({l,u-l,nh});
        if(b>r) z.push_back({max(a,r),b-max(a,r),e.h});
    }
    vector<Seg> m;
    for(auto e:z) if(e.w){
        if(!m.empty() && m.back().h==e.h && m.back().x+m.back().w==e.x) m.back().w+=e.w;
        else m.push_back(e);
    }
    q.swap(m);
}
static unsigned hsh(const string& x){ unsigned h=2166136261u; for(unsigned char c:x) h=(h^c)*16777619u; return h; }

int main(){
    JS j; j.ws(); if(j.p>=j.s.size()) return 0;
    int W=0,H=0; bool rotOK=false; vector<Item> it;
    j.ch('{');
    while(true){
        j.ws(); if(j.p<j.s.size() && j.s[j.p]=='}'){++j.p;break;}
        string key=j.str(); j.ch(':');
        if(key=="bin"){
            j.ch('{');
            while(true){ j.ws(); if(j.s[j.p]=='}'){++j.p;break;} string k=j.str(); j.ch(':');
                if(k=="W") W=(int)j.num(); else if(k=="H") H=(int)j.num(); else if(k=="allow_rotate") rotOK=j.boolean();
                j.ws(); if(j.p<j.s.size()&&j.s[j.p]==',') ++j.p;
            }
        } else if(key=="items"){
            j.ch('['); j.ws();
            while(j.p<j.s.size() && j.s[j.p]!=']'){
                Item a; a.w=a.h=a.lim=0; a.v=0; j.ch('{');
                while(true){ j.ws(); if(j.s[j.p]=='}'){++j.p;break;} string k=j.str(); j.ch(':');
                    if(k=="type") a.name=j.str(); else if(k=="w") a.w=(int)j.num(); else if(k=="h") a.h=(int)j.num(); else if(k=="v") a.v=j.num(); else if(k=="limit") a.lim=(int)j.num();
                    j.ws(); if(j.p<j.s.size()&&j.s[j.p]==',') ++j.p;
                }
                it.push_back(a); j.ws(); if(j.p<j.s.size()&&j.s[j.p]==',') ++j.p; j.ws();
            }
            if(j.p<j.s.size()&&j.s[j.p]==']') ++j.p;
        }
        j.ws(); if(j.p<j.s.size()&&j.s[j.p]==',') ++j.p;
    }
    int n=it.size(); Ans best;
    // Different score exponents make this a profile/shelf search over distinct item sequences,
    // rather than a residual-rectangle greedy construction.
    for(int trial=0;trial<16;++trial){
        vector<int> ord(n); iota(ord.begin(),ord.end(),0);
        vector<long double> score(n);
        for(int k=0;k<n;++k){
            long double area=(long double)it[k].w*it[k].h;
            long double d=it[k].v/max((long double)1,area);
            long double noise=1.0L + ((int)((hsh(it[k].name)+trial*1103515245u)%1001)-500)/2500.0L;
            if(trial==0) score[k]=d;
            else if(trial==1) score[k]=it[k].v;
            else if(trial==2) score[k]=sqrt(d*it[k].v);
            else if(trial==3) score[k]=d*sqrt((long double)min(it[k].w,it[k].h));
            else score[k]=powl(d, 0.55L+(trial%6)*0.18L)*powl((long double)it[k].v, (trial%3)*0.16L)*noise;
        }
        stable_sort(ord.begin(),ord.end(),[&](int a,int b){return score[a]>score[b];});
        vector<int> left(n); for(int k=0;k<n;++k) left[k]=it[k].lim;
        vector<Seg> profile; profile.push_back({0,W,0}); vector<Put> out; long long val=0;
        for(int step=0;step<1500;++step){
            int pick=-1,px=0,py=0,pr=0;
            // A sequence chooses the first still-available type that has a bottom-left skyline fit.
            for(int id:ord){
                if(left[id]<=0) continue;
                int x0,y0,x1,y1; bool f0=it[id].w<=W&&it[id].h<=H&&position(profile,W,it[id].w,it[id].h,x0,y0)&&y0+it[id].h<=H;
                bool f1=rotOK && it[id].w!=it[id].h && it[id].h<=W&&it[id].w<=H&&position(profile,W,it[id].h,it[id].w,x1,y1)&&y1+it[id].w<=H;
                if(!f0&&!f1) continue;
                pick=id;
                if(f1 && (!f0 || y1<y0 || (y1==y0 && x1<x0))){px=x1;py=y1;pr=1;} else {px=x0;py=y0;pr=0;}
                break;
            }
            if(pick<0) break;
            int ww=pr?it[pick].h:it[pick].w, hh=pr?it[pick].w:it[pick].h;
            out.push_back({pick,px,py,pr}); raiseProfile(profile,px,ww,py+hh); --left[pick]; val+=it[pick].v;
        }
        if(val>best.val){ best.val=val; best.a.swap(out); }
    }
    auto esc=[](const string& s){ string r; for(char c:s){ if(c=='"'||c=='\\') r+='\\'; r+=c; } return r; };
    cout << "{\"placements\":[";
    for(size_t k=0;k<best.a.size();++k){ if(k) cout<<','; auto p=best.a[k]; cout<<"{\"type\":\""<<esc(it[p.t].name)<<"\",\"x\":"<<p.x<<",\"y\":"<<p.y<<",\"rot\":"<<p.rot<<'}'; }
    cout << "]}\n";
}
