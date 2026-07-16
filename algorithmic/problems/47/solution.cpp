#include <bits/stdc++.h>
using namespace std;
struct Item { string id; long long w,h,v,lim; };
struct Pt { long long x,y; };
struct Pl { int t; long long x,y; int r; };
static long long numberAfter(const string&s, const string& key, size_t from=0){
    size_t p=s.find("\""+key+"\"",from); if(p==string::npos) return 0;
    p=s.find(':',p); ++p; while(p<s.size() && isspace((unsigned char)s[p])) ++p;
    long long z=0; bool neg=false; if(p<s.size()&&s[p]=='-') neg=true,++p;
    while(p<s.size()&&isdigit((unsigned char)s[p])) z=z*10+s[p++]-'0'; return neg?-z:z;
}
static string stringAfter(const string&s,const string& key,size_t from){
    size_t p=s.find("\""+key+"\"",from); if(p==string::npos) return "";
    p=s.find(':',p); p=s.find('"',p); size_t q=s.find('"',p+1); return s.substr(p+1,q-p-1);
}
static long long at(const vector<Pt>& a,long long x){
    int l=0,r=(int)a.size(); while(l<r){int m=(l+r)/2; if(a[m].x<=x)l=m+1;else r=m;}
    return a[max(0,l-1)].y;
}
static long long top(const vector<Pt>& a,int k,long long end){
    long long z=a[k].y;
    for(int j=k+1;j<(int)a.size()&&a[j].x<end;j++) z=max(z,a[j].y);
    return z;
}
static void put(vector<Pt>& a,long long x,long long w,long long nh){
    long long e=x+w, after=at(a,e); vector<Pt>b; b.reserve(a.size()+2);
    for(auto p:a) if(p.x<x||p.x>e) b.push_back(p);
    b.push_back({x,nh}); b.push_back({e,after});
    sort(b.begin(),b.end(),[](Pt p,Pt q){return p.x<q.x;});
    vector<Pt> c; for(auto p:b){
        if(!c.empty()&&c.back().x==p.x) c.back()=p;
        else if(c.empty()||c.back().y!=p.y) c.push_back(p);
    } a.swap(c);
}
int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    string s((istreambuf_iterator<char>(cin)),{});
    long long W=numberAfter(s,"W"), H=numberAfter(s,"H");
    size_t bp=s.find("\"allow_rotate\""); bool rotok=bp!=string::npos && s.find("true",bp)<s.find('}',bp);
    vector<Item> it;
    size_t p=0;
    while((p=s.find("\"type\"",p))!=string::npos){
        size_t e=s.find('}',p); Item q;
        q.id=stringAfter(s,"type",p); q.w=numberAfter(s,"w",p); q.h=numberAfter(s,"h",p);
        q.v=numberAfter(s,"v",p); q.lim=numberAfter(s,"limit",p); it.push_back(q); p=e+1;
    }
    auto start=chrono::steady_clock::now();
    auto expired=[&](){return chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now()-start).count()>805;};
    vector<Pl> best; long long bestV=-1; mt19937 rng(712367);
    // Each restart is an independent density-biased skyline construction, not a hybrid packer.
    for(int run=0; !expired(); ++run){
        vector<double> bias(it.size());
        for(int i=0;i<(int)it.size();i++){
            double d=(double)it[i].v/(double)(it[i].w*it[i].h);
            if(run==0) bias[i]=d;
            else if(run==1) bias[i]=(double)it[i].v;
            else if(run==2) bias[i]=d*sqrt((double)it[i].w*it[i].h);
            else bias[i]=d*(0.70+0.60*(double)(rng()%10000)/9999.0);
        }
        vector<Pt> sky={{0,0},{W,0}}; vector<long long> used(it.size()); vector<Pl> out; long long val=0;
        while(!expired()){
            int bt=-1,br=0,bk=-1; long long bx=0,by=0,bh=0; double bs=-1;
            for(int t=0;t<(int)it.size();t++) if(used[t]<it[t].lim){
                for(int rr=0;rr<=(rotok?1:0);rr++){
                    if(rr && it[t].w==it[t].h) continue;
                    long long w=rr?it[t].h:it[t].w, h=rr?it[t].w:it[t].h;
                    if(w>W||h>H) continue;
                    long long lx=0,ly=LLONG_MAX; int lk=-1;
                    for(int k=0;k<(int)sky.size();k++){
                        if((k&127)==0&&expired()) break;
                        long long x=sky[k].x; if(x+w>W) continue;
                        long long y=top(sky,k,x+w); if(y+h>H) continue;
                        if(y<ly || (y==ly && x<lx)){ly=y;lx=x;lk=k;}
                    }
                    if(lk<0) continue;
                    // Mild height penalty distinguishes the useful low-contour placement rule
                    // from a pure density order, while retaining density as the main objective.
                    double sc=bias[t]/(1.0+0.16*(double)ly/max(1LL,H));
                    if(sc>bs || (sc==bs && it[t].v>(bt<0?0:it[bt].v))){bs=sc;bt=t;br=rr;bx=lx;by=ly;bk=lk;bh=h;}
                }
            }
            if(bt<0) break;
            long long w=br?it[bt].h:it[bt].w;
            put(sky,bx,w,by+bh); used[bt]++; val+=it[bt].v; out.push_back({bt,bx,by,br});
        }
        if(val>bestV){bestV=val; best.swap(out);}
    }
    cout<<"{\"placements\":[";
    for(size_t i=0;i<best.size();i++){ if(i) cout<<','; auto &z=best[i];
        cout<<"{\"type\":\""<<it[z.t].id<<"\",\"x\":"<<z.x<<",\"y\":"<<z.y<<",\"rot\":"<<z.r<<'}'; }
    cout<<"]}";
}
