#include <bits/stdc++.h>
using namespace std;
struct Item { string id; int w,h,lim; long long v; };
struct Pl { int t,x,y,r; };
struct State { vector<pair<int,int>> s; vector<int> used; vector<Pl> p; long long val=0, area=0; };
int W,H; bool canrot; vector<Item> it;

static int heightAt(const vector<pair<int,int>>& s,int x){
    int k=0;
    for(int i=0;i<(int)s.size();++i) { if(s[i].first>x) break; k=i; }
    return s[k].second;
}
static int topOver(const vector<pair<int,int>>& s,int x,int e){
    int z=heightAt(s,x);
    for(auto q:s) if(q.first>x && q.first<e) z=max(z,q.second);
    return z;
}
static vector<pair<int,int>> raised(const vector<pair<int,int>>& s,int x,int e,int z){
    vector<pair<int,int>> a;
    for(auto q:s) if(q.first<x) a.push_back(q);
    a.push_back({x,z});
    a.push_back({e,heightAt(s,e)});
    for(auto q:s) if(q.first>e) a.push_back(q);
    vector<pair<int,int>> b;
    for(auto q:a) {
        if(!b.empty() && b.back().first==q.first) b.back().second=q.second;
        else b.push_back(q);
    }
    vector<pair<int,int>> c;
    for(auto q:b) if(c.empty() || c.back().second!=q.second || q.first==W) c.push_back(q);
    if(c.empty() || c.back().first!=W) c.push_back({W,heightAt(s,W)});
    return c;
}
static string signature(const State& a){
    string z;
    for(int x:a.used) z+=to_string(x)+",";
    z+="/";
    for(auto q:a.s) z+=to_string(q.first)+":"+to_string(q.second)+",";
    return z;
}
static long double key(const State& a, long double beta, long double bestden){
    // beta is deliberately small: it only breaks near-value ties in favor of filled area.
    return (long double)a.val + beta*bestden*(long double)a.area;
}
static State searchBeam(long double beta) {
    long double bd=0;
    for(auto &q:it) bd=max(bd,(long double)q.v/q.w/q.h);
    State root; root.s={{0,0},{W,0}}; root.used.assign(it.size(),0);
    vector<State> beam(1,root); State answer=root;
    const int BW=16, STEPS=260;
    for(int step=0;step<STEPS && !beam.empty();++step){
        vector<State> next;
        for(const State& a:beam){
            for(int ti=0;ti<(int)it.size();++ti){
                if(a.used[ti]>=it[ti].lim || it[ti].v==0) continue;
                for(int r=0;r<=1;r++){
                    if(r && (!canrot || it[ti].w==it[ti].h)) continue;
                    int ww=r?it[ti].h:it[ti].w, hh=r?it[ti].w:it[ti].h;
                    for(auto q:a.s){
                        int x=q.first;
                        if(x+ww>W) continue;
                        int y=topOver(a.s,x,x+ww);
                        if(y+hh>H) continue;
                        State b=a;
                        b.s=raised(a.s,x,x+ww,y+hh);
                        ++b.used[ti]; b.p.push_back({ti,x,y,r});
                        b.val+=it[ti].v; b.area+=(long long)ww*hh;
                        next.push_back(move(b));
                    }
                }
            }
        }
        if(next.empty()) break;
        sort(next.begin(),next.end(),[&](const State&a,const State&b){
            long double ka=key(a,beta,bd), kb=key(b,beta,bd);
            if(ka!=kb) return ka>kb;
            return a.val>b.val;
        });
        vector<State> keep; unordered_set<string> seen;
        for(auto &a:next){
            string h=signature(a);
            if(seen.insert(h).second) {
                if(a.val>answer.val) answer=a;
                keep.push_back(move(a));
                if((int)keep.size()==BW) break;
            }
        }
        beam.swap(keep);
    }
    return answer;
}
static string esc(const string& x){ string r; for(char c:x){ if(c=='"'||c=='\\') r+='\\'; r+=c; } return r; }
int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    string in((istreambuf_iterator<char>(cin)),{});
    smatch m;
    auto num=[&](const string& k)->long long { regex r("\\\""+k+"\\\"\\s*:\\s*(-?[0-9]+)"); if(regex_search(in,m,r)) return stoll(m[1]); return 0; };
    W=(int)num("W"); H=(int)num("H");
    regex br("\\\"allow_rotate\\\"\\s*:\\s*(true|false)"); canrot=regex_search(in,m,br)&&m[1]=="true";
    regex ob("\\{[^{}]*\\}");
    for(sregex_iterator a(in.begin(),in.end(),ob),e;a!=e;++a){
        string o=a->str();
        regex tr("\\\"type\\\"\\s*:\\s*\\\"([^\\\"]*)\\\"");
        smatch z; if(!regex_search(o,z,tr)) continue;
        auto get=[&](string k)->long long { regex r("\\\""+k+"\\\"\\s*:\\s*(-?[0-9]+)"); smatch q; return regex_search(o,q,r)?stoll(q[1]):0; };
        it.push_back({z[1],(int)get("w"),(int)get("h"),(int)get("limit"),get("v")});
    }
    State best; best.val=-1;
    // A nonzero width is the discriminating mechanism: alternatives survive early commitments.
    State q=searchBeam(0.025L); if(q.val>best.val) best=move(q);
    cout << "{\"placements\":[";
    for(int i=0;i<(int)best.p.size();++i){ auto q=best.p[i]; if(i) cout<<','; cout<<"{\"type\":\""<<esc(it[q.t].id)<<"\",\"x\":"<<q.x<<",\"y\":"<<q.y<<",\"rot\":"<<q.r<<'}'; }
    cout << "]}\n";
}
