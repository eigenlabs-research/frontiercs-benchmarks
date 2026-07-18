#include <bits/stdc++.h>
using namespace std;

struct Item{string type; int w,h; long long v; int lim;};
struct Rect{int x,y,w,h;};
struct Place{int id,x,y,rot,w,h;};

struct Parser{
    string s; size_t p=0;
    Parser(){ std::ostringstream ss; ss<<cin.rdbuf(); s=ss.str(); }
    void ws(){ while(p<s.size() && isspace((unsigned char)s[p])) p++; }
    void ch(char c){ ws(); if(p>=s.size()||s[p]!=c) exit(0); p++; }
    string str(){ ws(); ch('"'); string r; while(p<s.size()){ char c=s[p++]; if(c=='"') break; if(c=='\\'&&p<s.size()){ char e=s[p++]; if(e=='"'||e=='\\'||e=='/') r+=e; else if(e=='n') r+='\n'; else if(e=='t') r+='\t'; else r+=e; } else r+=c; } return r; }
    long long num(){ ws(); long long sign=1,v=0; if(s[p]=='-'){sign=-1;p++;} while(p<s.size()&&isdigit((unsigned char)s[p])) v=v*10+(s[p++]-'0'); return sign*v; }
    bool boolean(){ ws(); if(s.compare(p,4,"true")==0){p+=4; return true;} p+=5; return false; }
};

static inline bool inter(const Rect&a,const Rect&b){return a.x<b.x+b.w && b.x<a.x+a.w && a.y<b.y+b.h && b.y<a.y+a.h;}
static inline bool contains(const Rect&a,const Rect&b){return a.x<=b.x && a.y<=b.y && a.x+a.w>=b.x+b.w && a.y+a.h>=b.y+b.h;}

struct Solver{
    int W,H; bool allow; vector<Item> it;
    vector<Place> best;
    long long bestVal=0;

    void prune(vector<Rect>& fr){
        vector<Rect> v; v.reserve(fr.size());
        for(auto &r:fr) if(r.w>0&&r.h>0) v.push_back(r);
        sort(v.begin(), v.end(), [](const Rect&a,const Rect&b){
            if(a.x!=b.x) return a.x<b.x; if(a.y!=b.y) return a.y<b.y; if(a.w!=b.w) return a.w>b.w; return a.h>b.h;
        });
        vector<char> rem(v.size(),0);
        for(size_t i=0;i<v.size();++i) if(!rem[i]){
            for(size_t j=0;j<v.size();++j) if(i!=j&&!rem[j]){
                if(contains(v[i],v[j])) rem[j]=1;
            }
        }
        fr.clear();
        for(size_t i=0;i<v.size();++i) if(!rem[i]) fr.push_back(v[i]);
        if(fr.size()>4500){
            sort(fr.begin(), fr.end(), [](const Rect&a,const Rect&b){return 1LL*a.w*a.h>1LL*b.w*b.h;});
            fr.resize(4500);
        }
    }

    double candScore(int var, const Item&it, int rw,int rh,const Rect&f){
        double area=1.0*rw*rh, den=it.v/area;
        int waste=f.w*f.h-rw*rh;
        int ss=min(f.w-rw,f.h-rh), ls=max(f.w-rw,f.h-rh);
        if(var==0) return den*1e9 - waste*1e-3 - ss*1e-2;
        if(var==1) return (double)it.v - waste*0.01;
        if(var==2) return it.v/sqrt(area)*1e5 - ss*10 - ls;
        if(var==3) return it.v/(double)rw*1e6 - f.y*1e-3 - waste*1e-4;
        if(var==4) return it.v/(double)rh*1e6 - waste*1e-4;
        return den*1e9 + (double)it.v*1e-3 - (f.w-rw+f.h-rh);
    }

    vector<Place> runMaxRects(int var, chrono::steady_clock::time_point deadline){
        int n=it.size(); vector<int> rem(n); for(int i=0;i<n;i++) rem[i]=it[i].lim;
        vector<Rect> freeR={{0,0,W,H}}; vector<Place> out; out.reserve(5000);
        while(chrono::steady_clock::now()<deadline){
            int bf=-1, bi=-1, brot=0, brw=0, brh=0; double bs=-1e100; Rect bfr{};
            for(int fi=0; fi<(int)freeR.size(); ++fi){
                Rect f=freeR[fi];
                for(int i=0;i<n;i++) if(rem[i]>0){
                    for(int rr=0; rr<(allow?2:1); ++rr){
                        int rw= rr?it[i].h:it[i].w, rh= rr?it[i].w:it[i].h;
                        if(rw>f.w||rh>f.h) continue;
                        if(rr && it[i].w==it[i].h) continue;
                        double sc=candScore(var,it[i],rw,rh,f);
                        // Prefer low-left positions for equal choices to keep skyline compact.
                        sc -= f.y*1e-5 + f.x*1e-6;
                        if(sc>bs){bs=sc; bf=fi; bi=i; brot=rr; brw=rw; brh=rh; bfr=f;}
                    }
                }
            }
            if(bi<0) break;
            Rect used{bfr.x,bfr.y,brw,brh};
            out.push_back({bi,used.x,used.y,brot,brw,brh}); rem[bi]--;
            vector<Rect> nf; nf.reserve(freeR.size()+8);
            for(auto f:freeR){
                if(!inter(f,used)){ nf.push_back(f); continue; }
                int fx2=f.x+f.w, fy2=f.y+f.h, ux2=used.x+used.w, uy2=used.y+used.h;
                if(used.x>f.x) nf.push_back({f.x,f.y,used.x-f.x,f.h});
                if(ux2<fx2) nf.push_back({ux2,f.y,fx2-ux2,f.h});
                if(used.y>f.y) nf.push_back({f.x,f.y,f.w,used.y-f.y});
                if(uy2<fy2) nf.push_back({f.x,uy2,f.w,fy2-uy2});
            }
            freeR.swap(nf); prune(freeR);
            if(out.size()>25000) break;
        }
        return out;
    }

    vector<Place> runShelves(int var){
        int n=it.size(); vector<int> rem(n); for(int i=0;i<n;i++) rem[i]=it[i].lim;
        vector<Place> out; int y=0;
        while(y<H){
            int bestI=-1,brot=0,bh=0; double bs=-1e100;
            for(int i=0;i<n;i++) if(rem[i]) for(int rr=0; rr<(allow?2:1); ++rr){
                int rw=rr?it[i].h:it[i].w, rh=rr?it[i].w:it[i].h; if(rw>W||y+rh>H) continue;
                double sc = (var==0? it[i].v/(double)(rw*rh) : (var==1? it[i].v/(double)rh : it[i].v/(double)rw));
                if(sc>bs){bs=sc;bestI=i;brot=rr;bh=rh;}
            }
            if(bestI<0) break;
            int x=0;
            while(x<W){
                int ci=-1,cr=0,cw=0,ch=0; double cs=-1e100;
                for(int i=0;i<n;i++) if(rem[i]) for(int rr=0; rr<(allow?2:1); ++rr){
                    int rw=rr?it[i].h:it[i].w, rh=rr?it[i].w:it[i].h;
                    if(rw<=W-x && rh<=bh){
                        double sc=(var==0? it[i].v/(double)(rw*rh) : (var==1? it[i].v/(double)rw : (double)it[i].v));
                        if(sc>cs){cs=sc;ci=i;cr=rr;cw=rw;ch=rh;}
                    }
                }
                if(ci<0) break;
                out.push_back({ci,x,y,cr,cw,ch}); rem[ci]--; x+=cw;
            }
            y+=bh;
        }
        return out;
    }

    vector<Place> runDPShelves(int var, chrono::steady_clock::time_point deadline){
        struct Opt{int id,rot,w,h;};
        struct Pack{long long val=0; int usedW=0; vector<pair<int,int>> cnt;};
        int n=it.size(); vector<int> rem(n); for(int i=0;i<n;i++) rem[i]=it[i].lim;
        vector<int> heights;
        for(int i=0;i<n;i++){
            if(it[i].w<=W && it[i].h<=H) heights.push_back(it[i].h);
            if(allow && it[i].h<=W && it[i].w<=H && it[i].w!=it[i].h) heights.push_back(it[i].w);
        }
        sort(heights.begin(), heights.end()); heights.erase(unique(heights.begin(), heights.end()), heights.end());
        vector<Place> out; int y=0;
        while(y<H && chrono::steady_clock::now()<deadline){
            Pack bestPack; vector<Opt> bestOpts; int bestH=-1; double bestScore=-1e100;
            for(int sh: heights){
                if(y+sh>H) continue;
                vector<Opt> opts; opts.reserve(n);
                for(int i=0;i<n;i++) if(rem[i]>0){
                    int bw=INT_MAX, bh=0, br=0;
                    if(it[i].h<=sh && it[i].w<=W){ bw=it[i].w; bh=it[i].h; br=0; }
                    if(allow && it[i].w<=sh && it[i].h<=W && it[i].h<bw){ bw=it[i].h; bh=it[i].w; br=1; }
                    if(bw!=INT_MAX) opts.push_back({i,br,bw,bh});
                }
                if(opts.empty()) continue;
                struct Group{int opt,k,w; long long v;}; vector<Group> groups;
                for(int oi=0; oi<(int)opts.size(); ++oi){
                    int left=rem[opts[oi].id], p2=1;
                    while(left>0){ int k=min(p2,left); groups.push_back({oi,k,opts[oi].w*k,it[opts[oi].id].v*k}); left-=k; p2<<=1; }
                }
                const long long NEG=LLONG_MIN/4;
                vector<long long> dp(W+1,NEG); vector<int> preW(W+1,-1), preG(W+1,-1); dp[0]=0;
                for(int gi=0; gi<(int)groups.size(); ++gi){
                    auto &g=groups[gi];
                    for(int x=W; x>=g.w; --x) if(dp[x-g.w]!=NEG && dp[x-g.w]+g.v>dp[x]){
                        dp[x]=dp[x-g.w]+g.v; preW[x]=x-g.w; preG[x]=gi;
                    }
                }
                int bx=0; for(int x=1;x<=W;x++) if(dp[x]>dp[bx]) bx=x;
                if(dp[bx]<=0) continue;
                double sc;
                if(var==0) sc = (double)dp[bx]/sh;
                else if(var==1) sc = (double)dp[bx]/max(1,bx) + (double)dp[bx]/sh*0.002;
                else sc = (double)dp[bx] - (W-bx)*1000.0 - sh*10.0;
                if(sc>bestScore){
                    bestScore=sc; bestH=sh; bestPack.val=dp[bx]; bestPack.usedW=bx; bestPack.cnt.clear(); bestOpts=opts;
                    vector<int> cc(opts.size(),0); int cur=bx;
                    while(cur>0 && preG[cur]>=0){ auto &g=groups[preG[cur]]; cc[g.opt]+=g.k; cur=preW[cur]; }
                    for(int oi=0; oi<(int)opts.size(); ++oi) if(cc[oi]) bestPack.cnt.push_back({oi,cc[oi]});
                }
            }
            if(bestH<0 || bestPack.cnt.empty()) break;
            int x=0;
            for(auto pr: bestPack.cnt){
                Opt o=bestOpts[pr.first]; int c=min(pr.second, rem[o.id]);
                for(int k=0;k<c;k++){ out.push_back({o.id,x,y,o.rot,o.w,o.h}); x+=o.w; }
                rem[o.id]-=c;
            }
            y += bestH;
        }
        return out;
    }

    long long val(const vector<Place>&v){ long long s=0; for(auto&p:v) s+=it[p.id].v; return s; }
    void consider(vector<Place> v){ long long x=val(v); if(x>bestVal){bestVal=x; best.swap(v);} }
    void solve(){
        auto start=chrono::steady_clock::now(); auto deadline=start+chrono::milliseconds(880);
        for(int var=0; var<3 && chrono::steady_clock::now()<deadline; ++var) consider(runDPShelves(var,deadline));
        for(int var=0; var<6 && chrono::steady_clock::now()<deadline; ++var) consider(runMaxRects(var,deadline));
        for(int var=0; var<3; ++var) consider(runShelves(var));
    }
};

string esc(const string& s){ string r; for(char c:s){ if(c=='"'||c=='\\') {r+='\\'; r+=c;} else r+=c; } return r; }

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    Parser p; Solver sol;
    p.ch('{');
    for(int top=0; top<2; ++top){ if(top) p.ch(','); string k=p.str(); p.ch(':');
        if(k=="bin"){
            p.ch('{'); for(int i=0;i<3;i++){ if(i)p.ch(','); string bk=p.str(); p.ch(':'); if(bk=="W") sol.W=p.num(); else if(bk=="H") sol.H=p.num(); else sol.allow=p.boolean(); } p.ch('}');
        }else{
            p.ch('['); p.ws(); bool first=true; while(p.p<p.s.size() && p.s[p.p]!=']'){
                if(!first) p.ch(','); first=false; p.ch('{'); Item it; for(int f=0;f<5;f++){ if(f)p.ch(','); string ik=p.str(); p.ch(':'); if(ik=="type") it.type=p.str(); else if(ik=="w") it.w=p.num(); else if(ik=="h") it.h=p.num(); else if(ik=="v") it.v=p.num(); else if(ik=="limit") it.lim=p.num(); } p.ch('}'); sol.it.push_back(it); p.ws(); }
            p.ch(']');
        }
    }
    sol.solve();
    cout << "{\"placements\":[";
    for(size_t i=0;i<sol.best.size();++i){ auto &q=sol.best[i]; if(i) cout << ','; cout << "{\"type\":\"" << esc(sol.it[q.id].type) << "\",\"x\":" << q.x << ",\"y\":" << q.y << ",\"rot\":" << q.rot << "}"; }
    cout << "]}\n";
    return 0;
}
