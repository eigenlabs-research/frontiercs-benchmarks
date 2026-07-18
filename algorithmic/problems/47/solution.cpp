#include <bits/stdc++.h>
using namespace std;

struct Item{ string type; long long w,h,v,limit; };
struct Rect{ long long x,y,w,h; };
struct Place{ int idx,rot; long long x,y,w,h; };

string input_data; size_t posi=0;
void ws(){ while(posi<input_data.size() && isspace((unsigned char)input_data[posi])) posi++; }
void ch(char c){ ws(); if(posi>=input_data.size()||input_data[posi]!=c) exit(0); posi++; }
string str(){ ws(); ch('"'); string s; while(posi<input_data.size()){ char c=input_data[posi++]; if(c=='"') break; if(c=='\\' && posi<input_data.size()){ char e=input_data[posi++]; if(e=='n') s+='\n'; else if(e=='t') s+='\t'; else if(e=='r') s+='\r'; else s+=e; } else s+=c; } return s; }
long long num(){ ws(); long long sign=1,val=0; if(input_data[posi]=='-'){sign=-1;posi++;} while(posi<input_data.size()&&isdigit((unsigned char)input_data[posi])) val=val*10+(input_data[posi++]-'0'); return sign*val; }
bool boolean(){ ws(); if(input_data.compare(posi,4,"true")==0){posi+=4;return true;} posi+=5; return false; }

bool intersects(const Rect&a,const Rect&b){ return a.x < b.x+b.w && b.x < a.x+a.w && a.y < b.y+b.h && b.y < a.y+a.h; }
bool contains(const Rect&a,const Rect&b){ return a.x<=b.x && a.y<=b.y && a.x+a.w>=b.x+b.w && a.y+a.h>=b.y+b.h; }

struct Solver{
    long long W,H; bool allow; vector<Item> items; vector<long long> rem; vector<Rect> freeR; vector<Place> places;
    int mode;
    void prune(){
        vector<Rect> v; v.reserve(freeR.size());
        for(auto &r: freeR) if(r.w>0 && r.h>0) v.push_back(r);
        sort(v.begin(), v.end(), [](const Rect&a,const Rect&b){
            if(a.x!=b.x) return a.x<b.x; if(a.y!=b.y) return a.y<b.y; if(a.w!=b.w) return a.w>b.w; return a.h>b.h;
        });
        vector<char> dead(v.size(),0);
        for(size_t i=0;i<v.size();++i) if(!dead[i]){
            for(size_t j=0;j<v.size();++j) if(i!=j && !dead[j]){
                if(contains(v[i],v[j])) dead[j]=1;
            }
        }
        freeR.clear();
        for(size_t i=0;i<v.size();++i) if(!dead[i]) freeR.push_back(v[i]);
        if(freeR.size()>700){
            sort(freeR.begin(), freeR.end(), [](const Rect&a,const Rect&b){ return a.w*a.h > b.w*b.h; });
            freeR.resize(700);
        }
    }
    void addPlacement(const Place&p){
        Rect u{p.x,p.y,p.w,p.h}; vector<Rect> nf; nf.reserve(freeR.size()*2+4);
        for(auto r: freeR){
            if(!intersects(r,u)){ nf.push_back(r); continue; }
            long long rx2=r.x+r.w, ry2=r.y+r.h, ux2=u.x+u.w, uy2=u.y+u.h;
            if(u.x>r.x) nf.push_back({r.x,r.y,u.x-r.x,r.h});
            if(ux2<rx2) nf.push_back({ux2,r.y,rx2-ux2,r.h});
            if(u.y>r.y) nf.push_back({r.x,r.y,r.w,u.y-r.y});
            if(uy2<ry2) nf.push_back({r.x,uy2,r.w,ry2-uy2});
        }
        freeR.swap(nf); prune(); places.push_back(p); rem[p.idx]--;
    }
    double candScore(const Item&it,long long w,long long h,const Rect&r){
        double area=(double)w*h, dens=(double)it.v/area;
        double waste=(double)r.w*r.h-area;
        double shortfit=min(r.w-w,r.h-h), longfit=max(r.w-w,r.h-h);
        if(mode==0) return dens*1e12 - waste*1e-3 - shortfit;
        if(mode==1) return (double)it.v*1e6 - waste - shortfit*1000;
        if(mode==2) return (double)it.v/sqrt(area)*1e9 - waste*1e-2;
        if(mode==3) return dens*1e12 + (double)it.v*10 - (shortfit*1000+longfit);
        if(mode==4) return (double)it.v*1e6/ (1.0+0.15*waste/area) - shortfit*1000;
        return dens*1e12 - (double)(r.w-w)*(r.h-h)*10 - longfit;
    }
    void run(){
        rem.resize(items.size()); for(size_t i=0;i<items.size();++i) rem[i]=items[i].limit;
        freeR={{0,0,W,H}};
        int guard=0;
        while(++guard<200000){
            bool found=false; double best=-1e300; Place bp{}; long long bestY=LLONG_MAX,bestX=LLONG_MAX;
            for(int i=0;i<(int)items.size();++i) if(rem[i]>0){
                for(int rot=0;rot<=(allow?1:0);++rot){
                    long long w=rot?items[i].h:items[i].w, h=rot?items[i].w:items[i].h;
                    if(w>W||h>H) continue;
                    for(auto &r: freeR) if(w<=r.w && h<=r.h){
                        double sc=candScore(items[i],w,h,r);
                        if(!found || sc>best+1e-9 || (fabs(sc-best)<1e-9 && (r.y<bestY || (r.y==bestY && r.x<bestX)))){
                            found=true; best=sc; bp={i,rot,r.x,r.y,w,h}; bestY=r.y; bestX=r.x;
                        }
                    }
                }
            }
            if(!found) break;
            addPlacement(bp);
        }
    }
    long long value() const { long long s=0; for(auto&p:places) s+=items[p.idx].v; return s; }
};

vector<Place> shelfPack(long long W,long long H,bool allow,const vector<Item>&items,int mode){
    int n=items.size(); vector<long long> rem(n); for(int i=0;i<n;i++) rem[i]=items[i].limit;
    vector<Place> out; long long y=0;
    struct O{int i,rot; long long w,h;};
    vector<O> all;
    for(int i=0;i<n;i++) for(int r=0;r<=(allow?1:0);r++){
        long long w=r?items[i].h:items[i].w, h=r?items[i].w:items[i].h;
        if(w<=W && h<=H) all.push_back({i,r,w,h});
    }
    while(y < H){
        long long bestH=0,bestVal=-1; vector<Place> bestRow; vector<long long> bestUse(n,0);
        for(auto base: all) if(rem[base.i]>0 && y+base.h<=H){
            long long rowH=base.h, x=0, val=0; vector<Place> row; vector<long long> used(n,0);
            int guard=0;
            while(x < W && ++guard<10000){
                int bi=-1, brot=0; long long bw=0,bh=0; double bs=-1e300;
                for(auto o: all) if(o.h<=rowH && o.w<=W-x && used[o.i]<rem[o.i]){
                    double dens=(double)items[o.i].v/(double)(o.w*o.h);
                    double tight=(double)(W-x-o.w);
                    double sc;
                    if(mode==0) sc=dens*1e12 - tight;
                    else if(mode==1) sc=(double)items[o.i].v*1e6 - tight*10 - (rowH-o.h)*1000;
                    else if(mode==2) sc=(double)items[o.i].v/(double)o.w*1e9 - (rowH-o.h)*100;
                    else sc=dens*1e12 + (double)items[o.i].v*100 - (rowH-o.h)*10000 - tight;
                    if(sc>bs){ bs=sc; bi=o.i; brot=o.rot; bw=o.w; bh=o.h; }
                }
                if(bi<0) break;
                row.push_back({bi,brot,x,y,bw,bh}); used[bi]++; val += items[bi].v; x += bw;
            }
            if(row.empty()) continue;
            double obj;
            if(mode==1) obj=(double)val;
            else if(mode==2) obj=(double)val/(double)rowH;
            else obj=(double)val/(double)(rowH*W);
            double cur = (bestH==0?-1.0: (mode==1?(double)bestVal:(mode==2?(double)bestVal/(double)bestH:(double)bestVal/(double)(bestH*W))));
            if(bestH==0 || obj>cur){ bestH=rowH; bestVal=val; bestRow.swap(row); bestUse.swap(used); }
        }
        if(bestH==0 || bestRow.empty()) break;
        for(auto &p: bestRow){ out.push_back(p); rem[p.idx]--; }
        y += bestH;
    }
    return out;
}

vector<Place> columnPack(long long W,long long H,bool allow,const vector<Item>&items,int mode){
    int n=items.size(); vector<long long> rem(n); for(int i=0;i<n;i++) rem[i]=items[i].limit;
    vector<Place> out; long long x=0;
    struct O{int i,rot; long long w,h;};
    vector<O> all;
    for(int i=0;i<n;i++) for(int r=0;r<=(allow?1:0);r++){
        long long w=r?items[i].h:items[i].w, h=r?items[i].w:items[i].h;
        if(w<=W && h<=H) all.push_back({i,r,w,h});
    }
    while(x < W){
        long long bestW=0,bestVal=-1; vector<Place> bestCol; vector<long long> bestUse(n,0);
        for(auto base: all) if(rem[base.i]>0 && x+base.w<=W){
            long long colW=base.w, y=0, val=0; vector<Place> col; vector<long long> used(n,0);
            int guard=0;
            while(y < H && ++guard<10000){
                int bi=-1, brot=0; long long bw=0,bh=0; double bs=-1e300;
                for(auto o: all) if(o.w<=colW && o.h<=H-y && used[o.i]<rem[o.i]){
                    double dens=(double)items[o.i].v/(double)(o.w*o.h);
                    double tight=(double)(H-y-o.h);
                    double sc;
                    if(mode==0) sc=dens*1e12 - tight;
                    else if(mode==1) sc=(double)items[o.i].v*1e6 - tight*10 - (colW-o.w)*1000;
                    else if(mode==2) sc=(double)items[o.i].v/(double)o.h*1e9 - (colW-o.w)*100;
                    else sc=dens*1e12 + (double)items[o.i].v*100 - (colW-o.w)*10000 - tight;
                    if(sc>bs){ bs=sc; bi=o.i; brot=o.rot; bw=o.w; bh=o.h; }
                }
                if(bi<0) break;
                col.push_back({bi,brot,x,y,bw,bh}); used[bi]++; val += items[bi].v; y += bh;
            }
            if(col.empty()) continue;
            double obj;
            if(mode==1) obj=(double)val;
            else if(mode==2) obj=(double)val/(double)colW;
            else obj=(double)val/(double)(colW*H);
            double cur = (bestW==0?-1.0: (mode==1?(double)bestVal:(mode==2?(double)bestVal/(double)bestW:(double)bestVal/(double)(bestW*H))));
            if(bestW==0 || obj>cur){ bestW=colW; bestVal=val; bestCol.swap(col); bestUse.swap(used); }
        }
        if(bestW==0 || bestCol.empty()) break;
        for(auto &p: bestCol){ out.push_back(p); rem[p.idx]--; }
        x += bestW;
    }
    return out;
}

long long packValue(const vector<Place>&p,const vector<Item>&items){ long long s=0; for(auto &q:p) s+=items[q.idx].v; return s; }

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    input_data.assign((istreambuf_iterator<char>(cin)), istreambuf_iterator<char>());
    long long W=0,H=0; bool allow=false; vector<Item> items;
    ch('{');
    for(int top=0; top<2; ++top){ if(top) ch(','); string k=str(); ch(':');
        if(k=="bin"){
            ch('{'); for(int i=0;i<3;i++){ if(i) ch(','); string bk=str(); ch(':'); if(bk=="W") W=num(); else if(bk=="H") H=num(); else if(bk=="allow_rotate") allow=boolean(); } ch('}');
        }else if(k=="items"){
            ch('['); ws(); bool first=true; while(posi<input_data.size() && input_data[posi]!=']'){ if(!first) ch(','); first=false; ch('{'); Item it; for(int f=0;f<5;f++){ if(f) ch(','); string ik=str(); ch(':'); if(ik=="type") it.type=str(); else if(ik=="w") it.w=num(); else if(ik=="h") it.h=num(); else if(ik=="v") it.v=num(); else if(ik=="limit") it.limit=num(); } ch('}'); items.push_back(it); ws(); } ch(']');
        }
    }
    // Discard impossible and zero-value types; sort is only for internal search, type strings are preserved.
    vector<Item> filtered; for(auto &it:items) if(it.limit>0 && it.v>0 && ((it.w<=W&&it.h<=H) || (allow&&it.h<=W&&it.w<=H))) filtered.push_back(it);
    items.swap(filtered);
    vector<Place> bestP; long long bestV=-1;
    for(int mode=0; mode<6; ++mode){
        Solver s; s.W=W; s.H=H; s.allow=allow; s.items=items; s.mode=mode; s.run();
        if(s.value()>bestV){ bestV=s.value(); bestP=s.places; }
    }
    // Independent shelf/strip construction: often stronger when a few small high-density
    // pieces should tile long rows instead of being fragmented by maximal-rectangle splits.
    for(int mode=0; mode<4; ++mode){
        vector<Place> p=shelfPack(W,H,allow,items,mode); long long v=packValue(p,items);
        if(v>bestV){ bestV=v; bestP.swap(p); }
    }
    // Dual strip construction with vertical columns. This is the same shelf idea in the
    // orthogonal direction and helps no-rotation cases with tall/narrow profitable pieces.
    for(int mode=0; mode<4; ++mode){
        vector<Place> p=columnPack(W,H,allow,items,mode); long long v=packValue(p,items);
        if(v>bestV){ bestV=v; bestP.swap(p); }
    }
    cout << "{\"placements\":[";
    for(size_t i=0;i<bestP.size();++i){ auto&p=bestP[i]; if(i) cout << ','; cout << "{\"type\":\"" << items[p.idx].type << "\",\"x\":" << p.x << ",\"y\":" << p.y << ",\"rot\":" << p.rot << "}"; }
    cout << "]}\n";
    return 0;
}
