#include <bits/stdc++.h>
using namespace std;

struct Item { string id; int w,h,lim; long long v; };
struct Pl { int t,x,y,r; };
struct Rect { int x,y,w,h; };
struct Ans { long long value=-1; vector<Pl> p; };

// Minimal JSON reader: the input format contains only objects, arrays, strings, integers, and bools.
struct Json {
    string s; size_t p=0;
    explicit Json(string z):s(move(z)) {}
    void ws(){ while(p<s.size() && isspace((unsigned char)s[p])) ++p; }
    char get(){ ws(); return p<s.size()?s[p++]:0; }
    string str(){ ws(); string r; if(p>=s.size() || s[p++]!='"') return r; while(p<s.size() && s[p]!='"'){ if(s[p]=='\\' && p+1<s.size()) ++p; r+=s[p++]; } if(p<s.size()) ++p; return r; }
    long long num(){ ws(); long long sign=1,x=0; if(p<s.size() && s[p]=='-') sign=-1,++p; while(p<s.size() && isdigit((unsigned char)s[p])) x=x*10+s[p++]-'0'; return sign*x; }
    bool boolean(){ ws(); if(s.compare(p,4,"true")==0){p+=4;return true;} p+=5; return false; }
    void skip(){ ws(); if(p>=s.size())return; if(s[p]=='"'){str();return;} if(s[p]=='{'){get(); while(1){ws();if(s[p]=='}'){++p;return;}str();get();skip();ws();if(s[p]==',')++p;} } if(s[p]=='['){get();while(1){ws();if(s[p]==']'){++p;return;}skip();ws();if(s[p]==',')++p;} } if(s[p]=='t'||s[p]=='f') boolean(); else num(); }
};

static vector<Item> it;
static int W,H; static bool canRot;
static double maxDen,maxVal;

static inline bool contains(const Rect&a,const Rect&b){
    return b.x>=a.x && b.y>=a.y && b.x+b.w<=a.x+a.w && b.y+b.h<=a.y+a.h;
}
static inline bool intersects(const Rect&a,int x,int y,int w,int h){
    return a.x<x+w && x<a.x+a.w && a.y<y+h && y<a.y+a.h;
}

// Split every maximal empty rectangle hit by a placement.  Unlike a guillotine tree,
// the retained empty rectangles may overlap, which is precisely what retains useful L-shaped space.
static void splitFree(vector<Rect>& free, int x,int y,int w,int h){
    vector<Rect> add;
    for(int i=0;i<(int)free.size();){
        Rect f=free[i];
        if(!intersects(f,x,y,w,h)){ ++i; continue; }
        free[i]=free.back(); free.pop_back();
        if(x>f.x) add.push_back({f.x,f.y,x-f.x,f.h});
        if(x+w<f.x+f.w) add.push_back({x+w,f.y,f.x+f.w-x-w,f.h});
        if(y>f.y) add.push_back({f.x,f.y,f.w,y-f.y});
        if(y+h<f.y+f.h) add.push_back({f.x,y+h,f.w,f.y+f.h-y-h});
    }
    for(Rect q:add) if(q.w>0 && q.h>0) free.push_back(q);
    // Remove non-maximal rectangles.  This is also important for keeping the candidate scan bounded.
    vector<char> dead(free.size());
    for(int i=0;i<(int)free.size();++i) if(!dead[i])
        for(int j=0;j<(int)free.size();++j) if(i!=j && !dead[i] && contains(free[j],free[i]) &&
            (!contains(free[i],free[j]) || j<i)) dead[i]=1;
    int z=0; for(int i=0;i<(int)free.size();++i) if(!dead[i]) free[z++]=free[i];
    free.resize(z);
}

static Ans build(int mode){
    vector<Rect> free(1,{0,0,W,H});
    vector<int> used(it.size()); vector<Pl> out;
    long long value=0;
    // These variants change only the greedy objective, retaining the same maximal-rectangle mechanism.
    double denExp[] = {4.0,2.3,1.4,0.8,0.35,6.0};
    double valExp[] = {0.0,0.12,0.30,0.55,0.85,0.0};
    double fitWeight[] = {0.10,0.22,0.42,0.72,1.05,0.04};
    int family=mode%6;
    double de=denExp[family], ve=valExp[family], fw=fitWeight[family];
    int corner=(mode/6)%4;
    for(int steps=0; !free.empty() && steps<15000; ++steps){
        int bi=-1,bf=-1,br=0,bx=0,by=0; double best=-1;
        for(int fi=0;fi<(int)free.size();++fi){
            const Rect& f=free[fi];
            for(int k=0;k<(int)it.size();++k) if(used[k]<it[k].lim){
                for(int r=0;r<=(canRot && it[k].w!=it[k].h);++r){
                    int w=r?it[k].h:it[k].w, h=r?it[k].w:it[k].h;
                    if(w>f.w || h>f.h) continue;
                    double dn=(double)it[k].v/(double)(w*h)/maxDen;
                    double vn=(double)it[k].v/maxVal;
                    double itemScore=pow(max(1e-12,dn),de)*pow(max(1e-12,0.25+0.75*vn),ve);
                    int rw=f.w-w, rh=f.h-h;
                    // Reward consuming one dimension exactly; it avoids making thin unusable scraps.
                    double fit=1.0-(double)min(rw,rh)/(double)max(1,min(f.w,f.h));
                    double area=(double)(w*h)/(double)(f.w*f.h);
                    double score=itemScore*(1.0+fw*fit+0.10*area);
                    // stable mode-specific tiebreaks are enough to produce several different layouts.
                    score += 1e-12*(double)((k*31+r*7+fi*3+mode*19)%997);
                    if(score>best){
                        int px=f.x,py=f.y;
                        if(corner==1){px=f.x+rw; py=f.y;}
                        else if(corner==2){px=f.x; py=f.y+rh;}
                        else if(corner==3){px=f.x+rw; py=f.y+rh;}
                        best=score;bi=k;bf=fi;br=r;bx=px;by=py;
                    }
                }
            }
        }
        if(bi<0) break;
        int w=br?it[bi].h:it[bi].w, h=br?it[bi].w:it[bi].h;
        out.push_back({bi,bx,by,br}); ++used[bi]; value+=it[bi].v;
        splitFree(free,bx,by,w,h);
    }
    return {value,move(out)};
}

static string esc(const string&s){ string r; for(unsigned char c:s){ if(c=='"'||c=='\\')r+='\\',r+=c; else if(c=='\n')r+="\\n"; else if(c=='\r')r+="\\r"; else if(c=='\t')r+="\\t"; else r+=c; } return r; }

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    Json j(string((istreambuf_iterator<char>(cin)),{})); j.get();
    while(1){ j.ws(); if(j.p>=j.s.size()||j.s[j.p]=='}'){++j.p;break;} string key=j.str();j.get();
        if(key=="bin"){ j.get(); while(1){j.ws();if(j.s[j.p]=='}'){++j.p;break;}string k=j.str();j.get();if(k=="W")W=j.num();else if(k=="H")H=j.num();else if(k=="allow_rotate")canRot=j.boolean();else j.skip();j.ws();if(j.s[j.p]==',')++j.p;} }
        else if(key=="items"){ j.get();while(1){j.ws();if(j.s[j.p]==']'){++j.p;break;}j.get();Item q{};while(1){j.ws();if(j.s[j.p]=='}'){++j.p;break;}string k=j.str();j.get();if(k=="type")q.id=j.str();else if(k=="w")q.w=j.num();else if(k=="h")q.h=j.num();else if(k=="v")q.v=j.num();else if(k=="limit")q.lim=j.num();else j.skip();j.ws();if(j.s[j.p]==',')++j.p;}it.push_back(q);j.ws();if(j.s[j.p]==',')++j.p;} }
        else j.skip(); j.ws();if(j.p<j.s.size()&&j.s[j.p]==',')++j.p;
    }
    maxDen=1;maxVal=1; for(auto&q:it){maxDen=max(maxDen,(double)q.v/(q.w*q.h));maxVal=max(maxVal,(double)q.v);}
    Ans ans;
    // Twenty-four layouts keep the construction comfortably within the per-case time budget.
    for(int m=0;m<24;++m){ Ans q=build(m); if(q.value>ans.value) ans=move(q); }
    cout<<"{\"placements\":[";
    for(size_t i=0;i<ans.p.size();++i){if(i)cout<<',';auto&q=ans.p[i];cout<<"{\"type\":\""<<esc(it[q.t].id)<<"\",\"x\":"<<q.x<<",\"y\":"<<q.y<<",\"rot\":"<<q.r<<'}';}
    cout<<"]}\n";
}
