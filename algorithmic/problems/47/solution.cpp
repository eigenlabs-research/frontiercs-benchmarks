#include <bits/stdc++.h>
using namespace std;

struct Item { string id; int w,h; long long v,lim; };
struct Pl { int t,x,y,r; };
struct Tok {
    string s; size_t p=0;
    Tok(){ s.assign(istreambuf_iterator<char>(cin),{}); }
    void ws(){ while(p<s.size() && isspace((unsigned char)s[p])) ++p; }
    char ch(){ ws(); return p<s.size()?s[p]:0; }
    void take(char c){ ws(); if(p<s.size()&&s[p]==c) ++p; }
    string str(){ ws(); take('"'); string r; while(p<s.size()&&s[p]!='"'){ if(s[p]=='\\'&&p+1<s.size()) ++p; r+=s[p++]; } take('"'); return r; }
    long long num(){ ws(); int sg=1; if(s[p]=='-'){sg=-1;++p;} long long z=0; while(p<s.size()&&isdigit((unsigned char)s[p])) z=z*10+s[p++]-'0'; return sg*z; }
    bool boolean(){ ws(); if(s.compare(p,4,"true")==0){p+=4;return true;} p+=5;return false; }
};

struct Sky {
    int W,H; map<int,int> a;
    Sky(int w=0,int h=0):W(w),H(h){a[0]=0;a[W]=0;}
    int height(int x,int e) const {
        auto it=a.upper_bound(x); --it; int z=0;
        while(it!=a.end() && it->first<e){ z=max(z,it->second); ++it; }
        return z;
    }
    bool put(int x,int w,int h, int &y){
        if(x+w>W) return false;
        y=height(x,x+w); if(y+h>H) return false;
        int after=height(x+w,min(W,x+w+1));
        a[x]=y+h; a[x+w]=after;
        auto it=a.upper_bound(x); while(it!=a.end() && it->first<x+w) it=a.erase(it);
        // Coalesce equal-height adjacent segments; this keeps candidate positions bounded.
        for(auto q=a.begin(); q!=a.end(); ){
            auto n=next(q); if(n==a.end()) break;
            if(q->second==n->second && n->first!=W) a.erase(n); else q=n;
        }
        return true;
    }
    bool best(int w,int h,int style,int &bx,int &by){
        bool ok=false; long long best1=0,best2=0,best3=0;
        for(auto it=a.begin(); it!=a.end(); ++it){
            int x=it->first; if(x+w>W) continue;
            int y=height(x,x+w); if(y+h>H) continue;
            long long q1,q2,q3;
            if(style==0) q1=y,q2=y+h,q3=x;
            else if(style==1) q1=y+h,q2=y,q3=x;
            else q1=y,q2=(long long)(W-x-w)+(H-y-h),q3=x;
            if(!ok || tie(q1,q2,q3)<tie(best1,best2,best3)) ok=true,bx=x,by=y,best1=q1,best2=q2,best3=q3;
        }
        return ok;
    }
};

int main(){
    Tok z; vector<Item> it; int W=0,H=0; bool rot=false;
    z.take('{');
    while(z.ch()!='}' && z.ch()){
        string key=z.str(); z.take(':');
        if(key=="bin"){
            z.take('{'); while(z.ch()!='}') { string k=z.str(); z.take(':'); if(k=="W") W=z.num(); else if(k=="H") H=z.num(); else rot=z.boolean(); if(z.ch()==',')z.take(','); } z.take('}');
        } else {
            z.take('['); while(z.ch()!=']'){
                z.take('{'); Item q; while(z.ch()!='}') { string k=z.str(); z.take(':'); if(k=="type")q.id=z.str(); else if(k=="w")q.w=z.num(); else if(k=="h")q.h=z.num(); else if(k=="v")q.v=z.num(); else q.lim=z.num(); if(z.ch()==',')z.take(','); } z.take('}'); it.push_back(q); if(z.ch()==',')z.take(',');
            } z.take(']');
        }
        if(z.ch()==',') z.take(',');
    }
    vector<Pl> answer; long long answerVal=-1;
    long long area=1LL*W*H;
    // A bounded portfolio is intentional: different global priorities make very different skylines.
    for(int trial=0;trial<54;trial++){
        struct U{int t; double key; unsigned tie;}; vector<U> u;
        uint32_t seed=0x9e3779b9u*(trial+17);
        auto rnd=[&](){seed=seed*1664525u+1013904223u; return seed;};
        for(int k=0;k<(int)it.size();k++){
            long long geom=max(1LL,area/(1LL*it[k].w*it[k].h));
            int n=(int)min({it[k].lim,geom,(long long)420});
            double den=(double)it[k].v/(it[k].w*it[k].h);
            double val=(double)it[k].v;
            double shape=(double)min(it[k].w,it[k].h)/max(it[k].w,it[k].h);
            for(int c=0;c<n;c++){
                double noise=((rnd()&65535)/65535.0-.5);
                double key;
                switch(trial%6){
                    case 0:key=log1p(den);break;
                    case 1:key=log1p(val);break;
                    case 2:key=log1p(den)+noise*.45;break;
                    case 3:key=.65*log1p(den)+.35*log1p(val)+noise*.18;break;
                    case 4:key=log1p(den)+.30*shape+noise*.25;break;
                    default:key=log1p(val/sqrt((double)it[k].w*it[k].h))+noise*.35;
                }
                u.push_back({k,key,rnd()});
            }
        }
        sort(u.begin(),u.end(),[](const U&A,const U&B){return A.key!=B.key?A.key>B.key:A.tie<B.tie;});
        Sky sky(W,H); vector<long long> used(it.size()); vector<Pl> cur; long long score=0;
        for(auto e:u){
            int k=e.t; if(used[k]>=it[k].lim) continue;
            int bx,by,br=0; bool found=false;
            int x,y; int style=(trial/6)%3;
            if(sky.best(it[k].w,it[k].h,style,x,y)){ bx=x;by=y;found=true; }
            if(rot && it[k].w!=it[k].h){ int rx,ry; if(sky.best(it[k].h,it[k].w,style,rx,ry)){
                // Compare orientations using the same skyline criterion, favoring lower top.
                if(!found || make_pair(ry+it[k].w,ry)<make_pair(by+it[k].h,by)){bx=rx;by=ry;br=1;found=true;}
            }}
            if(found){ int ww=br?it[k].h:it[k].w, hh=br?it[k].w:it[k].h, yy; if(sky.put(bx,ww,hh,yy)){cur.push_back({k,bx,yy,br});used[k]++;score+=it[k].v;} }
        }
        if(score>answerVal){answerVal=score;answer.swap(cur);}
    }
    cout << "{\"placements\":[";
    for(size_t i=0;i<answer.size();i++){ if(i) cout<<','; auto&p=answer[i]; cout<<"{\"type\":\""<<it[p.t].id<<"\",\"x\":"<<p.x<<",\"y\":"<<p.y<<",\"rot\":"<<p.r<<"}"; }
    cout << "]}";
}
