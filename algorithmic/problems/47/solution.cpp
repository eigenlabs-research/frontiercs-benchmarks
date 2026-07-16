#include <bits/stdc++.h>
using namespace std;
struct Item { string id; int w,h,lim; long long v; };
struct P { int t,x,y,r,w,h; };
struct Seg { int l,r,h; };

static size_t afterKey(const string& s,const string& k){
    size_t p=s.find("\""+k+"\""); if(p==string::npos) return p;
    p=s.find(':',p); return p==string::npos?p:p+1;
}
static long long num(const string& s,const string& k){
    size_t p=afterKey(s,k); if(p==string::npos) return 0;
    while(p<s.size() && isspace((unsigned char)s[p])) ++p;
    bool neg=false; if(p<s.size()&&s[p]=='-') neg=true,++p;
    long long z=0; while(p<s.size()&&isdigit((unsigned char)s[p])) z=z*10+s[p++]-'0'; return neg?-z:z;
}
static string str(const string& s,const string& k){
    size_t p=afterKey(s,k); if(p==string::npos) return "";
    p=s.find('"',p); if(p==string::npos) return ""; ++p; string r;
    while(p<s.size()&&s[p]!='"'){ if(s[p]=='\\'&&p+1<s.size()) ++p; r+=s[p++]; } return r;
}
static int topAt(const vector<Seg>& a,int x,int w){
    int z=0, e=x+w;
    for(const auto& q:a) { if(q.r>x && q.l<e) z=max(z,q.h); if(q.l>=e) break; }
    return z;
}
static void raise(vector<Seg>& a,int x,int w,int nh){
    int e=x+w; vector<Seg> b; b.reserve(a.size()+3);
    for(auto q:a){
        if(q.r<=x||q.l>=e) b.push_back(q);
        else {
            if(q.l<x) b.push_back({q.l,x,q.h});
            b.push_back({max(q.l,x),min(q.r,e),nh});
            if(q.r>e) b.push_back({e,q.r,q.h});
        }
    }
    vector<Seg> c; c.reserve(b.size());
    for(auto q:b){ if(q.l>=q.r) continue; if(!c.empty()&&c.back().r==q.l&&c.back().h==q.h)c.back().r=q.r; else c.push_back(q); }
    a.swap(c);
}

static pair<long long,vector<P>> skyline(const vector<Item>& it,int W,int H,bool rot,double exponent,int flavor){
    vector<Seg> sky(1,{0,W,0}); vector<int> used(it.size()); vector<P> ans; long long value=0;
    int total=0; for(auto &q:it) total+=q.lim;
    // A contour position is always an existing skyline edge; this is the bottom-left skyline packing rule.
    for(int step=0;step<total;step++){
        bool found=false; P best{}; double bestU=-1; int bestY=INT_MAX, bestX=INT_MAX;
        vector<int> edges; edges.reserve(sky.size()*2);
        for(auto q:sky) { edges.push_back(q.l); edges.push_back(q.r); }
        sort(edges.begin(),edges.end()); edges.erase(unique(edges.begin(),edges.end()),edges.end());
        for(int t=0;t<(int)it.size();t++) if(used[t]<it[t].lim) {
            for(int r=0;r<=(rot?1:0);r++){
                int w=r?it[t].h:it[t].w, h=r?it[t].w:it[t].h;
                if(w>W||h>H) continue;
                vector<int> xs; xs.reserve(edges.size()*2);
                for(int e:edges) { if(e<=W-w) xs.push_back(e); if(e>=w) xs.push_back(e-w); }
                sort(xs.begin(),xs.end()); xs.erase(unique(xs.begin(),xs.end()),xs.end());
                for(int x:xs){
                    int y=topAt(sky,x,w); if(y+h>H) continue;
                    double area=(double)w*h;
                    double u;
                    if(flavor==0) u=it[t].v/pow(area,exponent);
                    else if(flavor==1) u=it[t].v/pow((double)h,exponent); // rewards thin row fillers
                    else u=it[t].v/pow((double)w,exponent);              // rewards narrow column fillers
                    // Bottom-most placement is primary; the policy only chooses among equal contour levels.
                    if(!found || y<bestY || (y==bestY && (u>bestU+1e-9 || (fabs(u-bestU)<1e-9 && x<bestX)))){
                        found=true; best={t,x,y,r,w,h}; bestY=y; bestX=x; bestU=u;
                    }
                }
            }
        }
        if(!found) break;
        raise(sky,best.x,best.w,best.y+best.h); used[best.t]++; value+=it[best.t].v; ans.push_back(best);
    }
    return {value,ans};
}
int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    string s((istreambuf_iterator<char>(cin)),{}); if(s.empty()) return 0;
    int W=(int)num(s,"W"), H=(int)num(s,"H");
    size_t bp=afterKey(s,"allow_rotate"); bool allow=false;
    if(bp!=string::npos){ while(bp<s.size()&&isspace((unsigned char)s[bp]))bp++; allow=s.compare(bp,4,"true")==0; }
    vector<Item> a; size_t ip=afterKey(s,"items"); if(ip!=string::npos){
        size_t p=s.find('[',ip); if(p!=string::npos) { ++p;
            while(p<s.size()) { while(p<s.size()&&(isspace((unsigned char)s[p])||s[p]==','))++p; if(p>=s.size()||s[p]==']')break;
                if(s[p]!='{'){++p;continue;} size_t st=p; int d=0; do { if(s[p]=='{')d++; else if(s[p]=='}')d--; ++p; } while(p<s.size()&&d);
                string o=s.substr(st,p-st); Item q; q.id=str(o,"type");q.w=num(o,"w");q.h=num(o,"h");q.v=num(o,"v");q.lim=num(o,"limit"); if(q.w>0&&q.h>0&&q.lim>0)a.push_back(q);
            }
        }
    }
    pair<long long,vector<P>> out={0,{}};
    // Four deliberately distinct skyline priorities; keeping this bounded matters on high-limit cases.
    const pair<double,int> policies[]={{0.65,0},{1.05,0},{0.85,1},{0.85,2}};
    for(auto policy:policies) {
        auto z=skyline(a,W,H,allow,policy.first,policy.second); if(z.first>out.first) out=move(z);
    }
    cout<<"{\"placements\":[";
    for(size_t i=0;i<out.second.size();i++){ auto q=out.second[i]; if(i)cout<<',';
        cout<<"{\"type\":\""<<a[q.t].id<<"\",\"x\":"<<q.x<<",\"y\":"<<q.y<<",\"rot\":"<<q.r<<"}";
    }
    cout<<"]}";
}
