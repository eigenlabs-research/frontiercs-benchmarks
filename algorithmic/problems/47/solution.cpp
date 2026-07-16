#include <bits/stdc++.h>
using namespace std;
struct Item{string id; int w,h; long long v,lim;};
struct R{int x,y,w,h;};
struct P{int t,x,y,rot,w,h;};

static long long numberAfter(const string&s,const string&key,size_t from=0){
    size_t p=s.find("\""+key+"\"",from); if(p==string::npos) return 0;
    p=s.find(':',p); p=s.find_first_of("-0123456789",p); size_t q=p;
    while(q<s.size() && (isdigit((unsigned char)s[q])||s[q]=='-')) ++q;
    return stoll(s.substr(p,q-p));
}
static string stringAfter(const string&s,const string&key,size_t from=0){
    size_t p=s.find("\""+key+"\"",from); p=s.find(':',p); p=s.find('"',p); size_t q=s.find('"',p+1); return s.substr(p+1,q-p-1);
}
static bool intersects(const R&a,const R&b){return a.x<b.x+b.w&&b.x<a.x+a.w&&a.y<b.y+b.h&&b.y<a.y+a.h;}
static bool contains(const R&a,const R&b){return a.x<=b.x&&a.y<=b.y&&a.x+a.w>=b.x+b.w&&a.y+a.h>=b.y+b.h;}

static void splitFree(vector<R>& f,const R& used){
    vector<R> z; z.reserve(f.size()*2);
    for(const R&r:f){
        if(!intersects(r,used)){z.push_back(r); continue;}
        int rx=r.x+r.w, ry=r.y+r.h, ux=used.x+used.w, uy=used.y+used.h;
        if(used.x>r.x) z.push_back({r.x,r.y,used.x-r.x,r.h});
        if(ux<rx) z.push_back({ux,r.y,rx-ux,r.h});
        if(used.y>r.y) z.push_back({r.x,r.y,r.w,used.y-r.y});
        if(uy<ry) z.push_back({r.x,uy,r.w,ry-uy});
    }
    // Fragmentation can be quadratic on adversarial instances.  Keeping the largest
    // valid empty rectangles is conservative (it can only omit future choices).
    if(z.size()>600){
        sort(z.begin(),z.end(),[](const R&a,const R&b){return 1LL*a.w*a.h>1LL*b.w*b.h;});
        z.resize(300);
    }
    vector<char> dead(z.size());
    for(size_t i=0;i<z.size();++i) if(z[i].w<=0||z[i].h<=0) dead[i]=1;
    for(size_t i=0;i<z.size();++i) if(!dead[i]) for(size_t j=0;j<z.size();++j)
        if(i!=j&&!dead[j]&&contains(z[j],z[i])) {dead[i]=1;break;}
    f.clear(); f.reserve(z.size()); for(size_t i=0;i<z.size();++i) if(!dead[i]) f.push_back(z[i]);
}

int main(){
 ios::sync_with_stdio(false);cin.tie(nullptr);
 string s((istreambuf_iterator<char>(cin)),{});
 size_t bp=s.find("\"bin\""); int W=(int)numberAfter(s,"W",bp),H=(int)numberAfter(s,"H",bp);
 size_t ap=s.find("\"allow_rotate\"",bp); bool rotate=s.find("true",ap)!=string::npos && s.find("true",ap)<s.find('}',ap);
 vector<Item>a; size_t pos=s.find('[',s.find("\"items\""));
 while(pos!=string::npos){
   pos=s.find('{',pos); if(pos==string::npos) break; size_t end=s.find('}',pos); if(end==string::npos)break;
   string o=s.substr(pos,end-pos+1); if(o.find("\"type\"")==string::npos)break;
   Item it; it.id=stringAfter(o,"type"); it.w=(int)numberAfter(o,"w");it.h=(int)numberAfter(o,"h");it.v=numberAfter(o,"v");it.lim=numberAfter(o,"limit"); a.push_back(it);pos=end+1;
 }
 vector<P> best; long long bestval=-1;
 // Different exponents deliberately cover profit-first through density-first policies.
 const double exps[]={1.35,1.0,.75,.5,0.0,1.15,.9,.6};
 for(int mode=0;mode<8;mode++){
   vector<long long> left(a.size());for(size_t i=0;i<a.size();i++)left[i]=a[i].lim;
   vector<R> free={{0,0,W,H}}; vector<P> out; long long val=0;
   while(!free.empty()){
     int bt=-1,bf=-1,br=0,bw=0,bh=0; double bs=-1e100, bfit=1e100;
     for(int t=0;t<(int)a.size();t++) if(left[t]){
       for(int rr=0;rr<(rotate?2:1);rr++){
         int iw=rr?a[t].h:a[t].w, ih=rr?a[t].w:a[t].h; if(iw>W||ih>H)continue;
         int choose=-1; double fit=1e100;
         for(int j=0;j<(int)free.size();j++) if(iw<=free[j].w&&ih<=free[j].h){
           int dw=free[j].w-iw,dh=free[j].h-ih;
           double q;
           if(mode%3==0) q=max(dw,dh)*1000000.0+min(dw,dh);
           else if(mode%3==1) q=(double)free[j].w*free[j].h-(double)iw*ih;
           else q=min(dw,dh)*1000000.0+max(dw,dh);
           if(q<fit){fit=q;choose=j;}
         }
         if(choose<0)continue;
         double area=(double)iw*ih;
         double pri=a[t].v/pow(max(1.0,area),exps[mode]);
         // stable, tiny geometry tie breaker only; priority is the mechanism.
         double sc=pri;
         if(sc>bs+1e-12 || (fabs(sc-bs)<1e-12 && fit<bfit)){
           bs=sc;bfit=fit;bt=t;bf=choose;br=rr;bw=iw;bh=ih;
         }
       }
     }
     if(bt<0)break;
     R u={free[bf].x,free[bf].y,bw,bh}; out.push_back({bt,u.x,u.y,br,bw,bh}); val+=a[bt].v;--left[bt]; splitFree(free,u);
   }
   if(val>bestval){bestval=val;best.swap(out);}
 }
 cout<<"{\"placements\":[";
 for(size_t i=0;i<best.size();i++){if(i)cout<<',';const P&p=best[i];cout<<"{\"type\":\""<<a[p.t].id<<"\",\"x\":"<<p.x<<",\"y\":"<<p.y<<",\"rot\":"<<p.rot<<'}';}
 cout<<"]}";
}
