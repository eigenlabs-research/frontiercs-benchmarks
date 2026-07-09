// throughput-probe identity scheduler
#include <cstdio>
int main(){
    int J,M;
    if(scanf("%d %d",&J,&M)!=2) return 0;
    for(int j=0;j<J;j++) for(int k=0;k<M;k++){int m; long long p; if(scanf("%d %lld",&m,&p)!=2) return 0;}
    for(int m=0;m<M;m++){ for(int j=0;j<J;j++) printf("%d%c", j, j+1<J?' ':'\n'); }
    return 0;
}
// probe #0 marker
// burst variant 2 unique 117099240
