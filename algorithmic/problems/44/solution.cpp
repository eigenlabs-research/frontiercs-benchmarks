// Valid identity tour: 0,1,...,N-1,0. Produces N+1 lines of output — for large N
// this is a big stdout, exactly the payload that triggered the Temporal size failure.
#include <cstdio>
int main(){
    int N;
    if(scanf("%d",&N)!=1) return 0;
    // skip N coordinate lines
    for(long long i=0;i<(long long)N;i++){ long long x,y; if(scanf("%lld %lld",&x,&y)!=2) break; }
    printf("%d\n", N+1);
    for(int i=0;i<N;i++) printf("%d\n", i);
    printf("0\n");
    return 0;
}
