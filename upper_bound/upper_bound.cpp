#include<bits/stdc++.h>
using namespace std;
#include "Int.h"

typedef unsigned long long ull;
const int N=105;
int f[N][N],n1,n2;
ull pow2[N],pow3[N];
inline int rand32(){  //31 bits
#ifdef _WIN32
	//assert(RAND_MAX==32767);
	return ((rand()&1)<<30)+(rand()<<15)+rand();
#else
	//assert(RAND_MAX==2147483647);
	return rand();
#endif
}
inline ull rand64(){  //64 bits
#ifdef _WIN32
	//assert(RAND_MAX==32767);
	return ((ull)rand()<<60)+((ull)rand()<<45)+((ull)rand()<<30)+((ull)rand()<<15)+rand();
#else
	//assert(RAND_MAX==2147483647);
	return ((ull)rand()<<62)+((ull)rand()<<31)+rand();
#endif
}
void pre(){
	pow2[0]=pow3[0]=1;
	for (int i=1;i<=n1;++i)pow2[i]=pow2[i-1]*2;
	for (int i=1;i<=n2;++i)pow3[i]=pow3[i-1]*3;
}
int calc(ull x){
	for (int i=0;i<=n1;++i)memset(f[i],0x3f,sizeof(int)*(n2+1));
	f[0][0]=0;
	for (int i=0;i<=n1;++i)
		for (int j=0;j<=n2;++j){
			ull v=x/pow2[i]/pow3[j];
			f[i+1][j]=min(f[i+1][j],f[i][j]+int(v%2));
			f[i][j+1]=min(f[i][j+1],f[i][j]+int(v%3));
		}
	return f[n1][n2];
}
int main(){
	srand(time(0));
	//n1=9; n2=8;
	//n1=11; n2=9;
	n1=n2=23;
	int T=10000; double s=0;
	pre();
	for (int i1=1;i1<=T;++i1)s+=calc(rand64()%(pow2[n1]*pow3[n2]));
	s/=T; s+=2*n1+3*n2;
	s/=(n1*log(2)+n2*log(3))/log(3);
	printf("%.6lf\n",s);
	return 0;
}

