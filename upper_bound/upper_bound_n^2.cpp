#include<bits/stdc++.h>
using namespace std;
#include "../C-Quadratic-Sieve_C++/utils.h"
#include "../C-Quadratic-Sieve_C++/statistics.h"
#include "Int1.h"

typedef unsigned long long ull;
const int N=5005;
int f[N][N],n1,n2;
//Int pow2[N],pow3[N];
Int rand_Int(int n){  //rand Int with n digits
	Int x=0;
	for (int i=0;i<n;++i)x=x*10+rand64()%10;
	return x;
}
/*void pre(){
	pow2[0]=pow3[0]=1;
	for (int i=1;i<=n1;++i)pow2[i]=pow2[i-1]*2;
	for (int i=1;i<=n2;++i)pow3[i]=pow3[i-1]*3;
}*/
int calc(Int x){
	for (int i=0;i<=n1;++i)memset(f[i],0x3f,sizeof(int)*(n2+1));
	f[0][0]=0;
	for (int j=0;j<=n2;++j){
		Int v=x;
		for (int i=0;i<=n1;++i){
			f[i+1][j]=min(f[i+1][j],f[i][j]+int(v.a[0]%2));
			f[i][j+1]=min(f[i][j+1],f[i][j]+int(v%3));
			v/=2;
		}
		x/=3;
	}
	return f[n1][n2];
}
double run(int T=100){
	int t1=clock();
	double s=0; //pre();
	for (int i1=1;i1<=T;++i1)s+=calc(rand_Int(n1));
	s/=T; s+=2*n1+3*n2;
	s/=(n1*log(2)+n2*log(3))/log(3);
	printf("n1=%d n2=%d T=%d ave=%.6lf\n",n1,n2,T,s);
	printf("time=%d\n",clock()-t1);
	return s;
}
int main(){
	srand(time(0));
	n1=9; n2=8;
	//n1=11; n2=9;
	//n1=n2=1000;
	//n1=1000; n2=800;
	vector<double> a;
	for (int i1=1;;++i1){
		int T=10000;
		double ans=run(T);
		a.push_back(ans);
		printf("--------i1=%d #samples=%d %.6lf--------\n",i1,i1*T,mean(a));
	}
	return 0;
}

