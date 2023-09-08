//Complexity: O(n^3/w) per sample.
#include<bits/stdc++.h>
using namespace std;
#include "../C-Quadratic-Sieve_C++/utils.h"
#include "../C-Quadratic-Sieve_C++/statistics.h"
#include "Int1.h"

typedef unsigned long long ull;
const int N=5005;
int f[N][N],n1,n2;  //base=2^n1*3^n2.
Int rand_Int(int n){  //rand Int with n digits
	Int x=0;
	for (int i=0;i<n;++i)x=x*10+rand64()%10;
	return x;
}
int calc(Int x){
	for (int i=0;i<=n1;++i)memset(f[i],0x3f,sizeof(int)*(n2+1));
	f[0][0]=0;
	for (int j=0;j<=n2;++j){
		Int v=x;
		//f[i][j]: already divided by i 2's, and j 3's.
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
	double s=0;
	for (int i1=1;i1<=T;++i1){
		// set x to be moderately larger than base, so it will only introduce
		// small error, even if we don't restrict the rand value in [0,base).
		Int x=rand_Int(max(n1,n2));
		s+=calc(x);
	}
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
	n1=1000; n2=100;
	vector<double> a;
	for (int i1=1;;++i1){
		int T=1;
		double ans=run(T);
		a.push_back(ans);
		printf("--------i1=%d #samples=%d %.6lf--------\n",i1,i1*T,mean(a));
	}
	return 0;
}

