#include<bits/stdc++.h>
using namespace std;
#include "../C-Quadratic-Sieve_C++/utils.h"
#include "../C-Quadratic-Sieve_C++/statistics.h"
#include "Int1.h"

typedef unsigned long long ull;
const int N=1005,N1=155;
int f[N][N][N1],n1,n2,n3,p,fp;  //base=2^n1*3^n2*p^n3.
inline void upd(int &x,int y){if (y<x)x=y;}
Int rand_Int(int n){  //rand Int with n digits
	Int x=0;
	for (int i=0;i<n;++i)x=x*10+rand64()%10;
	return x;
}
int calc(Int x){
	for (int i=0;i<=n1;++i)
		for (int j=0;j<=n2;++j)memset(f[i][j],0x3f,sizeof(int)*(n3+1));
	f[0][0][0]=0;
	//f[i][j][k]: already divided by i 2's, j 3's, and k p's.
	for (int i=0;i<=n1;++i){
		Int y=x;
		for (int j=0;j<=n2;++j){
			Int z=y;
			for (int k=0;k<=n3;++k){
				upd(f[i+1][j][k],f[i][j][k]+int(z.a[0]%2));
				upd(f[i][j+1][k],f[i][j][k]+int(z%3));
				upd(f[i][j][k+1],f[i][j][k]+int(z%p));
				z/=p;
			}
			y/=3;
		}
		x/=2;
	}
	return f[n1][n2][n3];
}
double run(int T=100){
	int t1=clock();
	double s=0;
	for (int i1=1;i1<=T;++i1){
		// set x to be moderately larger than base, so it will only introduce
		// small error, even if we don't restrict the rand value in [0,base).
		Int x=rand_Int(n1*log10(2)+n2*log10(3)+n3*log10(p)+15);
		s+=calc(x);
	}
	s/=T; s+=2*n1+3*n2+fp*n3;
	s/=(n1*log(2)+n2*log(3)+n3*log(p))/log(3);
	printf("n1=%d n2=%d n3=%d p=%d T=%d ave=%.6lf\n",n1,n2,n3,p,T,s);
	printf("time=%d\n",clock()-t1);
	return s;
}
int main(){
	srand(time(0));
	p=5; fp=5;
	//p=7; fp=6;
	//p=163; fp=15;
	//p=487; fp=18;
	//n1=9; n2=8; n3=1;
	n1=100; n2=100; n3=15;
	//n1=1000; n2=1000; n3=150;
	vector<double> a;
	for (int i1=1;;++i1){
		int T=10;
		double ans=run(T);
		a.push_back(ans);
		printf("--------i1=%d #samples=%d %.6lf--------\n",i1,i1*T,mean(a));
	}
	return 0;
}

