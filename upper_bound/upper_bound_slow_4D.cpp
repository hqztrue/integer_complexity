//Complexity: O(n^5/w) per sample.
#include<bits/stdc++.h>
using namespace std;
#include "../integer-complexity-128/utils.h"
#include "../integer-complexity-128/statistics.h"
#include "Int1.h"

typedef unsigned long long ull;
const int N=455,N1=65,N2=35;
int f[N][N][N1][N2],n1,n2,n3,n4,p,fp,q,fq;  //base=2^n1*3^n2*p^n3*q^n4.
inline void upd(int &x,int y){if (y<x)x=y;}
Int rand_Int(int n){  //rand Int with n digits
	Int x=0;
	for (int i=0;i<n;++i)x=x*10+rand64()%10;
	return x;
}
int calc(Int x){
	for (int i=0;i<=n1;++i)
		for (int j=0;j<=n2;++j)
			for (int k=0;k<=n3;++k)memset(f[i][j][k],0x3f,sizeof(int)*(n4+1));
	f[0][0][0][0]=0;
	//f[i][j][k][l]: already divided by i 2's, j 3's, k p's, and l q's.
	for (int i=0;i<=n1;++i){
		Int y=x;
		for (int j=0;j<=n2;++j){
			Int z=y;
			for (int k=0;k<=n3;++k){
				Int t=z;
				for (int l=0;l<=n4;++l){
					upd(f[i+1][j][k][l],f[i][j][k][l]+int(t.a[0]%2));
					upd(f[i][j+1][k][l],f[i][j][k][l]+int(t%3));
					upd(f[i][j][k+1][l],f[i][j][k][l]+int(t%p));
					upd(f[i][j][k][l+1],f[i][j][k][l]+int(t%q));
					t/=q;
				}
				z/=p;
			}
			y/=3;
		}
		x/=2;
	}
	return f[n1][n2][n3][n4];
}
double run_sample(int T=100){
	int t1=clock();
	double s=0;
	for (int i1=1;i1<=T;++i1){
		// set x to be moderately larger than base, so it will only introduce
		// small error, even if we don't restrict the rand value in [0,base).
		Int x=rand_Int(n1*log10(2)+n2*log10(3)+n3*log10(p)+n4*log10(q)+15);
		s+=calc(x);
	}
	s/=T; s+=2*n1+3*n2+fp*n3+fq*n4;
	s/=(n1*log(2)+n2*log(3)+n3*log(p)+n4*log(q))/log(3);
	printf("n1=%d n2=%d n3=%d n4=%d p=%d q=%d T=%d ave=%.6lf\n",n1,n2,n3,n4,p,q,T,s);
	printf("time=%d\n",clock()-t1);
	return s;
}
void run(){
	p=5; fp=5;
	q=7; fq=6;
	n1=100; n2=100; n3=15; n4=7;
	n1=400; n2=400; n3=60; n4=20;
	n1=400; n2=400; n3=40; n4=20;
	n1=400; n2=440; n3=40; n4=20;
	vector<double> a;
	for (int i1=1;;++i1){
		int T=1;
		double ans=run_sample(T);
		a.push_back(ans);
		double ave=mean(a),mu=stddev(a);
		printf("--------i1=%d #samples=%d mean=%.6lf stddev=%.6lf--------\n",i1,i1*T,ave,mu);
	}
}
int main(){
	srand(time(0));
	run();
	return 0;
}

