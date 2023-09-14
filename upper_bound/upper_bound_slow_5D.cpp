//Complexity: O(n^6/w) per sample.
#include<bits/stdc++.h>
using namespace std;
#include "../integer-complexity-128/utils.h"
#include "../integer-complexity-128/statistics.h"
#include "Int1.h"

typedef unsigned long long ull;
//const int N=125,N1=35,N2=35,N3=25;
const int N1=255,N2=215,N3=30,N4=20,N5=6;
int f[N1][N2][N3][N4][N5],n1,n2,n3,n4,n5,p,fp,q,fq,r,fr;  //base=2^n1*3^n2*p^n3*q^n4*r^n5.
inline void upd(int &x,int y){if (y<x)x=y;}
Int rand_Int(int n){  //rand Int with n digits
	Int x=0;
	for (int i=0;i<n;++i)x=x*10+rand64()%10;
	return x;
}
int calc(Int x){
	for (int i=0;i<=n1;++i)
		for (int j=0;j<=n2;++j)
			for (int k=0;k<=n3;++k)
				for (int l=0;l<=n4;++l)memset(f[i][j][k][l],0x3f,sizeof(int)*(n5+1));
	f[0][0][0][0][0]=0;
	//f[i][j][k][l][l1]: already divided by i 2's, j 3's, k p's, l q's, and l1 r's.
	for (int i=0;i<=n1;++i){
		Int y=x;
		for (int j=0;j<=n2;++j){
			Int z=y;
			for (int k=0;k<=n3;++k){
				Int t=z;
				for (int l=0;l<=n4;++l){
					Int t1=t;
					for (int l1=0;l1<=n5;++l1){
						upd(f[i+1][j][k][l][l1],f[i][j][k][l][l1]+int(t1.a[0]%2));
						upd(f[i][j+1][k][l][l1],f[i][j][k][l][l1]+int(t1%3));
						upd(f[i][j][k+1][l][l1],f[i][j][k][l][l1]+int(t1%p));
						upd(f[i][j][k][l+1][l1],f[i][j][k][l][l1]+int(t1%q));
						upd(f[i][j][k][l][l1+1],f[i][j][k][l][l1]+int(t1%r));
						t1/=r;
					}
					t/=q;
				}
				z/=p;
			}
			y/=3;
		}
		x/=2;
	}
	return f[n1][n2][n3][n4][n5];
}
double run_sample(int T=100){
	int t1=clock();
	double s=0;
	for (int i1=1;i1<=T;++i1){
		// set x to be moderately larger than base, so it will only introduce
		// small error, even if we don't restrict the rand value in [0,base).
		Int x=rand_Int(n1*log10(2)+n2*log10(3)+n3*log10(p)+n4*log10(q)+n5*log10(r)+15);
		s+=calc(x);
	}
	s/=T; s+=2*n1+3*n2+fp*n3+fq*n4+fr*n5;
	s/=(n1*log(2)+n2*log(3)+n3*log(p)+n4*log(q)+n5*log(r))/log(3);
	printf("n1=%d n2=%d n3=%d n4=%d n5=%d p=%d q=%d r=%d T=%d ave=%.6lf\n",n1,n2,n3,n4,n5,p,q,r,T,s);
	printf("time=%d\n",clock()-t1);
	return s;
}
void run(){
	p=5; fp=5;
	q=7; fq=6;
	r=11; fr=8;
	//n1=100; n2=100; n3=14; n4=8; n5=2;
	n1=205; n2=200; n3=28; n4=16; n5=5;
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

