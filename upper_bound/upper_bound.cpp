#include<bits/stdc++.h>
using namespace std;
#include "../C-Quadratic-Sieve_C++/utils.h"
#include "../C-Quadratic-Sieve_C++/statistics.h"
#include "Int.h"

const int N=100005;
int _f[2][N],n1,n2;  //base=2^n1*3^n2.
inline void upd(int &x,int y){if (y<x)x=y;}
int calc(Int &x){
	int *f=_f[0],*f1=_f[1];
	memset(f,0x3f,sizeof(int)*(n1+1));
	f[0]=0;
	for (int i=0;i<=n2;++i){
		memset(f1,0x3f,sizeof(int)*(n1+1));
		Int v=x;
		v.init_div2();
		for (int j=0;j<=n1;++j){
			upd(f[j+1],f[j]+v.mod2);
			upd(f1[j],f[j]+v.mod3);
			v.div2();
		}
		x/=3; swap(f,f1);
	}
	return f[n1];
}
double run(int T=1){
	int t1=clock();
	double s=0; Int x;
	for (int i1=1;i1<=T;++i1){
		x.init_rand();
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
	//n1=9; n2=8;
	//n1=11; n2=9;
	n1=n2=10000;
	//n1=1000; n2=800;
	vector<double> a;
	for (int i1=1;;++i1){
		int T=1;
		double ans=run(T);
		a.push_back(ans);
		printf("--------i1=%d #samples=%d %.6lf--------\n",i1,i1*T,mean(a));
	}
	return 0;
}

