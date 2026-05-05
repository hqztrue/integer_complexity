//Complexity: O(n^3/w) per sample.
#include<bits/stdc++.h>
using namespace std;
#include "../integer-complexity-128/utils.h"
#include "../integer-complexity-128/statistics.h"
#include "Int1.h"

typedef unsigned long long ull;
const int N=5005;
int f[N][N],pre[N][N],n1,n2;  //base=2^n1*3^n2.
int ub[N],lb[N],diffl,diffr;
inline void upd(int &x,int y,int &pre,int v){if (y<x)x=y,pre=v;}
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
			upd(f[i+1][j],f[i][j]+int(v.a[0]%2),pre[i+1][j],0);
			upd(f[i][j+1],f[i][j]+int(v%3),pre[i][j+1],1);
			v/=2;
		}
		x/=3;
	}
	return f[n1][n2];
}
void recover_path(){
	int v1=n1,v2=n2;
	while (v1||v2){
		lb[v1]=min(lb[v1],v2);
		ub[v1]=max(ub[v1],v2);
		if (pre[v1][v2]==0)--v1;
		else --v2;
	}
	for (int i=0;i<=n1;++i)
		diffl=max(diffl,i-lb[i]),diffr=max(diffr,ub[i]-i);
	for (int i=0;i<=n1;++i)printf("%d %d %d\n",i,lb[i],ub[i]);
	printf("diff: %d %d\n",diffl,diffr);
}
double run_sample(int T=100){
	int t1=clock();
	double s=0;
	for (int i1=1;i1<=T;++i1){
		// set x to be moderately larger than base, so it will only introduce
		// small error, even if we don't restrict the rand value in [0,base).
		Int x=rand_Int(max(n1,n2));
		s+=calc(x);
		recover_path();
	}
	s/=T; s+=2*n1+3*n2;
	s/=(n1*log(2)+n2*log(3))/log(3);
	printf("n1=%d n2=%d T=%d ave=%.6lf\n",n1,n2,T,s);
	printf("time=%d\n",clock()-t1);
	return s;
}
void run(){
	//n1=9; n2=8;
	//n1=11; n2=9;
	//n1=100; n2=100;
	//n1=1000; n2=1000;
	n1=5000; n2=5000;
	
	for (int i=0;i<N;++i)ub[i]=0,lb[i]=N;
	vector<double> a;
	for (int i1=1;;++i1){
		int T=10;
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

