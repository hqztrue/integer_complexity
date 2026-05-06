//Complexity: O(n^2) per sample.
#include<bits/stdc++.h>
using namespace std;
#include "../integer-complexity-128/utils.h"
#include "../integer-complexity-128/statistics.h"
#include "Int.h"
#include "consts.h"

const int N=200005,inf=1e9;
int _f[2][N],n1,n2;  //base=2^n1*3^n2.
inline void upd(int &x,int y){if (y<x)x=y;}
int calc(Int x,int n1,int n2){
	int *f=_f[0],*f1=_f[1]; Int v;
	memset(f,0x3f,sizeof(int)*(n1+1));
	f[0]=0;
	for (int i=0;i<=n2;++i){
		memset(f1,0x3f,sizeof(int)*(n1+1));
		v=x; v.init_div2();
		//f[j]: already divided by i 3's, and j 2's.
		for (int j=0;j<=n1;++j){
			upd(f[j+1],f[j]+v.mod2);
			upd(f1[j],f[j]+v.mod3);
			v.div2();
		}
		x/=3; swap(f,f1);
	}
	return f1[n1]+2*n1+3*n2;
}
int calc_all(Int x){
	int ans=inf;
	double logx=x.log2();
	int n1=(int)ceil(logx),n2=(int)ceil(logx/log2(3));
	int *f=_f[0],*f1=_f[1]; Int v;
	memset(f,0x3f,sizeof(int)*(n1+1));
	f[0]=0;
	for (int i=0;i<=n2;++i){
		n1=max((int)ceil(logx-i*log2(3)),0);
		memset(f1,0x3f,sizeof(int)*(n1+1));
		v=x; v.init_div2();
		//f[j]: already divided by i 3's, and j 2's.
		for (int j=0;j<=n1;++j){
			upd(f[j+1],f[j]+v.mod2);
			upd(f1[j],f[j]+v.mod3);
			v.div2();
		}
		//printf("calc %.5lf %d %d %d\n",logx,n1,i,f[n1]+2*n1+3*i);
		ans=min(ans,f[n1]+2*n1+3*i);
		x/=3; swap(f,f1);
	}
	return ans;
}
double run_sample(int T=1){
	int t1=clock();
	double s=0; Int x;
	for (int i1=1;i1<=T;++i1){
		x.init_rand();
		s+=calc(x,n1,n2);
	}
	s/=T; s/=(n1*log(2)+n2*log(3))/log(3);
	printf("n1=%d n2=%d T=%d ave=%.6lf\n",n1,n2,T,s);
	printf("time=%d\n",clock()-t1);
	return s;
}
void run(){
	//n1=9; n2=8;
	//n1=11; n2=9;
	//n1=n2=10000;
	n1=11500; n2=10000;
	//n1=5600; n2=5000;
	//n1=115000; n2=100000;
	Int::N=(n1*log(2)+n2*log(3))/log(2)/32+5;
	vector<double> a;
	for (int i1=1;;++i1){
		int T=1;
		double ans=run_sample(T);
		a.push_back(ans);
		double ave=mean(a),mu=stddev(a);
		printf("--------i1=%d #samples=%d mean=%.6lf stddev=%.6lf--------\n",i1,i1*T,ave,mu);
	}
}
int heuristic_solve(Int x){
	int n=ceil(x.bit_length()/log2(6));
	return calc(x,n,n);
}
void test_heuristic(){
	Int n; Int::N=500;
	for (auto x:{x1_str,x2_str,x3_str}){
		n.read(x);
		int ans=calc_all(n);
		//int ans=heuristic_solve(n);
		printf("%d %.6lf\n",ans,ans/(n.bit_length()/log2(3)));
	}
}

#include "upper_bound_test_heuristics.cpp"

int main(){
	srand(time(0));
	test_heuristic();
	//test_heuristic_conjectures();
	//run();
	return 0;
}

