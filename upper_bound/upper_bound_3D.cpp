//Complexity: O(n^3) per sample.
#include<bits/stdc++.h>
using namespace std;
#include "../integer-complexity-128/utils.h"
#include "../integer-complexity-128/statistics.h"
#include "Int_3D.h"
#include "consts.h"

const int N=10005,N1=5005,p=5,fp=5;
int _f[2][N1][N];
int n1,n2,n3;  //base=2^n1*3^n2*p^n3.
inline void upd(int &x,int y){if (y<x)x=y;}
int calc(Int x,int n1,int n2,int n3){
	auto *f=_f[0],*f1=_f[1]; Int y,z;
	for (int i=0;i<=n3;++i)memset(f[i],0x3f,sizeof(int)*(n1+1));
	f[0][0]=0;
	for (int j=0;j<=n2;++j){
		for (int i=0;i<=n3;++i)memset(f1[i],0x3f,sizeof(int)*(n1+1));
		y=x;
		for (int k=0;k<=n3;++k){
			z=y; z.init_div2();
			//f[k][i]: already divided by i 2's, j 3's, and k p's.
			for (int i=0;i<=n1;++i){
				auto f0=f[k][i];
				upd(f[k][i+1],f0+z.mod2);
				upd(f1[k][i],f0+z.mod3);
				upd(f[k+1][i],f0+z.mod5);
				z.div2();
			}
			y/=p;
		}
		x/=3; swap(f,f1);
	}
	return f1[n3][n1]+2*n1+3*n2+fp*n3;
}
int calc_heuristic(Int x,int n1,int n2,int n3){
	auto *f=_f[0],*f1=_f[1]; Int y,z;
	int len=x.bit_length(),ans=1<<30,id1,id2,id3;
	for (int i=0;i<=n3;++i)memset(f[i],0x3f,sizeof(int)*(n1+1));
	f[0][0]=0;
	for (int j=0;j<=n2;++j){
		for (int i=0;i<=n3;++i)memset(f1[i],0x3f,sizeof(int)*(n1+1));
		y=x;
		for (int k=0;k<=n3;++k){
			z=y; z.init_div2();
			//f[k][i]: already divided by i 2's, j 3's, and k p's.
			int _n1=min(n1,max((int)ceil(len-j*log2(3)-k*log2(p)),0));
			for (int i=0;i<=_n1;++i){
				auto f0=f[k][i];
				upd(f[k][i+1],f0+z.mod2);
				upd(f1[k][i],f0+z.mod3);
				upd(f[k+1][i],f0+z.mod5);
				z.div2();
			}
			if (_n1<n1){
				int v=f[k][_n1]+2*_n1+3*j+fp*k;
				if (v<ans)ans=v,id1=_n1,id2=j,id3=k;
			}
			y/=p;
		}
		x/=3; swap(f,f1);
	}
	printf("id=%d %d %d\n",id1,id2,id3);
	printf("%d %.6lf\n",ans,ans/(len/log2(3)));
	return ans;
}
double run_sample(int T=1){
	int t1=clock();
	double s=0; Int x;
	for (int i1=1;i1<=T;++i1){
		x.init_rand();
		s+=calc(x,n1,n2,n3);
	}
	s/=T; s/=(n1*log(2)+n2*log(3)+n3*log(p))/log(3);
	printf("n1=%d n2=%d n3=%d p=%d T=%d ave=%.6lf\n",n1,n2,n3,p,T,s);
	printf("time=%d\n",clock()-t1);
	return s;
}
void run(){
	//n1=100; n2=100; n3=15;
	n1=1000; n2=1100; n3=150;
	//n1=2000; n2=2100; n3=300;
	Int::N=(n1*log(2)+n2*log(3)+n3*log(p))/log(2)/32+5;
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
	double r1=0.6,r2=0.35,r3=0.13;
	//double r1=0.8,r2=0.8,r3=0;
	int len=x.bit_length();
	int n1=ceil(len*r1),n2=ceil(len*r2),n3=ceil(len*r3);
	printf("n1=%d n2=%d n3=%d\n",n1,n2,n3);
	return calc_heuristic(x,n1,n2,n3);
}
void test_heuristic(){
	//const char _x1[]="1110010111001111101011100001111110010111010111000010000111110100001101001010110000000100011001001000001000011110000001101101000111000010110000100111100111100011110110111010110001110110010001111011100011001001101100111011111101101101100100111110100101010000000110101110101001010001001100111001100010111010110101010110011110111111110111";
	Int n; Int::N=500;
	for (auto x:{x1_str,x2_str,x3_str}){
		n.read(x);
		int ans=heuristic_solve(n);
	}
}
int main(){
	srand(time(0));
	test_heuristic();
	//run();
	return 0;
}

