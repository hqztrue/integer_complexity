//Complexity: O(n^3) per sample.
#include<bits/stdc++.h>
using namespace std;
#include "../C-Quadratic-Sieve_C++/utils.h"
#include "../C-Quadratic-Sieve_C++/statistics.h"
#include "Int_3D.h"

const int N=10005,N1=2005,p=5,fp=5;
typedef int arr[N1][N];
arr _f[2];
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
				upd(f[k][i+1],f[k][i]+z.mod2);
				upd(f1[k][i],f[k][i]+z.mod3);
				upd(f[k+1][i],f[k][i]+z.mod5);
				z.div2();
			}
			y/=p;
		}
		x/=3; swap(f,f1);
	}
	return f1[n3][n1]+2*n1+3*n2+fp*n3;
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
/*int heuristic_solve(Int x){
	int n=ceil(x.bit_length()/log2(6));
	return calc(x,n,n);
}*/
void run(){
	n1=100; n2=100; n3=1;
	//n1=1000; n2=1000; n3=150;
	Int::N=(n1*log(2)+n2*log(3)+n3*log(p))/log(2)/32+5;
	vector<double> a;
	for (int i1=1;;++i1){
		int T=10;
		double ans=run_sample(T);
		a.push_back(ans);
		printf("--------i1=%d #samples=%d %.6lf--------\n",i1,i1*T,mean(a));
	}
}
/*void test_heuristic(){
	const char x1[]="1111001110001101101100011111100111011101001111011010110000000101011110001101001010010110100101010011100110111100101101101101001010101111101000011111001010101100101010111000110010110001111010001111101000111011001111001101100100010100111011010110011010101000100001010001101011001011110001000100100001110101011110111111010111111100010101110011010011011101110101010010001110101011110110001000000011010000010000010010001111010110011010011100100010111100010010101010000100100010101111101101000100001110101110100011000001100011111100011011100110001010111111011000010111001010010101010101010001000010101101111110111000010011100110101110101010000101100110101000000101011010000110010000011011100111101000011010010100100111111001100001011101101001111101001001010011001110111110000010011011011001111010000011100101010000101100001110110101100101010100000000100010100011011001100011101100010011100011110110011101011100010000111011100000111111100111011101110010100111010111111001000110100110001111000011000000100000010110001010000000101000101000001011100101001110110011010011100101011000001010000000010110010111010101010101010010001010101110101001111000011111011010000010111100001111101011010110101111111101111101010010010000100100111100101110001000001111101111010100111100110100111001110010000010101100111001011011111001010111001010110110101000100011000010111110110111001011011100110010110010010000111101001101100111100010111000010111000110001000010001111011001100100001100111101000000000000011010111011110011100001000100101100110011001110011001101110100000110111111000001101010000010100110001011111101111111011011100101000100111011010010010001100110000001101010110001100001011101010010101111100011111101000001111011010100010100101010010011101011011001010110110110011100001101101111010111011101001110110001101000011001001010001000101100001110000111010110010001101100111011110001111001110100010110110101001100110001110110101000011011110010101001101010010100101101111010110100010000010110000100101001001000011101110100011110000110000100101101110110110110100110000000010111100111110000110010110001010010011110010000110100111001101111001011000110010100100001001011010100111000111010000000001000011111010001010101101101000010011011100101110101011000111110011010110000101001011110100111010111110100011000111100111001110110111000100110011000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001110010111001111101011100001111110010111010111000010000111110100001101001010110000000100011001001000001000011110000001101101000111000010110000100111100111100011110110111010110001110110010001111011100011001001101100111011111101101101100100111110100101010000000110101110101001010001001100111001100010111010110101010110011110111111110111";
	const char x2[]="11100111101101100100111001001101111000101111110001000010010100010100101001110010100111110110111001001010101011111010101111101000111101111101001111000001001110001100000010110101010101111110100110101011100011011100111001001011101110110001101001011000111110000100001000011001011100010011010010011000001110110000101011100011001101101100101111011111001111011001011100001011111010011110000000100011001010111001001011000111011000001000000000001011100011100101010110100101110011110101110001001010110100110001111101110111110111001000001001000110001000010010111001110101010010010101110101001010001100011010010010000110000001010011101000011111100111110101001100010101010110011100010011001010100100001101111100101010000011111101010010011100111111101110001010101011001100000011011001101110010000000011110101011010001101100000111111011011000000011001001111110111000010100100000100011101101011110010111001000011110000111101101001001100111100000011101111100000110111000000101100001000100001111101011110101001110010110101111111001000101100010010100111101100010101010000011000001100000101111100110100011100011001110101001101011110000011000100011001111010110110100011111100110111101111110101010001100001010000110001100101110110100001011111011001011110111110001000000010111000011000111000001111000001011010111001110110011001010110000101111011110100001000100001110101110011001001000101100010111110001101101110010010001011111001100001001011000101111100010000010101011100100111111100100100001111010000101101100000001110010000011111000110011001100001101111011010010100010010010101010101111100010110101100000010000110110110111010000011100110001000101111001001111111100010111001101100101101010110000010110101000111000100000100010100100010010111101101001010001001100111011000111111000010101110100001000001011100010000110111111101110000100100000000011001110001011001010101111000000011001000100111111101010110111101111100011001110101100100110000011101101010101100100010010101011000001111111110010000111000110001100011001010101100100101100100101010011011010110100011001100111110111011001111110110001110011101010010000101011111111011111011111001111110010001001110001101111101011101001010101111010110001001010001101110001100100111011010101000101101011111111011111000000111000100011110010101100010001001111011000111000110111011010101000010110101011000010010110100101010110001011001011111010001000111100011011000000010110010001011111010101101011100100101111000000111001111011011111011001100001100111000010111101001000001000101101001000000001111000000101000010111010000110000111101001100001010100010100001110101100110100010000100110110010000010100111000101010111111100110111100001001100010111001001011000010101110100011111000110100100101011100110000101001100110110111101011010100000000111100000100101010111110101101111100111000001000100011110111001000001001011111100001111011011110001101010100111001011000010110000111000111111000000001110000101100111111111110110000100101011011111111100111010111001010111000111011110111100111101001010100001011011011000101010011000001101000011001100011001000001010111001011001110110010101101110000100001001100110011010111010110111101001001110001011011010111111110001111000000001011001110011101101101001101001111010111001010000100101001011000100101110110011101100000101000011101001011111111000001011001010111110100000110001000000000101110111000111111010000110110001010010010010001100110101110101001000000100011110000011100011100011001010100101110010110101101011101000100101000101100110110001110001101100110011011111011100010110100011011110010000001001101100001000001011011010111000110101101000010110111010010000001001010111101110101000111001001010100001001100011111101000010100110000001001111010001111101011101111001100001001010000101101001101001110000000010001100111111111000111001101000101110111010110111011111001111110000011101111000100001111110001000110011000011011100100101001001111011000101001110011001101010100010010001101100111000011100100000101010001000000010001001001100110110100001101110001001000001100110010001111110001110110100100011110000100001110110000010010100101000111101100100101011101111011110011001001000010001010100001010110110000000100110101101110100001101011000001110111111100011101010011001100010100101011001111011100100101101011011100010111100101101111011010111111101010111110110111100110100010000000110111100101111100111001001010101111011010100010101110110111100001001100100110101011110111111110101011000011010001000100111110111100101110110100110001011111001110001001000110011111011001000101011100110011011011100000101000010100101011001110110010000101001110100100111101010100000010001000000000111100011100111100001100000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000110011101110011100011111000111101100101001100100000110001110010001010110000111001100011011001101001110001101000011101111001010110100010001110010110110111011100111000100101000110110001001101001101001100110001011110001011001110000100010010010110000100101111101010111100011000110001000100111011010110011000011111011000011001101011101111";
	const char x3[]="110111000111001001001100110101000100010000010001101011110101001100110111101100100100101111110001101100100000010111100000100011000010110010010101001001110110010111000000000001110010100111110000000010010000110101111100110110010110111101101111000111000110010010011001000111101110000010110000010101101010101100100110000010000111111110011010000001001010111111100110010001001111011001110111111011110001100010001011011100111101011011001100110111101010110100101101101111000010111011000000100011010000111101001001100001010111000001011000010101001000110010110001000110101011010011100100110011101011011000001111001110111100110101010010100110100010101100000100110100001010011111100011001001010100000111001011100001001100100011101011111011111011001110010110100111011111100101111111010000000000001110110011010010011001000101001000110110010010110010011001110000010010010010101110100011011001010101001110111100001101000100001001011010000111001010001011111010011010101110110010010101000101000110000110100000110000011010100001000011000111100000011010010000000011111000101100101110010100010100000000010010110100000111000100001111100001111110010001110011111011010111001010110000010000101101111110101101110001101101101101101110111100010111111001101010000100101000001101111100110111011101001100001111111100001111110101001100100011000001011101001110111111100011101101100010101010011111111111100010111010111010010000000010001011100010000011110110101111010111101010010000011000001011011010000011001000011010000100000100100001100110100000110011111001111011010101011001110011101111000110110011100111010011001011101011100101000110011001010010011110010001011011110101101011111011101001001010011001100100001111000100000100011000010001011100001101111110011100000110011010100110101110111010100110001110110001110010101000100010111111001101111001100100100010001000010010101001111110110111001111001110111110011010000101101011111100011011001100011000011000111011010111110101010101101110010000000000110011111101111100100100101011000100011000110011000111100111010000100000001100011000001011010010001000111010001010101011001010011010001001110000010110010010100101101111000000100011100001001100001010010011110110001100010111100111101111011111110001101000001110100111010100001011010101001011000010001101000110100001101011101010101010101110110101111111001101001100110010111010110000111010010111010101011111110000100010000011001111101010100011010010110101011110110001101000000010010101001110000101010100111001011011100111000011000000011001100110011011011111001000101010001111000101111111000110000101001110010111011111110111110000110100001101101000111001101110001001011001111110111010000010011011000110001001101000111011010000011100000000011100100001010111110011111010010100001000001101101010101101000100001111001001111110010010111010110111010010000110011100101110101000110101001111101100111011000110010101111110110111000001010011111111101100001100000101110110001001110011100111101011110010111011101100001111111101010011110110000111111101000011100100110000001111001001000011110000000001011111010111111101110100110110000001110011111010000010000010100110101100010111111001000010101000111101110110101100110001100001011100100000011110000100010110110000010101111101101011101010000101110010000110101110010110101110111011001010110100010010010001010010101000111100010001110111001010100011111111010101011100000110010100000001100110110001010111010001111011100111110010101010000001110101110101011001111010110110101100111101111101011101001011001100000111100000110100000101001011111011010100101011001111010001101000010101110110011110010100111110111100101010110010010101101011001001101001010000010100010001101100000001010000010101110010110101101011010010111001001000000110110010011111111101111000010001100000110010000010101011001100100110010010110011001001111011000101011111011110010101111010101101011101100100111011101101011100101000000110001011110111011101010101111011010001110111101010100010010000001010110011001010100110100011100111001000011111111011001111101110001101100101011001001011000110010101010000100110110100101101011101000101110001101110110001110111101101111010010111010110101011101001110100110111101001101101011010000001100111101110011010110111110011010110110000100000011011111100011000100110001100111101001011110000001101010011010011010101101110010000100001010010010011101111000110001001101111000101101101101101011011011011111010011101110110000001100111010001101100000110101010100100010001000111011010101011001111000011011001000101100001011100110101010000010011011101000101011000100101101111001111011011011001110010100001010100011011111111100001011100011110101011000100010000111101010100000101000110101001000100010000100000110100101111001110001100011101100001111111000000011001101011000010110101111101111000111001001110000001001000010100111011101110111110001101001010111110100110101010001010001110101010101011010101011011001000100000010110110110001101001110001110111111011000000000100110100000001001101001011010100100111001010110110001110110011000101010010011111001100011101101110110001110111110101011010101110000100001110001001110011100101110101011101101101011000110100110100110101011011000111101011000100111110011101111000010011001101011101001110100100110011011000111001001010000100010100001001100000000100100010011001011110010110100100000100011001111000011001111110011001001001000010000000001011011000101000010111010000010011001010110010100110001011110010110101010111110100001100011101100101001011011100101001100000110000011000101000011111100010010111000000100101101101100101000000011001100110101111001110111111001011001010000110101110011100010111001010100001111100111100110000111001000010101010100001001110010010111110110000001111010111111011110110111001001000100011101010011000000100100001110110011000001100100010100111110101111001011101011011001010101000001010100110100101011010111001100101110011011111011111001110011010110011110101100011000100100111010010000011111111011011101111000001100101100000111011000111100001000110001000111011101011100011001010101101001001011100100010100101010101001010110110001011011100110111011000111100100010001101010011100110110001000010100000100101100000101011000101001001111111011111010000000110111100111110000001000100101100110101001011101110100110010011100010000111010010101011001001100100101011111111010110011000110100101000011101011111110101101001101100100100111001011000100111101111011001101010100110100111110001000001110110110011010000001111010001010111001101111100010111010000100010011100111110001010001010111101111001110011110000011110000101101100010110000110001001001010001010111110111001110010110100101101111011000110110100110111011110100000011110100111111100011110001000011100011000011101010000100100100110010101000111010000010010100101101110001110011100011001000100101001010001011010101111011001011010100011011111100001111110000010100010100011010010000011001110010011100110110001100000111011100010010110010101100010111110100110110011111101000011101111100100001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001100011011011000011101101000101000111011010000010011001111101011110111110011000100111011111111110000101011000111110101000101111101100011101000101100111011001111100000101111100100100001001010111100001100011101011001000101001001001110110011010101010001100000000010101001100010000101001000100111111100000001000000111000000101000110010010";
	Int n; Int::N=500;
	for (auto x:{x1,x2,x3}){
		n.read(x);
		int ans=heuristic_solve(n);
		printf("%d %.6lf\n",ans,ans/(n.bit_length()/log2(3)));
	}
}*/
int main(){
	srand(time(0));
	//test_heuristic();
	run();
	return 0;
}

