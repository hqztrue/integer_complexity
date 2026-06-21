#include<bits/stdc++.h>
using namespace std;
//const int f0[]={0,1,2,3,4,5,5,6,6};
void Minus(int &x,int y,int &cnt){
	assert(x>=y);
	x-=y;
	cnt+=y;
}
void Div(int &x,int y,int &cnt){
	assert(x>0);
	assert(x%y==0);
	x/=y;
	cnt+=y;
}
double Ratio(int v,int n){
	assert(n>1);
	return v/log(n)*log(3);
}
int reduce2(int &n,int a,int b){
	int cnt=0;
	for (int i=0;i<a;++i){
		Minus(n,1,cnt);
		Div(n,2,cnt);
	}
	Div(n,2,cnt);
	for (int i=0;i<b-1;++i){
		Minus(n,1,cnt);
		Div(n,3,cnt);
	}
	return cnt;
}
int reduce3(int &n,int a,int b){
	int cnt=0;
	for (int i=0;i<b;++i){
		Minus(n,2,cnt);
		Div(n,3,cnt);
	}
	bool sign=0;
	if (n%3==1)Minus(n,1,cnt),sign=1;
	Div(n,3,cnt);
	for (int i=0;i<a-2;++i){
		if (i%2==sign)Minus(n,1,cnt);
		Div(n,2,cnt);
	}
	return cnt;
}
int ipow(int x,int y){
	int v=1;
	for (int i=1;i<=y;++i)v*=x;
	return v;
}
int f(int n){
	if (n<=5)return n;
	if (n%2==0)return f(n/2)+2;
	if (n%4==1)return f(n/4)+5;
	if (n%8==3)return f(n/8)+8;
	if (n%3==0)return f(n/3)+3;
	if (n%3==1)return f(n/3)+4;
	if (n%9==2)return f(n/9)+8;
	if (n%9==5)return f(n/9)+9;
	int a=0,b=0;
	while ((n+1)%ipow(2,a+1)==0)++a;
	while ((n+1)%ipow(3,b+1)==0)++b;
	assert(a>=3);
	assert(b>=2);
	int n0=n,n1=n;
	int cnt0=reduce2(n0,a,b);
	int cnt1=reduce3(n1,a,b);
	//printf("%d %d %d %d\n",n,a,b,n1);
	double ratio0=Ratio(cnt0,ipow(2,a+1)*ipow(3,b-1));
	double ratio1=Ratio(cnt1,ipow(2,a-2)*ipow(3,b+1));
	//printf("%d %.6lf %.6lf\n",n,ratio0,ratio1);
	//printf("%d %d %d\n",n0,cnt0,f[n0]);
	return ratio0<ratio1?f(n0)+cnt0:f(n1)+cnt1;
}
int main()
{
	//freopen("1.in","r",stdin);
	//freopen("1.out","w",stdout);
	//for (int i=2;i<9;++i)printf("%.6lf\n",Ratio(f0[i],i));
	int T=1000000000;
	double ma=0;
	for (int i=2;i<=T;++i){
		int ans=f(i);
		double r=Ratio(ans,i);
		if (r>ma){
			ma=r;
			printf("%d %d %.6lf %.6lf\n",i,ans,r,ma);
		}
	}
	//system("pause");for (;;);
	return 0;
}
