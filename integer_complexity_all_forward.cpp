#include<bits/stdc++.h>
using namespace std;
const int N=1000005,inf=~0u>>2;
int f[N],n;
inline void upd(int &x,int y){if (y<x)x=y;}
vector<int> min_plus_conv(const vector<int> &a,const vector<int> &b){  //O(n^2) naive implementation
	int n=a.size(),m=b.size();
	vector<int> c(n+m-1);
	for (int i=0;i<n+m-1;++i)c[i]=inf;
	for (int i=0;i<n;++i)
		for (int j=0;j<m;++j)upd(c[i+j],a[i]+b[j]);
	return c;
}
void func(int l,int r){
	if (l==r){
		for (int i=2;i<=min(l,n/l);++i)upd(f[l*i],f[l]+f[i]);
		return;
	}
	int m=(l+r)/2;
	func(l,m);
	vector<int> a(f+l,f+m+1),b(f+1,f+r-l+1);
	vector<int> c=min_plus_conv(a,b);
	for (int i=m+1;i<=r;++i)upd(f[i],c[i-l-1]);
	func(m+1,r);
}
void calc(int n){
	f[1]=1; for (int i=2;i<=n;++i)f[i]=inf;
	func(1,n);
	for (int i=1;i<=n;++i)printf("%d %d\n",i,f[i]);
}
int main()
{
	//freopen("1.in","r",stdin);
	freopen("2.txt","w",stdout);
	int t1=clock();
	n=10000;
	calc(n);
	//printf("time=%d\n",clock()-t1);
	//system("pause");for (;;);
	return 0;
}

