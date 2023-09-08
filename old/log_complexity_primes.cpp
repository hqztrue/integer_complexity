#include<bits/stdc++.h>
using namespace std;
typedef unsigned char uchar;
typedef unsigned int uint;
typedef unsigned long long ull;
const int N=10000005;
const uint L=36,inf=127;
uchar *a;
uint calc_g(uchar n){
	uint res=1;
	while (n>=5||n==3)res*=3,n-=3;
	return res<<n/2;
}
void calc_all(uint n){  //computes f(i) for all i<=n.
	a=(uchar*)malloc(n+1); a[1]=1;
	for (uint i=2;i<=n;++i)a[i]=inf;
	uint g[inf+1];
	for (uchar i=0;i<=inf;++i)g[i]=calc_g(i);
	for (uint i=2;i<=n;++i){
		if (a[i-1]+1<a[i])a[i]=a[i-1]+1;
		uchar t=a[i-1],k=1;
		while (g[k]+g[t-k]>=i&&k<t/2)++k;
		for (uint j=6,limit=g[k];j<=limit;++j)
			if (a[j]+a[i-j]<a[i])a[i]=a[j]+a[i-j];
		for (uint j=2,ij=i*2;j<=i&&ij<=n;++j,ij+=i)
			if (a[i]+a[j]<a[ij])a[ij]=a[i]+a[j]; 
	}
}
double log_comp(int x){
	return a[x]/log(x)*log(3);
}
pair<double,int> v[N];
int p[N],p1,g[N];
int get_prime(int p[],int n){
	int len=0;
	memset(g+1,0,sizeof(int)*n);
	for (int i=2;i<=n;++i){
		if (!g[i])g[i]=p[++len]=i;
		for (int j=1;j<=len&&p[j]<=g[i]&&i*p[j]<=n;++j)g[i*p[j]]=p[j];
	}
	return len;
}
int main()
{
	//freopen("1.in","r",stdin);
	//freopen("1.txt","w",stdout);
	int n=1000000;
	int t1=clock();
	calc_all(n);
	p1=get_prime(p,n);
	for (int i=1;i<=n;++i)v[i].first=log_comp(i),v[i].second=i;
	sort(v+1,v+1+n);
	for (int i=1;i<=10000;++i){
		int x=v[i].second;
		if (g[x]==x)printf("%.8lf %d\n",v[i].first,x);
	}
	printf("time=%d\n",clock()-t1);
	//system("pause");for (;;);
	return 0;
}

