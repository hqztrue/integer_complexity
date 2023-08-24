#include<bits/stdc++.h>
//#include<omp.h>
using namespace std;
typedef unsigned char uchar;
typedef unsigned int uint;
typedef unsigned long long ull;
const uint N=1e8,N1=2e7,inf=127,L=60;  //This L is enough for our application.
uchar f[N1+1][L],*a;  //f[i][j] stores the integer complexity for floor(n/i)-j.
int p[N],p1,n0;
int get_prime(int p[],int n){
	int len=0,*f=new int[n+1];
	memset(f+1,0,sizeof(int)*n);
	for (int i=2;i<=n;++i){
		if (!f[i])f[i]=p[++len]=i;
		for (int j=1;j<=len&&p[j]<=f[i]&&i*p[j]<=n;++j)f[i*p[j]]=p[j];
	}
	delete[] f;return len;
}
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
const vector<ull>& get_factors(const vector<pair<uint,int>> &p){
	static vector<ull> a;
	a.clear(); a.push_back(1);
	for (int i=0;i<p.size();++i){
		int n=a.size(),t=p[i].second; uint p0=p[i].first;
		for (ull m=p0;t;--t,m*=p0)
			for (int j=0;j<n;++j)a.push_back(a[j]*m);
	}
	return a;
}
uchar calc_single(ull n){  //computes a single f(n).
	int t1=clock();
	/*int n0=1e8;
	calc_all(n0);
	printf("time=%d\n",clock()-t1);
	p1=get_prime(p,sqrt(n));
	printf("time=%d\n",clock()-t1);*/
	if (n<=n0)return a[n];
	memset(f,0x7f,sizeof(f));
	vector<pair<uint,int>> factors[L];
	for (int i=n/n0;i;--i){
		uchar *f0=f[i]; ull cur=n/i;
		int s=sqrt(cur);
		int pend=upper_bound(p,p+p1,s)-p;
		for (int j=0;j<L;++j)factors[j].clear();
		for (int j=1;j<pend;++j){
			uint p0=p[j];
			if (cur%p0>=L)continue;
			for (ull t=cur/p0*p0;cur-t<L;t-=p0){
				ull x=t; int cnt=0;
				while (x%p0==0)x/=p0,++cnt;
				factors[cur-t].emplace_back(p0,cnt);
			}
		}
		for (int j=L-1;j>=0;--j){
			ull t=cur-j;
			const vector<ull>& g=get_factors(factors[j]);
			for (auto &k:g)
				if (k<=t/k){
					ull r=t/k; uchar v=a[k];
					if (r<=n0)v+=a[r];
					else {
						int id=n/r,len=n/id-r;
						if (len>=L)continue;
						v+=f[id][len];
					}
					if (v<f0[j])f0[j]=v;
				}
		}
		for (int j=L-1;j>=0;--j)  //naive (min,+)-conv
			for (int k=L-1;k>j;--k)
				if (f0[k]+a[k-j]<f0[j])f0[j]=f0[k]+a[k-j];
	}
	//printf("%d\n",f[1][0]);
	printf("time=%d\n",clock()-t1);
	//delete []a;
	return f[1][0];
}
int main()
{
	//freopen("1.in","r",stdin);
	//freopen("1.out","w",stdout);
	int t1=clock();
	ull n;
	//n=1LL<<42;
	//calc_all(n);
	n0=1e8;
	calc_all(n0);
	ull N=1LL<<52,N0=1LL<<50;
	assert(N0<=N*N1);
	p1=get_prime(p,sqrt(N));
	
	// test the conjecture f(2^i3^j5^k)=2i+3j+5k.
	ull x=1;
	for (int i=0;;++i,x*=2){
		if (i<12)continue;
		ull y=1;
		for (int j=0;;++j,y*=3){
			ull z=1;
			if (x*y>N0)break;
			for (int k=0;k<=5;++k,z*=5){
				n=x*y*z;
				if (n>N0)break;
				if (n==1)continue;
				printf("i=%d j=%d k=%d n=%I64d\n",i,j,k,n);
				uchar ans=calc_single(n);
				uchar v=i*2+j*3+k*5;
				if (ans!=v){
					printf("err: %I64d %d %d\n",n,ans,v);
					exit(0);
				}
			}
		}
	}
	//printf("time=%d\n",clock()-t1);
	//system("pause");for (;;);
	return 0;
}

