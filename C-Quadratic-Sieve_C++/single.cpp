#include "avl.c"            // the trees.
#include "cint.c"           // the integers.
#include "fac_headers.h"    // factor headers.
#include "fac_utils.c"      // utilities and front-end.
#include "fac_quadratic.c"  // quadratic sieve source.
#include "fac_lanczos.c"    // quadratic sieve Lanczos.
#include "fac_testing.c"    // quadratic sieve tests.
//#include<omp.h>
#include<bits/stdc++.h>
using namespace std;
#include "factorization.h"  // pollard-rho for <62 bits
#include "statistics.h"

const uchar inf=127;
const ushort inf1=10000;
const ull U=1ull<<62;
unordered_map<ull,pair<uchar,ushort>> H;
unordered_map<u128,pair<uchar,ushort>> H1;
uint g[inf+1];
uchar *a;
uint n0; u128 n,N0;
void print(u128 x){
	if (x<0)putchar('-'),x=-x;
	if (x>9)print(x/10);
	putchar(x%10+'0');
}
void println(u128 x){print(x); putchar('\n');}
u128 u128_from_str(const char *s){u128 x=0; while (*s)x=x*10+*s++-'0'; return x;}
inline u128 _rand128(){static u128 x=u128_from_str("199609092119960909211996090921996090921");x+=(x<<17)+(x>>29)+1;return x;}

u128 calc_g(int n){
	u128 res=1;
	while (n>=5||n==3)res*=3,n-=3;
	return res<<n/2;
}
int defect(u128 n,int cur){
	const double eps=1e-3;
	for (int i=1;;++i)
		if (calc_g(i)*(1+eps)>=n)return max(cur-i,0);
}
void calc_all(uint n){  //computes f(i) for all i<=n.
	a=(uchar*)malloc(n+1); a[1]=1;
	for (uint i=2;i<=n;++i)a[i]=inf;
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
template<class T>
const vector<pair<T,uchar>> prime_factors(T x){
	vector<pair<T,uchar>> ans;
	u64 *pf = factorize(x);
	int m=pf[64];
	for (int i=0;i<m;++i){
		int cnt=0;
		while (x%pf[i]==0)x/=pf[i],++cnt;
		ans.emplace_back(pf[i],cnt);
	}
	free(pf);
	return ans;
}
/*static inline void fac_display_verbose(fac_cint ** ans) {
	for(int i = 0; i < 100; ++i)
		putchar(' ');
	putchar('\r');
	char * str = fac_answer_to_string(ans);
	puts(str);
	free(str);
}*/
template<class T>
const vector<pair<T,uchar>> prime_factors128(T x){  // factorize 128 bits
	vector<pair<T,uchar>> a;
	cint N;
	fac_params config = {0};
	static char s[105],*p;
	//sprintf(s,"%I64d",x);
	for (p=s;x;++p)*p=x%10+'0',x/=10;
	reverse(s,p); *p=0;
	const int bits = 64 + 4 * (int) strlen(s);
	cint_init_by_string(&N, bits, s, 10); // init the number as a cint.
	fac_cint ** ans = c_factor(&N, &config); // execute the routine.
	
	//fac_display_verbose(ans); // print answer.
	for(char i = 0; ans[i]; ++i){
		char * s = cint_to_string(&ans[i]->cint, 10);
		T y=0;
		for (char *p=s;*p;++p)y=y*10+*p-'0';
		a.emplace_back(y,ans[i]->power);
		free(s);
	}
	//for (auto &t:a)printf("%I64d %d\n",t.first,t.second);
	return a;
}
template<class T>
vector<T> get_factors_from_primes(const vector<pair<T,uchar>> &p){
	vector<T> a; a.push_back(1);
	for (int i=0;i<p.size();++i){
		int n=a.size(),t=p[i].second; T p0=p[i].first;
		for (T m=p0;t;--t,m*=p0)
			for (int j=0;j<n;++j)a.push_back(a[j]*m);
	}
	return a;
}
bool is_prime_slow(u128 x){
	auto v=prime_factors128(x);
	return v.size()==1&&v[0].second==1;
}
template<class T>
vector<T> get_factors(T x){
	if (sizeof(T)==sizeof(ull))return get_factors_from_primes(prime_factors(x));
	else return get_factors_from_primes(prime_factors128(x));
}
ushort dfs(ull x,int t=10){  //computes f(x) for a single x, overall using <=t 1's for addition.
	if (x<=n0)return a[x];
	auto it=H.find(x);
	if (it!=H.end()&&it->second.first>=t)return it->second.second;
	ushort ans=inf1;
	for (int i=1;i<=g[t]&&i<x;++i)
		if (a[i]<=t)ans=min(ans,ushort(dfs(x-i,t-a[i])+a[i]));
	for (auto &g:get_factors(x))
		if (g>1&&g<=x/g)ans=min(ans,ushort(dfs(g,t)+dfs(x/g,t)));
	H[x]=make_pair(t,ans);
	return ans;
}
ushort dfs128(u128 x,int t=10){
	if (x<=n0)return a[x];
	if (x<U)return dfs(x,t);
	auto it=H1.find(x);
	if (it!=H1.end()&&it->second.first>=t)return it->second.second;
	ushort ans=inf1;
	for (int i=1;i<=g[t]&&i<x;++i)
		if (a[i]<=t)ans=min(ans,ushort(dfs128(x-i,t-a[i])+a[i]));
	for (auto &g:get_factors(x))
		if (g>1&&g<=x/g)ans=min(ans,ushort(dfs128(g,t)+dfs128(x/g,t)));
	H1[x]=make_pair(t,ans);
	return ans;
}

void init(uint _n0=1e9,u128 _N0=1e23){
	int t1=clock();
	n0=_n0; N0=_N0;
	calc_all(n0);
	printf("time=%d\n",clock()-t1);
}
void check1(){  // test the conjecture f(p^i)=i*f(p).
	int t1=clock();
	vector<int> primes={577,811,109};
	//vector<int> primes={433,163,487,2};
	//vector<int> primes={2};
	for (auto mul:primes){
		printf("mul=%I64d\n",mul);
		u128 x=mul;
		for (int i=1;;++i,x*=mul){
			//if (i<10)continue;
			n=x;
			//printf("i=%d n=%I64d\n",i,n);
			int v=i*a[mul],d=defect(n,v-1);  // v: the conjectured complexity. d: defect.
			printf("def=%d\n",d);
			int ans=dfs128(n,d);
			printf("i=%d ans=%d\n",i,ans);
			//printf("time=%d\n",clock()-t1);
			//printf("size=%d\n",H.size());
			if (ans<v){
				printf("err: %d^%d=",mul,i);
				print(n);
				printf(" ans=%d tgt=%d\n",ans,v);
				//exit(0);
				break;
			}
			if (n>N0/mul)break;
		}
	}
	printf("time=%d\n",clock()-t1);
}
void check2(){  // test the conjecture f(2^i3^j5^k)=2i+3j+5k (k<=5)
	int t1=clock();
	u128 x=2;
	for (int i=1;;++i,x*=2){
		//if (i<40)continue;
		if (x>N0)break;
		u128 y=1;
		for (int j=0;;++j,y*=3){
			if (x*y>N0)break;
			u128 z=1;
			for (int k=0;k<=5;++k,z*=5){
				n=x*y*z;
				if (n>N0)break;
				if (n==1)continue;
				if (i<60||i==60&&j<3||i==60&&j==3&&k<5)continue;
				//if (i<58)continue;
				printf("i=%d j=%d k=%d ",i,j,k); println(n);
				int v=i*2+j*3+k*5,d=defect(n,v-1);
				int ans=dfs128(n,d);
				if (ans<v){
					printf("err: %d %d",ans,v); println(n);
					exit(0);
				}
			}
		}
	}
	printf("time=%d\n",clock()-t1);
}
void check3(){  // test the conjecture f(2^i+1)=2i+1.
	int t1=clock();
	vector<int> primes={2};
	for (auto mul:primes){
		printf("mul=%I64d\n",mul);
		u128 x=mul;
		for (int i=1;;++i,x*=mul){
			if (i==3||i==9)continue;
			if (i<=76)continue;
			n=x+1;
			//printf("i=%d n=%I64d\n",i,n);
			int v=i*2+1,d=defect(n,v-1);
			printf("def=%d\n",d);
			int ans=dfs128(n,d);
			printf("i=%d ans=%d\n",i,ans);
			//printf("time=%d\n",clock()-t1);
			//printf("size=%d\n",H.size());
			if (ans<v){
				printf("err: %d^%d+1=",mul,i);
				print(n);
				printf(" ans=%d tgt=%d\n",ans,v);
				//exit(0);
				break;
			}
			if (n>N0/mul)break;
		}
	}
	printf("time=%d\n",clock()-t1);
}
void check4(){  // test the conjecture f(3p)=min{f(3p-1)+1,f(p)+3} for prime p.
	int t1=clock();
	u128 U=1e20,u=_rand128()%U+1;
	//u=U;
	const int d=7;
	for (u128 p=u;;){
		//++p;
		p=_rand128()%U+1+rand();
		if (!is_prime_slow(p))continue;
		printf("p="); println(p);
		auto v1=dfs128(p*3-1,d),v2=dfs128(p,d),v3=dfs128(p*3,d);
		//println(v1); println(v2); println(v3);
		if (min(v1+1,v2+3)!=v3){
			printf("err: p=",p); println(p);
			printf("3p="); println(p*3);
			printf(" f(3p-1)=%d f(p)=%d f(3p)=%d\n",v1,v2,v3);
			//exit(0);
			break;
		}
	}
	printf("time=%d\n",clock()-t1);
}
void run_sample(int num_samples=1e4){
	//freopen("data.txt","w",stdout);
	int t1=clock(),t2=t1;
	vector<double> a;
	for (int i=0;i<1000*rand();++i)_rand128();
	for (int i=0;i<num_samples;++i){
		//if (i%100==0)printf("i=%d\n",i);
		if (clock()-t2>6e4){
			printf("i=%d\n",i);
			t2=clock();
		}
		//n=N0-i;
		n=_rand128()%N0+2;
		//int d=6;
		ushort ans=10000;
		for (int d=0;d<=1000&&d<=defect(n,ans-1);++d){
			//printf("def=%d\n",d);
			ans=min(ans,dfs128(n,d));
		}
		double x=ans/log(n)*log(3);
		//printf("%d, ",i); print(n);
		//printf(", %d, %.6lf, %.6lf\n",ans,x,1-1.*n/calc_g(ans));
		a.push_back(x);
		//fflush(stdout);
	}
	double ave=mean(a),mu=stddev(a);
	printf("N0=%.2lf ",log10(N0)); println(N0);
	printf("#samples=%d\n",num_samples);
	printf("mean=%.6lf stddev=%.6lf\n",ave,mu);
	printf("time=%d\n",clock()-t1);
	//fclose(stdout);
}
void test(){
	//prime_factors(123456789012345678);
	
	/*//n=162879576091729561ull;
	n=155104303499468569ull;
	printf("%d\n",defect(n,119));
	cout<<(int)dfs(n,2)<<endl;
	printf("time=%d\n",clock()-t1);
	exit(0);*/
}
int main()
{
	srand(time(0));
	//u128 N0=1; //N0<<=120;
	//for (int i=1;i<=23;++i)N0*=10;
	//init(1e7,1e16);
	init(1e8,1e12);
	//init(1e9,1e20);
	//init(2e9,1e25);
	
	//check1();
	//check2();
	//check3();
	//check4();
	//run_sample(1e4);
	for (u128 i=1;i<=1e14;i*=10){
		N0=i;
		run_sample(1e6);
	}
    return 0;
}

