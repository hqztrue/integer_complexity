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
#include "allocator.h"
#include "utils.h"
#include "c_rho-128-bits.h"

const uchar inf=127;
const ushort inf1=10000;
const ull U=1ull<<62;
const u128 zero=0,inf128=~zero>>1;
unordered_map<ull,pair<uchar,ushort>> H;
unordered_map<u128,pair<uchar,ushort>> H1;
//unordered_map<ull,pair<uchar,ushort>,std::hash<ull>,std::equal_to<ull>,myallocator<ull>> H;
//unordered_map<u128,pair<uchar,ushort>,std::hash<u128>,std::equal_to<u128>,myallocator<u128>> H1;
u128 g[inf1],g1;
uchar *a;
uint n0; u128 n,N0;
uint log3_floor(u128 n){
	uint ans=0;
	while (n>=3)n/=3,++ans;
	return ans;
}
u128 calc_g(int n){  // compute the maximum integer with complexity n
	u128 res=1;
	while (n>=5||n==3)res*=3,n-=3;
	return res<<n/2;
}
int defect_approx(u128 n,int cur){
	const double eps=1e-3;
	for (int i=1;;++i)
		if (g[i]*(1+eps)>=n)return max(cur-i,0);
}
int complexity_LB_naive(u128 n){
	//for (int i=1;;++i)
	//	if (g[i]>=n)return i;
	return lower_bound(g+1,g+g1,n)-g;
}
unordered_map<u128,int> M_lb;
int complexity_LB(u128 n){
	if (n==1)return 1;
	int l=upper_bound(g+1,g+g1,n)-g-1;
	auto it=M_lb.find(n);
	if (it!=M_lb.end())return l+it->second;
	return l+2;
}
int _init=[](){ // precompute
	for (int i=0;;++i){
		g[i]=calc_g(i);
		if (i&&g[i]<g[i-1]){g1=i; break;}
	}
	for (int i=0;i<300;++i)g[g1++]=inf128;
	
	for (u128 i=1,i1=0;i1<=10;i*=2,++i1)
		for (u128 j=1,j1=0;;j*=3,++j1){
			if (j<=inf128/i)M_lb[i*j]=1;
			if (j>inf128/3)break;
		}
	for (u128 i=1,i1=0;i1<=2;i*=2,++i1)
		for (u128 j=1,j1=0;;j*=3,++j1){
			for (u128 k=1,k1=0;k1<=2;k*=2,++k1)
				for (u128 l=1,l1=0;;l*=3,++l1){
					if (i1+k1<=2&&1.*i*(1.*k*l+1)*j<=1.*inf128)M_lb[i*(k*l+1)*j]=1;
					if (l>inf128/3)break;
				}
			if (j>inf128/3)break;
		}
	for (u128 i=1,i1=0;i1<=2;i*=2,++i1)
		for (u128 j=1,j1=0;;j*=3,++j1){
			if (j<=inf128/i)M_lb[i*j]=0;
			if (j>inf128/3)break;
		}
	return 0;
}();
void calc_all(uint n){  // computes f(i) for all i<=n.
	a=(uchar*)malloc(n+1); a[0]=0; a[1]=1;
	for (uint i=2;i<=n;++i)a[i]=inf;
	uint g32[inf+1];
	for (uchar i=0;i<=inf;++i)g32[i]=calc_g(i);
	for (uint i=2;i<=n;++i){
		if (a[i-1]+1<a[i])a[i]=a[i-1]+1;
		uchar t=a[i-1],k=1;
		while (k<t/2&&g32[k+1]+g32[t-k-1]>=i)++k;
		for (uint j=6;j<=g32[k];++j)
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
const vector<pair<T,uchar>> prime_factors128_old(T x){  // factorize 128 bits
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
		if (a.size()&&a[a.size()-1].first==y)a[a.size()-1].second+=ans[i]->power;
		else a.emplace_back(y,ans[i]->power);
		free(s);
	}
	//for (auto &t:a)printf("%I64d %d\n",t.first,t.second);
	return a;
}
template<class T>
const vector<pair<T,uchar>> prime_factors128(T x){  // factorize 128 bits
	vector<pair<T,uchar>> a;
	// allocate memory for 128 factors.
    positive_number *factors = (positive_number*)calloc(128, sizeof(positive_number)), n = x;
	// generate a random number of ~ n_bits bits.
	//printf("n=%40s\n", to_string_128_bits(n));
	// fill the "factors" array with the prime factors.
	factor(n, factors);
	// iterate over the factors (zero terminated array).
	vector<u128> v;
	for (int j = 0; factors[j]; ++j)v.push_back(factors[j]);
	sort(v.begin(),v.end());
	for (int i=0;i<v.size();++i)
		if (i==0||v[i]!=v[i-1])a.push_back({v[i],1});
		else ++a[a.size()-1].second;
    // release memory.
    free(factors);
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
ushort dfs(ull x,int t){
	if (x<=n0)return a[x];
	auto it=H.find(x);
	if (it!=H.end()&&it->second.first>=t)return it->second.second;
	ushort ans=inf1,k=1;
	while (k<t/2&&g[k+1]+g[t-k-1]>=x)++k;
	assert(g[k]<=n0);
	for (int i=1;i<=g[k];++i){
		ushort lb=complexity_LB(x-i);
		if (a[i]+lb<=t)ans=min(ans,ushort(dfs(x-i,t-a[i])+a[i]));
	}
	auto fac=get_factors(x);
	for (auto &g:fac)
		if (g>1&&g<=x/g){
			ushort lb1=complexity_LB(g),lb2=complexity_LB(x/g);
			if (lb1+lb2>t)continue;
			//ushort v1=dfs(g,t-lb2);
			ushort v1=inf1;
			for (int i=lb1;i<=t-lb2&&i<v1;++i){
				v1=min(v1,dfs(g,i));
			}
			if (v1+lb2>t)continue;
			ushort v2=dfs(x/g,t-v1);
			ans=min(ans,ushort(v1+v2));
		}
	H[x]=make_pair(t,ans);
	return ans;
}
ushort dfs128(u128 x,int t){  //decides whether f(x)<=t. If true, return the optimal f(x).
	if (x<=n0)return a[x];
	if (x<U)return dfs(x,t);
	auto it=H1.find(x);
	if (it!=H1.end()&&it->second.first>=t)return it->second.second;
	ushort ans=inf1,k=1;
	while (k<t/2&&g[k+1]+g[t-k-1]>=x)++k;
	assert(g[k]<=n0);
	for (int i=1;i<=g[k];++i){
		ushort lb=complexity_LB(x-i);
		if (a[i]+lb<=t)ans=min(ans,ushort(dfs128(x-i,t-a[i])+a[i]));
	}
	auto fac=get_factors(x);
	for (auto &g:fac)
		if (g>1&&g<=x/g){
			ushort lb1=complexity_LB(g),lb2=complexity_LB(x/g);
			if (lb1+lb2>t)continue;
			/*if (t-complexity_LB(x)>=3){
				if (dfs128(x/g,lb2)>lb2)++lb2;
				if (lb1+lb2>t)continue;
			}*/
			//ushort v1=dfs128(g,t-lb2);
			ushort v1=inf1;
			for (int i=lb1;i<=t-lb2&&i<v1;++i){
				v1=min(v1,dfs128(g,i));
			}
			if (v1+lb2>t)continue;
			ushort v2=dfs128(x/g,t-v1);
			ans=min(ans,ushort(v1+v2));
		}
	H1[x]=make_pair(t,ans);
	return ans;
}
ushort calc_single(u128 n){  //compute f(n)
	int lb=complexity_LB(n);
	int ans=inf1;
	for (int t=lb;t<ans;++t){
		int res=dfs128(n,t);
		ans=min(ans,res);
	}
	return ans;
}

void init(uint _n0=1e9,u128 _N0=1e23){
	int t1=clock();
	n0=_n0; N0=_N0;
	calc_all(n0);
	printf("init time=%d\n",clock()-t1);
}
void clear(){
	H.clear(); H1.clear();
	//unordered_map<ull,pair<uchar,ushort>>().swap(H);
	//unordered_map<u128,pair<uchar,ushort>>().swap(H1);
}
void check1(){  // test the conjecture f(p^i)=i*f(p). (In particular, f(2^i)=2i.)
	int t1=clock();
	//vector<int> primes={733,379,739,541};  // conjecture fails
	//vector<int> primes={577,811,109};
	//vector<int> primes={433,163,487,2};
	//vector<int> primes={2};
	vector<int> primes={811};
	for (auto mul:primes){
		printf("mul=%I64d\n",mul);
		u128 x=mul;
		for (int i=1;;++i,x*=mul){
			//if (i<80)continue;
			n=x;
			//printf("i=%d n=%I64d\n",i,n);
			int v=i*a[mul],lb=complexity_LB(n);  // v: the conjectured complexity. d: defect.
			//printf("def=%d\n",d);
			int ans=v;
			for (int t=lb;t<v;++t){
				int res=dfs128(n,t);
				ans=min(ans,res);
				if (ans<v){
					printf("improve: %d %d\n",ans,v);
					break;
				}
			}
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
				int v=i*2+j*3+k*5,d=defect_approx(n,v-1);
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
			int v=i*2+1,d=defect_approx(n,v-1);
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
	u128 U=1e20,u=rand128()%U+1;
	//u=U;
	const int d=7;
	for (u128 p=u;;){
		//++p;
		p=rand128()%U+1;
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
void run_sample(int num_samples=1e4,bool debug=1){
	//freopen("data.txt","w",stdout);
	int t1=clock(),t2=t1;
	vector<double> a;
	clear();
	for (int i=0;i<num_samples;++i){
		//if (i%100==0)printf("i=%d\n",i);
		if (clock()-t2>6e4){
			if (debug)printf("i=%d\n",i);
			t2=clock();
		}
		//n=N0-i;
		n=rand128()%N0+2;
		ushort ans=calc_single(n);
		double x=ans/log(n)*log(3);
		//printf("%d, ",i); print(n);
		//printf(", %d, %.6lf, %.6lf\n",ans,x,1-1.*n/calc_g(ans));
		a.push_back(x);
		//fflush(stdout);
	}
	double ave=mean(a),mu=stddev(a);
	if (debug){
		printf("N0=%.2lf ",log10(N0)); println(N0);
		printf("#samples=%d\n",num_samples);
		printf("mean=%.6lf stddev=%.6lf\n",ave,mu);
		printf("time=%d\n",clock()-t1);
	}
	else {
		printf("%.6lf,",ave);
	}
	//fclose(stdout);
}
void verify(int T=1e4){
	init(1e8,1e23);
	u128 U=1e8-1; n0=1e4;
	int t1=clock();
	for (int i=1;i<=T;++i){
		u128 n=rand128()%U+2;
		assert(calc_single(n)==a[n]);
		//assert(dfs128(n,a[n])==a[n]);
	}
	printf("verify time=%d\n",clock()-t1);
}
void factorize_test(u128 N=1e18,int T=1000){
	int t1=clock();
	for (int i1=1;i1<=T;++i1){
		u128 n=rand128()%N+1;
		//println(n);
		auto a1=prime_factors128(n);
		auto a2=prime_factors128_old(n);
		if (a1!=a2){
			prln(a1);
			prln(a2);
			exit(0);
		}
	}
	printf("fac time=%d\n",clock()-t1);
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
	
	//verify();
	//factorize_test(1e18,1e3);
	//return 0;
	
	//init(1e6,1e18);
	//init(1e7,1e30);
	init(1e8,1e27);
	//init(1e9,1e36);
	//init(1e9,1e20);
	//init(2e9,1e35);
	
	check1();
	//check2();
	//check3();
	//check4();
	//run_sample(1e4);
	/*while (1){
		N0=1e19;
		run_sample(1000,0);
		return 0;
	}*/
	/*for (u128 i=1;i<=1e30;i*=10){
		if (i<=1e10)continue;
		N0=i;
		//run_sample(1e6);
		run_sample(1e5);
	}*/
    return 0;
}

