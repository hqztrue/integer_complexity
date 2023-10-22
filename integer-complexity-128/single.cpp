/*
	Author: hqztrue.
	paper: "Improved Algorithms for Integer Complexity".
	link: https://arxiv.org/abs/2308.10301
*/
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
#include "utils.h"
#include "factorization.h"  // pollard-rho for <62 bits
#include "statistics.h"
#include "allocator.h"
#include "c_rho-128-bits.h"

//f(n): integer complexity of n.
//f'(n): integer complexity of n, where subtractions are allowed.
const uchar inf=127;
const ushort inf1=10000;
const uint inf32=~0u>>1;
const ull U=1ull<<62;
const u128 zero=0,inf128=~zero>>1;
unordered_map<ull,pair<ushort,ushort>> H;  //memoization for dfs()
unordered_map<u128,pair<ushort,ushort>> H1;  //memoization for dfs128()
//unordered_map<ull,pair<ushort,ushort>,std::hash<ull>,std::equal_to<ull>,myallocator<ull>> H;
//unordered_map<u128,pair<ushort,ushort>,std::hash<u128>,std::equal_to<u128>,myallocator<u128>> H1;
u128 g[inf1],g1;
uchar *a;  //record f(n) (n<=n0).
string s[10000005];  //record the formulas for f(n) (n<=n0).
uint n0;  //threshold for precomputing f(n) (n<=n0).
u128 N0;  //we only compute f(n) for n<=N0.
ull CNT; int T0=clock();  //timer for the current progress
uint log3_floor(u128 n){
	uint ans=0;
	while (n>=3)n/=3,++ans;
	return ans;
}
u128 calc_g(int n){  //computes the maximum integer with complexity n.
	u128 res=1;
	while (n>=5||n==3)res*=3,n-=3;
	return res<<n/2;
}
int defect_approx(u128 n,int cur){  //approximately computes the defect.
	const double eps=1e-3;
	for (int i=1;;++i)
		if (g[i]*(1+eps)>=n)return max(cur-i,0);
}
int complexity_LB_naive(u128 n){  //naive lower bound for f(n).
	//for (int i=1;;++i)
	//	if (g[i]>=n)return i;
	return lower_bound(g+1,g+g1,n)-g;
}
unordered_map<u128,int> M_lb;  //record the numbers with low defect.
int complexity_LB(u128 n){  //better lower bound for f(n). (improvable)
	//return complexity_LB_naive(n);
	if (n==1)return 1;
	int l=upper_bound(g+1,g+g1,n)-g-1;
	auto it=M_lb.find(n);
	if (it!=M_lb.end())return l+it->second;
	return l+2;
}
int complexity_UB(u128 n){
	if (n==1)return 1;
	const double c=4.125;
	return floor(c*log((double)n)/log(3));
}
int _init=[](){ //precomputation
	for (int i=0;;++i){
		g[i]=calc_g(i);
		if (i&&g[i]<g[i-1]){g1=i; break;}
	}
	for (int i=0;i<300;++i)g[g1++]=inf128;
	
	//precompute the complexity lower bounds.
	for (u128 i=1,i1=0;i1<=10;i*=2,++i1)  //defect 1
		for (u128 j=1,j1=0;;j*=3,++j1){
			if (j<=inf128/i)M_lb[i*j]=1;
			if (j>inf128/3)break;
		}
	for (u128 i=1,i1=0;i1<=2;i*=2,++i1)  //defect 1
		for (u128 j=1,j1=0;;j*=3,++j1){
			for (u128 k=1,k1=0;k1<=2;k*=2,++k1)
				for (u128 l=1,l1=0;;l*=3,++l1){
					if (i1+k1<=2&&1.*i*(1.*k*l+1)*j<=1.*inf128)M_lb[i*(k*l+1)*j]=1;
					if (l>inf128/3)break;
				}
			if (j>inf128/3)break;
		}
	for (u128 i=1,i1=0;i1<=2;i*=2,++i1)  //defect 0
		for (u128 j=1,j1=0;;j*=3,++j1){
			if (j<=inf128/i)M_lb[i*j]=0;
			if (j>inf128/3)break;
		}
	return 0;
}();

void calc_all(uint n){  //computes f(i) for all i<=n. From Martin N. Fuller.
	a=(uchar*)malloc(n+1); a[0]=0; a[1]=1;
	for (uint i=2;i<=n;++i)a[i]=inf;
	uint g32[inf+1];
	for (uchar i=0;i<=inf;++i)g32[i]=min(calc_g(i),inf32);
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
u128 max_sub(u128 n,ushort f){
	return g[f-complexity_LB_naive(n)];
}
void calc_all_subtraction(uint n1){  //computes f'(i) for all i<=n. From Janis Iraids.
	uint g32[inf+1];
	for (uchar i=0;i<=inf;++i)g32[i]=min(calc_g(i),(u128)inf32);
	auto max_sub=[&](uint n,uchar f){
		return g32[f-complexity_LB_naive(n)];
	};
	uint n=n1+max_sub(n1,complexity_UB(n1));
	a=(uchar*)malloc(n+1); a[0]=0; a[1]=1;
	for (uint i=2;i<=n;++i)a[i]=inf;
	bool flag=1;
	while (flag){
		flag=0;
		for (uint i=2;i<=n;++i){
			if (a[i-1]+1<a[i])a[i]=a[i-1]+1,flag=1;
			uchar t=a[i-1],k=1;
			while (k<t/2&&g32[k+1]+g32[t-k-1]>=i)++k;
			for (uint j=6;j<=g32[k];++j)
				if (a[j]+a[i-j]<a[i])a[i]=a[j]+a[i-j],flag=1;
			for (uint j=2,ij=i*2;j<=i&&ij<=n;++j,ij+=i)
				if (a[i]+a[j]<a[ij])a[ij]=a[i]+a[j],flag=1;
		}
		for (uint i=n1;i;--i){
			uint d=max_sub(i,a[i]);
			for (uint j=1;j<=d;++j)
				if (a[i+j]+a[j]<a[i])a[i]=a[i+j]+a[j],d=max_sub(i,a[i]),flag=1;
		}
	}
	/*ll s=0;
	for (int i=1;i<=n1;++i)s+=a[i];
	printf("s=%I64d\n",s);*/
}
void calc_all_print(uint n){  //computes f(i) for all i<=n, and track the formulas.
	a=(uchar*)malloc(n+1); a[0]=0; a[1]=1; s[1]="1";
	for (uint i=2;i<=n;++i)a[i]=inf;
	uint g32[inf+1];
	for (uchar i=0;i<=inf;++i)g32[i]=min(calc_g(i),(u128)inf32);
	for (uint i=2;i<=n;++i){
		if (a[i-1]+1<a[i])a[i]=a[i-1]+1,s[i]=s[i-1]+"+1";
		uchar t=a[i-1],k=1;
		while (k<t/2&&g32[k+1]+g32[t-k-1]>=i)++k;
		for (uint j=6;j<=g32[k];++j)
			if (a[j]+a[i-j]<a[i])a[i]=a[j]+a[i-j],s[i]=s[j]+"+"+s[i-j];
		for (uint j=2,ij=i*2;j<=i&&ij<=n;++j,ij+=i)
			if (a[i]+a[j]<a[ij])a[ij]=a[i]+a[j],s[ij]="("+s[i]+")*("+s[j]+")";
	}
}

template<class T>
const vector<pair<T,uchar>> prime_factors(T x){  //factorize 64 bits, by (mostly) pollard-rho.
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
const vector<pair<T,uchar>> prime_factors128_old(T x){  //factorize 128 bits, by C-Quadratic-Sieve.
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
const vector<pair<T,uchar>> prime_factors128(T x){  //factorize 128 bits, by C-RHO-128.
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
vector<T> get_factors_from_primes(const vector<pair<T,uchar>> &p){  //get all factors from prime factorization.
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
vector<T> get_factors(T x){  //get all factors of x.
	if (sizeof(T)==sizeof(ull))return get_factors_from_primes(prime_factors(x));
	else return get_factors_from_primes(prime_factors128(x));
}
void clear_hash(){  //clear the hash tables to save memory.
	H.clear(); H1.clear(); P.clear(); P1.clear(); M_primes.clear();
}
ushort dfs128(u128 x,ushort t);
//similar to dfs128, decides whether f(x)<=t, for x with 62 bits.
//If true, return the optimal f(x). Otherwise return an upper bound on f(x).
ushort dfs(ull x,ushort t){
	if (x<=n0)return a[x];
	auto it=H.find(x);
	if (it!=H.end()&&it->second.first>=t)return it->second.second;
	++CNT;
	if (clock()-T0>60000){
		printf("time=%d CNT=%I64d\n",clock()-T0,CNT);
		T0=clock();
		/*if (H.size()>60000000)H.clear();
		if (H1.size()>25000000)H1.clear();
		if (P.size()>20000000)P.clear();
		if (P1.size()>5000000)P1.clear();
		if (M_primes.size()>10000000)M_primes.clear();*/
	}
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
			ushort v1=inf1;
			for (int i=lb1;i<=t-lb2&&i<v1;++i){
				v1=min(v1,dfs(g,i));
			}
			if (v1+lb2>t)continue;
			ushort v2=dfs(x/g,t-v1);
			ans=min(ans,ushort(v1+v2));
		}
	/*//subtraction
	u128 d=max_sub(x,t);
	for (int i=1;i<=d;++i){
		ushort lb=complexity_LB(x+i);
		if (a[i]+lb<=t)ans=min(ans,ushort(dfs128(x+i,t-a[i])+a[i]));
	}*/
	H[x]=make_pair(t,ans);
	return ans;
}
//decides whether f(x)<=t, for x with 128 bits.
//If true, return the optimal f(x). Otherwise return an upper bound on f(x).
ushort dfs128(u128 x,ushort t){
	if (x<=n0)return a[x];
	if (x<U)return dfs(x,t);
	auto it=H1.find(x);
	if (it!=H1.end()&&it->second.first>=t)return it->second.second;
	ushort ans=inf1,k=1;
	//addition
	while (k<t/2&&g[k+1]+g[t-k-1]>=x)++k;
	assert(g[k]<=n0);
	for (int i=1;i<=g[k];++i){
		ushort lb=complexity_LB(x-i);
		if (a[i]+lb<=t)ans=min(ans,ushort(dfs128(x-i,t-a[i])+a[i]));
	}
	//multiplication
	auto fac=get_factors(x);
	for (auto &g:fac)
		if (g>1&&g<=x/g){
			ushort lb1=complexity_LB(g),lb2=complexity_LB(x/g);
			if (lb1+lb2>t)continue;
			/*if (t-complexity_LB(x)>=3){
				if (dfs128(x/g,lb2)>lb2)++lb2;
				if (lb1+lb2>t)continue;
			}*/
			//search within the smaller part first.
			//ushort v1=dfs128(g,t-lb2);
			ushort v1=inf1;
			for (int i=lb1;i<=t-lb2&&i<v1;++i){
				v1=min(v1,dfs128(g,i));
			}
			if (v1+lb2>t)continue;
			ushort v2=dfs128(x/g,t-v1);
			ans=min(ans,ushort(v1+v2));
		}
	/*//subtraction
	u128 d=max_sub(x,t);
	for (int i=1;i<=d;++i){
		ushort lb=complexity_LB(x+i);
		if (a[i]+lb<=t)ans=min(ans,ushort(dfs128(x+i,t-a[i])+a[i]));
	}*/
	H1[x]=make_pair(t,ans);
	return ans;
}
ushort dfs128_lazy(u128 x,ushort t){  //without performing new computation, only use the stored hash tables.
	if (x<=n0)return a[x];
	if (x<U){
		auto it=H.find(x);
		return it!=H.end()?it->second.second:inf1;
	}
	else {
		auto it=H1.find(x);
		return it!=H1.end()?it->second.second:inf1;
	}
}
ushort calc_single(u128 n){  //computes f(n).
	int lb=complexity_LB(n);
	int ans=inf1;
	for (int t=lb;t<ans;++t){  //iteratively increase the depth.
		int res=dfs128(n,t);
		ans=min(ans,res);
	}
	return ans;
}

void dfs128_print(u128 x,ushort t){  //print the formula for f(x).
	if (x<=n0){
		printf("(");
		cout<<s[x];
		//print(x);
		printf(")");
		return;
	}
	//if (x<U)return dfs(x,t);
	ushort ans=inf1,k=1;
	//addition
	while (k<t/2&&g[k+1]+g[t-k-1]>=x)++k;
	assert(g[k]<=n0);
	for (int i=1;i<=g[k];++i){
		ushort lb=complexity_LB(x-i);
		if (a[i]+lb<=t){
			auto v=dfs128_lazy(x-i,t-a[i]);
			ans=min(ans,ushort(v+a[i]));
			if (v+a[i]==t){
				printf("(");
				dfs128_print(x-i,v);
				printf("+%s)",s[i].c_str());
				return;
			}
		}
	}
	//multiplication
	auto fac=get_factors(x);
	for (auto &g:fac)
		if (g>1&&g<=x/g){
			ushort lb1=complexity_LB(g),lb2=complexity_LB(x/g);
			if (lb1+lb2>t)continue;
			ushort v1=inf1;
			for (int i=lb1;i<=t-lb2&&i<v1;++i){
				v1=min(v1,dfs128_lazy(g,i));
			}
			if (v1+lb2>t)continue;
			ushort v2=dfs128_lazy(x/g,t-v1);
			ans=min(ans,ushort(v1+v2));
			if (v1+v2==t){
				printf("(");
				dfs128_print(x/g,v2);
				printf("*");
				dfs128_print(g,v1);
				printf(")");
				return;
			}
		}
}
void print_expr(u128 x,ushort t=0){  //t: f(x). d: dfs depth.
	if (!t)t=calc_single(x);
	ushort d=complexity_LB(x);
	for (;dfs128(x,d)>t;)++d;
	printf("%d %d\n",d,t);
	print(x); cout<<" "<<t<<":"<<endl;
	dfs128_print(x,t);
	cout<<endl;
}

void init(uint _n0=1e9,u128 _N0=1e23,bool print=0){
	int t1=clock();
	n0=_n0; N0=_N0;
	if (!print){
		//calc_all_subtraction(n0);
		calc_all(n0);
	}
	else calc_all_print(n0);
	printf("init time=%d\n",clock()-t1);
}
void clear(){
	H.clear(); H1.clear();
	//unordered_map<ull,pair<ushort,ushort>>().swap(H);
	//unordered_map<u128,pair<ushort,ushort>>().swap(H1);
}
void check1(){  // test the conjecture f(p^i)=i*f(p). (In particular, f(2^i)=2i.)
	int t1=clock();
	//vector<int> primes={733,379,739,541};  // conjecture fails
	//vector<int> primes={577,811};  // also fails
	//vector<int> primes={109,433,163,487,2};
	//vector<int> primes={2};
	vector<int> primes={811};
	for (auto mul:primes){
		printf("mul=%I64d\n",mul);
		u128 x=mul;
		for (int i=1;;++i,x*=mul){
			//if (i<80)continue;
			u128 n=x;
			//printf("i=%d n=%I64d\n",i,n);
			int v=i*a[mul],lb=complexity_LB(n);  // v: the conjectured complexity.
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
				println(n);
				printf("ans=%d tgt=%d\n",ans,v);
				//exit(0);
				break;
			}
			printf("time=%d\n",clock()-t1); t1=clock();
			//printf("hash size: H=%d H1=%d M_lb=%d P=%d P1=%d M_primes=%d\n",H.size(),H1.size(),M_lb.size(),P.size(),P1.size(),M_primes.size());
			if (n>N0/mul)break;
		}
	}
	printf("time=%d\n",clock()-t1);
}
void check2(){  // test the conjecture f(2^i3^j5^k)=2i+3j+5k (k<=5).
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
				u128 n=x*y*z;
				if (n>N0)break;
				if (n==1)continue;
				//if (n<1e35)continue;
				//if (i<60||i==60&&j<3||i==60&&j==3&&k<5)continue;
				int v=i*2+j*3+k*5,lb=complexity_LB(n);
				int ans=v;
				for (int t=lb;t<v;++t){
					int res=dfs128(n,t);
					ans=min(ans,res);
					if (ans<v){
						printf("improve: %d %d\n",ans,v);
						printf("n="); println(n);
						for (;;);
						//break;
					}
				}
				printf("i=%d j=%d k=%d ans=%d\n",i,j,k,ans);
				printf("n="); println(n);
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
			//if (i<=76)continue;
			u128 n=x+1;
			//printf("i=%d n=%I64d\n",i,n);
			int v=i*2+1,lb=complexity_LB(n);
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
void check4(){  // test the conjecture f(3p)=min{f(3p-1)+1,f(p)+3} for prime p, by sampling a few integers.
	int t1=clock();
	u128 U=1e22;
	for (u128 p;;){
		p=rand128()%U+1; //++p;
		if (!is_prime128(p))continue;
		printf("p="); println(p);
		auto v1=calc_single(p*3-1),v2=calc_single(p),v3=calc_single(p*3);
		//println(v1); println(v2); println(v3);
		if (min(v1+1,v2+3)!=v3){
			printf("err: p=",p); println(p);
			printf("3p="); println(p*3);
			printf(" f(3p-1)=%d f(p)=%d f(3p)=%d\n",v1,v2,v3);
			for (;;); //exit(0); //break;
		}
		clear_hash();
	}
	printf("time=%d\n",clock()-t1);
}
void check5(){  // test the conjecture f(3^k*n)=f(n)+3k for n=73(3^21+1)+6=763605783898.
	int t1=clock(),t0=t1;
	u128 x=763605783898ull;
	int fn=calc_single(x);  //fn=82
	for (int i=0;;++i,x*=3){
		u128 n=x;
		int v=fn+i*3,lb=complexity_LB(n);  // v: the conjectured complexity.
		int ans=v;
		for (int t=lb;t<v;++t){
			int res=dfs128(n,t);
			ans=min(ans,res);
			if (ans<v){
				printf("improve: %d %d\n",ans,v);
				break;
			}
		}
		//int ans=calc_single(n);
		printf("i=%d ans=%d\n",i,ans);
		//printf("time=%d\n",clock()-t1);;
		if (ans<v){
			printf("err: n=");
			println(n);
			printf("ans=%d tgt=%d\n",ans,v);
			//exit(0);
			break;
		}
		printf("time=%d\n",clock()-t1); t1=clock();
		if (n>N0/3)break;
	}
	printf("total time=%d\n",clock()-t0);
}
double run_sample(int num_samples=1e4,bool debug=1){  //computes the average integer complexity, by sampling.
	//freopen("data.txt","w",stdout);
	int t1=clock(),t2=t1;
	vector<double> a;
	clear();
	for (int i=0;i<num_samples;++i){
		//if (i%100==0)printf("i=%d\n",i);
		if (clock()-t2>6e4){  //track the speed
			if (debug)printf("i=%d\n",i);
			t2=clock();
		}
		//u128 n=N0-i;
		u128 n=rand128()%(N0-1)+2;
		ushort ans=calc_single(n);
		double x=ans/log(n)*log(3);
		//printf("%d, ",i); print(n);
		//printf(", %d, %.6lf, %.6lf\n",ans,x,1-1.*n/calc_g(ans));
		a.push_back(x);
		clear_hash();
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
	return ave;
}
void verify(int T=1e4){  //verify the correctness, by checking f(n) for a few samples.
	u128 U=1e8; init(U,1e23); n0=1e4;
	int t1=clock();
	for (int i=1;i<=T;++i){
		u128 n=rand128()%(U-1)+2;
		assert(calc_single(n)==a[n]);
	}
	printf("verify time=%d\n",clock()-t1);
}
void factorize_test(u128 N=1e18,int T=1000){  //verify the correctness for different factorization algorithms.
	int t1=clock();
	ull s=0;
	for (int i1=1;i1<=T;++i1){
		u128 n=rand128()%N+1;
		//println(n);
		auto a1=prime_factors128(n);
		s+=a1[0].second;
		/*auto a2=prime_factors128_old(n);
		if (a1!=a2){
			prln(a1);
			prln(a2);
			exit(0);
		}*/
	}
	println(s);
	printf("fac time=%d\n",clock()-t1);
}
void test(){
	//prime_factors(123456789012345678);
	
	/*u128 n;
	n=162879576091729561ull;
	n=155104303499468569ull;
	printf("%d\n",defect(n,119));
	cout<<(int)dfs(n,2)<<endl;
	printf("time=%d\n",clock()-t1);
	exit(0);*/
}
int main()
{
	srand(time(0));
	//srand(1);
	//u128 N0=pow128(10,23);
	
	//verify();
	//factorize_test(1e30,1e2);
	//calc_all_subtraction(1e7);
	//return 0;
	
	//init(1e6,1e18);
	//init(1e7,1e30);
	//init(1e7,1e38,1);
	init(1e8,1e38);
	//init(1e9,1e38);
	//init(2e9,1e38);
	
	//print_expr(u128_from_str("1234567890"));
	//print_expr(u128_from_str("151770612880318395249730891"),179);
	//print_expr(u128_from_str("1361788799550131972374553991985921"),227);
	
	check1();
	//check2();
	//check3();
	//check4();
	//check5();
	//run_sample(1e4);
	
	/*vector<double> a;
	for (int T=1;;++T){
		N0=pow128(10,18);
		int num_samples=10;
		double res=run_sample(num_samples,1);
		a.push_back(res);
		printf("--------T=%d #samples=%d %.6lf--------\n",T,T*num_samples,mean(a));
		//return 0;
	}*/
	/*for (u128 i=1;i<=1e30;i*=10){
		if (i<=1e10)continue;
		N0=i;
		//run_sample(1e6);
		run_sample(1e5);
	}*/
    return 0;
}

