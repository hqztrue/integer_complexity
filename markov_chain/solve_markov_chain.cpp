#include<bits/stdc++.h>
using namespace std;
typedef long long ll;
typedef unsigned int uint;
vector<int> mod={2,3,5,6,7,11,13,14,55};
const int M=55;
//const int B=2310;  //2*3*5*7*11
//const int B=997920;  //2**5*3**4*5*7*11
//const int B=1081080;  //2**3*3**3*5*7*11*13
//const int B=299376000;  //2**7*3**5*5**3*7*11
//const int B=1945944000;  //2**6*3**5*5**3*7*11*13
//const int B=1077753600;  //2**8*3**7*5**2*7*11
const int B=1167566400;  //2**6*3**6*5**2*7*11*13
int *a,*cnt;
float *v,*v1,l[M+1],c=3.4;
ll sizei,sized;
int f[]={0,1,2,3,4,5,5,6,6,6,7,8,7,8,8,8,8,9,8,9,9,9,10,11,
	9,10,10,9,10,11,10,11,10,11,11,11,10,11,11,11,11,
	12,11,12,12,11,12,13,11,12,12,12,12,13,11,12,12,
	12,13,14,12,13,13,12,12,13,13,14,13,14,13,14,12,
	13,13,13,13,14,13,14};
void prepare(){
	sizei=sizeof(int)*B; sized=sizeof(float)*B;
	a=(int*)malloc(sizei);
	cnt=(int*)malloc(sizei);
	v=(float*)calloc(1,sized);
	v1=(float*)malloc(sized);
	for (int i=2;i<=M;++i)l[i]=log(i)/log(3);
	assert(a&&cnt&&v&&v1);
	/*mod.clear();
	for (int i=2;i<=M;++i)
		if (B%i==0)mod.push_back(i);*/
	auto mod1=mod;
	mod.clear();
	for (int p:mod1)
		if(B%p==0)mod.push_back(p);
}
void inita(){
	for (int i=0;i<B;++i)a[i]=2; //mod[rand()%mod.size()];
}
double compute(int i,int p){
	double v1=0;
	for (uint j=i/p,add=B/p;j<B;j+=add){
		v1+=v[j];
	}
	v1=v1/p+f[p]+f[i%p]-c*l[p];
	return v1;
}
void iteration(int T=1){
	while (T--){
		for (int i=0;i<B;++i){
			v1[i]=compute(i,a[i]);
		}
		memcpy(v,v1,sized);
		double c1=0;
		for (int i=0;i<B;++i)c1+=v[i];
		c1/=B;
		c+=c1;
		for (int i=0;i<B;++i)v[i]-=c1;
	}
}
void modify(){
	for (int i=0;i<B;++i){
		int pbest=0;
		double vbest=1e9;
		for (int p:mod){
			double v0=compute(i,p);
			if (v0<vbest)vbest=v0,pbest=p;
		}
		a[i]=pbest;
	}
}
void print_mod(){
	for (int i=0;i<B;++i)cnt[i]=0;
	for (int i=0;i<B;++i)++cnt[a[i]];
	puts("mod:");
	for (int i=0;i<B;++i)
		if (cnt[i])printf("%d %d\n", i, cnt[i]);
}
int main()
{
	//freopen("1.in","r",stdin);
	//freopen("1.out","w",stdout);
	int T=100000000;
	prepare();
	inita();
	printf("B=%d\n", B);
	for (int i1=1;i1<=T;++i1){
		if (i1%1==0){
			printf("i=%d c=%.10lf\n", i1, c);
			//for (int i=0;i<B;++i)printf("%.6lf%c", v[i], i+1==B?'\n':' ');
		}
		iteration(1);
		modify();
		print_mod();
	}
	//system("pause");for (;;);
	return 0;
}
