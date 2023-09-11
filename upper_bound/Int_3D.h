#include "../C-Quadratic-Sieve_C++/utils.h"

struct Int{
	const static int N0=300015;
	const uint mask=0x55555555u;
	const static int rem[];
	static int N;
	uint a[N0/32+1];
	int mod2,mod3,mod5,sum[4],cur_w,cur_b;
	Int& operator =(const Int &A){
		memcpy(a,A.a,sizeof(uint)*N);
		return *this;
	}
	Int& operator /=(int x){
		ull d=0;
		for (int i=N-1;i>=0;--i)
			d=(d<<32)+a[i],a[i]=d/x,d%=x;
		return *this;
	}
	int at(int x)const{return (a[x/32]>>x%32)&1;}
	void init_div2(){
		mod2=a[0]&1;
		cur_w=cur_b=0;
		mod3=0;
		for (int i=0;i<N;++i){
			//binary bits at even locations contribute 1 to mod 3,
			//and odd locations contribute -1.
			mod3+=__builtin_popcount(a[i]&mask);
			mod3-=__builtin_popcount((a[i]>>1)&mask);
		}
		mod3=(mod3%3+3)%3;
		
		for (int i=0;i<4;++i)sum[i]=0;
		for (int i=0;i<N*32;++i)sum[i%4]+=at(i);
		mod5=0;
		for (int i=0;i<4;++i)mod5+=rem[i]*sum[i];
		mod5%=5;
	}
	void div2(){  //update mod2, mod3 and mod5 after dividing by 2.
		sum[0]-=mod2;
		int t=sum[0]; for (int i=0;i<3;++i)sum[i]=sum[i+1]; sum[3]=t;
		int m=(mod5-mod2+5)*3%5;
		mod5=0;
		for (int i=0;i<4;++i)mod5+=rem[i]*sum[i];
		mod5%=5;
		assert(mod5==m);
		
		mod3=(3-mod3+mod2)%3;
		if (++cur_b==32)cur_b=0,++cur_w;
		mod2=(a[cur_w]>>cur_b)&1;
	}
	int bit_length(){
		for (int i=N-1;i>=0;--i)
			if (a[i])return (i+1)*32-__builtin_clz(a[i]);
		return 1;
	}
	void clear(){memset(a,0,sizeof(uint)*N);}
	void read(const char *s){
		int l=strlen(s); clear();
		for (int i=0;i<l;++i)a[i/32]|=(s[l-1-i]-'0')<<i%32;
	}
	// set N0 to be slightly larger than log_2(base), so it will only introduce
	// small error, even if we don't restrict the rand value in [0,base).
	void init_rand(){
		for (int i=0;i<N;++i)a[i]=rand32();
	}
};
const int Int::rem[]={1,2,4,3};
int Int::N;

