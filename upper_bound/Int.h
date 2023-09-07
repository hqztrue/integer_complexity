#include "../C-Quadratic-Sieve_C++/utils.h"
struct Int{
	const static int N0=100005,N=N0/32+1;
	const uint mask=0x55555555u;
	uint a[N];
	int mod2,mod3,cur_w,cur_b;
	Int& operator /=(int x){
		ull d=0;
		for (int i=N-1;i>=0;--i)
			d=(d<<32)+a[i],a[i]=d/x,d%=x;
		return *this;
	}
	void init_div2(){
		mod2=a[0]&1;
		cur_w=cur_b=0;
		mod3=0;
		for (int i=0;i<N;++i){
			mod3+=__builtin_popcount(a[i]&mask);
			mod3-=__builtin_popcount((a[i]>>1)&mask);
		}
		mod3=(mod3%3+3)%3;
	}
	void div2(){
		if (++cur_b==32)cur_b=0,++cur_w;
		mod3=(3-mod3+mod2)%3;
		mod2=(a[cur_w]>>cur_b)&1;
	}
	// set N0 to be slightly larger than log_2(base), so will only introduce
	// small error, even if we don't restrict the rand value in [0,base).
	void init_rand(){
		for (int i=0;i<N;++i)a[i]=rand32();
	}
};

