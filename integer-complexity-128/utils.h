#ifndef __UTILS__
#define __UTILS__
typedef unsigned char uchar;
typedef unsigned short ushort;
typedef unsigned int uint;
typedef long long ll;
typedef unsigned long long ull;
typedef unsigned __int128 u128;

// provided for convenience, take a string and return a 128-bit unsigned integer.
__uint128_t from_string_128_bits(const char *str) {
    __uint128_t res = 0;
    for (; *str; res = res * 10 + *str++ - '0');
    return res;
}

// provided to print 128-bit unsigned integers.
static char *to_string_128_bits(__uint128_t num) {
    static char s[40];
    __uint128_t mask = -1;
    size_t a, b, c = 1, d;
    strcpy(s, "0");
    for (mask -= mask / 2; mask; mask >>= 1) {
        for (a = (num & mask) != 0, b = c; b;) {
            d = ((s[--b] - '0') << 1) + a;
            s[b] = "0123456789"[d % 10];
            a = d / 10;
        }
        for (; a; ++c, memmove(s + 1, s, c), *s = "0123456789"[a % 10], a /= 10);
    }
    return s;
}

void print(u128 x){
	if (x<0)putchar('-'),x=-x;
	if (x>9)print(x/10);
	putchar(x%10+'0');
}
void println(u128 x){print(x); putchar('\n');}
u128 u128_from_str(const char *s){u128 x=0; while (*s)x=x*10+*s++-'0'; return x;}

//inline u128 _rand128(){static u128 x=u128_from_str("199609092119960909211996090921996090921");x+=(x<<17)+(x>>29)+1;return x;}
inline uint rand32(){  //32 bits
#ifdef _WIN32
	//assert(RAND_MAX==32767);
	return ((uint)rand()<<30)+(rand()<<15)+rand();
#else
	//assert(RAND_MAX==2147483647);
	return ((uint)rand()<<31)+rand();
#endif
}
inline ull rand64(){  //64 bits
#ifdef _WIN32
	//assert(RAND_MAX==32767);
	return ((ull)rand()<<60)+((ull)rand()<<45)+((ull)rand()<<30)+((ull)rand()<<15)+rand();
#else
	//assert(RAND_MAX==2147483647);
	return ((ull)rand()<<62)+((ull)rand()<<31)+rand();
#endif
}
inline u128 rand128(){  //128 bits
	return ((u128)rand64()<<64)+rand64();
}

//print functions
string to_string(char c){return string(1,c);}
string to_string(const char *s){return (string)s;}
string to_string(const string &s){return s;}
string to_string(u128 x){return to_string(to_string_128_bits(x));}
template<class T> string to_string(const complex<T> &c){stringstream ss; ss<<c; return ss.str();}
string to_string(const vector<bool> &v){string res="{"; for (int i=0;i<v.size();++i)res+=char('0'+v[i]); return res+="}";}
template<size_t sz> string to_string(const bitset<sz> &b){string res=""; for (int i=0;i<sz;++i)res+=char('0'+b[i]); return res;}
template<class T,class U> string to_string(const pair<T,U> &p){return "("+to_string(p.first)+", "+to_string(p.second)+")";}
#if __cplusplus>199711L
string to_string(bool b){
	//return b?"true":"false";
	return to_string((int)b);
}
template<class T> string to_string(T v){  //containers with begin(),end()
	bool fst=1; string res="{";
	for (const auto &x:v){
		if (!fst)res+=", ";
		fst=0; res+=to_string(x);
	}
	return res+="}";
}
void _DBG(){cerr<<"]"<<endl;}
template<class T,class ...U> void _DBG(const T &t,const U&... u){
	cerr<<to_string(t); if (sizeof...(u))cerr<<", "; _DBG(u...);
}
#define dbg(...) cerr<<__FUNCTION__<<"() @ "<<"L"<<__LINE__<<" -> ["<<#__VA_ARGS__<<"]: [",_DBG(__VA_ARGS__)
#endif
template<class T> void pr(const T &x){cout<<to_string(x);}
#if __cplusplus>199711L
template<class T,class ...U> void pr(const T &t,const U&... u){pr(t); pr(u...);}
void prln(){pr("\n");}
template<class T,class ...U> void prln(const T &t,const U&... u){
	pr(t); if (sizeof...(u))pr(" "); prln(u...);
}
#endif

u128 pow128(u128 x,int y){
	u128 ans=1;
	for (int i=1;i<=y;++i)ans*=x;
	return ans;
}

#endif //__UTILS__

