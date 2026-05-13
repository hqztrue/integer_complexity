namespace test_heuristics{
const int LEN_MAX=100000;
char x_str[LEN_MAX+2];
void test_heuristic_conjectures(){
	const int LEN_MIN=10;
	for (int len=LEN_MIN;len<LEN_MAX;++len){
		printf("len=%d\n", len);
		
		//2^n
		x_str[0]='1';
		for (int i=1;i<=len;++i)x_str[i]='0';
		x_str[len+1]=0;
		
		//2^n+1
		/*x_str[0]='1';
		for (int i=1;i<len;++i)x_str[i]='0';
		x_str[len]='1';
		x_str[len+1]=0;*/
		
		//2^n-1
		/*for (int i=0;i<len;++i)x_str[i]='1';
		x_str[len]=0;*/
		
		Int n; Int::N=10000;
		for (auto x:{x_str}){
			n.read(x);
			int ans=calc_all(n);
			//printf("len %d %.6lf\n", n.bit_length(), n.log2());
			printf("%d %.6lf\n",ans,ans/(n.log2()/log2(3)));
			if (ans<2*len){  //<=2*len
				printf("counter example: %d\n",len);
				for (;;);
			}
		}
	}
}
} using namespace test_heuristics;
