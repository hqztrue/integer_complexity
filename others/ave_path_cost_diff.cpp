#include<bits/stdc++.h>
using namespace std;
const int N=10005;
int a[N][N],f[N][N],pre[N][N],n;
int ub[N],lb[N],diffl,diffr;
mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());
//inline void upd(int &x,int y,int &pre,int v){if (y<x||y==x&&rng()%2)x=y,pre=v;}
inline void upd(int &x,int y,int &pre,int v){if (y<x)x=y,pre=v;}
int run(){
	for (int i=0;i<n;++i)
		for (int j=0;j<n;++j)a[i][j]=rng()%2;  //rand()
	memset(f,0x5,sizeof(f));
	f[0][0]=0;
	for (int i=0;i<n;++i)
		for (int j=0;j<n;++j){
			if (i)upd(f[i][j],f[i-1][j]+a[i][j],pre[i][j],0);
			if (j)upd(f[i][j],f[i][j-1]+a[i][j],pre[i][j],1);
		}
	return f[n-1][n-1];
}
void recover_path(){
	int v1=n-1,v2=n-1;
	while (v1||v2){
		lb[v1]=min(lb[v1],v2);
		ub[v1]=max(ub[v1],v2);
		if (pre[v1][v2]==0)--v1;
		else --v2;
	}
	for (int i=0;i<n;++i)
		diffl=max(diffl,i-lb[i]),diffr=max(diffr,ub[i]-i);
	//for (int i=0;i<n;++i)printf("%d %d %d\n",i,lb[i],ub[i]);
	printf("diff: %d %d\n",diffl,diffr);
}
int main()
{
	//freopen("1.in","r",stdin);
	//freopen("1.out","w",stdout);
	srand(time(0));
	//n=1000;
	n=5000;
	for (int i=0;i<n;++i)ub[i]=0,lb[i]=n;
	int T=1000;
	for (int i1=1;i1<=T;++i1){
		int ans=run();
		recover_path();
		printf("ans=%d\n",ans);
	}
	//system("pause");for (;;);
	return 0;
}
