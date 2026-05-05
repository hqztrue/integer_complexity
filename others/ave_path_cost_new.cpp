#include<bits/stdc++.h>
using namespace std;
const int N=12005;
int a[N][N],b[N][N],f[N][N];
int main()
{
	//freopen("1.in","r",stdin);
	//freopen("1.out","w",stdout);
	//srand(1);
	srand(time(0));
	//int n=1000,m=900; //2^n*3^m
	int n=11500,m=10000;
	for (int i=0;i<n;++i)
		for (int j=0;j<m;++j){
			a[i][j]=rand()%2;
			b[i][j]=rand()%3;
		}
	memset(f,0x5,sizeof(f));
	f[0][0]=0;
	for (int i=0;i<n;++i)
		for (int j=0;j<m;++j){
			if (i)f[i][j]=min(f[i][j],f[i-1][j]+a[i][j]);
			if (j)f[i][j]=min(f[i][j],f[i][j-1]+b[i][j]);
		}
	int ans=f[n-1][m-1];
	printf("%d\n", ans);
	double v=1.*(2*n+3*m+ans)/(log(2)/log(3)*n+m);
	printf("%.10lf\n", v);
	//system("pause");for (;;);
	return 0;
}
