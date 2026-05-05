#include<bits/stdc++.h>
using namespace std;
const int N=10005;
int a[N][N],f[N][N];
int main()
{
	//freopen("1.in","r",stdin);
	//freopen("1.out","w",stdout);
	srand(1);
	int n=10000;
	for (int i=0;i<n;++i)
		for (int j=0;j<n;++j)a[i][j]=rand()%2;
	memset(f,0x5,sizeof(f));
	f[0][0]=0;
	for (int i=0;i<n;++i)
		for (int j=0;j<n;++j){
			if (i)f[i][j]=min(f[i][j],f[i-1][j]+a[i][j]);
			if (j)f[i][j]=min(f[i][j],f[i][j-1]+a[i][j]);
		}
	printf("%d\n", f[n-1][n-1]);
	//system("pause");for (;;);
	return 0;
}
