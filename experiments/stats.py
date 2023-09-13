import math,scipy
from scipy.stats import norm

def conf_interval(a,n,sigma):
	#d=5
	d=norm.ppf(1-a/2)
	#print('d=',d)
	return 2*d*sigma/math.sqrt(n)

a=1e-5

print('1e24',conf_interval(a,20000,0.02))
print('1e25',conf_interval(a,2000,0.017))
print('D(b,r) (2,3)',conf_interval(a,10**5,0.0011))
print('D(b,r) (2,3,5)',conf_interval(a,10**5,0.0026))
