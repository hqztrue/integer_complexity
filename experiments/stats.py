import math,scipy
from scipy.stats import norm

a=1e-5

n=20000; sigma=0.02
#n=2000; sigma=0.017
#n-=1

d=norm.ppf(1-a/2)
#d=5
print('d=',d)
print('confidence interval=',2*d*sigma/math.sqrt(n))
