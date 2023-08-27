import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

#x=[12,14,16,18,20,22]
#y=[3.350673,3.342282,3.335129,3.329165,3.324267,3.320562]

x=[4,5,6,7,8,9,10,11,12,13,14,16,18,20,22]
y=[3.400887,3.395414,3.388457,3.380310,3.373065,3.365958,3.360919,3.355967,
   3.350673,3.345934,3.342282,3.335129,3.329165,3.324267,3.320562]

x=x[:]
y=y[:]

x=np.array(x); y=np.array(y)

def f1(x,a,b):
	return a/x+b

def f2(x,a,b,c):
	return a/pow(x,b)+c

def f3(x,a,b,c):
	return a/(x+b)+c

def f4(x,a,b,c,d):
	return a/(x*x+b*x+c)+d

def f5(x,a,b,c):
	return a/(x*x+b)+c

f=f5

p_est,err_est=curve_fit(f,x,y,maxfev=100000)
print('params=',p_est)
#print('err=',err_est)
p = plt.plot(x, y, "rx")
p = plt.plot(x, f(x, *p_est), "k--")
plt.show()

