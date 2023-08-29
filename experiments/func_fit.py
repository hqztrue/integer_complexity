import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

#err
#x=[4,5,6,7,8,9,10,11,12,13,14,16,18,20,22]
#y=[3.400887,3.395414,3.388457,3.380310,3.373065,3.365958,3.360919,3.355967,
#   3.350673,3.345934,3.342282,3.335129,3.329165,3.324267,3.320562]

#x=list(range(1,23))
x=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,20,22]

y=[3.238373,3.349894,3.393001,3.400376,3.395626,3.388161,3.380561,3.373263,
   3.366459,3.360493,3.355246,3.350317,3.346067,3.342115,3.338113,3.335232,
   3.331408,3.329342]

y+=[3.324267,3.320562]

I=3
y=y[I:]
x=x[I:I+len(y)]

x=np.array(x); y=np.array(y)

def f1(x,a,b):
	return a/x+b

def f2(x,a,b,c):
	return a/(x+b)+c

def f3(x,a,b,c,d):
	return a/(pow(x,b)+c)+d

def f4(x,a,b,c,d):
	return a/(x*x+b*x+c)+d

def f5(x,a,b,c):
	return a/(x*x+b)+c

f=f5

p_est,err_est=curve_fit(f,x,y,maxfev=100000)
print('params=',p_est)
#print('err=',err_est)
p=plt.plot(x,y,"rx")
p=plt.plot(x,f(x,*p_est),"k--")
plt.title('average integer complexity')
plt.xlabel("log_10 n")
plt.ylabel("average f_log")
plt.show()

