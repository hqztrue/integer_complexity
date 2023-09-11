import numpy as np
import scipy, math
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

#err
#x=[4,5,6,7,8,9,10,11,12,13,14,16,18,20,22]
#y=[3.400887,3.395414,3.388457,3.380310,3.373065,3.365958,3.360919,3.355967,
#   3.350673,3.345934,3.342282,3.335129,3.329165,3.324267,3.320562]

#x=list(range(1,23))
data=[(1,3.238373),(2,3.349894),(3,3.393001),(4,3.400376),(5,3.395626),
	  (6,3.388161),(7,3.380561),(8,3.373263),(9,3.366459),(10,3.360493),
	  (11,3.355246),(12,3.350321),(13,3.346080),(14,3.342084),(15,3.338514),
	  (16,3.335244),(17,3.332222),(18,3.329551),(19,3.327147),(20,3.324807),
	  (21,3.322698),(22,3.320632),(23,3.318700),(24,3.316791),(25,3.315385)]
data=data[4:]

x,y=zip(*data)
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

def f6(x,a,b,c):
	return a/(x**1.5+b)+c

f=f3

popt,pcov=curve_fit(f,x,y,maxfev=100000)
print('params=',popt)
#perr = np.sqrt(np.diag(pcov))
perr=f(x,*popt)-y
mse=math.sqrt(sum(x**2 for x in perr)/len(perr))
print('err=',mse)
p=plt.plot(x,y,"rx")
p=plt.plot(x,f(x,*popt),"k--")
plt.title('average integer complexity')
plt.xlabel("log_10 n")
plt.ylabel("average f_log")
plt.show()

