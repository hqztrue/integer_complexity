def greedy(x): #Steinerberger
	ans=0
	while x>5:
		m=x%6
		if m==0 or m==3: x//=3; ans+=3
		elif m==1: x//=3; ans+=3+1
		elif m==2 or m==4: x//=2; ans+=2
		elif m==5: x//=2; ans+=2+1
	return ans+x

a=(3, 3, 2, 3, 2, 5, 3, 7, 2, 3, 2, 11, 2, 6, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2, 2, 2, 2, 3, 2, 2,
3, 3, 2, 3, 2, 7, 6, 6, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2, 3, 7, 2, 3, 2, 2, 3, 55, 2, 3, 2, 2,
3, 3, 2, 3, 2, 5, 3, 3, 2, 3, 2, 2, 2, 2, 2, 3, 2, 7, 3, 3, 2, 3, 2, 2, 6, 6, 2, 3, 2, 2,
3, 7, 2, 3, 2, 5, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 3, 3, 3, 3, 55, 3, 2, 2, 3, 5, 2, 3, 2, 7,
3, 2, 2, 3, 2, 5, 3, 3, 2, 3, 2, 2, 3, 7, 2, 3, 2, 2, 3, 3, 2, 3, 2, 11, 6, 6, 2, 3, 2, 2,
3, 3, 2, 3, 2, 5, 3, 3, 2, 3, 2, 7, 3, 3, 2, 3, 3, 3, 3, 3, 2, 3, 2, 2, 3, 7, 2, 3, 2, 2,
3, 3, 2, 3, 2, 5, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 7, 3, 3, 2, 3, 2, 11,
3, 3, 2, 3, 2, 5, 2, 7, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2, 2, 2, 2, 3, 2, 2, 3, 3, 2, 3, 7, 7,
3, 3, 2, 3, 2, 5, 3, 3, 2, 3, 2, 2, 6, 6, 2, 3, 3, 3, 3, 7, 2, 3, 2, 2, 3, 2, 2, 3, 2, 2,
3, 3, 2, 3, 2, 5, 3, 3, 2, 3, 2, 2, 2, 2, 2, 3, 2, 7, 2, 2, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2,
3, 3, 2, 3, 2, 5, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 3, 3, 3, 3, 2, 3, 7, 7, 3, 3, 2, 3, 2, 7,
3, 3, 2, 3, 2, 5, 3, 3, 2, 3, 2, 11, 3, 7, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2,
3, 3, 2, 3, 7, 5, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2, 2, 55, 55, 3, 2, 2,
3, 3, 2, 3, 2, 5, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 11, 3, 3, 2, 3, 2, 7, 3, 3, 2, 3, 2, 2,
3, 3, 2, 3, 2, 5, 2, 7, 2, 3, 2, 2, 6, 6, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2, 3, 5, 2, 3, 2, 2,
3, 3, 2, 3, 2, 7, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2, 3, 7, 2, 3, 2, 11, 3, 3, 2, 3, 2, 2,
3, 3, 2, 3, 2, 5, 3, 3, 2, 3, 2, 2, 3, 2, 2, 3, 3, 7, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2,
3, 3, 2, 3, 2, 5, 3, 2, 2, 3, 2, 2, 3, 3, 2, 3, 3, 3, 3, 3, 2, 3, 7, 7, 3, 3, 2, 3, 2, 7,
3, 3, 2, 3, 2, 5, 3, 3, 2, 3, 2, 2, 3, 7, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2,
3, 3, 2, 3, 2, 5, 3, 3, 2, 3, 2, 2, 3, 2, 2, 3, 3, 3, 3, 3, 2, 3, 2, 2, 3, 7, 2, 3, 2, 2,
3, 3, 2, 3, 2, 5, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 7, 3, 3, 2, 3, 2, 2,
3, 3, 2, 3, 2, 5, 3, 7, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2, 2, 2, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2,
3, 3, 2, 3, 2, 5, 3, 3, 2, 3, 2, 11, 3, 3, 2, 3, 2, 2, 3, 7, 2, 3, 2, 2, 2, 2, 2, 3, 2, 2,
3, 3, 2, 3, 3, 5, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 7, 3, 3, 2, 3, 2, 2, 3, 55, 2, 3, 2, 2,
3, 3, 2, 3, 2, 5, 3, 6, 2, 3, 2, 2, 3, 3, 2, 3, 3, 11, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 7,
3, 3, 2, 3, 2, 5, 6, 6, 2, 3, 2, 2, 3, 7, 2, 3, 2, 2, 3, 3, 2, 2, 2, 2, 3, 5, 2, 7, 2, 2,
3, 2, 2, 3, 2, 5, 3, 3, 2, 3, 2, 7, 3, 3, 2, 3, 3, 3, 3, 2, 2, 3, 2, 11, 3, 7, 2, 3, 2, 2,
3, 3, 2, 3, 2, 5, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 3, 3, 3, 3, 2, 3, 2, 7, 3, 3, 2, 3, 2, 2,
3, 3, 2, 3, 2, 5, 2, 7, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2, 2, 2, 2, 3, 2, 11,
3, 3, 2, 3, 2, 5, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2, 3, 7, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2,
3, 3, 2, 3, 2, 5, 3, 3, 2, 3, 2, 2, 2, 2, 2, 3, 2, 7, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2,
3, 7, 2, 3, 2, 5, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 7,
3, 3, 2, 3, 2, 5, 3, 3, 2, 3, 2, 2, 3, 7, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2,
3, 3, 2, 3, 2, 5, 3, 3, 2, 3, 2, 11, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2, 3, 7, 2, 3, 2, 2,
3, 3, 2, 3, 2, 5, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 7, 3, 55, 2, 3, 2, 2,
3, 3, 2, 3, 2, 5, 3, 7, 2, 3, 2, 2, 3, 3, 2, 3, 2, 11, 3, 3, 2, 3, 2, 2, 2, 2, 2, 3, 2, 2,
3, 3, 2, 3, 2, 5, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2, 3, 7, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2,
3, 3, 2, 3, 2, 5, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 7, 2, 2, 2, 3, 2, 11, 3, 3, 2, 3, 2, 2,
3, 2, 2, 3, 2, 5, 2, 2, 2, 3, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 5, 3, 7, 7, 3, 3, 2, 3, 2, 7,
3, 3, 2, 3, 2, 5, 3, 3, 2, 3, 2, 2, 3, 7, 2, 3, 3, 3, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 11,
3, 3, 2, 3, 7, 5, 3, 3, 2, 3, 55, 55, 2, 2, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2, 3, 7, 2, 3, 2, 2,
3, 3, 2, 3, 2, 5, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 7, 3, 2, 2, 3, 2, 2,
3, 3, 2, 3, 2, 5, 3, 7, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2, 3, 2, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2,
3, 3, 2, 3, 2, 5, 6, 6, 2, 3, 2, 2, 3, 3, 2, 3, 3, 3, 3, 7, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2,
3, 3, 2, 3, 2, 5, 2, 2, 2, 3, 2, 11, 3, 3, 2, 3, 2, 7, 2, 2, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2,
3, 2, 2, 3, 2, 5, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 3, 3, 3, 3, 2, 3, 7, 7, 3, 55, 2, 3, 2, 7,
3, 3, 2, 3, 2, 5, 3, 3, 2, 3, 2, 2, 3, 7, 2, 3, 2, 11, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2,
3, 3, 2, 3, 7, 5, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 3, 3, 3, 3, 55, 3, 2, 2, 3, 7, 2, 3, 2, 2,
3, 6, 2, 3, 2, 5, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 7, 3, 3, 2, 3, 2, 2,
3, 3, 2, 3, 2, 5, 3, 7, 2, 3, 2, 2, 3, 3, 2, 3, 3, 3, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 7, 7,
3, 3, 2, 3, 2, 5, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 3, 3, 3, 7, 2, 3, 2, 2, 3, 2, 2, 3, 2, 11,
3, 3, 2, 3, 2, 5, 2, 2, 2, 3, 2, 2, 2, 2, 2, 3, 2, 7, 3, 3, 2, 3, 2, 2, 7, 7, 2, 3, 2, 2,
3, 3, 2, 3, 2, 5, 3, 3, 2, 3, 2, 2, 3, 6, 2, 3, 3, 3, 3, 3, 2, 3, 7, 7, 3, 2, 2, 3, 2, 7,
3, 3, 2, 3, 2, 55, 3, 3, 2, 3, 2, 2, 3, 7, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 3, 3,
3, 3, 2, 3, 2, 5, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2, 3, 7, 2, 3, 2, 2,
3, 3, 2, 3, 2, 5, 3, 3, 2, 3, 2, 11, 3, 3, 2, 3, 2, 2, 2, 2, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2,
3, 3, 2, 3, 2, 5, 3, 7, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2, 2, 55, 2, 3, 2, 14,
3, 3, 2, 3, 2, 5, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 3, 11, 3, 7, 2, 3, 2, 2, 3, 2, 2, 3, 2, 2,
3, 3, 2, 3, 2, 5, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 7, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2,
3, 3, 2, 3, 2, 5, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 3, 3, 3, 3, 2, 3, 2, 11, 3, 3, 2, 3, 2, 7,
3, 3, 2, 3, 2, 5, 3, 3, 2, 3, 2, 2, 3, 7, 2, 3, 3, 3, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2,
3, 3, 2, 3, 2, 5, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 3, 3, 3, 3, 2, 3, 2, 2, 3, 7, 2, 3, 2, 11,
3, 3, 2, 3, 2, 5, 3, 3, 2, 3, 55, 55, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 7, 3, 3, 2, 3, 2, 2,
6, 3, 2, 3, 2, 5, 3, 7, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2,
3, 3, 2, 3, 2, 5, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2, 3, 7, 2, 3, 2, 2, 2, 2, 2, 3, 2, 2,
3, 3, 2, 3, 2, 5, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 7, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2,
3, 3, 2, 3, 2, 5, 2, 2, 2, 3, 2, 11, 3, 3, 2, 3, 3, 3, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 7,
3, 3, 2, 3, 2, 5, 3, 3, 2, 3, 2, 2, 3, 7, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2, 3, 55, 2, 3, 2, 2,
3, 3, 2, 3, 2, 5, 3, 3, 2, 3, 2, 2, 2, 2, 2, 3, 3, 11, 3, 3, 2, 3, 2, 2, 3, 7, 2, 3, 2, 2,
3, 3, 2, 3, 2, 5, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 7, 7, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2,
3, 2, 2, 3, 2, 5, 3, 7, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 11, 3, 3, 2, 3, 2, 2,
3, 3, 2, 3, 2, 5, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 3, 3, 3, 7, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2,
3, 3, 2, 3, 2, 5, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 7, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 11,
3, 2, 2, 3, 2, 5, 3, 3, 2, 3, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 7,
3, 3, 2, 3, 2, 5, 3, 3, 2, 3, 2, 2, 3, 7, 2, 3, 3, 3, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2,
3, 3, 2, 3, 2, 5, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2, 6, 6, 2, 3, 2, 2, 2, 7, 2, 3, 2, 2,
3, 3, 2, 3, 2, 5, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2, 2, 2, 2, 3, 2, 7, 3, 3, 2, 3, 2, 2)
B=2310

f=[10**10 for i in range(100)]
f[2]=2; f[3]=3; f[5]=5; f[6]=5; f[7]=6; f[11]=8; f[14]=8; f[55]=12

def greedy_new(x): #Shriver
	ans=0
	while x>1:
		b=x%B
		d=a[b]
		ans+=b%d+f[d]
		x//=d
	return ans

x1=10**1000+31415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679
x2=10**2000+14142135623730950488016887242096980785696718753769480731766797379907324784621070388503875343276415727
x3=10**3000+27182818284590452353602874713526624977572470936999595749669676277240766303535475945713821785251664274

'''x1=10**1000+314159265358979323846264338327950288419716939937510
x2=10**2000+141421356237309504880168872420969807856967187537694
x3=10**3000+271828182845904523536028747135266249775724709369995'''

_x1=31415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679
#print(bin(_x1)[2:])

for x in [x1,x2,x3]:
	#print(bin(x)[2:])
	print('old: ',greedy(x),'new: ',greedy_new(x))

