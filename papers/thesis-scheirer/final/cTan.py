import scipy as sp
from scipy.optimize import fsolve
import pylab as plt
import matplotlib

def f(x):
	return (1*(x**4)) - 4*(x**3) - 5*(x**2)

def df(x):
	dx=0.000001
	return (f(x+dx)-f(x-dx))/(2.0*dx)
	
def g(x):
	x1=x[0]
	x2=x[1]
	y1=f(x[0])
	y2=f(x[1])
	dy1=df(x[0])
	dy2=df(x[1])
	
	out=[dy1-dy2]
	out.append(dy1-((y2-y1)/(x2-x1)))
	return out

sol=fsolve(g,[-2,4])
print(sol)
print(df(sol[0]))
print(f(sol[0]))

x = plt.linspace(-2,5,400)

plt.figure()
plt.title('Generic function of x vs x')
plt.ylabel('f(x)')
plt.xlabel('x')
plt.gca().title.set_fontsize(20)
plt.gca().xaxis.label.set_fontsize(20)
plt.gca().yaxis.label.set_fontsize(20)
plt.xticks(fontsize = 12)
plt.yticks(fontsize = 12)
plt.plot(x,f(x),color='#f36118',linewidth='2',label='$x^4 - 4 x^3 - 5 x^2$')
plt.plot(x, df(sol[0])*(x-sol[0])+f(sol[0]),color='#00cdd8',linewidth='2')
plt.plot(sol[0],f(sol[0]),'ko')
plt.plot(sol[1],f(sol[1]),'ro')
plt.legend(loc='best')
##plt.show()
plt.savefig('generic_common_tangent.png')
