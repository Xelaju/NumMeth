from numpy import *
from scipy.linalg import *
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
# given constants 

radius0= 43164000
G= 6.67300e-11
M= 5.9722e24
m=10
con=G * M / radius0


#equations relating r, r' and r"
# Mach keine Matrix multiplikation, ein Grauss fuer den Rechner O(n^2), besser ist folgendes:
def f(y):
	dydt = hstack((y[2:4],-M*G/norm(y[0:2])**3 * y[0:2]))
	return dydt


#explicit Euler Method
def eE(f, y0, T, N, flag=False):
	""" f     -> is the function that returns the derivative of the argument 
		y0    -> is the starting value
		T     -> is the time in Days to be evaluated 
		Steps -> is the time step taken for one evaluation point in seconds
		flag  -> False returns all the points evaluated 
		 	  -> True  returns only the last point evaluated
	"""

	T = T *60*60*24    #changing the units of days to seconds for the interval
	t,h = linspace(0, T, N+1, retstep=True)
	y = zeros(( 4 , size(t) ))
	y[:, 0]= y0 # r sieht fuer mich zu sehr nach Radius aus ^

	for i in xrange(N):
		y[:, i+1]= y[:, i]+ h* f(y[:,i])

	if flag==True:
		return y[:,-1]
	else:
		return t, y # Auch die Zeit zurueckgeben

#implicit Euler Method
def iE(f, y0, T, N, flag=False):
	""" f     -> is the function that returns the derivative of the argument 
		y0    -> is the starting value
		T     -> is the time in Days to be evaluated 
		Steps -> is the time step taken for one evaluation point in seconds
		flag  -> False returns all the points evaluated 
		 	  -> True  returns only the last point evaluated
	"""
	T= T *60*60*24    #changing the units of days to seconds for the interval
	t,h =linspace(0, T, N+1, retstep=True)
	r=zeros(( 4 , size(t) ))
	r[:, 0]= y0
	for i in xrange(N):
		func= lambda alph: alph - h* f(alph) - r[:,i]
		r[:, i+1]= fsolve(func, r[:,i])

	if flag==True:
		return r[:,-1]
	else:
		return t, r


# implicit mid-point method
def iM(f, y0, T, N, flag=False):
	""" f     -> is the function that returns the derivative of the argument 
		y0    -> is the starting value
		T     -> is the time in Days to be evaluated 
		Steps -> is the time step taken for one evaluation point in seconds
		flag  -> False returns all the points evaluated 
		 	  -> True  returns only the last point evaluated
	"""
	T= T *60*60*24    #changing the units of days to seconds for the interval
	t,h =linspace(0, T, N+1, retstep=True)
	r=zeros(( 4 , size(t) ))
	r[:, 0]= y0
	for i in xrange(N):
		func= lambda alph: alph - h * f(alph/2 + r[:,i]/2) - r[:,i]
		r[:, i+1]= fsolve(func, r[:,i])

	if flag==True:
		return r[:,-1]
	else:
		return t, r

def SV(y0, T, N, G=6.67300*1e-11, M=5.9722*1e24):
	v0 = y0[2:4].copy()
	y0 = y0[0:2].copy()

	d = size(y0)
	t, h = linspace(T[0],T[1],N+1, retstep=True)

	y = zeros((d, size(t)))
	v = zeros((d, size(t)))
	y[:,0] = y0
	v[:,0] = v0

	f = lambda r: - M*G/norm(r)**3 * r

	for i in xrange(int(N)):
		y[:,i+1] = y[:,i] + h * v[:,i] + h**2 * f(y[:,0]) / 2
		v[:,i+1] = v[:,i] + h * (f(y[:,i]) + f(y[:,i+1])) / 2

	return t, vstack((y,v))

def Energy(y, m=10, G=6.67300*1e-11, M=5.9722*1e24):
	E_kin = m * sum(abs(y[2:4,::])**2, axis=0) / 2
	E_pot = - G * m * M / sum(abs(y[0:2,::])**2, axis=0)**(1./2)
	E_tot = E_kin + E_pot
	return (E_kin, E_pot, E_tot)


# starting value
y0 = array([radius0,0,0,sqrt(con)])
T = [0,5]  # Days to be evaluated
N = 100 # Number of intervals to be evaluated


# implementing the functions 

# explicit Euler
t_EE, y_EE = eE(f, y0, T[1], N)
# implicit Euler
t_IE, y_IE = iE(f,y0,T[1],N)
# implicit Midpoint 
t_IM, y_IM = iM(f,y0,T[1],N)

t_SV, y_SV = SV(y0, T, N)


plt.figure()
plt.plot(t_EE, Energy(y_EE)[2],'g-', label='total Energy EE')
plt.plot(t_IE, Energy(y_IE)[2],'b-', label='total Energy IE')
plt.plot(t_IM, Energy(y_IM)[2],'r-', label='total Energy IM')
plt.plot(t_SV, Energy(y_SV)[2],'y-', label='total Energy SV')
plt.legend()
plt.show()




