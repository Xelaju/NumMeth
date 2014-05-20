from numpy import *
from scipy.linalg import *
from scipy.optimize import fsolve
from matplotlib.pyplot import *
# given constants 

radius0= 43164000
G= 6.67300e-11
M= 5.9722e24
m=10
con=G * M / radius0


#equations relating r, r' and r"
def f(x):
	const= -G*M / (norm (x[:2]))**3
	## !! why dont I have to seperate the force in r direction to the x,y direction with a cosine and sine 
	A=array([[0,0, 1.,0], [0,0,0, 1.], [const, 0,0,0], [0, const , 0, 0]])
	return dot( A, x)


#explicit Euler Method
def eE(f, y0, T, N, flag=False):
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
		r[:, i+1]= r[:, i]+ h* f(r[:,i])

	if flag==True:
		return r[:,-1]
	else:
		return r

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
		return r


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
		func= lambda alph: alph - h* f(alph/2 + r[:,i]/2) - r[:,i]
		r[:, i+1]= fsolve(func, r[:,i])

	if flag==True:
		return r[:,-1]
	else:
		return r


# starting value
y0=array([radius0,0,0,sqrt(con)])
T=13  # Days to be evaluated
N= 86400# Number of intervals to be evaluated


# implementing the functions 

# explicit Euler
projection= eE(f, y0, T, N)
print projection
# implicit Euler
projiE= iE(f,y0,T,N)
print projiE
# implicit Midpoint 
projiM= iM(f,y0,T,N)
print projiM

figure()
plot(projection[0,:],projection[1,:], label='eE')
plot(projiE[0,:],projiE[1,:], label='iE')
plot(projiM[0,:],projiM[1,:], label='iM')
legend()
show()




