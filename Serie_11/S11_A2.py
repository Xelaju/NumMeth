from numpy import *
import math

def mcquad(d, k):
	'''
	Input:	d ... dimension of the array
			k ... number of samples (10^k)

	Output:	I 		... Integral
			var2	square of sampple variance
	'''

	N = 10**k
	fx = []
	Fx = []

	x = random.rand(d,N)
	for n in xrange(N):
		fx.append(sum(x[:,n])**2)
	I = sum(fx)/N

	return I

M = 100; k = 3; I = []; var = []
for d in xrange(2,10):
	fx = []
	for m in xrange(M):
		fx.append(mcquad(d,k))

	I_c = sum(fx)/M
	I.append(I_c)
	fx = asarray(fx)
	var.append((sum(fx**2)/M - I_c**2)/(M-1))

I = asarray(I)
var = asarray(var)
I_lower = I - var
I_upper = I + var

print 'Vertrauensintervalle und Varianz: '
for i in xrange(len(I_lower)):
	print '[{},{}] \t {}'.format(I_lower[i], I_upper[i], var[i])