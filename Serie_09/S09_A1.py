import numpy as np
import matplotlib.pyplot as plt

# |x| must be even
def clenshaw_curtis(f,a,b):
	n = 50
	x = np.linspace(a,b,n)
	a = np.zeros((n/2))
	if n%2 == 0:
		c = np.fft.fft(f(x))
		c = np.fft.fftshift(c)
		c_real = c.real
		a[1:] = 2*c_real[n/2+1:]
		a[0] = c_real[n/2]
		n_a = len(a)
		f_int = 0
		for k in xrange(n_a/2):
			f_int += 2./(1 - (2*k)**2) * a[2*k]
		return f_int
	else: raise TypeError, 'odd_length'

f = lambda x: np.abs(x)

print clenshaw_curtis(f,-1,1)

