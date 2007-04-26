import numpy as np
import matplotlib.pyplot as plt

def evaliptrig(y,N):
	n = len(y)
	if n % 2 == 0:
		c = np.fft.fft(y)*1./n
		a = np.zeros(N, dtype = complex)
		a[:n/2] = c[:n/2]
		a[N-n/2:] = c[n/2:]
		v = np.fft.ifft(a) * N
		return v
	else: raise TypeError, 'odd_length'

alpha = [0.5,0.9,0.95,0.99]
Linf = [[],[],[],[]]; L2 = [[],[],[],[]]; x = [[],[],[],[]]

N = 4096;
for i, a in enumerate(alpha):

	f = lambda t: 1./np.sqrt(1-a*np.sin(2*np.pi*t))

	n = 2
	while n < 1 + 2**7:
		t = np.linspace(0,1,n+1); y = f(t[:-1])
		v = np.real(evaliptrig(y,N))
		t = np.linspace(0,1,N+1); fv = f(t[:-1])
		d = np.abs(v-fv); Linf[i].append(d.max())
		L2[i].append(np.linalg.norm(d)/np.sqrt(N))
		x[i].append(n)
		n *= 2
	print x

#############################
# Plot
############################
plt.figure()
plt.semilogy(x[0],Linf[0],color='blue', linewidth=1., linestyle='--')
plt.semilogy(x[0],L2[0],color='blue', linewidth=1., linestyle='-')

plt.semilogy(x[1],Linf[1],color='green', linewidth=1., linestyle='--')
plt.semilogy(x[1],L2[1],color='green', linewidth=1., linestyle='-')

plt.semilogy(x[2],Linf[2],color='black', linewidth=1., linestyle='--')
plt.semilogy(x[2],L2[2],color='black', linewidth=1., linestyle='-')

plt.semilogy(x[3],Linf[3],color='red', linewidth=1., linestyle='--')
plt.semilogy(x[3],L2[3],color='red', linewidth=1., linestyle='-')

plt.show()
