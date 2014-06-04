from numpy import linspace, array, hstack, vstack, size, zeros, shape
from numpy.linalg import norm
import matplotlib.pyplot as plt

def SV(y0, T, N):
	v0 = y0[2:4].copy()
	y0 = y0[0:2].copy()

	d = size(y0)
	t, h = linspace(T[0],T[1],N+1, retstep=True)

	y = zeros((d, size(t)))
	v = zeros((d, size(t)))
	y[:,0] = y0
	v[:,0] = v0

	f = lambda r: 48 * ((1/norm(r))**14 - (1/norm(r))**8 / 2) * r

	for i in xrange(int(N)):
		y[:,i+1] = y[:,i] + h * v[:,i] + h**2 * f(y[:,i]) / 2
		v[:,i+1] = v[:,i] + h * (f(y[:,i]) + f(y[:,i+1])) / 2

	return t, vstack((y,v))


if __name__ == '__main__':
	B = array(range(3,62,3)) * 1./20 # [0.15,0.3,0.45,...,3]
	v0 = array([1, 0])
	T = [0,15]
	Z = 0.02 # Zeitschritt
	N = int(T[1]/Z)
	print N

	result = zeros((size(B), 4, N+1))
	for i, b in enumerate(B):
		r0 = array([-10, b])
		y0 = hstack((r0,v0))
		result[i,:,:] = SV(y0, T, N)[1] # Ich hatte das Problem, dass Python glaubte die letze Zeile gehoere nicht mehr
										# zur Schleife, gab aber auch keine Warnung aus. Dieses Problem taucht 
										# haeufiger auf. Gibt es eine Methode indentation Probleme schnell zu beheben, bzw.
										# zu vermeiden?


	plt.figure()
	for i in xrange(size(B)):
		plt.plot(result[i,0,:], result[i,1,:])
	plt.show()