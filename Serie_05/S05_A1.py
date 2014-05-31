import numpy as np
from matplotlib.pyplot import *

# Erstellung der Matrix
def get_Ab(n,m):
	N = n + 1 				# N = n + 1
	X = np.linspace(-5,5,m)
	A = np.vander(X,N) # Sehr zu empfehlender Trick

	F = lambda x: 1 / (1 + x**2)
	return A, F(X)

# QR-Zerlegung
def qr_solve(A, b):
	Q, R = np.linalg.qr(A)
	c = np.linalg.solve(R, np.dot(Q.T,b))
	return c

def svd_solve(A, b):
	U, s, VT = np.linalg.svd(A)
	r = np.sum(s > 1e-10)
	return np.dot((VT.T)[:r, :], np.dot(U[:, :r].T, b) / s[:r])

def cheating(A,b):
	return np.linalg.lstsq(A,b)[0]

for n in xrange(2,14,2):
	m = 20
	A,b = get_Ab(n, m)
	print qr_solve(A,b)
	print svd_solve(A,b)
	print cheating(A,b)
	m = 40
	A,b = get_Ab(n,m)
	print qr_solve(A,b)
	print svd_solve(A,b)
	print cheating(A,b)