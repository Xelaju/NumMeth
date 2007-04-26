from numpy import *
from scipy import sparse as sp
from scipy.sparse import linalg as spl
from numpy import linalg as la
from matplotlib import pyplot as pp

T = 1.
m = 1.
L = 10.
N = 513.
h = L/N

x = linspace(0,L,N+1)
F = ( eye(N-1,N-1,1) + eye(N-1,N-1,-1) + -2 * eye(N-1) )/h**2

def Fm(F,L,x):
	dm = 0.99
	m1 = 1-dm
	m2 = 1+dm
	mx = (m1*(L-x)+m2*x)/2.
	mx = mx[1:-1]
	MX = diag(mx)
	A = dot(MX,F)
	return A


A1 = Fm(F,L,x)


w1, y1 = la.eigh(A1)

n = arange(0,20)

# plot
P1 = pp.figure()
for i in xrange(10):
	P1 = pp.subplot(5,2,i+1)
	P1 = pp.plot(x[1:N],y1.T[i])
	P1 = pp.grid(True)
P1 = pp.savefig('S06_E1g_Eigenschwingungen.eps')

P2 = pp.figure()
P2 = pp.plot(n,w1[0:20])
P2 = pp.show()

