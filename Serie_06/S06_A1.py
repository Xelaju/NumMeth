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
A1 = T*F/m

w1, y1 = la.eigh(A1)

data = (array([ones(N-1),-2*ones(N-1),ones(N-1)]))/h**2
diags = array([-1,0,1])
A2 = sp.spdiags(data,diags,N-1,N-1)

w2, y2 = spl.eigsh(A2, k=50)
#print y2.T[1]

w = sqrt(-w2)

n = arange(0,20)

# plot
P1 = pp.figure()
for i in xrange(10):
	P1 = pp.subplot(5,2,i+1)
	P1 = pp.plot(x[1:N],y2.T[i])
	P1 = pp.grid(True)
P1 = pp.savefig('S06_E1_Eigenschwingungen.eps')

P2 = pp.figure()
P2 = pp.plot(n,w2[0:20])
P2 = pp.show()

