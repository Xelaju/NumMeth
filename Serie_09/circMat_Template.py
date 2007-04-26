from numpy import *
from numpy.fft import fft, ifft
import time


def createA(n):
    # Matrix
    A = diagflat(-n * ones(n - 1), -1) + diagflat((n + 1) * ones(n), 0)
    A[0, n - 1] = -n
    return A


def createb(n):
    # rechte Seite
    return sin(linspace(1, 2 * n - 1, n) / (2. * n))


def createu(n):
    # der Vektor u, der A erzeugt
    u = zeros(n)
    u[0] = n + 1
    u[1] = -n
    return u

n = 1000


A = createA(n)
b = createb(n)
u = createu(n)


t0 = time.clock()
d = fft(u)
d_inv = 1./d
D_inv = diag(d_inv)
x_fft = ifft(dot(D_inv,fft(b)))
print time.clock() - t0

t0 = time.clock()
x_solve = linalg.solve(A,b)
print time.clock() - t0

