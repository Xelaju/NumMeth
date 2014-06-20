# -*- encoding: utf-8 -*-
from numpy import *
from matplotlib.pyplot import *

from numpy.fft import fft
from scipy.special import erf


def cc2(f, a, b, N):
    """
    Clenshaw-Curtis quadrature

    Keyword Arguments:
    func -- f(x)
    a    -- 
    b    -- 
    N    -- Anz. Punkte
    """

    bma = 0.5 * (b - a)
    x = cos( pi * arange(N+1).astype(double)/N ) # Chebychev - Abszissen x in [-1,1]
    x *= bma 
    x += 0.5 * (b + a)

    fx = f(x)
    vx = hstack( (fx, fx[-2:0:-1]) )    # Für die Indexes siehe IndexesSlices.md
                                        # Durch Spiegelung erhält man 2pi-Periodizität
                                        # Achtung! Ränder nicht doppelt zählen
    g = real(fft(vx)) * 0.5 / N         # Die Länge des Arrays beträgt 2N
    A = zeros(N+1)

    A[1:N] = g[1:N] + flipud(g[N+1:])   # flipud() == g[-1:N:-1] kehrt das array um

    A[[0,N]] = g[[0,N]]                 # Kurzschreibweise für A[0] = g[0] 
                                        # und A[N] = g[N]
    w = zeros_like(x)
    w[::2] = 2./(1. - r_[:N+1:2]**2)    # r_[:N+1:2] == arange(0,N+1,2) == [0,2,4,...,N]
                                        # Vorsicht! Jeder zweite Index von w ist 0
    return dot(w,A) * bma               #@all Warum muss hier ein weiteres Mal skaliert werden?


ns = arange(2, 10).astype(double)

f = lambda x: abs(x)
F = 1
errors = [abs(cc2(f, -1, 1, 2 ** i) - F) for i in ns]

figure()
semilogy(2 ** ns, errors, '.-')
xlabel('$N$')
ylabel('Abs. Fehler')
grid(True)
savefig('f1.pdf')

f = lambda x: exp(-x ** 2)
F = 1.4936482656248540508
errors = [abs(cc2(f, -1, 1, 2 ** i) - F) for i in ns]

figure()
semilogy(2 ** ns, errors, '.-')
grid(True)
xlabel('$N$')
ylabel('Abs. Fehler')
savefig('f2.pdf')

# A rather wiggly function
f = lambda x: tanh(20 * sin(12 * x)) + 0.02 * exp(3 * x) * sin(300 * x)
F = 0.0000160929413069564119442643376710

errors = [abs(cc2(f, -1, 1, 2 ** i) - F) for i in ns]

figure()
semilogy(2 ** ns, errors, '.-')
grid(True)
xlabel('$N$')
ylabel('Abs. Fehler')
savefig('f3.pdf')

show()
