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

    #################
    # Dein Code ... #
    #################

    #x = a + 0.5 * (b - a) * (cos(pi * (2 * arange(2 * N) + 1) / 2 * (2 * N)) + 1.) # Funktioniert relativ schlecht
    #x = a + 0.5 * (b - a) * (cos(pi * (2 * arange(2 * N) + 1) / 2 * (2 * N+1)) + 1.) # Funktioniert merkwuerdigerweise exzellent fuer das erste Beispiel
    x = a + 0.5 * (b - a) * (cos(pi * arange(N * 2) / (N * 2 - 1)) + 1.) # Funktioniert mehr oder minder gut

    N = N.astype(int)
    a = np.zeros((N))
    if N%2 == 0:
        z = np.fft.fft(f(x))/(2*N)        # Wir erhalten 2 * N Fourierkoeffizienten
        temp = z[N+1:]          
        a[1:] = (z[1:N] + temp[::-1]).real
        a[0] = z[0].real            # Wir erhalten N Fourierkoeffizienten fuer a

        f_int = 0.0
        for k in xrange(N/2):
            f_int += 2./(1 - (2*k)**2) * a[2*k]
        return f_int
    else: raise TypeError, 'odd_length'

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
