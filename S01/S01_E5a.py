from numpy import *
from scipy.special import gamma
from matplotlib.pyplot import *

def sinApprox(x, tol, N):
    t = lambda n, x: ((-1)**n / gamma(2 * n + 2)) * x**(2 * n+1)  
    S = 0
    n = 0
    while ( all( abs(t(n,x)) > abs(tol * S) ) ) and (n <= N):   # Was bewirkt die Funktion all() ?
        S += t(n,x)
        n += 1
    return S

if __name__ == '__main__':
    x1 = linspace(pi, 6*pi, 101)
    x2 = linspace(6*pi, 30*pi, 101)

    semilogy(x1, abs( sin(x1) - sinApprox(x1, 1e-7, 50) ) ) # Da bisher alle Funktionen arrays als Argument akzeptierten, nahm ich an, dass die plot() bzw. semilogy() Funktion über die gegebene
                                                            # Funktion itteriert. Sie übergibt aber einfach das gesamte array. Es wäre nett, wenn du die funktionsweise der plot() Funktion erläutern könntest.
    semilogy(x2, abs( sin(x2) - sinApprox(x2, 1e-7, 50) ) )
    grid(True)
    show()




    
        

         
