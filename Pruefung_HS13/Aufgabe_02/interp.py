import numpy as np
import warnings
warnings.simplefilter("ignore", np.RankWarning)
import matplotlib.pylab as plt
from numpy.linalg import solve, cond, norm

# Stuetzstellenvektoren (nicht veraendern!)
np.random.seed(3)
x02 = np.random.random( 2+1)
x04 = np.random.random( 4+1)
x06 = np.random.random( 6+1)
x08 = np.random.random( 8+1)
x10 = np.random.random(10+1)
x12 = np.random.random(12+1)
x14 = np.random.random(14+1)
x16 = np.random.random(16+1)
x18 = np.random.random(18+1)
x20 = np.random.random(20+1)


# Unteraufgabe b) ######################################################
def interp_monom(x):

    """
    Erstellt die Matrix A des Gleichungssystems in der Monombasis.

    Input: x ... Vektor mit Stuetzstellen

    Output: A ... Matrix

    """

    # Groesse von x
    n = x.size
    A = np.vander(x,n)

    return A


def unteraufgabe_b():
    A10 = interp_monom(x10)
    A20 = interp_monom(x20)
    return A10,A20

# Unteraufgabe c) ######################################################
def unteraufgabe_c():
    f = lambda x: np.sin(10*x*np.cos(x))

    [A10,A20] = unteraufgabe_b()
    b10 = np.zeros_like(x10)
    b20 = np.zeros_like(x20)
    alpha10 = np.zeros_like(x10)
    alpha20 = np.zeros_like(x20)

    b10 = f(x10)
    b20 = f(x20)
    alpha10 = solve(A10,b10)
    alpha20 = solve(A20,b20)

    x = np.linspace(0.0, 1.0, 100)
    pi10 = np.zeros_like(x)
    pi20 = np.zeros_like(x)

    pi10 = np.polyval(alpha10,x)
    pi20 = np.polyval(alpha20,x)

    plt.figure()
    plt.plot(x,f(x),"-b",label=r"$f(x)$")
    plt.plot(x,pi10 ,"-g",label=r"$p_{10}(x)$")
    plt.plot(x10,b10,"dg")
    plt.plot(x,pi20 ,"-r",label=r"$p_{20}(x)$")
    plt.plot(x20,b20,"or")
    plt.grid(True)
    plt.xlabel(r"$x$")
    plt.ylabel(r"$y$")
    plt.legend()
    plt.savefig("interpolation.eps")

# Unteraufgabe d) ######################################################
def unteraufgabe_d():
    n = np.arange(2,21,2)
    xs = [x02,x04,x06,x08,x10,x12,x14,x16,x18,x20]
    f = lambda x: np.sin(10*x*np.cos(x))

    residuals = np.zeros_like(n,dtype=np.floating)
    condition = np.ones_like(n,dtype=np.floating)

    for i, x in enumerate(xs):
        b = f(x)
        A = interp_monom(x)
        alpha = solve(A,b)
        residuals[i] = norm(np.dot(A,alpha) - b)
        condition[i] = cond(A)

    plt.figure()
    plt.plot(n,residuals,"-o")
    plt.grid(True)
    plt.xlabel(r"$n$")
    plt.ylabel(r"$\|A \alpha - b\|_2$")
    plt.savefig("residuals.eps")

    plt.figure()
    plt.semilogy(n,condition,"-o")
    plt.grid(True)
    plt.xlabel(r"$n$")
    plt.ylabel(r"$\log(\mathrm{cond}(A))$")
    plt.savefig("condition.eps")


# Unteraufgabe e) ######################################################
def lagrange_poly(xs):

    """
    Berechnet die Koeffizienten der Lagrange-Polynome

    Input: xs ... Stuetzstellen

    Output: l ... Koeffizienten der Lagrange-Polynome, wobei
                  l[i,:] die Koeffizienten des i-ten Polynom sind
    """
    n = xs.size
    l = np.zeros((n,n))
    eps = 1.e-15
    x = np.linspace(0.,1.,15) + eps


    for k in xrange(n):
        y = np.zeros_like(x)
        for i, xn in enumerate(x):
            temp = 1.
            for j in xrange(n):
                if j != k:
                    temp *= (xn - xs[j])/(xs[k] - xs[j])
            y[i] = temp
        l[k,:] = np.polyfit(x,y,n-1)

    return l


def unteraufgabe_e():
    n = 20
    xs = np.linspace(0.0, 1.0, n+1)

    l = lagrange_poly(xs)
    print l[-1]



# Unteraufgabe f) ######################################################
def divdiff(x,y):

    """
    Berechnet die dividierten Differenzen

    Input: x ... Stuetzstellen
           y ... Datenpunkte, d.h. f(x)

    Output: y ... dividierte Differenzen

    """

    n = y.size
    T = np.zeros((n,n))
    T[:,0] = y

    for level in xrange(1,n):
        for i in xrange(n-level):
            T[i, level] = (T[i+1,level-1] - T[i,level-1])/(x[i+level]-x[i])

    # Rueckgabe
    return T[0,:]


def evalNewton(xs,dd,xx):

    """
    Evaluiert das Newton-Polynom

    Input: xs  ... Stuetzstellen
           dd ... dividierte Differenzen
           xx ... Auswertungspunkte

    Output: yy ... Newton-Polynom ausgewertet an xx

    """

    n = xs.size
    yy = np.zeros_like(xx)

    for k in xrange(len(yy)):
        yy[k] = dd[-1]
        for i in xrange(2, len(dd)+1):
            yy[k] = (xx[k] - xs[-i]) * yy[k] + dd[-i]

    # Rueckgabe
    return yy


def unteraufgabe_f():
    n = 20
    x = np.linspace(0.0, 1.0, n+1)
    f = lambda x: np.sin(10*x*np.cos(x))

    X = np.linspace(0.0, 1.0, 1000)
    Y = X

    dd = divdiff(x,f(x))
    Y = evalNewton(x,dd,X)

    plt.figure()
    plt.plot(X,f(X),"-b",label=r"$f(x)$")
    plt.plot(x,f(x),"or")
    plt.plot(X,Y,"-r",label=r"$p_{20}(x)$")
    plt.grid(True)
    plt.xlabel(r"$x$")
    plt.ylabel(r"$y$")
    plt.legend()
    plt.savefig("divdiff.eps")


# Unteraufgabe g) ######################################################
def sum_abs_lagrange_poly(l,x):

    """
    Berechnet \sum_{i=0}^n |l_i(x)|.

    Input: l ... Koeffizienten der Lagrange-Polynome wie produziert durch
                 die Funktion lagrange_poly
           x ... Auswertungspunkte

    Output: sum_alp ... \sum_{i=0}^n |l_i(x)|

    """

    n = l.shape[0]
    sum_alp = np.zeros_like(x)

    ####################################################################
    #
    # Berechnen Sie hier \sum_{i=0}^n |l_i(x)|
    #
    ####################################################################

    return sum_alp

def unteraufgabe_g():
    # Sampling punkte
    x = np.linspace(0.0,1.0,1000)
    N = np.arange(2,16)

    LU = np.ones_like(N,dtype=np.floating)
    LT = np.ones_like(N,dtype=np.floating)

    # Approximiere Lebesgue-Konstante
    for i,n in enumerate(N):
        ################################################################
        #
        # xU = np.linspace(0.0,1.0,n)
        #
        # LU[i] = ...
        #
        # j = np.arange(n+1)
        # xT = 0.5*(np.cos((2.0*j+1.0)/(2.0*(n+1.0))*np.pi) + 1.0)
        #
        # LT[i] = ...
        #
        ################################################################
        continue

    # Plot
    plt.figure()
    plt.semilogy(N,LU,"-ob",label=r"Aequidistante Punkte")
    plt.semilogy(N,LT,"-og",label=r"Chebyshev Punkte")
    plt.grid(True)
    plt.xlim(N.min(),N.max())
    plt.xlabel(r"$n$")
    plt.ylabel(r"$\Lambda^{(n)}$")
    plt.legend(loc="upper left")
    plt.savefig("lebesgue.eps")

########################################################################
unteraufgabe_b()
unteraufgabe_c()
unteraufgabe_d()
unteraufgabe_e()
unteraufgabe_f()
unteraufgabe_g()

