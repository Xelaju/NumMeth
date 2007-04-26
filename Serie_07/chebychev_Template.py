from numpy import *
from matplotlib.pyplot import *

def interp_chebychev_barycentric(x, y, xx):
    """
    Purpose: berechnet das Interpolations-Polynom p(xx) mittels
             baryzentrischer Interpolationsformel mit Stueztstellen x,
             Stuetzwerten y mit der vereinfachten Formel fuer Chebychev Abszissen

    Keyword Arguments:
    x  -- Chebychev Abszissen
    y  -- Stuetzwerte
    xx -- Auswertungspunkte
    """ 
    
    ####################################################################################
    # Implementiere hier die skalierungsinvariante baryzentrische Interpolationsformel #
    # fuer Chebychev Abszissen                                                         #
    ####################################################################################
    p = zeros(len(xx))
    eps = 1.e-5
    for i, s in enumerate(xx):
        P = 0.; Q = 0.
        for k, c in enumerate(x):
            q = (-1)**k / (s+eps - c)
            Q += q
            P += q * y[k]
        p[i] = P/Q
    return p

def chebpts(a, b, n):
    """
    Gibt die Chebychev-Abszissen auf dem Intervall a, b zurueck
    
    Keyword Arguments:
    a -- 
    b -- 
    n -- Anz. Knoten
    """
    return a + 0.5 * (b - a) * (cos(pi * arange(n) / (n - 1)) + 1.)

if __name__ == '__main__':

    # ----------------------------------------------------------------------
    # Runge
    f = lambda x: 1. / (1 + x ** 2)
    a = -5
    b = 5
    n = 20
    x = chebpts(a, b, n)
    xx = linspace(a, b, 1000)
    y = f(x)
    pxx = interp_chebychev_barycentric(x, y, xx)

    figure()
    semilogy(xx, abs(pxx - f(xx)))
    grid(True)
    title(r'$1/(1+x^2)$ mit %d Knoten' % n)
    ylabel('Abs. Error')
    xlabel('x')
    grid(True)
    savefig('ci-runge.pdf')
    show()

    # ----------------------------------------------------------------------
    # A rather wiggly function
    f = lambda x: tanh(20 * sin(12 * x)) + 0.02 * exp(3 * x) * sin(300 * x)
    ns = array([50, 150, 300,  600, 1000, 2000, 5000])
    a = -1
    b = 1

    figure()
    for n in ns:
        xx = linspace(a, b, 5 * n)
        x = chebpts(a, b, n)
        y = f(x)

        pxx = interp_chebychev_barycentric(x, y, xx)

        semilogy(xx, abs(pxx - f(xx)), label='n = %d' % n, lw=0.5)
    xlabel('x')
    ylabel('Abs. Error')

    ax = gca()
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    grid(True)
    savefig('ci-rather-wiggly-function.pdf')
    show()

    # ----------------------------------------------------------------------
    f = lambda x: sin(10 ** 5 * x)
    a = -1
    b = 1
    n = 10 ** 6
    x = chebpts(a, b, n)
    xx = linspace(a, b, 1000)
    y = f(x)
    pxx = interp_chebychev_barycentric(x, y, xx)

    figure()
    semilogy(xx, abs(f(xx) - pxx))
    grid(True)
    title(r'$\sin(10^5 x)$ mit %.e Knoten' % n)
    ylabel('Abs. Error')
    xlabel('x')
    savefig('ci-sin10p5x.pdf')
    show()

    # ----------------------------------------------------------------------
    f = lambda x: sin(exp(10 * x))
    a = -1
    b = 1
    n = 1000
    x = chebpts(a, b, n)
    xx = linspace(a, b, 1000)
    y = f(x)
    pxx = interp_chebychev_barycentric(x, y, xx)
    semilogy(xx, abs(pxx - f(xx)), lw=0.5)
    grid(True)
    title(r'$\sin(e^{10x})$ mit %d Knoten' % n)
    ylabel('Abs. Error')
    xlabel('x')
    ylim(ymin=1e-6)
    grid(True)
    savefig('ci-sinexp10x.pdf')
    show()
