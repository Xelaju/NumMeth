from numpy.linalg import cond, solve
from numpy import *
from matplotlib.pyplot import *
from matplotlib.mlab import find
from functools import partial

def bary_eval(yy, xs, lam, x):
    """
    Baryzentrische Interpolations Formel
    
    Eq (5.9) `Approximation Theory and Approximation Practice` p. 35
    Keyword Arguments:
    yy  -- f_i
    xs  -- x_i
    lam -- lambda_i
    x   -- evaluation points
    """
    # Funktionswerte an der Stelle x
    y = zeros_like(x)

    xs = xs.reshape(-1, 1)
    x = x.reshape(1, -1)

    # finde Evaluationspunkte welche sehr nahe an den Stuetzstellen liegen
    lxx = product(xs - x, axis=0)
    tol = 1e-20
    # array of boolean's
    idx = abs(lxx) < tol

    idx_n = find(idx)
    # finde die zugehoerigen Funktionswerte
    yidx = zeros(len(idx_n), dtype=int)
    for i, ii in enumerate(idx_n):
        d = abs(squeeze(x)[ii] - xs)
        yidx_local = find(d == min(d))
        yidx[i] = yidx_local
    # an den Stuetzstellen settz man direkt die
    # bekannten Funktionswerte
    y[idx_n] = yy[yidx]

    # x, mit genuegend Abstand zu Stuetzstellen
    idx = logical_not(idx)  # invertiere idx
    xr = squeeze(x)[idx].reshape(1, -1)
    tmp = (squeeze(lam) * squeeze(yy)).reshape(-1, 1)
    y2 = (lxx[idx]).reshape(1, -1) * sum(tmp / (xs - xr), axis=0)

    y[idx] = y2
    return y

def chebyshev_pts(n, a, b):
    """
    Gibt die n-Chebychev Abszissen auf dem Intervall [a,b]
    zurueck

    Keyword Arguments:
    n -- Anzahl Punkte
    a -- Intervall Anfang
    b -- Intervall Ende
    """
    assert b > a

    ii = arange(n)
    cc = sort(cos(ii * pi / (n - 1)))
    return (cc + 1) * (b - a) / 2 + a

if __name__ == '__main__':


    #######################################################################
    # Beispiele mit exakter Loesung f, u = f'', f(0) = 0, f(1) = 1        #
    # Ueberpruefe damit Deine Implementation.                             #
    #######################################################################
    ex0 = lambda a: {'u': lambda x: x ** a,
                     'f': lambda x: (x ** (a + 2) - x) / ((a + 1) * (a + 2))}

    ex1 = lambda a: {'u': lambda x: exp(x) * sin(a * x),
                     'f': lambda x: (2 * a - 2 * a * x + 2 * a * exp(1) * x * cos(a) - 2 * a * exp(x) * cos(a * x) - exp(1) * x * sin(a) 
                                     + a ** 2 * exp(1) * x * sin(a) + exp(x) * sin(a * x) - a ** 2 * exp(x) * sin(a * x)) / (1. + a ** 2) ** 2}    
    ## Code Beispiel:
    # 
    u = ex0(a=10)['u']
    f = ex0(a=10)['f']
    # 
    xx = linspace(0,1,1000)
    f_exact = f(xx)
    u_exact = u(xx)
    plot(xx,u_exact)
    show()

    
