from numpy import *

def interval(f,a,b,tol):
    if a > b:
        t = b; b = a; a = t
    fa = f(a); fb = f(b)
    if fa * fb > 0: return
    v = 1
    if fa > 0: v = -1
    x = 0.5 * (a + b)                # Attention (a + b) / 2 = 0
    while (b - a > tol) and (a < x) and (x < b):
        if v * f(x) > 0: b = x
        else: a = x
        x = 0.5 * (a + b)
    return x

if __name__ == '__main__':
    f = lambda x: x * exp(x) - 1
    x = interval(f,0,1,1e-9)
    print 'x_interval = ', x
    from scipy.optimize import fsolve
    x = fsolve(f,1)
    print 'x_fsolve = ', x



