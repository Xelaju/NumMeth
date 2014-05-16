from scipy.linalg import solve, norm, lu_factor, lu_solve
from numpy import *
from matplotlib.pyplot import *

def newton(y,F,DF,tol=1e-14,maxit=10000):
    yn = y.copy()
    for i in xrange(maxit):
        s = solve(DF(yn), F(yn))
        yn -= s
        if norm(s) < tol*norm(yn): return yn
    print 'Failure of convergence!'

def dampnewton(x,F,DF,q=0.5,tol=1e-10):
    cvg = []
    lup = lu_factor(DF(x))
    s = lu_solve(lup,F(x))
    xn = x-s
    lam = 1
    st = lu_solve(lup,F(xn)) # simplified Newton                                                                                                                                     
    while norm(st) > tol*norm(xn):
        while norm(st) > (1-lam*0.5)*norm(s):
            lam *= 0.5
            if lam < 1e-10:
                cvg = -1
                print 'Failure of convergence'
                return x, cvg
            xn = x-lam*s
            st = lu_solve(lup,F(xn)) # simplified Newton                                                                                                                             
        cvg += [[lam, norm(xn), norm(F(xn))]]
        x = xn
        lup = lu_factor(DF(x))
        s = lu_solve(lup,F(x))
        lam = min(lam/q, 1.) # Wozu dieser Test?                                                                                                                                     
        xn = x-lam*s
        st = lu_solve(lup,F(xn)) # simplified Newton                                                                                                                                 
    x = xn
    return x

X = linspace(-2,2,200)
Y = zeros(200)
itr = 0
for x in X:
    y = array([0])
    F = lambda y: array([exp(x**2-y[0]**2)+x**2+y[0]-6./5.])
    DF = lambda y: array([[-2*y[0]*exp(x**2-y[0]**2)+1]])
    #Y[itr] = newton(y,F,DF) Das Newton-Verfahren konvergiert leider nicht
    Y[itr] = dampnewton(y,F,DF)
    itr += 1

figure()
plot(X,Y)
show()


