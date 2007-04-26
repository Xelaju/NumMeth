from scipy.linalg import lu_solve, lu_factor, norm, solve
from numpy import *
from matplotlib.pyplot import *

def newton(x,x0,F,DF,tol=1e-14,maxit=10000):
    xn = x.copy()
    res = zeros(4); res[1:3] = xn; res[3] = norm(xn-x0)
    k = 1
    for i in xrange(maxit):
        s = solve(DF(xn), F(xn))
        xn -= s
        res1 = zeros(4); res1[0] = k; res1[1:3] = xn; res1[3] = norm(xn-x0)
        res = vstack((res,res1))
        k += 1
        if norm(s) < tol*norm(xn): break
    logdiff = diff(log(res[:,3]))
    rates = logdiff[1:]/logdiff[:-1]
    return xn, rates        


if __name__=='__main__':
    F = lambda x: array([ exp(x[0]*x[1]) + x[0]**2 + x[1] - 6./5., x[0]**2 + x[1]**2 + x[0] - 11./20.])
    DF = lambda x: array([[ x[1]*exp(x[0]*x[1]) + 2*x[0], x[0]*exp(x[0]*x[1]) + 1.],[2*x[0] - 1, 2*x[1]]])
    
    # Nullstelle (mit Newton-Ver. ermittelt) 
    x0 = array([-0.60891357, -0.88777127])
    
    # Berechnung der Nullstellen und Konvergenz nach Startwert
    x = array([[3./5.,1./2.],[2./5.,1./4.],[-23./5.,41./5.]])
    Y = zeros(2)
    y = vstack((Y,Y,Y))
    rates = []
    itr = 0
    for i in x:
        y[itr], rates1 = newton(i,x0,F,DF)
        itr += 1
    
    # Ausgabe der Werte
    title = '{0:>15s}\t{1:>15s}'.format('Startwert','Nullstellen')
    print(''.join(len(title)*['-']))
    print(title)
    print(''.join(len(title)*['-']))
    for i in xrange(3):
        sw = x[i]
        nst = y[i]
        row = '{}\t{}'
        print(row.format(sw, nst))
    print(''.join(len(title)*['-']))

    print rates
