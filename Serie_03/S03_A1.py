from scipy.linalg import lu_solve, lu_factor, norm
from numpy import array, arctan 

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
    return x, array(cvg)
    
if __name__=='__main__':
    print '____________________ arctan ___________________'
    F = lambda x: array([arctan(x[0])])
    DF = lambda x: array([[1./(1.+x[0]**2)]])
    x = array([20.])
    print 'Startwert: ', x
    x, cvg = dampnewton(x,F,DF)
    print 'Nullstelle: ', x
    print cvg

    print '_________ Newton Fixpunkt von arctan __________'
    F = lambda x: array([2.*x[0]-(1.+x[0]**2.)*arctan(x[0])])
    DF = lambda x: array([[1.-2.*x[0]*arctan(x[0])]])
    x = array([1.5])
    x, cvg = dampnewton(x,F,DF)
    print x
    print cvg

    print '###############################################'
    print 'Trotz eines Startwerts von x0 = 20, konvergiert'
    print 'das gedaempfte Newtonverfahren. Das allg. Newton-'
    print 'Verfahren divergiert bereits bei ||x|| > ', x

