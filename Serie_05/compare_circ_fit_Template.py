from numpy import array, ones, transpose, dot, sum
from numpy import sin, cos, pi, sqrt, linspace, max
from numpy.linalg import solve, qr, inv, norm
from numpy import random
from scipy.optimize import leastsq
from matplotlib import pyplot as plt
from numpy import * # Darf ich in der Pruefung mehr verwenden als zunaechst eingebunden ist?


def circ_linear_fit(x, y):
    ######################################################
    # Implementiere hier ein lineares Ausgleichsproblem  #
    #                                                    #
    # Input: x, y: Datenpunkte                           #
    # Output: m1, m2, r: Mittelpunkt und Radius          #
    ######################################################
    m = shape(x)[0]
    A = ones((m,3))
    b = zeros(m)
    for i in xrange(m):
        A[i,0] = 2*x[i]
        A[i,1] = 2*y[i]
        b[i] = x[i]**2 + y[i]**2

    U, s, VT = linalg.svd(A)
    r = sum(s > 1e-10)

    return dot((VT.T)[:r, :], dot(U[:, :r].T, b) / s[:r]) 



def GradHess(x, y, z):
    # Gradient and Hessian for Phi
    N = x.size
    xm = x - z[0]
    ym = y - z[1]
    R = sqrt(xm**2 + ym**2)
    # Gradient vector of Phi (Fehler ?)
    GradPhi = array([[dot(z[2]/R - 1, xm)],
                    [dot(z[2]/R - 1, ym)],
                    [N*z[2] - sum(R)]])
    SumX = sum(xm/R)
    SumY = sum(ym/R)
    MixedSum = z[2] * sum(xm*ym/R**3)
    InvSum = z[2] * sum(1/R)
    # Hessian matrix of Phi
    HessPhi = array([[N - InvSum + z[2]*sum(xm**2/R**3), MixedSum, SumX],
                     [MixedSum, N - InvSum + z[2]*sum(ym**2/R**3), SumY],
                     [SumX, SumY, N]])
    return GradPhi, HessPhi

def Jacobi(x,y,z):
    N = x.size
    xm = x - z[0]
    ym = y - z[1]
    R = sqrt(xm**2 + ym**2)
    Jacobi = array([-xm/R,
                    -ym/R,
                    -ones((N,1))]).T
    return Jacobi

def Fz(z,x,y):
    N = x.size
    xm = x - z[0]
    ym = y - z[1]
    R = sqrt(xm**2 + ym**2)
    F = array([R-z[2]]).T
    return F

def func(param, xdata, ydata):
    return (ydata - dot(xdata,param))



def circ_newton_fit(x, y, m1, m2, r):
    ############################################################
    # Implementiere hier ein nichtlineares Ausgleichsproblem   #
    # und berechne die Loesung mit dem Newton Algorithmus.     #
    #                                                          #
    # Input: x, y: Datenpunkte                                 #
    #        m1, m2, r: Startwerte fuer Mittelpunkt und Radius #
    # Output: m1e, m2e, re: Mittelpunkt und Radius             #
    ############################################################
    maxit = 10000
    tol = 1.e-5
    z = array([[m1],[m2],[r]])

    for i in xrange(maxit):
        GradPhi, HessPhi = GradHess(x,y,z)
        HessPhi_inv = inv(HessPhi)
        s = dot(HessPhi_inv,GradPhi)
        z -= s
        if norm(s) < tol * norm(z): return z
    print 'Failure of convergence'



def circ_gaussnewton_fit(x, y, m1, m2, r):
    ##############################################################
    # Implementiere hier ein nichtlineares Ausgleichsproblem     #
    # und berechne die Loesung mit dem Gauss-Newton Algorithmus. #
    #                                                            #
    # Input: x, y: Datenpunkte                                   #
    #        m1, m2, r: Startwerte fuer Mittelpunkt und Radius   #
    # Output: m1e, m2e, re: Mittelpunkt und Radius               #
    ##############################################################
    z = array([[m1],[m2],[r]])
    maxit = 10000
    tol = 1.e-5
    for i in xrange(maxit):
        J = Jacobi(x,y,z)
        F = Fz(z,x,y)
        s = solve(dot(J.T,J),dot(J.T,F))
        z -= s
        if norm(s) < tol * norm(z): return z
    print 'Failure of convergence'



def circ_scipy_least_squares_fit(x, y, m1, m2, r):
    ##############################################################################
    # Benutze hier die Funktion leastsq aus scipy.optimize                       #
    # fuer das nichtlineare Ausgleichsproblem.                                   #
    #                                                                            #
    # Input: x,y: Datenpunkte                                                    #
    #        m1, m2, r: Startwerte fuer Mittelpunkt und Radius                   #
    # Output: m1e, m2e, re: Mittelpunkt und Radius                               #
    ##############################################################################

    z = array([[m1],[m2],[r]])
    maxit = 10000
    tol = 1.e-5
    for i in xrange(maxit):
        J = Jacobi(x,y,z)
        F = Fz(z,x,y)
        F = squeeze(F)
        z = squeeze(z)
        s = leastsq(func, z, args=(J,F))
        z -= s[0]
        if norm(s[0]) < tol * norm(z): return z
    print 'Failure of convergence'




def draw_circle(m1, m2, r, s, lab):
    # simple routine for plotting circles
    theta = linspace(0.0, 2*pi, num=500)
    plt.plot(m1 + r*cos(theta), m2 + r*sin(theta), s, label=lab)


if __name__ == "__main__":
    # Data
    random.seed(0)
    N = 30
    sig = 0.2
    t = linspace(0, 2*pi, N)
    x = cos(t) + random.normal(0, sig, N)
    y = sin(t) + random.normal(0, sig, N)

    plt.figure()
    plt.plot(x, y, '+', label='Data')
    plt.axis('equal')

    print('Lineares Ausgleichsproblem:')
    m1_i, m2_i, r_i = circ_linear_fit(x, y)
    print('(m_1, m_2) = (%.8f, %.8f)  und  r = %.8f' % (m1_i, m2_i, r_i))
    draw_circle(m1_i, m2_i, r_i, 'r-', 'Linear')

    print('Nicht-lineares Ausgleichsproblem mit Newton:')
    # we take solution of linear problem as inital guess
    m1_nl, m2_nl, r_nl = circ_newton_fit(x, y, m1_i, m2_i, r_i)
    print('(m_1, m_2) = (%.8f, %.8f)  und  r = %.8f' % (m1_nl, m2_nl, r_nl))
    draw_circle(m1_nl, m2_nl, r_nl, 'g-', 'Newton')

    print('Nicht-lineares Ausgleichsproblem mit Gauss-Newton:')
    # we take solution of linear problem as inital guess
    m1_nl, m2_nl, r_nl = circ_gaussnewton_fit(x, y, m1_i, m2_i, r_i)
    print('(m_1, m_2) = (%.8f, %.8f)  und  r = %.8f' % (m1_nl, m2_nl, r_nl))
    draw_circle(m1_nl, m2_nl, r_nl, 'g-', 'Gauss-Newton')
    
    print('Nicht-lineares Ausgleichsproblem mit Scipy Least Squares:')
    m1_lsq, m2_lsq, r_lsq = circ_scipy_least_squares_fit(x, y, m1_i, m2_i, r_i)
    print('(m_1, m_2) = (%.8f, %.8f)  und  r = %.8f' % (m1_lsq, m2_lsq, r_lsq))
    draw_circle(m1_nl, m2_nl, r_nl, 'c-', 'Least Squares aus scipy.optimize')



    plt.legend()
    plt.savefig('CircleFit.pdf')
