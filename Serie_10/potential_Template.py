def trapez(func, a, b, N):
    """
    Numerical quadrature based on trapezoidal rule
    func: handle to y = f(x)
    a,b: bounds of integration interval
    N+1: number of equidistant quadrature points
    """

    from numpy import linspace, sum

    # quadrature nodes
    x = linspace(a, b, N + 1)
    h = x[1] - x[0]

    # quadrature weights: internal nodes: w=1, boundary nodes: w=0.5
    I = sum(func(x[1:-1])) + 0.5 * (func(x[0]) + func(x[-1]))

    return I * h


def simpson(func, a, b, N):
    """
    Numerical quadrature based on Simpson rule
    func: handle to y = f(x)
    a,b: bounds of integration interval
    N+1: number of equidistant quadrature points
    """

    from numpy import linspace, sum

    # quadrature nodes
    x = linspace(a, b, N + 1)
    h = x[1] - x[0]
    xm = linspace(a + h / 2.0, b - h / 2.0, N)

    # quadrature weights: internal nodes: w=1/3, boundary nodes: w=1/6
    I = 2.0 / 3.0 * sum(func(xm)) + 1.0 / 3.0 * sum(func(x[1:-1])) \
        + 1.0 / 6.0 * (func(x[0]) + func(x[-1]))

    return I * h

def trapez2d(func2dim,a,b,Nx,c,d,Ny):
    y_quadr = lambda x: trapez(lambda y: func2dim(x, y), c, d, Ny)
    return trapez(np.vectorize(y_quadr), a, b, Nx)

def simpson2d(func2dim,a,b,Nx,c,d,Ny):
    y_quadr = lambda x: simpson(lambda y: func2dim(x, y), c, d, Ny)
    return simpson(np.vectorize(y_quadr), a, b, Nx)

import numpy as np

# potential
f = lambda x, y: 1. / np.sqrt((x - xp) ** 2 + (y - yp) ** 2)

# point xp = yp = 2.
xp = 2.
yp = 2.
potxpyp = 1.449394876268660
potxpyp_trapez2d = trapez2d(f, -1., 1., 512, -1., 1., 512)
potxpyp_simpson2d = simpson2d(f, -1., 1., 32, -1., 1., 32)
print '!- (xp,yp) = ( 2, 2) ---------------'
print 'trapez2d        : ', potxpyp_trapez2d
print 'error trapez2d  : ', abs(potxpyp_trapez2d - potxpyp)
print 'simpson2d       : ', potxpyp_simpson2d
print 'error simpson2d : ', abs(potxpyp_simpson2d - potxpyp)

# point xp = yp = 10.
xp = 10.
yp = 10.
potxpyp = 0.283080070385743
potxpyp_trapez2d = trapez2d(f, -1., 1., 512, -1., 1., 512)
potxpyp_simpson2d = simpson2d(f, -1., 1., 32, -1., 1., 32)
print '!- (xp,yp) = (10,10) ---------------'
print 'trapez2d        : ', potxpyp_trapez2d
print 'error trapez2d  : ', abs(potxpyp_trapez2d - potxpyp)
print 'simpson2d       : ', potxpyp_simpson2d
print 'error simpson2d : ', abs(potxpyp_simpson2d - potxpyp)

# point xp = yp = 20.
xp = 20.
yp = 20.
potxpyp = 0.141450870624223
potxpyp_trapez2d = trapez2d(f, -1., 1., 512, -1., 1., 512)
potxpyp_simpson2d = simpson2d(f, -1., 1., 32, -1., 1., 32)
print '!- (xp,yp) = (20,20) ---------------'
print 'trapez2d        : ', potxpyp_trapez2d
print 'error trapez2d  : ', abs(potxpyp_trapez2d - potxpyp)
print 'simpson2d       : ', potxpyp_simpson2d
print 'error simpson2d : ', abs(potxpyp_simpson2d - potxpyp)
