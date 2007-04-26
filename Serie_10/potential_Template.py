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


def trapez2d(func2dim, a, b, Nx, c, d, Ny):
    """
    2 dimensional numerical quadrature based on trapezoidal rule
    func2d: handle to z = f(x,y)
    a,b   : bounds of integration interval in x-direction
    Nx+1  : number of equidistant quadrature points in x-direction
    c,d   : bounds of integration interval in y-direction
    Ny+1  : number of equidistant quadrature points in y-direction
    """

    from numpy import linspace

    # remove the following 2 lines and implement!
    I2d = 0.
    hy = 0.

    return I2d * hy


def simpson2d(func2dim, a, b, Nx, c, d, Ny):

    from numpy import linspace

    """
    2 dimensional numerical quadrature based on Simpson rule
    func2d: handle to z = f(x,y)
    a,b   : bounds of integration interval in x-direction
    Nx+1  : number of equidistant quadrature points in x-direction
    c,d   : bounds of integration interval in y-direction
    Ny+1  : number of equidistant quadrature points in y-direction
    """
    # remove the following 2 lines and implement!
    I2d = 0.
    hy = 0.

    return I2d * hy

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
print 'error trapez2d  : ', potxpyp_trapez2d - potxpyp
print 'simpson2d       : ', potxpyp_simpson2d
print 'error simpson2d : ', potxpyp_simpson2d - potxpyp

# point xp = yp = 10.
xp = 10.
yp = 10.
potxpyp = 0.283080070385743
potxpyp_trapez2d = trapez2d(f, -1., 1., 512, -1., 1., 512)
potxpyp_simpson2d = simpson2d(f, -1., 1., 32, -1., 1., 32)
print '!- (xp,yp) = (10,10) ---------------'
print 'trapez2d        : ', potxpyp_trapez2d
print 'error trapez2d  : ', potxpyp_trapez2d - potxpyp
print 'simpson2d       : ', potxpyp_simpson2d
print 'error simpson2d : ', potxpyp_simpson2d - potxpyp

# point xp = yp = 20.
xp = 20.
yp = 20.
potxpyp = 0.141450870624223
potxpyp_trapez2d = trapez2d(f, -1., 1., 512, -1., 1., 512)
potxpyp_simpson2d = simpson2d(f, -1., 1., 32, -1., 1., 32)
print '!- (xp,yp) = (20,20) ---------------'
print 'trapez2d        : ', potxpyp_trapez2d
print 'error trapez2d  : ', potxpyp_trapez2d - potxpyp
print 'simpson2d       : ', potxpyp_simpson2d
print 'error simpson2d : ', potxpyp_simpson2d - potxpyp
