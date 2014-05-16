# Provides:
# 1. Compute for any interpolation points:
#      diffmat():   pseudospectral differentiation matrix
# 2. Additionally compute for the Chebyshev extreme (Gauss-Chebyshev-Lobatto)
# points on [a,b] (specified using kwargs):
#      chebx():     the points themselves
#      cheb():      the differentiation matrix
#      interpwts(): row vector of interpolation weights
# 3. Compute for Chebyshev extreme points on any [a,b]; no need to give [a,b]
# because these are scale-invariant:
#      coefft():    map collocation values to Chebyshev expansion coefficients
#      coefftmat(): build explicit matrix equivalent to coefft()
#      wdegree():   generate invertible weight matrix  W  which expands in
#                   Chebyshev polynomials and weights by a power of the degree

# References:
#   * Based on code by Greg von Winckel (with major additions); see
#       http://www.scientificpython.net/1/post/2012/04/pseudospectral-differentiation.html
#   * See also L. N. Trefethen, (2000), "Spectral Methods in MATLAB", SIAM Press.
#   * See also L. N. Trefethen, (2013), "Approximation Theory and Approximation
#     Practice", SIAM Press.
#   * For coefft() see p. 2533 in E. Bueler (2007), 'Error bounds for
#     approximate eigenvalues of periodic-coefficient linear delay differential
#     equations', SIAM J. Num. Analysis 45 (6), 2510--2536.  

# Copyright (C) 2013 Ed Bueler

import scipy as sp
from scipy.linalg import solve
from operator import mul

def diffmat(x):
  """Compute the differentiation matrix for  x  is an ordered array
  of grid points.  Uses barycentric formulas for stability.
  """
  n = sp.size(x)
  e = sp.ones((n,1))
  Xdiff = sp.outer(x,e)-sp.outer(e,x)+sp.identity(n)
  xprod = -reduce(mul,Xdiff) # product of rows
  W = sp.outer(1/xprod,e)
  D = W/sp.multiply(W.T,Xdiff)
  d = 1-sum(D)
  for k in range(0,n):  # Set diagonal elements
    D[k,k] = d[k]
  return -D.T

def chebx(N,**kwargs):
  """Compute the Chebyshev extreme points for  N  degree poly interpolation
  in the interval  [-1,1]  (or on  [a,b]  if given by kwargs).  Gives vector of
  N+1  points.  Points are in decreasing order, as in Trefethen (2000,2013).
  """
  x = sp.cos(sp.pi * sp.arange(N+1) / N)
  a = kwargs.get('a',-1.0)
  b = kwargs.get('b',1.0)
  xscale = (b - a) / 2.0
  return xscale * (x + 1.0) + a

def cheb(N,**kwargs):
  """Computes the same results as 'cheb.m' in Trefethen (2000), namely the
  (N+1) x (N+1)  differentiation matrix and the  N+1  points  x  in the interval
  [-1,1]  (or on  [a,b]  if given by kwargs).
  """
  x = chebx(N,**kwargs)
  D = diffmat(x)
  return D, x

def interpwts(N,xx,**kwargs):
  """Given  xx  in  [-1,1],  computes  c_j  so that
    p(xx) = \sum_{j=0}^N c_j f_j
  if  f_j  are the  N+1  collocation values of  p():  f_j = p(x_j).
  Thus we interpolate at  xx  to get  F = p(xx)  by
    c = interpwts(N,xx)
    F = sp.dot(c, f)
  Input  xx  must be a scalar.  See formula (5.13) in Trefethen (2013).
  """
  xc = chebx(N,**kwargs)
  c = sp.zeros((1,N+1))          # a scipy row vector, not a list
  for j in range(N+1):
    if xx == xc[j]:
      c[0,j] = 1.0
      return c
  c[0,0] = 0.5 / (xx - xc[0])
  sign = -1.0
  for j in range(1,N):           # j = 1,2,3,...,N-1
    c[0,j] = sign / (xx - xc[j])
    sign *= -1.0
  c[0,N] = sign * 0.5 / (xx - xc[N])
  return c / sp.sum(c)

def coefft(v):
  """Computes Chebyshev coefficients of polynomial represented by collocation
  values in v.  This method is scale invariant; the coefficients of the
  expansion on [a,b] depend on function values in the same way as for [-1,1].
  """
  N = len(v)-1
  if N==0:
    return v
  vrv = sp.concatenate((v,sp.flipud(v[1:-1])))
  U = sp.fft(vrv) / N
  a = sp.ones((N+1,))
  a[[0,-1]] = 0.5
  return a * U[0:N+1]

def coefftmat(N):
  """Computes real (N+1,N+1) matrix whose action is equivalent to coefft().
  """
  A = sp.zeros((N+1,N+1))
  for j in range(N+1):
    col = sp.zeros(N+1)
    col[j] = 1.0
    A[:,j] = sp.real(coefft(col))
  return A

def wdegree(N,omega):
  """Generate invertible weight matrix  W  whose action   g = W f
  is to expand  f  in Chebyshev polynomials and weight by the  omega-th  power
  of the Chebyshev degree:
    f_j = f(x_j) --> c_j = coeffs in {T_j(x)}
                 --> (1 + j^omega) c_j = weighted coeffs in {T_j(x)}
                 --> g_j = g(x_j)
  The first map is constructed using coefft() to get an explicit matrix  W0.
  The last map  W0^{-1}  is applied using solve().  Returned matrix is
  (N+1) by (N+1).
  """
  W0 = coefftmat(N)
  # now we want  W = W0^{-1} (1 + j^omega) W0
  # so we solve  W0 W = (1 + j^omega) W0  for  W
  B = sp.dot(sp.diag(1.0 + sp.arange(N+1)**omega), W0)
  return solve(W0, B)

if __name__ == "__main__":
    print "pseudospectral.py:  module providing 1D tools for spectral method"

