from numpy import sqrt, pi, arange, array
from numpy.random import rand

# formulate calculation of pi as integration on [0,1]
func = lambda x, y: sqrt(x**2 + y**2)

# implement Monte-Carlo-method


def mcquad(k, func):
    # Monte-Carlo-Quadrature for function func on [0,1]^2
    # M repetitions with 10**k samples
    M = 10**k
    dim = 2

    # generate M samples in [0,1]
    x = rand(dim,M)
    # evaluate function at sample points
    fx = func(x[0,:],x[1,:])
    # evaluate integral (expectation value)
    I = sum(fx)/M

    return I

# we draw 10**1, ... 10**k samples
k = 7
results = []
for m in xrange(k):
    results.append(mcquad(m, func))

print('Monte-Carlo to compute pi/4 = %f:' % (0.25 * pi))
print(array(results))
