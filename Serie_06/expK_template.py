from scipy import zeros, dot, random, mat, linalg, diag, sqrt, sum, hstack, ones, exp
from scipy.linalg import norm, eig, expm, expm2, expm3
import time
from S06_E2 import expD, expM, expA, expL
from numpy import *


def delta(i, j):
    if i == j:
        return 1
    else:
        return 0


def ConstructMatrix(N, case='minij'):
    H = mat(zeros([N, N], dtype=complex))
    for i in xrange(N):
        for j in xrange(N):
            if case == 'sqrt':
                H[i, j] = 1L * \
                    sqrt((i + 1) ** 2 + (j + 1) ** 2) + (i + 1) * delta(i, j)
            elif case == 'dvr':
                if i != j:
                    t = i - j
                    H[i, j] = 2. / t ** 2
                    H[i, j] *= (-1) ** t
            elif case == 'minij':
                H[i, j] = float(1 + min(i, j))
            else:
                return None
    return H


if __name__ == "__main__":
    from scipy import linspace, diag

    k = 50
    dt = 0.01
    N = 2 ** 9

    print('Construct the matrix ...')
    case = 'minij'
    case = 'sqrt'
    # case = 'dvr'

    if case == 'dvr':
        eps = 0.5  # physiscal parameter
        a = -10.0
        b = 10.0  # comp. interval
        N = 2 ** 9
        h = (b - a) / N
        x = linspace(a, b, N)
        # potential
        potential = lambda x: x ** 2
        t0 = time.clock()
        A = ConstructMatrix(N, case='dvr')
        A *= (eps / h) ** 2
        A += diag(potential(x))
        A *= dt / eps
        print(time.clock() - t0)
    else:
        t0 = time.clock()
        A = ConstructMatrix(N, case)
        print(time.clock() - t0)

    # initial value
    v = mat(ones(N, dtype=complex)).T
    v = v / norm(v)
    t = 1.e-2
    print 'expM braucht: '
    t0 = time.clock()
    exp_M = expM(A,v,t)
    print(time.clock() - t0)
    print '\n'

    print 'expD braucht: '
    t0 = time.clock()
    exp_D = expD(A,v,t)
    print(time.clock() - t0)
    print 'Der Fehler bezueglich expM ist: '
    print norm(exp_M - exp_D)
    print '\n'
    
    print 'expA braucht: '
    t0 = time.clock()
    exp_A = expA(A,v,t)
    print(time.clock() - t0)
    print 'Der Fehler bezueglich expM ist: '
    print norm(exp_M - exp_A)
    print '\n'

    print 'expL braucht: '
    t0 = time.clock()
    exp_L = expL(A,v,t)
    print(time.clock() - t0)
    print 'Der Fehler bezueglich expM ist: '
    print norm(exp_M - exp_L)
    print '\n'

    #expA(A,v,t)[0:2,0]
    #expL(A,v,t)[0:2,0]
