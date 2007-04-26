from numpy import *
from numpy.linalg import qr, norm


def GR(a):
    """
    Gram-Schmidt Verfahren

    Keyword Arguments:
    a -- (m x n) Matrix

    Returns: q, r
    q -- Matrix Q
    r -- Matrix R
    """

    #
    # implementiere das Gram-Schmidt Verfahren hier  #
    #
    a_sh = a.shape
    q = zeros((a_sh[0],a_sh[0]))
    r = zeros(a_sh)
    
    for j in xrange(a_sh[1]):
        v = a[:,j]
        for i in xrange(j):
            r[i,j] = dot(transpose(q[i,:]),v)
            v = v - dot(r[i,j],q[i,:])
        r[j,j] = norm(v)
        q[:,j] = v / r[j,j]
                
    return q, r


def GRmod(a):
    """
    modifiziertes Gram-Schmidt Verfahren
    Keyword Arguments:
    a -- (m x n) Matrix

    Returns: q, r
    q -- Matrix Q
    r -- Matrix R
    """

    #
    # implementiere das modifizierte Gram-Schmidt Verfahren hier #
    #
    a_sh = a.shape
    q = zeros((a_sh[0],a_sh[0]))
    r = zeros(a_sh)
    v = zeros(a_sh)
    for i in xrange(a_sh[0]):
        v[i,:] = a[i,:]
    for i in xrange(a_sh[0]):
        r[i,i] = norm(v[i,:])
        q[i,:] = v[i,:]/r[i,i]
        for j in xrange(i+1,a_sh[1]):
            r[i,j] = dot(transpose(q[i,:]),v[:,j])
            v[:,j] = v[:,j] - dot(r[i,j],q[i,:])

    return q, r


# main

# Matrix Definition
n = 50
m = 50
Z = zeros((m, n))
for i in xrange(m):
    for j in xrange(n):
        Z[i, j] = 1 + min(i, j)

# numpy QR-implementation (als Vergleich)
q0, r0 = qr(Z, mode='full')

## Berechne hier die Guete deines oben implementierten Gram-Schmidt Verfahren
q2, r2 = GR(Z)

## Guete des modifizierten Gram-Schmidt Verfahren
q3, r3 = GRmod(Z)

## print statement
print("numpys qr liefert:         %.10e" % (dot(q0, q0) - eye(n)).max())
print("Gram-Schmidt liefert:      %.10e" % (dot(q2, q2) - eye(n)).max())
print("mod. Gram-Schmidt liefert: %.10e" % (dot(q3, q3) - eye(n)).max())

