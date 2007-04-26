from numpy import *

# Startwerte
x0 = array([-2, -1, 1, 2, 1.39174520027])
# Toleranz
tol = 1e-12
# maximale Anzahl Iterationen
maxit = 1000

for x0_it in xrange(x0.size):
    x = x0[x0_it]
    print 'Startwert: ', x
    for i in xrange(5):
        fn = arctan(x)
        dfn = 1. / (1 + x**2)
        s = fn / dfn
        x -= s
        print(x)

G = lambda x: 2*x - (1 + x**2)*arctan(x)
DG = lambda x: 1 - 2*x*arctan(x)

xs = 1.5
x = 0
for i in xrange(maxit):
    gn = G(xs)
    dgn = DG(xs)
    s = gn / dgn
    xs -= s
    if abs(s) < tol * abs(xs): x = xs; break
    
print 'Fixpunkt: '
print(x)



    
