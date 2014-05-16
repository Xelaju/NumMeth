from numpy import *

def mcquad(k):
    k = 10**k
    X = random.rand(3,k)
    r = lambda x: sqrt((sqrt(x[0]**2 + x[1]**2) - 0.6)**2 + x[2]**2)
    f = lambda r: 1 + cos(pi * (r/0.3)**2)
    rx = r(X)
    rx_bool = (rx < 0.3).astype(int)
    s = sum(rx_bool)
    I = s * 8./k
    fx = f(rx) * rx_bool
    
    
    var = abs(sum(fx**2)/s - I**2)/(s-1)

    return I, var

for k in xrange(1,7):
    I, var = mcquad(k)
    print "k\tValue\tVarianz\t\t\tVertrauensintervall"
    print "{}\t{}\t{}\t[{},{}]".format(k,I,var,I-var,I+var)

print "Vergleichswert: ",2 * pi**2 * 0.6 * 0.3**2
    
