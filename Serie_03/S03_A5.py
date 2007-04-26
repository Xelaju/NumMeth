from numpy import *
from scipy.optimize import fsolve, bisect
from scipy.integrate import quad
from matplotlib.pyplot import *

#E = 0.2
#b = 2.6
#eps = 1.0

#F = lambda r: E * r**12 - E * b**2 * r**10 + 4 * eps * r**6 - 4 * eps
#DF = lambda r: 12 * E * r**11 - 10 * E * b**2 * r**9 + 24 * eps * r**5
#D2F = lambda r: 112 * E * r**10 - 90 * E * b**2 * e**8 + 120 * eps * r**4 

def bisectNew(x1,x2,F,DF,D2F): 
    x0 = [x1,x2]
    try:
        xd2f = bisect(D2F,x1,x2)
        x0.insert(1,xd2f)
    except:
        print 'D2F hat keine Nullstelle zwischen 0.01 und 100'
    try:
        xdf_1 = bisect(DF,x1,x0[1])
        x0.insert(1,xdf_1)
    except:
        print 'DF hat keine Nullstelle zwischen 0.01 und ', x0[1]
    if x0[1] != x2:
        try:
            xdf_2 = bisect(DF,x0[1],x0[2])
            x0.append(xdf_2)
        except:
            print 'DF hat keine Nullstelle zwischen ', x0[1] ,' und 100'

    nst = []
    for i in xrange(len(x0)-1):
        try:
            nst.append(bisect(F,x0[i],x0[i+1]))
        except:
            print 'F hat keine Nullstelle zwischen ', x0[i] ,' und ', x0[i+1]

    max_value = max(nst)
    return max_value



