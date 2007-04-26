from numpy import *
from scipy.optimize import fsolve
from scipy.integrate import quad
from matplotlib.pyplot import *
from S03_E5 import bisectNew

EPS = linspace(0.01,2,20)
num = 0
for eps in EPS:
    B = linspace(0,3,500)
    Y1 = zeros(500)
    E = 0.2
    itr = 0
    for b in B:
        F = lambda r: E * r**12 - E * b**2 * r**10 + 4 * eps * r**6 - 4 * eps
        DF = lambda r: 12 * E * r**11 - 10 * E * b**2 * r**9 + 24 * eps * r**5
        D2F = lambda r: 112 * E * r**10 - 90 * E * b**2 * e**8 + 120 * eps * r**4
        Y1[itr] = bisectNew(0.01,100,F,DF,D2F)
        itr += 1

    Y2 = zeros(500)
    itr = 0
    for b in B:
        r0 = Y1[itr]
        U = lambda r: 4 * eps * (1/r**12-1/r**6)
        dthetadr = lambda r: b / ( r**2 * sqrt(1-b**2/r**2 - (U(r))/E ))
        Y2[itr] = pi - 2*quad(dthetadr, r0, np.inf)[0]
        itr += 1
    
    figure()
    plot(B,Y1)
    plot(B,Y2)
    savefig('plotNEW_S03_E4_{}.eps'.format(num))
    num += 1



