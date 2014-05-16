import numpy as np
from scipy.optimize import fsolve

# constants
m = 1.  # Masse
L = 1.  # Stablaenge
k = 9.5 # Federkonstante
g = 10. # g

# U
U = lambda theta: 0.5*k*theta**2 + m*g*L*cos(theta)
# dU/dtheta
dUdtheta = lambda theta: k*theta - m*g*L*np.sin(theta)
# dU2/dtheta2
dU2dtheta2 = lambda theta: k - m*g*L*np.cos(theta)

# neue Gleichgewichtslage
theta0 = 1. # waehle einen guten Startwert fuer fsolve... Tipp: NICHT NULL!!!
# -- implementiere die Nullstellen suche HIER -- #

def newton(x, tol, maxit):
    for i in xrange(maxit):
        s = dUdtheta(x) / dU2dtheta2(x)
        x -= s
        if abs(s) < tol * abs(x): return x

# Tipp: theta = fsolve(funktion,startwert)
theta = newton(theta0, 1e-7, 200)

# ausgabe
print 'Theta       = ',theta,' rad, bzw ',theta*(360./(2.*np.pi)),' Grad'
print 'dU/dtheta   = ',dUdtheta(theta)
print 'dU2/dtheta2 = ',dU2dtheta2(theta)

