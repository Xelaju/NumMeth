from scipy.linalg import solve, norm
from numpy import *
from matplotlib.pyplot import *

def newton(x,F,DF,tol=1e-14,maxit=10000):
    xn = x.copy() # Diese Methode ist unschoen, habe aber bisher keine andere Loesung gesehen
    for i in xrange(maxit):
        s = solve(DF(xn), F(xn))
        xn -= s
        if norm(s) < tol*norm(xn): return xn
    print 'Failure of convergence!'

F = lambda x: array([exp(x[0]**2-x[1]**2)+x[0]**2+x[1]-6./5.,(x[0]+1./4.)**2+(x[1]+3./4.)**2-1])
DF = lambda x: array([[2*x[0]*exp(x[0]**2-x[1]**2)+2*x[0],-2*x[1]*exp(x[0]**2-x[1]**2)+1],[2*(x[0]+1./4.),2*(x[1]+3./4.)]])

# Nullstellen bestimmen mit Newton verfahren
X = array([[-1.2,-1.2],[0.19,0.13],[-0.1,0.23],[0.74,-0.62]])
A = zeros(2)
NST = vstack([A,A,A,A])
it = 0
for i in X:
    NST[it] = newton(i,F,DF)
    it += 1

# Ergebnisse in Tabelle ausgeben
title = '{0:>15s}\t{1:>15s}'.format('Nullstellen (geschaetzt)','Nullstellen (Newton-Verfahren)')
print(''.join(len(title)*['-'])) # separator '-------'
print(title)
print(''.join(len(title)*['-']))
for i in xrange(4):
    x = X[i]
    nst = NST[i]
    row = '{}\t\t\t{}'
    print(row.format(x, nst))
print(''.join(len(title)*['-']))

# Graph formatieren und plotten
x = linspace(-8,8,100)
y = linspace(-20,1,100)
xx, yy = meshgrid(x,y)
z1 = exp(xx**2-yy**2)+xx**2+yy-6./5.
z2 = (xx+1./4.)**2+(yy+3./4.)**2-1

figure()
C1 = contour(x,y,z1,colors='k')
C2 = contour(x,y,z2,colors='r')

clabel(C1, inline=1, fontsize=10)
clabel(C2, inline=1, fontsize=10)
grid(True)
savefig('plot_S03_E3.eps')
