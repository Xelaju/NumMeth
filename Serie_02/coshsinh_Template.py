from numpy import *
from matplotlib.pyplot import *

# Plot of functions cosh(x)-1 and sinh(x)
x = linspace(-3.0, 3.0, 1000)

figure()
plot(x, cosh(x)-1, label=r'$\cosh(x)-1$')
plot(x, sinh(x),   label=r'$\sinh(x)$')
grid(True)
legend(loc='lower right')
xlabel(r'$x$')
ylabel(r'$f(x)$')
xlim(-5, 5)
ylim(-6, 6)
savefig('coshsinh.pdf')
show()

# Newton Methode fuer Nullstellensuche von 1.) cosh-1 und  2.) sinh

# Maximale Anzahl Iterationen
maxit = 100000
# Toleranz
tol = 1e-5

# 1.) cosh-1
xch = 10.0  # Startwert
xch0 = xch

##################################################################
# Implementiere hier ein Newton-Verfahren, mit maximal maxit     #
# Iterationsschritten. Verwende einen geeigneten Fehlerindikator #
# und die Toleranz tol fuer die Abbruchbedingung. Beachte, dass  #
# die Loesung sehr klein ist (Stichwort: Ausloeschung).          #
##################################################################
itrch = 0
for i in xrange(maxit):
    itrch += 1
    fn = cosh(xch0) - 1
    dfn = sinh(xch0)
    s = fn / dfn
    xch0 -= s
    if abs(s) < tol: xch = xch0; break

# 2. sinh
xsh = 10.0  # Startwert
xsh0 = xsh

##################################################################
# Implementiere hier ein Newton-Verfahren, mit maximal maxit     #
# Iterationsschritten. Verwende einen geeigneten Fehlerindikator #
# und die Toleranz tol fuer die Abbruchbedingung. Beachte, dass  #
# die Loesung sehr klein ist (Stichwort: Ausloeschung).          #
##################################################################
itrsh = 0
for i in xrange(maxit):
    itrsh += 1
    fn = sinh(xsh0)
    dfn = cosh(xsh0)
    s = fn / dfn
    xsh0 -= s
    print(xsh0)
    if abs(s) < tol: xsh = xsh0; break


print('cosh(x)-1:')
print('{:>10}\t{:>20}'.format('Nullstelle', 'Iteration'))
print('{xch:10.6f}\t{iter:10d}'.format(xch= xch, iter= itrch) )

print('sinh:')
print('{:>10}\t{:>20}'.format('Nullstelle', 'Iteration'))
print('{xsh:10.6f}\t{iter:10d}'.format(xsh= xsh, iter= itrsh) )
