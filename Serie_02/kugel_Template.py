from numpy import *
from matplotlib.pyplot import *

# Dichte Werte
rho = array([0.0001, 0.05, 0.4, 0.6, 0.95, 0.9999])
# Toleranz
tol = 1e-5
# Maximale Anzahl Iterationen
maxit = 1000

print('Fixpunktiteration:')
print('{:>10}\t{:>20}\t{:>20}'.format('Dichte', 'Eintauchtiefe', 'Iteration'))
# Schleife ueber Dichte Werte
for rho_it in xrange(rho.size):
    ######################################################################
    # Implementiere hier eine Fixpunktiteration mit maximal maxit        #
    # Iterationsschritten. Verwende einen geeigneten Fehlerindikator und #
    # die Toleranz tol fuer die Abbruchbedingung.                        #
    ######################################################################
    h_curr = 1
    iter = 0
    for i in xrange(maxit):
        h_next = sqrt((h_curr**3 + 4*rho[rho_it])/3)
        iter += 1
        if abs(h_next - h_curr) <= (tol * abs(h_next)): break
        else: h_curr = h_next
        
    ## aligned output
    print('{rho:10.6f}\t{h:10.6f}\t{iter:10d}'.format( rho= rho[rho_it], h= h_next, iter= iter) )
    pass


print('\nSekantenverfahren:')
print('{:>10}\t{:>20}\t{:>20}'.format('Dichte', 'Eintauchtiefe', 'Iteration'))
# Schleife ueber Dichte Werte
for rho_it in xrange(rho.size):
    ######################################################################
    # Implementiere hier ein Newton-Verfahren mit maximal maxit          #
    # Iterationsschritten. Verwende einen geeigneten Fehlerindikator und #
    # die Toleranz tol fuer die Abbruchbedingung.                        #
    ######################################################################
    iter = 0
    F = lambda h: h**3 - 3 * h**2 + 4 * rho[rho_it]
    h0 = 0
    h1 = 2
    f0 = F(h0)
    for i in xrange(maxit):
        fn = F(h1)
        s = fn * (h1 - h0)/(fn - f0)
        h0 = h1; h1 -= s
        iter += 1
        if abs(s) < tol:
            h = h1
            break
        f0 = fn

    ## aligned output
    print('{rho:10.6f}\t{h:10.6f}\t{iter:10d}'.format(rho=rho[rho_it], h=h, iter=iter))

    pass


print('\nNewtonverfahren:')
print('{:>10}\t{:>20}\t{:>20}'.format('Dichte', 'Eintauchtiefe', 'Iteration'))
# Schleife ueber Dichte Werte
for rho_it in xrange(rho.size):
    ######################################################################
    # Implementiere hier ein Sekanten-Verfahren mit maximal maxit        #
    # Iterationsschritten. Verwende einen geeigneten Fehlerindikator und #
    # die Toleranz tol fuer die Abbruchbedingung.                        #
    ######################################################################
    iter = 0
    F = lambda h: h**3 - 3 * h**2 + 4 * rho[rho_it]
    DF = lambda h: 3 * h**2 - 6 * h
    h0 = 1
    for i in xrange(maxit):
        fn = F(h0)
        dfn = DF(h0)
        s = fn / dfn
        h0 -= s
        iter += 1
        if abs(s) < tol: h = h0; break

    ## aligned output
    print('{rho:10.6f}\t{h:10.6f}\t{iter:10d}'.format(rho=rho[rho_it], h=h, iter=iter))

    pass

F0 = lambda h: h**3 - 3 * h**2 + 4 * rho[0]
F1 = lambda h: h**3 - 3 * h**2 + 4 * rho[1]
F2 = lambda h: h**3 - 3 * h**2 + 4 * rho[2]
F3 = lambda h: h**3 - 3 * h**2 + 4 * rho[3]
F4 = lambda h: h**3 - 3 * h**2 + 4 * rho[4]
F5 = lambda h: h**3 - 3 * h**2 + 4 * rho[5]
x = linspace(0, 2, 100)
figure() 
subplot(3,2,1)
plot(x, F0(x))
grid(True)
subplot(3,2,2)
plot(x, F1(x))
grid(True)
subplot(3,2,3)
plot(x, F2(x))
grid(True)
subplot(3,2,4)
plot(x, F3(x))
grid(True)
subplot(3,2,5)
plot(x, F4(x))
grid(True)
subplot(3,2,6)
plot(x, F5(x))
grid(True)
savefig('Kugel_plot.eps')
show()
