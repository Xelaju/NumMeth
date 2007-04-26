import numpy as np
import matplotlib.pylab as plt


def simpson(func, a, b, N):
    """
    Zusammengesetzte Simpson-Regel

    Input: func ... Funktions handle zu y = f(x)
           a,b  ... Integrations-Intervallgrenzen
           N    ... Anzahl aequidistanter Teilintervalle

    Output: I ... Integral I
    """
    I = 0
    #########################################################
    #                                                       #
    # Implementiere hier die zusammengesetzte Simpson-Regel #
    #                                                       #
    #########################################################

    # quadrature nodes
    x = np.linspace(a, b, N + 1)
    h = x[1] - x[0]
    xm = np.linspace(a + h / 2.0, b - h / 2.0, N)

    # quadrature weights: internal nodes: w=1/3, boundary nodes: w=1/6
    I = 2.0 / 3.0 * sum(func(xm)) + 1.0 / 3.0 * sum(func(x[1:-1])) \
        + 1.0 / 6.0 * (func(x[0]) + func(x[-1]))

    return I * h


def mc(func, a, b, N):
    """
    Monte-Carlo Quadratur

    Input: func ... Funktions handle zu y = f(x)
           a,b  ... Integrations-Intervallgrenzen
           N    ... Anzahl Samples

    Output: I ... Integral
    """
    I = 0
    ################################################
    #                                              #
    # Implementiere hier die Monte-Carlo Quadratur #
    #                                              #
    ################################################

    # Ausgabe
    x = np.random.rand(N)
    x = a + (b - a) * x
    fx = func(x)
    I = sum(fx) / N * (b - a)

    #var = (sum(fx**2)/N - I**2)/(N-1)

    return I

# Konstanten
a = 0.0  # untere Intervallsgrenze
b = 2.0  # obere  Intervallsgrenze
# exakter Wert des Integrals I
Ifexakt = (np.exp(b) - b - np.exp(a) + a) / (np.exp(1.) - 1.)

#############################################################
#                                                           #
# Exakter Wert des Integrals von phi(x)                     #
# phi(x) ist die approx. durch eine Taylor-Reihe von f(x)   #
# Iphiexakt = None                                          #
#                                                           #
#############################################################

##########################################################################
#                                                                        #
# Funktionsdefinitionen                                                  #
# f = lambda x: (np.exp(x) - 1.)/(np.exp(1.) - 1.)                       #
# phi = lambda x: .../(np.exp(1.) - 1.)                                  #
# fminusphi = lambda x: ...                                              #
#                                                                        #
##########################################################################


phi = lambda x: (x + 1./2. * x**2 + 1./6. * x**3) / (np.exp(1.) - 1.)
f = lambda x: (np.exp(x) - 1.)/(np.exp(1.) - 1.)

# Groesse der Konvergenzstudie N = [Nmin,...,Nmax]
Nmin = 2
Nmax = 500
N = range(Nmin, Nmax + 1)

# Speicher Reservierung fuer die numerischen Quadratur Resultate I
res_simp = np.zeros((Nmax - Nmin + 1, 1))  # Resultate der Simpson-Regel
res_mc = np.zeros((Nmax - Nmin + 1, 1))  # Resultate der Monte-Carlo Quadratur
# Resultate der MC Quad. + Varianzkontrolle
res_mcvc = np.zeros((Nmax - Nmin + 1, 1))

# Speicher Reservierung fuer die Fehler approx. I - exaktem I
err_simp = np.zeros((Nmax - Nmin + 1, 1)) + 1.e-17  # Fehler der Simpson-Regel
# Fehler der Monte-Carlo Quadratur
err_mc = np.zeros((Nmax - Nmin + 1, 1)) + 1.e-17
# Fehler der MC Quad. + Varianzkontrolle
err_mcvc = np.zeros((Nmax - Nmin + 1, 1)) + 1.e-17


# Konvergenzstudie (integriere fuer verschiedene N und berechne die Fehler)
i = 0
for n in xrange(Nmin,Nmax + 1):
    res_simp[i] = simpson(f, a, b, n)
    res_mc[i] = mc(f, a, b, n)
    res_mcvc[i] = mc(phi, a, b, n)
    i += 1

err_simp = abs(res_simp - Ifexakt)
err_mc = abs(res_mc - Ifexakt)
err_mcvc = abs(res_mcvc - Ifexakt)

# Plot der Konvergenzstudie
plt.loglog(N, err_simp, 'bo', label='Simpson ')
plt.loglog(N, err_mc, 'ro', label='Monte-Carlo')
plt.loglog(N, err_mcvc, 'go', label='Monte-Carlo mit Varianzkontrolle')
plt.xlabel('N')
plt.ylabel('error')
plt.legend(loc=3)
plt.savefig('konvergenzstudie.pdf')


#############################################################
#                                                           #
# Berechnung der Konvergenzraten (verwenden Sie np.polyfit) #
#                                                           #
#############################################################

p_simp = np.polyfit(np.log(N[-20:]),np.log(err_simp[-20:]),1)
p_mc = np.polyfit(np.log(N[-20:]),np.log(err_mc[-20:]),1)
p_mcvc = np.polyfit(np.log(N[-20:]),np.log(err_mcvc[-20:]),1)

# Ausgabe der Konvergenzraten
print'Simpson                        Quadratur Konvergenzrate: {}'.format(-p_simp[0][0])
print'Monte-Carlo                    Quadratur Konvergenzrate: {}'.format(-p_mc[0][0])
print'Monte-Carlo + Varianzkontrolle Quadratur Konvergenzrate: {}'.format(-p_mcvc[0][0])
