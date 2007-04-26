from numpy import size, linspace, zeros, zeros_like


def klassRK(f, y0, T, N):
    r"""
    Wendet auf eine System erster Ordnung der Form y'=f(t,y) das klassische
    Runge-Kutta-Verfahren an

    INPUT:
    f   - lambda function mit Argumenten (t,y)
    y0  - Vektor der Anfangswerte zum Zeitpunkt t0
    T   - T=[t0 tEnd] zu approximierendes Intervall
    N   - Anzahl der Iterationsschritte

    OUTPUT
    t   - Diskretisierung des Zeitintervalls
    y   - Matrix der approximierten Funktionswerte
    """
    # Anzahl der Gleichungen = Dimension der Loesung
    d = size(y0)
    # Zeitwerte und Schrittweite
    t, h = linspace(T[0], T[1], N + 1, retstep=True)

    y = zeros((d, size(t)))      # Loesungsmatrix
    y[:, 0] = y0                 # Anfangswert

    # Loesung berechnen
    for i in xrange(N):

        #################################################
        #                                               #
        # Implementiere die Runge-Kutta Iteration hier. #
        #                                               #
        # Hinweis: Benutze die Funtkion rkStep.         #
        # y[:,i+1] =                                    #
        #                                               #
        #################################################
        y[:,i+1] = rkStep(f, y[:,i], t[i], h) 

    return t, y


def rkStep(f, y0, t0, h):
    # Berechnet einen RK-Schritt der Schrittweite h

    ###########################################################
    #                                                         #
    # Implementiere einen einzelnen Runge-Kutta Schritt hier. #
    #                                                         #
    ###########################################################
    K1 = f(t0, y0)
    K2 = f(t0 + h/2, y0 + h * K1 / 2)
    K3 = f(t0 + h/2, y0 + h * K2 / 2)
    K4 = f(t0 + h, y0 + h * K3)

    return y0 + h * (K1/6 + K2/3 + K3/3 + K4/6)

if __name__ == '__main__':

    from numpy import array, dot, exp, cos, shape

    A = array([[0,1],[-101,-2]])
    f = lambda t, y: dot(A,y)

    y0 = array([1,-1])
    T = array([0,3])

    t, y = klassRK(f, y0, T, 10000)

    print shape(y)

    print 'Approximative Loesung mit klass. Runge-Kutta: y(3) = ', y[0,-1]
    print 'Refenrenzloesung: y(3) = ', exp(-3) * cos(10 * 3)
