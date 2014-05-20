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
        pass

    return t, y


def rkStep(f, y0, t0, h):
    # Berechnet einen RK-Schritt der Schrittweite h
    y1 = zeros_like(y0)

    ###########################################################
    #                                                         #
    # Implementiere einen einzelnen Runge-Kutta Schritt hier. #
    #                                                         #
    ###########################################################

    return y1
