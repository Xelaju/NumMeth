from numpy import *

# kleine Änderung an der Datei

def VorRek(y0,N):
    S = ones(N+1)
    S[0] = y0
    for x in xrange(N):
        S[x+1] = 1 - (x+1) * S[x]
    return S

def RueckRek(y0,N):
    S = ones(N+1)
    S[N] = y0
    for x in xrange(N):
        S[(N-1) - x] = (1 - S[N-x]) * 1./(N-x)
    return S 

if __name__ == '__main__':
    y0 = 1 - 1 / exp(1)
    V21 = VorRek(y0, 21)
    V41 = VorRek(y0, 41)
    V121 = VorRek(y0, 121)
    R21 = RueckRek(0, 21)
    R41 = RueckRek(0, 41)
    R121 = RueckRek(0, 121)
    print 'Startwert: y0 = ', y0
    print 'VorRek y0,...,yN fuer N = 21 ist: ', V21
    print 'VorRek y0,...,yN fuer N = 41 ist: ', V41
    print 'VorRek y0,...,yN fuer N = 121 ist: ', V121
    print 'RuekRek y0,...,yN fuer N = 21 ist: ', R21
    print 'RuekRek y0,...,yN fuer N = 41 ist: ', R41
    print 'RuekRek y0,...,yN fuer N = 121 ist: ', R121



#hopefully someone finds this, SOS!!