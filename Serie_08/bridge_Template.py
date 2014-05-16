from numpy import *
from numpy.linalg import eigvals
from matplotlib.pyplot import *
import matplotlib.patches as mpatches
import string

def bridge_read(FileIn):

    # oeffne Datei
    f = open(FileIn, 'r')

    # lese header
    f.readline()

    # lese Daten
    t = []
    x = []
    for line in f.readlines():
        tmp = map(float, line.split())
        t.append(tmp[0])
        x.append(tmp[1:])

    # schliesse Datei
    f.close()

    # konvertiere zu numpy arrays
    t = array(t)
    x = array(x)
    x = transpose(x)

    # return
    return t, x

def get_matrix():
    # Anfangs-Gelenk Koordinaten
    theta = 0.25*pi
    r = array([1.0, #  0,x_1
                  1.0, #  1,y_1
                  2.0, #  2,x_2
                  1.0, #  3,y_2
                  3.0, #  4,x_3
                  1.0, #  5,y_3
                  4.0, #  6,x_4
                  1.0, #  7,y_4
                  1.0, #  8,x_5
                  0.0, #  9,y_5
                  2.0, # 10,x_6
                  0.0, # 11,y_6
                  3.0, # 12,x_7
                  0.0, # 13,y_7
                  4.0, # 14,x_8
                  0.0  # 15,x_8
                 ])

    # Gleichgewichtsmatrix E
    c = cos(theta)
    s = sin(theta)
    #            (  1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 ,10 ,11 ,12 ,13 ,14 ,15 ,16 ,17 ,18 )
    E = array([[-1., 0., 0., c , 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.], #  1
                  [ 0., 0., 0., s , 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.], #  2
                  [ 1.,-1., 0., 0., 0., c , 0.,-c , 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.], #  3
                  [ 0., 0., 0., 0., 0., s , 1., s , 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.], #  4
                  [ 0., 1.,-1., 0., 0., 0., 0., 0., c , 0.,-c , 0., 0., 0., 0., 0., 0., 0.], #  5
                  [ 0., 0., 0., 0., 0., 0., 0., 0., s , 1., s , 0., 0., 0., 0., 0., 0., 0.], #  6
                  [ 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.,-c , 0., 0., 0., 0., 0.], #  7
                  [ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., s , 0., 0., 0., 0., 0.], #  8
                  [ 0., 0., 0., 0., 0.,-c , 0., 0., 0., 0., 0., 0., 0., 1.,-1., 0., 0., 0.], #  9
                  [ 0., 0., 0., 0.,-1.,-s , 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.], # 10
                  [ 0., 0., 0., 0., 0., 0., 0., 0.,-c , 0., 0., 0., 0., 0., 1.,-1., 0., 0.], # 11
                  [ 0., 0., 0., 0., 0., 0.,-1., 0.,-s , 0., 0., 0., 0., 0., 0., 0., 0., 0.], # 12
                  [ 0., 0., 0., 0., 0., 0., 0., c , 0., 0., 0., 0., 0., 0., 0., 1.,-1., 0.], # 13
                  [ 0., 0., 0., 0., 0., 0., 0.,-s , 0.,-1., 0., 0., 0., 0., 0., 0., 0., 0.], # 14
                  [ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., c , 0., 0., 0., 0., 0., 1.,-1.], # 15
                  [ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,-s ,-1., 0., 0., 0., 0., 0., 0.]  # 16
                 ])

    # aeussere Kraefte p
    m = 1.  # Gelenk Masse
    g = 10. # Schwerkraft
    p = array([  0.0, #  1
                  -m*g,  #  2
                    0.0, #  3
                  -m*g,  #  4
                    0.0, #  5
                  -m*g,  #  6
                    0.0, #  7
                  -m*g,  #  8
                    0.0, #  9
                  -m*g,  # 10
                    0.0, # 11
                  -m*g,  # 12
                    0.0, # 13
                  -m*g,  # 14
                    0.0, # 15
                  -m*g,  # 16
                 ])

    # Diagonal Matrix D^-1 = Di
    Di = diag([1./linalg.norm(array([r[ 0],r[ 1]]) - array([r[ 2],r[ 3]])), # Stab  1
                  1./linalg.norm(array([r[ 2],r[ 3]]) - array([r[ 4],r[ 5]])), # Stab  2
                  1./linalg.norm(array([r[ 4],r[ 5]]) - array([r[ 6],r[ 7]])), # Stab  3
                  1./linalg.norm(array([   0.,   0.]) - array([r[ 0],r[ 1]])), # Stab  4
                  1./linalg.norm(array([r[ 8],r[ 9]]) - array([r[ 0],r[ 1]])), # Stab  5
                  1./linalg.norm(array([r[ 8],r[ 9]]) - array([r[ 2],r[ 3]])), # Stab  6
                  1./linalg.norm(array([r[10],r[11]]) - array([r[ 2],r[ 3]])), # Stab  7
                  1./linalg.norm(array([r[ 2],r[ 3]]) - array([r[12],r[13]])), # Stab  8
                  1./linalg.norm(array([r[10],r[11]]) - array([r[ 4],r[ 5]])), # Stab  9
                  1./linalg.norm(array([r[12],r[13]]) - array([r[ 4],r[ 5]])), # Stab 10
                  1./linalg.norm(array([r[ 4],r[ 5]]) - array([r[14],r[15]])), # Stab 11
                  1./linalg.norm(array([r[14],r[15]]) - array([r[ 6],r[ 7]])), # Stab 12
                  1./linalg.norm(array([   5.,   0.]) - array([r[ 6],r[ 7]])), # Stab 13
                  1./linalg.norm(array([   0.,   0.]) - array([r[ 8],r[ 9]])), # Stab 14
                  1./linalg.norm(array([r[ 8],r[ 9]]) - array([r[10],r[11]])), # Stab 15
                  1./linalg.norm(array([r[10],r[11]]) - array([r[12],r[13]])), # Stab 16
                  1./linalg.norm(array([r[12],r[13]]) - array([r[14],r[15]])), # Stab 17
                  1./linalg.norm(array([   5.,   0.]) - array([r[14],r[15]]))  # Stab 18
                 ])

    # Matrix A austellen
    eta = 2500. # Elastizitaetsmodule
    # Stelle hier die Matrix A auf...
    # A ...
    A = eta*dot(E, dot(Di, transpose(E)))
    return A

if __name__ == '__main__':

    A = get_matrix()
    #
    # Berechne die Eigenwerte der Matrix #
    # ew = ...                           #
    #

    # Zeitintervall der Messung [0,T) mit T = 20 s
    T = 20.
    # lese die Messdaten
    t, x = bridge_read('bridge.dat')
    messdaten = x[1, :]

    #
    # Berechne das Energiespektrum der Messdaten #
    # power_2 = ...                              #
    #

    # Frequenzen omega_k = k/T
    freq = r_[-len(t) / 2:len(t) / 2] / T

    # plotte betraege der Fourier Koeffs. als Funktion der Eigenfrequenzen
    figure()
    #
    # Plotte das Energiespektrum         #
    # semilogy((2.*pi*freq)**2, power_2) #
    #

    grid(True)
    xlabel(r'$\lambda$')
    ylabel(r'$\log_{10} ( \mid c_k\mid^2 )$')
    xlim([0, 10000])

    # Plotte bei jedem Eigenwert v eine vertikale Linie
    # vlines(sorted(ew), 10**-3, 10**5)
    grid(True)
    savefig('power-spektrum.pdf')
    show()
