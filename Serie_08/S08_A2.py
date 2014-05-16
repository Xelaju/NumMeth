from numpy import *
import numpy.linalg as lg
from bridge_Template import get_matrix, bridge_read
from matplotlib.pyplot import *
import matplotlib.patches as mpatches
import string

if __name__ == '__main__':

	A = get_matrix()
	# Eigenwete von A
	ew = lg.eig(A)[0]

	# Zeitintervall der Messung [0,T) mit T = 20 s
	T = 20.
	# lese die Messdaten
	t, x = bridge_read('bridge.dat')
	messdaten = x[1, :]

	# Energiespektrum der Messdaten
	power_2 = fft.fft(messdaten)**2
	# Frequenzen omega_k = k/T
	freq = r_[-len(t) / 2:len(t) / 2] / T

	figure()
	semilogy((2.*pi*freq)**2, power_2)

	grid(True)
	xlabel(r'$\lambda$')
	ylabel(r'$\log_{10} ( \mid c_k\mid^2 )$')
	xlim([0, 10000])

	vlines(sorted(ew), 10**-3, 10**5)
	grid(True)
	savefig('power-spektrum.pdf')
	show()
