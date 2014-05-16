import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':

	s = lambda t, a: np.cos(150*t - a/(1+t)**2) * (np.exp(-12*t) - np.exp(-13*t) + 0.02 * np.exp(-10*(t-1./2.)**2))

	N = 2e6
	t = np.linspace(0,1,N) 

	# ------------------------------- #
	# power spectrum
	data_0 = s(t,a=0)
	data_60 = s(t,a=60)
	coeffs_0 = np.fft.fft(data_0)
	coeffs_60 = np.fft.fft(data_60)
	power_spec_0 = np.abs(coeffs_0)**2
	power_spec_60 = np.abs(coeffs_60)**2
	freqs = np.fft.fftfreq(data_0.size, 1. / N)

	# -------------------------------- #
	# auto-correlation
	auto_corr_0 = np.fft.ifft(coeffs_0 * np.conjugate(coeffs_0))
	auto_corr_60 = np.fft.ifft(coeffs_60 * np.conjugate(coeffs_60))

	# -------------------------------- #
	# cumulative auto-correlation
	cum_auto_corr_0 = np.cumsum(np.abs(auto_corr_0))
	cum_auto_corr_0 /= cum_auto_corr_0[-1]
	cum_auto_corr_60 = np.cumsum(np.abs(auto_corr_60))
	cum_auto_corr_60 /= cum_auto_corr_60[-1]

	d = np.linspace(0,10,N)

	#############################
	# Plot 
	#############################

	# -------------------------------- #
	# Plot der Signale
	plt.figure()
	plt.plot(t,data_0,label=r'$s(t,0)$')
	plt.plot(t,data_60,label=r'$s(t,60)$')

	plt.grid(True)
	plt.xlabel(r'$t$')
	plt.ylabel(r'$s(t,a) = \cos(150 t - \frac{a}{(1+t)^2}) (e^{-12t} - e^{-13t} + 0.02 e^{-10(t-\frac{1}{2})^2)}$')
	plt.title('Signals')
	plt.legend(loc='lower right')
	#plt.xlim(-1e4, 1e4)
	#plt.ylim()
	plt.show()

	# --------------------------------------- #
	plt.semilogy(freqs, power_spec_0,label=r'$E_{t,0}$')
	plt.semilogy(freqs, power_spec_60,label=r'$E_{t,60}$')

	plt.grid(True)
	plt.xlabel(r'$\lambda$')
	plt.ylabel(r'$E_{t,a} = \log_{10} ( \mid c_k\mid^2 )$')
	plt.title('Power Spectrum')
	plt.legend(loc='upper right')
	#plt.xlim(-1e4, 1e4)
	#plt.ylim()
	plt.show()

	# -------------------------------- #
	# Plot der auto-correlation
	plt.plot(t,auto_corr_0,label=r'$R(t)_0$')
	plt.plot(t,auto_corr_60,label=r'$R(t)_{60}$')

	plt.grid(True)
	plt.xlabel(r'$t$')
	plt.ylabel(r'$R(t)_a = (s(\tau,a) * s(\tau,a))(-t)$')
	plt.title('auto-correlation')
	plt.legend(loc='lower right')
	#plt.xlim(-1e4, 1e4)
	#plt.ylim()
	plt.show()

	# -------------------------------- #
	# Plot der auto-correlation
	plt.plot(d,cum_auto_corr_0,label=r'$R(t)_0$')
	plt.plot(d,cum_auto_corr_60,label=r'$R(t)_{60}$')

	plt.grid(True)
	plt.xlabel(r'$d$')
	plt.ylabel(r'$c_d = \frac{\sum_{j=0}^{d-1}|R_j|}{c_{N-1}}$')
	plt.title('cumulative-auto-correlation')
	plt.legend(loc='lower right')
	#plt.xlim(-1e4, 1e4)
	#plt.ylim()
	plt.show()


