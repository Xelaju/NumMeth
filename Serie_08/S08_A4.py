import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':

	f = lambda x: np.exp(-x/300)*(np.sin(2*np.pi*x/4.8)+np.sin(2*np.pi*x/4.5))

	N = 2e6
	x = np.arange(N)
	#x = np.linspace(0,N-1,N)
	data = f(x)

	# ------------------------------------- #
	# power - spectrum
	power_spec = np.abs(np.fft.fft(data))**2
	freqs = np.fft.fftfreq(data.size, 1. / N)

	# ------------------------------------- #
	# add noise
	data_noise = data + np.random.uniform(-3./2.,3./2.,N)

	# ------------------------------------- #
	# power - spec with noise
	power_spec_noise = np.abs(np.fft.fft(data_noise))**2

	# ------------------------------------- #
	# data noise time exp(-x/200)
	data_noise_e = data_noise * np.exp(-x/200)
	power_spec_noise_e = np.abs(np.fft.fft(data_noise_e))**2

	#########################################
	# Plot
	#########################################

	plt.figure()
	plt.semilogy(freqs, power_spec)

	plt.grid(True)
	plt.xlabel(r'$\lambda$')
	plt.ylabel(r'$\log_{10} ( \mid c_k\mid^2 )$')
	plt.title('Power Spectrum')
	#plt.xlim(-1e4, 1e4)
	#plt.ylim()
	plt.show()

	# --------------------------------------- #
	plt.plot(x,data)

	plt.grid(True)
	plt.xlabel(r'$x$')
	plt.ylabel(r'$e^{-\frac{x}{300}}(\sin{\frac{2\pi x}{4.8}}+\sin{\frac{2\pi x}{4.5}})$')
	plt.xlim(0,1e3)
	plt.show()

	# --------------------------------------- #
	plt.plot(x,data_noise)

	plt.grid(True)
	plt.xlabel(r'$x$')
	plt.ylabel(r'$e^{-\frac{x}{300}}(\sin{\frac{2\pi x}{4.8}}+\sin{\frac{2\pi x}{4.5}}) + Noise$')
	plt.xlim(0,1e3)

	plt.show()

	# --------------------------------------- #
	plt.semilogy(freqs, power_spec_noise)

	plt.grid(True)
	plt.xlabel(r'$\lambda$')
	plt.ylabel(r'$\log_{10} ( \mid c_k\mid^2 )$')
	plt.title('Power Spectrum with Noise')
	#plt.xlim(-1e4, 1e4)
	#plt.ylim()
	plt.show()

	# --------------------------------------- #
	plt.semilogy(freqs, power_spec_noise_e)

	plt.grid(True)
	plt.xlabel(r'$\lambda$')
	plt.ylabel(r'$\log_{10} ( \mid c_k\mid^2 )$')
	plt.title(r'Power Spectrum with Noise and mult. by $e^{-\frac{x}{200}}$')
	#plt.xlim(-1e4, 1e4)
	#plt.ylim()
	plt.show()
