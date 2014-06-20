import numpy as np
import matplotlib.pyplot as plt


def plot_coeffs(freqs, c):
	plt.semilogy(freqs, np.abs(c)**2)
	plt.xlim(-1e3, 1e3)
	plt.ylim(1e3, 1e8)
	plt.show()


if __name__ == '__main__':

 	f = lambda x: np.sin((x/200)**2) * np.exp(-((x-1000)/400)**2)

	s_rate = 2e5
	x = np.linspace(0,2000,s_rate)
	data = f(x)
	data += np.random.uniform(-1,1,s_rate)

	plt.figure()
	plt.plot(x,data)
	plt.show()

	coeffs = np.fft.fft(data)

	freqs = np.fft.fftfreq(data.size, 1. / s_rate)

	plot_coeffs(freqs, coeffs)

	peak = np.argmax(np.abs(coeffs)) # Return: index_array : ndarray of ints
	filtered_coeffs = np.zeros_like(coeffs)
	span = 3
	print peak
	if peak == 0:
		filtered_coeffs[0: span + 1] = coeffs[0 : span + 1]
		filtered_coeffs[-span : 0] = coeffs[-span : 0]
	else:
		filtered_coeffs[peak - span : peak + span + 1] = coeffs[peak - span : peak + span + 1] 
		filtered_coeffs[-peak - span : -peak + span + 1] = coeffs[-peak - span : -peak + span + 1]

	plot_coeffs(freqs,filtered_coeffs)

	filtered_data = np.real(np.fft.ifft(filtered_coeffs))
	plt.plot(x, filtered_data)
	plt.show()