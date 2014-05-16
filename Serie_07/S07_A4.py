from numpy import *
from matplotlib.pyplot import *
from S07_E1c import interp_barycentric, barycentric_weights
from chebychev_Template import chebpts, interp_chebychev_barycentric



if __name__ == '__main__':

	# -----------------------------------------------------
	# Aequidistante Stuetzstellen
	f = lambda x: 1./(1.+x**4)
	ns = array([5,10,19])
	a = -5
	b = 5

	figure()
	for n in ns:
		maxdiff = 0.
		x = linspace(a,b,n)
		xx = linspace(a,b,1000)
		y = f(x)
		lambdas = barycentric_weights(x)
		pxx = interp_barycentric(x,y,lambdas,xx)

		semilogy(xx, abs(pxx - f(xx)), label='n = %d' % n, lw=0.5)
		diff = max(abs(pxx - f(xx)))
		if diff > maxdiff: maxdiff = diff
		print 'Der maximale Fehler bei %d Aequdistanten Stuetzstellen ist: ' %n , maxdiff
	xlabel('x')
	ylabel('Abs. Error')

	ax = gca()
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
	ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

	grid(True)
	savefig('plot-aequidistant.pdf')
	show()

	#-----------------------------------------------------------
	# Chebychev-Abzissen
	ns = array([5,10,19])
	a = -5
	b = 5

	figure()
	for n in ns:
		maxdiff = 0.
		x = chebpts(a,b,n)
		xx = linspace(a,b,1000)
		y = f(x)
		lambdas = barycentric_weights(x)
		pxx = interp_chebychev_barycentric(x,y,xx)

		semilogy(xx, abs(pxx - f(xx)), label='n = %d' % n, lw=0.5)
		diff = max(abs(pxx - f(xx)))
		if diff > maxdiff: maxdiff = diff
		print 'Der maximale Fehler mit %d Chebychev-Abzissen ist: ' %n , maxdiff
	xlabel('x')
	ylabel('Abs. Error')

	ax = gca()
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
	ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

	grid(True)
	savefig('plot-chebychev.pdf')
	show()

