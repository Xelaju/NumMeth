from numpy import *

def trapez(g, a, b, n):
	x = linspace(a,b,n); h = x[1] - x[0]
	I = sum(g(x[1:-1])) + 0.5 * (g(x[0]) + g(x[-1]))
	return I * h

def simpson(g, a, b, n):
	# ensure even number
	if n % 2 == 1: n = n + 1
	x = linspace(a, b, n+1); h = x[1] - x[0]
	I = h * sum(g(x[0:-2:2]) + 4 * g(x[1:-1:2]) + g(x[2::2]))/3.0
	return I

if __name__ == '__main__':
	import matplotlib.pyplot as plt
	from scipy import integrate

	g = lambda x: x**(3./5.) - x
	f = lambda u: 5 * (u**7 - u**9)
	left = 0.0; right = 1.0
	exact_g, e = integrate.quad(g, left, right)
	exact_f, e = integrate.quad(f, left, right)
	print exact_f
	print exact_g

	N = linspace(2,101,100)
	res_gs = array(N); res_gt = array(N)
	res_fs = array(N); res_ft = array(N)
	for i in xrange(size(N)):
		res_gs[i] = simpson(g, left, right, N[i]); res_gt[i] = trapez(g, left, right, N[i])
		res_fs[i] = simpson(f, left, right, N[i]); res_ft[i] = trapez(f, left, right, N[i])
	err_gs = abs(res_gs - exact_g); err_gt = abs(res_gt - exact_g)
	err_fs = abs(res_fs - exact_f); err_ft = abs(res_ft - exact_f)

	plt.loglog(N,err_gs,'o'); plt.loglog(N,err_gt,'o')
	plt.loglog(N,err_fs,'o'); plt.loglog(N,err_ft,'o')
	plt.show()

	p_gs = polyfit(log(N[-20:]),log(err_gs[-20:]),1)
	p_gt = polyfit(log(N[-20:]),log(err_gt[-20:]),1)

	p_fs = polyfit(log(N[-20:]),log(err_fs[-20:]),1)
	p_ft = polyfit(log(N[-20:]),log(err_ft[-20:]),1)

	print 'convergence order of simpson function g ', -p_gs[0]
	print 'convergence order of trapez function g ', -p_gt[0]

	print 'convergence order of simpson function f ', -p_fs[0]
	print 'convergence order of trapez function f ', -p_ft[0]