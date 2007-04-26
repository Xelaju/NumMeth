from numpy import *
from scipy.linalg import solve
from S07_E1c import barycentric_weights, interp_barycentric, barycentric_weights_add
from spectral_method_Template import chebyshev_pts, bary_eval
from matplotlib.pyplot import *

from scipy.interpolate import BarycentricInterpolator as bi
from pseudospectral import diffmat

def A_lag(x,lambdas):
	# lagrange koeffizienten unter-ueber der Diagonalen d.h i != j
	n = size(x)
	A = zeros((n,n))

	for i in xrange(n):
		for j in xrange(n):
			# Erstellung der Koeff. unter/ueber der Diagonalen
			if j != i:
				R = 0.
				for k in xrange(n):
					if k != i:
						R += (lambdas[k]/lambdas[i])/(x[i]-x[k])
				R -= 1./(x[i]-x[j])
				# Im gesamten Schritt wird eine Zeile von A bis auf den Diagonaleintrag erstellt
				A[i,j] = -2 * (lambdas[j]/lambdas[i])/(x[i]-x[j]) * R
				# Erstellung der Hauptdiagonalen
				A[i,i] -= A[i,j]
	return A

def l_check(x,l,j,i):
	n = size(x)
	j = j-1
	i = i-1
	P = 0.
	for k in xrange(n):
		if k != i:
			P += l[k]/(l[i]*(x[i]-x[k]))

	R = -2 * l[j]/(l[i]*(x[i]-x[j])) * (P-1./(x[i]-x[j]))

	return R




if __name__ == '__main__':


	ex0 = lambda a: {'f': lambda x: x ** a,
                    'u': lambda x: (x ** (a + 2) - x) / ((a + 1) * (a + 2))}

   	ex1 = lambda a: {'f': lambda x: exp(x) * sin(a * x),
                    'u': lambda x: (2 * a - 2 * a * x + 2 * a * exp(1) * x * cos(a) - 2 * a * exp(x) * cos(a * x) - exp(1) * x * sin(a) 
                                    + a ** 2 * exp(1) * x * sin(a) + exp(x) * sin(a * x) - a ** 2 * exp(x) * sin(a * x)) / (1. + a ** 2) ** 2}

	#####################
	# Chebychev-Abzissen
	#####################
	x_cheb = chebyshev_pts(10,0,1)

	##################
	# Referenzloesung aus scipy
	##################
	L = bi(x_cheb)

	#####################
	# Funktionen
	#####################
	#f = ex0(a=10)['f']
	#u = ex0(a=10)['u']
	f = ex1(a=10)['f']
	u = ex1(a=10)['u']
	f_cheb = f(x_cheb)
	u_cheb = u(x_cheb)
	L.set_yi(f_cheb)

	#######################
	# Auswertung der Baryzentrischen Gewichte
	########################
	xx = linspace(0,10,1000)
	lambdas = barycentric_weights(x_cheb)
	lambdas_check = L.wi
	print 'B-Gewichte Referenzwerte: \n', lambdas_check
	print 'B-Gewichte berechnete Werte: \n', lambdas

	######################
	# Interpolation von f zum Test
	######################
	f_int_reference = L(xx)
	f_int_own = interp_barycentric(x_cheb,f_cheb,lambdas,xx)
	f_int = bary_eval(f_cheb,x_cheb,lambdas,xx)
	f_exact = f(xx)

	###################
	# Interpolation von u
	##################
	I_ref = np.identity(10)
	D_ref = diffmat(x_cheb) + I_ref
	D2_ref = dot(D_ref,D_ref) + I_ref

	t = lambda x: x**4/4
	

	A_ref = D_ref
	A_int = A_lag(x_cheb,lambdas)

	print A_int
	print A_ref


	y_ref = solve(A_ref,t(x_cheb))
	u_ref = interp_barycentric(x_cheb,y_ref,lambdas,xx)
	y_int = solve(A_int,f_cheb)
	u_int = interp_barycentric(x_cheb,y_int,lambdas,xx)
	u_exact = u(xx)
	
	##################
	# Terminal Ausgabe der Werte
	###############
	

	#################
	# Plot
	################
	#plot(xx,f_int_reference)
	#plot(xx,f_int_own)
	#plot(xx,f_exact,label=r'f(x)_exact')
	#plot(xx,u_exact,label=r'u(x)_exact')
	#plot(xx,u_int,label=r'u(x)_int')
	plot(xx,u_ref,label=r'u(x)_ref')
	legend(loc='lower right')
	grid(True)
	savefig('S07_E2c.eps')
	show()
	