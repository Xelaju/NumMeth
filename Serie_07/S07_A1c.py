from numpy import *

def barycentric_weights(x):
	lambdas = zeros(size(x))
	for j, k in enumerate(x):
		P = 1.
		for i in x:
			if i != k:
				P *= (i - k) # Problem mit den Vorzeichen
		lambdas[j] = 1./P
	return lambdas

def interp_barycentric(x,y,lambdas,xx):
	p = zeros(len(xx))
	eps = 1.e-5
	for i, s in enumerate(xx):
		P = 0.; Q = 0.
		for k, c in enumerate(x):
			q = lambdas[k]/(s+eps - c)
			Q += q
			P += q * y[k]
			p[i] = P/Q
	return p

def barycentric_weights_add(old_x,old_lambdas,new_x):
	# Test falls Duplikate dabei sind
	itr = 0
	for xx in new_x:
		for x in old_x:
			if xx == x: new_x = delete(new_x,itr); itr -= 1
		itr += 1

	n_old = size(old_x)
	n_new = size(new_x)
	n_tot = n_old + n_new

	new_lambdas = zeros(n_tot)

	for k, l in enumerate(old_lambdas):
		L = l
		for i in new_x:
			L /= (i - old_x[k])
		new_lambdas[k] = L

	
	for i, nx in enumerate(new_x):
		P = 1.
		for j in xrange(n_tot):
			if j < n_old:
				P *= (nx - old_x[j])
			if j >= n_old:
				s = j - n_old
				if nx != new_x[s]:					
					P *= (nx - new_x[s])
		new_lambdas[i + n_old] = 1./P

	return new_lambdas


