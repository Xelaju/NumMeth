# -*- encoding: utf-8 -*-
from numpy import *
from numpy.linalg import norm
from matplotlib.pyplot import *
from time import time

def init_pos_vel(d, n, box_width):
    x = linspace(0, box_width, n)
    tx = meshgrid(*[x]*d)
    pos = zip(*[xs.flatten() for xs in tx])
    vel = zeros_like(pos)
    return array(pos), array(vel)


epsilon = 1.
sigma = 1.
LJ = lambda r: epsilon* ((sigma/r)**12 - (sigma/r)**6)
DLJ = lambda x: 6*epsilon*(2*x*(sigma/norm(x))**14 - x*(sigma/norm(x))**8)

def compute_energy(q, p, V):
    """
    Keyword Arguments:
    q -- positions
    p -- momentum
    V -- Potential
    """
    Vpot = 0.
    Vkin = 0.
    for i in xrange(shape(q)[0]):
        Vkin += 0.5 * dot(p[i,:],p[i,:])
        for j in xrange(shape(q)[0]):
            if j == i:
                continue
            Vpot += 0.5 * V(norm(q[i,:] - q[j,:]))

    return Vpot+Vkin, Vpot, Vkin


def force_wrapper(f):    
    def F(q, DLJ):
        force = zeros_like(q)
        n = len(q)
        for i in xrange(n):
            force[i,:] = f(q, DLJ, i)
        return force
    return F

#@  force_wrapper erzeugt aus eval_DV die Funktion,
#@  die die Berechnung für jedes Teilchen durchführt
#@  (mit @force_wrapper wird sie als sogenannter Decorator verwendet; die
#@  Definition ist praktisch äquivalent zu
#@  eval_DV = force_wrapper(lambda x, DV, i: <eigentliche Funktion hier>)
#@  Sehr praktisch, um Code lesbar zu halten.

# Mehr infos unter: http://www.artima.com/weblogs/viewpost.jsp?thread=240808
#                   http://www.artima.com/weblogs/viewpost.jsp?thread=240845
#                   http://www.artima.com/weblogs/viewpost.jsp?thread=241209

@force_wrapper
def eval_DV(q, DLJ, i=None):
    """
    Berechnet die Kraft auf das i-te Teilchen
    Keyword Arguments:
    i  -- 
    x  -- Positionen
    DLJ -- -1.0*Grad U
    """
    DV = 0
    for j in xrange(len(q)):
        if j==i:
            continue
        DV += DLJ(q[i,:] - q[j,:])
    return DV


# --------- splitting parts --------------
F = lambda q: eval_DV(q, DLJ)

def PhiV(dt, u): # PhiA
	t = u[2] + dt
	q = u[0].copy()
	p = u[1] + dt * F(q) # das (-) ist durch DLJ absobiert
	return array([q,p,t])

def PhiT(dt, u): # PhiB
	t = u[2].copy()
	p = u[1].copy()
	q = u[0] + dt * p 
	return array([q,p,t])



if __name__ == '__main__':
    # dimension
    d = 2
    # Anz. Teilchen in jeder Dimension, 10 x ... x 10
    N = 3
    box_width = 3
    T = 60
    dt = 0.05

    # Anfangsbedingungen
    # u, p sind N x d Arrays, d.h. 1 Teilchen pro Zeile
    q, p = init_pos_vel(d, N, box_width)

    #-------------------------------------------
    from SplittingParameters import *
    S = SplittingParameters()
    y0 = q
    y0p = p
    tspan = [0., T]
    N = 10000
    u0 = array([y0,y0p,0.])
    ys = u0.copy()
    Methods = ['LT', 'S2', 'Y4', 'KL6']
    
    figure(figsize=(8,8))
    for method in Methods:
		a,b = S.build(method)
		print method ,' integrate ...'
		starttime = time()
		tt, uu = S.intsplit(PhiV, PhiT, a, b, tspan, N, ys, getall=True)
		t2 = time() - starttime
		print method ,' with', N ,' equal_timesteps in ', t2 ,' seconds'
		data = []
		for i in xrange(N+1):
			E, Vpot, Vkin = compute_energy(uu[:,0][i], uu[:,1][i], LJ)
			data.append([E, Vpot, Vkin])

		data = array(data)

		plot(tt, data[:,2], label='T')
		plot(tt, data[:,1], label='U')
		plot(tt, data[:,0], label='E')

		grid(True)
		xlabel('$t$')
		ylabel('Energy_method_%s' % method)
		legend()
		savefig('energy_method_%s.pdf' % method)
		show()




