from numpy import *
from numpy.linalg import norm
from matplotlib.pyplot import *

def init_pos_vel(d, n, box_width):
    x = linspace(0, box_width, n)
    tx = meshgrid(*[x]*d)
    pos = zip(*[xs.flatten() for xs in tx])
    vel = zeros_like(pos)
    return array(pos), array(vel)

def eval_DV(i, x, DV):
    """
    Berechnet die Kraft auf das i-te Teilchen
    Keyword Arguments:
    i  --
    x  -- Positionen
    DV -- -1.0*Grad U
    """
    val = 0
    for j in xrange(len(x)):
        if j==i:
            continue
        val += DV(x[i,:]-x[j,:])
    return val


def compute_energy(x, v, V):
    """
    Keyword Arguments:
    x -- positions
    v -- velocity
    V -- Potential
    """
    m = 1.
    Vpot = 0.

    for i in xrange(shape(x)[0]):
        for j in xrange(shape(x)[0]):
            if j == i:
                continue
            Vpot += V(norm(x[i,:] - x[j,:]))

    Vpot /= 2.
    Vkin = sum(m * sum(abs(v)**2, axis=1) / 2.)

    return Vpot+Vkin, Vpot, Vkin

# Als Beispiel
class StoermerVerlet:
    def init(self, x, v, DV, dt):
        return x, v

    def finalize(self, x, v, t):
        return x, v, t

    def run(self, x, v, DV, dt):
        """
        Keyword Arguments:
        x  -- pos
        v  -- vel
        V  --
        dt --
        """
        n = len(x)
        vk2 = zeros_like(v)
        xn = zeros_like(x)
        vn = zeros_like(v)
        # update positions
        for i in xrange(n):
            # Forces acting on x[i]
            Fij = eval_DV(i, x, DV)
            vk2[i,:] = v[i,:] + 0.5*dt*Fij
            xn[i,:] = x[i,:] + dt*vk2[i,:]

        # update velocities
        for i in xrange(n):
            Fij = eval_DV(i, xn, DV)
            vn[i,:] = vk2[i,:] + 0.5*dt*Fij
        return xn, vn

class LeapFrog:
    def run(self, x, v, DV, dt):
        """
        Keyword Arguments:
        x  --
        v  --
        DV --
        dt --
        """
        n = len(x)
        xn = zeros_like(x)
        vn = zeros_like(v)

        for i in xrange(n):
            vn[i,:] = v[i,:] + dt * eval_DV(i, x, DV)
            xn[i,:] = x[i,:] + dt * vn[i,:]

        return xn, vn

    # Achtung hier muessen DV, dt uebergeben werden, somit auch in der run() Funktion
    def init(self, x, v, DV, dt):
        vn = zeros_like(v)
        for i in xrange(len(x)):
            Fij = eval_DV(i, x, DV)
            vn[i] = v[i] - 0.5*dt*Fij
        return x, vn

    def finalize(self, x, v, t):
        x = array(x); v = array(v); t = array(t)
        v = 0.5*(v[1:,:,:] + v[:-1,:,:]) # Fehler im Template ? x und v sind Listen und keine Arrays
        x = x[:-1,:,:]
        t = t[:-1]
        return x, v, t

class VelocityVerlet:
    def init(self, x, v, DV, dt):
        return x, v

    def finalize(self, x, v, t):
        return x, v, t

    def run(self, x, v, DV, dt):
        """
        Keyword Arguments:
        x  -- pos
        v  -- vel
        V  --
        dt --
        """

        n = len(x)
        xn = zeros_like(x)
        vn = zeros_like(v)

        # Ausser dass ich Zeilen gespart habe, ist dies identisch mit dem StoermerVerlet Verfahren
        for i in xrange(n):
            xn[i,:] = x[i,:] + dt * v[i,:] + dt**2 * eval_DV(i, x, DV) * 0.5
        for i in xrange(n):
            vn[i,:] = v[i,:] + dt * (eval_DV(i, x, DV) + eval_DV(i, xn, DV)) * 0.5

        return xn, vn

def run(x, v, V,  T=1, dt=1e-3, P=StoermerVerlet()):
    """
    Keyword Arguments:
    x  -- initial positions
    v  -- initial velocities
    V  -- V(x), scalar potential
    P  -- time Propagator
    T  -- (default 1)
    dt -- (default 1e-3)
    """
    x, v = P.init(x, v, V, dt)
    xhist = [x]
    vhist = [v]
    thist = [0]
    nsteps = int(T/float(dt))

    for i in xrange(nsteps):
        x,v = P.run(x, v, V, dt)
        xhist.append(x)
        vhist.append(v)
        thist.append((i+1)*dt)
    xhist, vhist, thist = P.finalize(xhist, vhist, thist)

    return array(xhist), array(vhist), array(thist)

epsilon = 1.
sigma = 1.
LJ = lambda r: epsilon* ((sigma/r)**12 - (sigma/r)**6)
DLJ = lambda x: 6*epsilon*(2*x*(sigma/norm(x))**14 - x*(sigma/norm(x))**8)


def small_perturbation(d=2, N=5, box_width=5, T=50, dt=1e-2, Integrator=StoermerVerlet()):
    """
    Keyword Arguments:
    d         -- (default 2)
    N         -- (default 5)
    box_width -- (default 5)
    T         -- (default 50)
    dt        -- (default 1e-2)
    """
    phi = 2*pi*random.rand(N*N)
    vref =  1e-6 * array([cos(phi), sin(phi)]).T
    x, v = init_pos_vel(d, N, box_width)

    v = vref
    xhist, vhist, thist = run(x, v, DLJ, T, dt, Integrator)

    data = []
    for i,ti in enumerate(thist):
        E, Vpot, Vkin = compute_energy(xhist[i,:,:], vhist[i,:,:], LJ)
        data.append([E, Vpot, Vkin])

    data = array(data)
    figure()
    subplot(121)
    title('v1')
    plot(thist, data[:,2], label='T')
    plot(thist, data[:,1], label='U')
    plot(thist, data[:,0], label='E')

    grid(True)
    xlabel('$t$')
    ylabel('Energy')
    title('Original Problem')
    legend()

    # Perturbed experiment
    x, v = init_pos_vel(d, N, box_width)
    v = vref
    v[0,:] += 1e-10
    xhistp, vhistp, thistp = run(x, v, DLJ, T, dt, Integrator)

    datap = []
    for i,ti in enumerate(thistp):
        E, Vpot, Vkin = compute_energy(xhistp[i,:,:], vhistp[i,:,:], LJ)
        datap.append([E, Vpot, Vkin])

    datap = array(datap)
    subplot(122)
    title('Perturbed Problem')
    plot(thist, datap[:,2], label='T')
    plot(thist, datap[:,1], label='U')
    plot(thist, datap[:,0], label='E')

    grid(True)
    xlabel('$t$')
    ylabel('Energy')
    title('Perturbed Problem')
    legend()
    grid(True)
    tight_layout()
    savefig('perturbation-energies.pdf')
    print('--- `perturbation-energies.pdf` written to disk. ---')

    figure()
    colors = ['r', 'g', 'b', 'm', 'k']*100
    for i in xrange(xhist.shape[1]):
        plot(xhist[0,i,0], xhist[0,i,1],'o', color=colors[i])
        plot(xhist[-1,i,0], xhist[-1,i,1],'s', color=colors[i])
        plot(xhistp[-1,i,0], xhistp[-1,i,1],'d', color=colors[i])
        plot(xhist[:,i,0], xhist[:,i,1], color=colors[i], label='original')
        plot(xhistp[:,i,0], xhistp[:,i,1], '--', color=colors[i], label='perturbed')
    axis('equal')
    title('Trajektorien')
    xlabel('$x$')
    ylabel('$y$')
    grid(True)
    legend( [Line2D([],[], linestyle='--', color='k'),
             Line2D([],[], linestyle='-', color='k'),
             Line2D([],[], linestyle='', marker='o', color='k'),
             Line2D([],[], linestyle='', marker='s', color='k'),
             Line2D([],[], linestyle='', marker='d', color='k')],
            ["perturbed",
             "original",
             "$x(t=0)$",
             r'$x(t=t_{end})$',
             r'$x_{\epsilon}(t=t_{end})$'],
            numpoints=1)
    savefig('Trajektorien.pdf')
    print('--- `Trajektorien.pdf` written to disk. ---')

    figure()
    di=  [sqrt(sum((xhistp[i]-xhist[i])**2)) for i in xrange(len(thist))]
    vi = [sqrt(sum((vhistp[i]-vhist[i])**2)) for i in xrange(len(thist))]
    semilogy(thist, di, label=r'$|| x - x_{\epsilon}||$')
    semilogy(thist, vi, label=r'$|| v- v_{\epsilon}||$')
    xlabel('$t$')
    ylabel('Abweichung')
    grid(True)
    title('Phasenraumdistanz')
    legend(loc='lower left')
    savefig('phasenraumdistanz.pdf')
    print('--- `phasenraumdistanz.pdf` written to disk. ---')

    figure()
    semilogy(thist, abs(datap[:,0]-data[:,0]))
    ylabel(r'$|E-E_{\epsilon}|$')
    xlabel(r'$t$')
    grid(True)
    savefig('energie-differenz.pdf')
    print('--- `energie-differenz.pdf` written to disk. ---')

    show()
    savez('data.npz', xhist, vhist, thist, xhistp, vhistp, thistp, data, datap, x, v, vref)

if __name__ == '__main__':
    # dimension
    d = 2
    # Anz. Teilchen in jeder Dimension, 10 x ... x 10
    N = 3
    box_width = 3
    T = 60
    dt = 0.05

    # Anfangsbedingungen
    # x, v sind N x d Arrays, d.h. 1 Teilchen pro Zeile
    x, v = init_pos_vel(d, N, box_width)

    ######################################################
    #                                                    #
    # Hinweis ersetze P durch deinen eigenen Integrator: #
    # P = VelocityVerlet()                               #
    # oder                                               #
    # P = LeapFrog()                                     #    
    #                                                    #
    ######################################################
    xhist, vhist, thist = run(x, v, DLJ, T=T, dt=dt, P=LeapFrog())
    data = []
    for i,ti in enumerate(thist):
        E, Vpot, Vkin = compute_energy(xhist[i,:,:], vhist[i,:,:], LJ)
        data.append([E, Vpot, Vkin])


    data = array(data)
    figure(figsize=(8,8))
    plot(thist, data[:,2], label='T')
    plot(thist, data[:,1], label='U')
    plot(thist, data[:,0], label='E')

    grid(True)
    xlabel('$t$')
    ylabel('Energy')
    legend()
    savefig('energy_%d.pdf' % d)
    show()
