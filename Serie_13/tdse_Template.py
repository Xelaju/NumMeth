from numpy import *
from numpy.fft import  fft, ifft
from numpy.linalg import norm
from matplotlib.pyplot import *
from scipy.linalg import expm


def simulation(potential, initialvalue):
    scale = pi

    # Time interval, step size
    tend = 5.0*sqrt(ieps)
    nl = 1
    n = int(tend * 10**nl)
    print("Doing %d timesteps" % n)
    h = tend/n
    ih = 1j*h

    # place for norms
    normsol = zeros(n, dtype=floating)
    cinE = zeros(n)
    potE = zeros(n)

    # Space discretization
    l = 12
    N = 2**l

    # Laplacian
    x0 = arange(-N/2,N/2)
    x = zeros(N)
    x[:N/2] = x0[N/2:]
    x[N/2:] = x0[:N/2]
    A = eps * x**2 * 0.5

    #  Potential
    x = scale * arange(-1, 1, 2.0/N)
    V = ieps * potential(x) # Potential am Punkt x

    v = initialvalue(x) # v(x,0) Zum Zeitpunkt t = 0
    u = fft(v)          # u(x,0) Im Frequenzbereich
    normsol[0] = sqrt(2*pi) * norm(u)/N
    print("norm: %f" % normsol[0])

    # Exponentials
    eD = exp(- ih * A / 2)
    eV = exp(- ih * V)


    cinE[0] = real(eps * dot(conj(u), (A*u))/N**2 * (2.0*pi))
    # Erklaerung: 
    # - eps hebt multiplikation mit ieps auf
    # - 2.0 * pi um w anstatt f (frequenz) zu erhalten
    # - V*v ist potE im Zeitbereich
    # - N ist fuer die Skalierung weil fuer u und fft(V*v) erhaelt man N**2
    # - conj(u) gibt den Realteil (Was ist mir der Skalierung?)
    potE[0] = real(eps * dot(conj(u), fft(V*v)/N**2) * (2.0*pi))
    print("E kin: %f" % cinE[0])
    print("E pot: %f" % potE[0])
    print("E tot: %f" % (cinE[0]+potE[0]))

    plot(x, v, V, 0, h)

    # Time stepping
    for k in range(1, n):
        print("Timestep: %d" % k)
        # Propagate

        u = fft(v)
        dv = ifft(eD * fft(eV * ifft(eD * u)))
        v = dv

        # Compute norm
        u = fft(v)
        normsol[k] = sqrt(2*pi) * norm(u)/N
        # Compute energies
        cinE[k] = real(eps * dot(conj(u), (A*u))/N**2 * (2.0*pi))
        potE[k] = real(eps * dot(conj(u), fft(V*v)/N**2) * (2.0*pi))

        print("E kin: %f" % cinE[k])
        print("E pot: %f" % potE[k])
        print("E tot: %f" % (cinE[k]+potE[k]))

        plot(x, v, V, k, h)

    fig = figure()
    ax = fig.gca()
    ax.plot(h*arange(n), normsol)
    ax.grid(True)
    ax.set_xlabel(r"Time $t$")
    ax.set_ylabel(r"$\| u \|$")
    fig.savefig("norm.png")
    close(fig)

    fig = figure()
    ax = fig.gca()
    ax.plot(h*arange(n), cinE, label=r"$E_{kin}$")
    ax.plot(h*arange(n), potE, label=r"$E_{pot}$")
    ax.plot(h*arange(n), cinE+potE, label=r"$E_{tot}$")
    ax.grid(True)
    ax.legend(loc="best")
    ax.set_xlabel(r"Time $t$")
    ax.set_ylabel(r"Energy $E$")
    fig.savefig("energies.png")
    close(fig)


def plot(x, v, V, k, h, ymin=-3, ymax=3):
    fig = figure()
    ax = fig.gca()
    ax.fill_between(x, V*eps, ymin, color="k", alpha=0.2)
    ax.plot(x, V*eps,"k-")
    ax.plot(x, real(v), "b-")
    ax.plot(x, imag(v), "g-")
    ax.plot(x, abs(v), "r-")
    ax.set_title("Time t=%f" % k*h)
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(ymin, ymax)
    fig.savefig("solution_at_timestep %04d.png" % k)
    close(fig)


if __name__ == "__main__":
    # Parameter
    eps = 0.01
    ieps = 1.0/eps

    # Potentials
    def morse(x, c = -1.0, De = 8.0):
        V = De*(1-exp(-(x-c)/sqrt(2*De)))**2
        return V

    def harmonic(x):
        V = 0.5*(x**2)
        return V

    #  Initial value
    g0 = lambda x: (ieps/pi)**(0.25) * exp(-(0.5*ieps)*(x+0.5)**2) * exp(-1j*x*ieps)
    g1 = lambda x: (ieps/pi)**(0.25) * exp(-(0.5*ieps)*x**2)

    simulation(morse, g1)
