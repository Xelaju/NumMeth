from numpy import array, cos, sin, shape, min, double, zeros
import matplotlib.pyplot as plt
from ode45 import ode45
import scipy.optimize


def PendulumODE(y, l, g):
    """PendulumODE return the right-hand side of the math. pendulum"""
    dydt = array([y[1], -g * sin(y[0]) / l])
    return dydt


def IntegratePendulum(phi0, tEnd=1.8, l=0.6, g=9.81, flag=False):
    """IntegratePendulum solve the mathematical pendulum with ode45
    flag == false   -->  y = complete solution phi,phi'
    flag == true    -->  y = phi(tEnd)
    (order of the 2 outputs reversed wrt the usual to use fzero)
    """
    if shape(phi0) == (1,):
        phi0 = phi0[0]
    y0 = array([phi0, 0])
    tspan = (0, tEnd)
    t, y = ode45(lambda t, y: PendulumODE(y, l, g), tspan, y0)
    if flag:
        return y[-1][0]
    else:
        return (y, t)


def IntegratePendulumEE(phi0, N, tEnd=1.8, l=0.6, g=9.81, flag=False):
    """IntegratePendulum solve the mathematical pendulum with Explicit Euler
    flag == false   -->  y = complete solution phi,phi'
    flag == true    -->  y = phi(tEnd)
    """
    if shape(phi0) == (1,):
        phi0 = phi0[0]
    y0 = array([phi0, 0])

    # stepsize
    h = double(tEnd) / N
    # memory allocation
    y = zeros((N + 1, 2))
    y[0, :] = y0
    t = zeros(N + 1)

    #################################
    #                               #
    # Implement explicit Euler here #
    #                               #
    #################################

    d2phi = - g/l * sin(phi)
    dphi = g/l * cos(phi) - g/l

    for i in xrange(N+1):
        y[i+1,:] = y[]

    if flag:
        return y[-1][0]
    else:
        return (y, t)


def IntegratePendulumIE(phi0, N, tEnd=1.8, l=0.6, g=9.81, flag=False):
    """IntegratePendulum solve the mathematical pendulum with Implicit Euler
    flag == false   -->  y = complete solution phi,phi'
    flag == true    -->  y = phi(tEnd)
    """
    if shape(phi0) == (1,):
        phi0 = phi0[0]
    y0 = array([phi0, 0])

    # stepsize
    h = double(tEnd) / N
    # memory allocation
    y = zeros((N + 1, 2))
    y[0, :] = y0
    t = zeros(N + 1)

    #################################
    #                               #
    # Implement implicit Euler here #
    #                               #
    #################################

    if flag:
        return y[-1][0]
    else:
        return (y, t)


def IntegratePendulumIM(phi0, N, tEnd=1.8, l=0.6, g=9.81, flag=False):
    """IntegratePendulum solve the mathematical pendulum with Implicit Midpoint
    flag == false   -->  y = complete solution phi,phi'
    flag == true    -->  y = phi(tEnd)
    """
    if shape(phi0) == (1,):
        phi0 = phi0[0]
    y0 = array([phi0, 0])

    # stepsize
    h = double(tEnd) / N
    # memory allocation
    y = zeros((N + 1, 2))
    y[0, :] = y0
    t = zeros(N + 1)

    ####################################
    #                                  #
    # Implement implicit midpoint here #
    #                                  #
    ####################################

    if flag:
        return y[-1][0]
    else:
        return (y, t)


if __name__ == '__main__':
    from time import time

    starttime = time()
    tEnd = 4  # 2*1.8
    N = 500
    l = 0.6
    g = 9.81
    phiT = 1.48551280723

    # ODE45
    yODE45, tODE45 = IntegratePendulum(phiT, tEnd, l, g, False)
    # Explicit Euler
    y_EE, t_EE = IntegratePendulumEE(phiT, N, tEnd, l, g, False)
    # Implicit Euler
    y_IE, t_IE = IntegratePendulumIE(phiT, N, tEnd, l, g, False)
    # Implicit Midpoint
    y_IM, t_IM = IntegratePendulumIM(phiT, N, tEnd, l, g, False)

    plt.figure()
    plt.plot(yODE45.transpose()[0], yODE45.transpose()[1], 'g-', label='ODE45')
    plt.plot(y_EE.transpose()[0], y_EE.transpose()[1], 'r--', label='EE')
    plt.plot(y_IE.transpose()[0], y_IE.transpose()[1], 'b--', label='IE')
    plt.plot(y_IM.transpose()[0], y_IM.transpose()[1], 'c--', label='IM')
    #ax = axis;
    # plot([ax(1) ax(2)],[0 0],'k-');
    # plot([0 0],[ax(3) ax(4)],'k-');
    plt.xlabel('alpha')
    plt.ylabel('p')
    plt.legend()
    plt.axis('equal')
    plt.savefig('solutions.pdf')

    # Tracking total energy

    # Explicit Euler
    y = y_EE

    E_kin = 0.5 * (y[:, 1] ** 2)
    E_pot = -g / l * cos(y[:, 0])
    E_pot = E_pot - min(E_pot) + min(E_kin)
    E_tot = E_kin + E_pot

    plt.figure()
    plt.plot(t_EE, E_kin.transpose(), 'b-', label='kinetic energy')
    plt.plot(t_EE, E_pot.transpose(), 'r-', label='potential energy')
    plt.plot(t_EE, E_tot.transpose(), 'g-', label='total energy')
    plt.xlabel('time t')
    plt.ylabel('energy')
    plt.legend()
    plt.savefig('energieEE.pdf')

    # Implicit Euler
    y = y_IE

    E_kin = 0.5 * (y[:, 1] ** 2)
    E_pot = -g / l * cos(y[:, 0])
    E_pot = E_pot - min(E_pot) + min(E_kin)
    E_tot = E_kin + E_pot

    plt.figure()
    plt.plot(t_EE, E_kin.transpose(), 'b-', label='kinetic energy')
    plt.plot(t_EE, E_pot.transpose(), 'r-', label='potential energy')
    plt.plot(t_EE, E_tot.transpose(), 'g-', label='total energy')
    plt.xlabel('time t')
    plt.ylabel('energy')
    plt.legend()
    plt.savefig('energieIE.pdf')

    # Implicit Midpoint
    y = y_IM

    E_kin = 0.5 * (y[:, 1] ** 2)
    E_pot = -g / l * (cos(y[:, 0]))
    E_pot = E_pot - min(E_pot) + min(E_kin)
    E_tot = E_kin + E_pot

    plt.figure()
    plt.plot(t_EE, E_kin.transpose(), 'b-', label='kinetic energy')
    plt.plot(t_EE, E_pot.transpose(), 'r-', label='potential energy')
    plt.plot(t_EE, E_tot.transpose(), 'g-', label='total energy')
    plt.xlabel('time t')
    plt.ylabel('energy')
    plt.legend()
    plt.savefig('energieIM.pdf')
