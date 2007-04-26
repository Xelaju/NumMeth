from numpy import array, sin, pi, shape
import matplotlib.pyplot as plt
from ode45 import ode45
from scipy.optimize import *


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
    if shape(phi0)==(1,): phi0 = phi0[0]
    y0 = array([phi0 , 0])
    tspan= (0, tEnd)
    t ,y = ode45(lambda t,y: PendulumODE(y,l,g), tspan, y0)
    if flag: 
        return y[-1][0]
    else:
        return (y,t)


if __name__ == '__main__':
    from time import time

    starttime = time()

    tEnd = 1.8
    l = 0.6
    g = 9.81
    phiT=0
    # solve ODE for extreme values to check the period length:
    a1 = 0.8
    y1, t1 = IntegratePendulum(a1 * pi * 0.5, tEnd, l, g, False)
    a2 = 0.99
    y2, t2 = IntegratePendulum(a2 * pi * 0.5, tEnd, l, g, False)

    plt.figure()
    ax1 = plt.subplot(211)
    plt.plot(t1, y1.transpose()[0], 'b', label=r'$f \, | \, f_0 = %1.2f \, \frac{\pi}{2}$' % a1)
    plt.plot(t1, y1.transpose()[1], 'b--', label=r'$\frac{df}{dt} \, | \, f_0 = %1.2f \, \frac{\pi}{2}$' % a1)
    plt.plot(t2, y2.transpose()[0], 'r', label=r'$f \, | \, f_0 = %1.2f \, \frac{\pi}{2}$' % a2)
    plt.plot(t2, y2.transpose()[1], 'r--', label=r'$\frac{df}{dt} \, | \, f_0 = %1.2f \, \frac{\pi}{2}$' % a2)
    plt.title('Solutions to the pendulum ODE for different initial amplitudes')
    plt.xlabel('t')
    plt.legend(loc='upper left')
    plt.grid(True)


    # zero finding:
    def f(phi):
        return IntegratePendulum(phi, 0.45 , 0.6 , 9.81 ,True)
    
    phiT=fsolve(f,0.9*pi*0.5)
    

    #compute complete solution with period = tEnd for plotting
    y3, t3 = IntegratePendulum(phiT, tEnd, l, g, False)
    y3 = y3.transpose()

    ax1 = plt.subplot(212)
    plt.plot(t3, y3[0], '-', label=r'$f$')
    plt.plot(t3, y3[1], 'b--', label=r'$\frac{df}{dt}$')
    plt.title(r'Solution to the pendulum ODE for $f_0 = %1.3f $' % phiT)
    plt.xlabel('t')
    plt.legend(loc='upper left')
    plt.grid(True)
    plt.savefig('Pendulum.pdf')
    plt.show()

    t = time() - starttime
    print('time: %f sec' % round(t, 3))
