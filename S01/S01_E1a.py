from numpy import *
from matplotlib.pyplot import *

s = linspace(0, 1, 200)
x = exp(-s)
y = s

plot (s, x, label=r'$\psi_'+str(x)+'(x)$')
plot (s,y)

grid(True)
xlim(0,1)
savefig('plot_S01_E1a.png')
