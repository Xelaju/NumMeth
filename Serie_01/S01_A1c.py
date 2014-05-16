from numpy import *
from matplotlib.pyplot import *

x = linspace(0, 1, 200)
psi = x + 1 - x * exp(x)
f = x

plot(x, f)
plot(x, psi)

grid(True)
savefig("plot_S01_E1c.png")
