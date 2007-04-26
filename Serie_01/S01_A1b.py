from numpy import *
from matplotlib.pyplot import *

x = linspace(0, 1, 200) # range of function

psi = (1 + x) / (1 + exp(x))
f = x

plot(x, f)
plot(x, psi)

grid(True)
xlim(0, 1)

savefig("plot_S01_E1b.png")
#haha
