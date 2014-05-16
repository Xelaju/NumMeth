from numpy import *
from numpy.linalg import cond, norm
from matplotlib.pyplot import *

e = linspace(1.e-5,1,500)

Al = []; ATAl = []; B1l = []; Balphal = []

for eps in e:
    # A
    A = array([[1+eps, 1],[1-eps,1],[eps,eps]])
    Al.append(cond(A))
    # ATA
    AT = transpose(A)
    ATA = dot(AT,A)
    ATAl.append(cond(ATA))
    # B1
    alpha = eps*norm(A)/sqrt(2)
    I = eye(shape(A)[0])
    alphaI = -alpha*I
    B1 = vstack([hstack([I,A]),hstack([AT,zeros((2,2))])])
    B1l.append(cond(B1))
    # Balpha
    Balpha = vstack([hstack([alphaI,A]),hstack([AT,zeros((2,2))])])
    Balphal.append(cond(Balpha))
    
figure()
semilogy(e,Al,label=r"$A$")
semilogy(e,ATAl,label=r"$A^TA$")
semilogy(e,B1l,label=r"$B_1$")
semilogy(e,Balphal,label=r"$B_\alpha$") # Ich bekomme eine Fehlermeldung, wenn ich $B_{\alpha}$ schreibe 
xlabel(r"$\varepsilon$")
ylabel(r"Konditionszahl")
title(r"Konditionszahl der Normalengleichung")
legend(loc="upper right")
grid(True)
savefig("S04_E5.eps")





