from numpy import *
from numpy.linalg import det
from numpy.linalg import solve

l1 = 1./10.
l2 = 1./100.

A = array([[l1*exp(-l1*1),l2*exp(-l2*1)],[l1*exp(-l1*3),l2*exp(-l2*3)],[l1*exp(-l1*5),l2*exp(-l2*5)],[l1*exp(-l1*7),l2*exp(-l2*7)]])
print A
b = array([24.04, 20.64, 17.84, 15.53])
print b

ATA = dot(transpose(A),A)
print ATA
ATb = dot(transpose(A),b)
print ATb

x = solve(ATA,ATb)
print x
