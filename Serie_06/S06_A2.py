from scipy import zeros, dot, random, mat, linalg, diag, sqrt, sum, hstack, ones, exp
from scipy.linalg import norm, eig, expm, expm2, expm3, inv
import time
from numpy import *

def expD(A,y0,t):
	EW, S = eig(-1.j * A)
	Sinv = inv(S)
	D = diag(EW)
	E = dot(Sinv.T,dot(expm(D*t),S.T))
	return dot(E,y0)

def expM(A,y0,t):
	B = -1.j * A
	return dot(expm(B*t),y0)

def arnoldi(A,v0,k):
	V = mat(zeros((v0.shape[0],k+1), dtype=complex))
	V[:,0] = v0.copy()/norm(v0)
	H = mat(zeros((k+1,k), dtype=complex))

	for m in xrange(k):
		w = dot(A,V[:,m])
		for j in xrange(m+1):
			H[j,m] = dot(V[:,j].H, w)[0,0]	# Warum muss ich hier den ersten Eintrag nehmen, sollte dies nicht ohnehin ein Skalar sein?
											# ohne erhaelt man: TypeError: can't convert complex to float
			w -= dot(H[j,m],V[:,j])
		H[m+1,m] = norm(w)
		V[:,m+1] = w.copy()/H[m+1,m]
	return V, H

def lanscoz(A,v0,k):
	V = mat(zeros((v0.shape[0],k+1), dtype=complex))
	V[:,0] = v0.copy()/norm(v0) # Erstellung von q1
	alpha = zeros(k,complex)
	beta = zeros(k+1,complex)

	for m in xrange(k):
		w = dot(A, V[:,m])  			
		if m > 0: 						
			w -= dot(beta[m],V[:,m-1])	
		alpha[m] = dot(V[:,m].H,w)[0,0] 
		w -= dot(alpha[m],V[:,m])		
		beta[m+1] = norm(w)				
		V[:,m+1] = w.copy()/beta[m+1]	
	
	T = diag(alpha) + diag(beta[1:-1],1) + diag(beta[1:-1],-1) 
	return V, T	


def expA(A,v0,t,k=100):
	V, H = arnoldi(A,v0,k)
	#print 'H: \n', H[0:9,0:9]
	B = -1.j * H[:-1,:]
	E = expm(B * t)
	y = dot(V[:,:-1],dot(E,dot(V[:,:-1].H,v0)))
	return y

def expL(A,v0,t,k=100):
	V, T = lanscoz(A,v0,k) 	# Leider funktioniert das Lanscoz Verfahren sehr schlecht, wobei 
							# die ersten Werte sehr gut sind.
	#print 'T: \n', T[0:9,0:9]
	B = -1.j * T
	E = expm(B * t)
	y = dot(V[:,:-1],dot(E,dot(V[:,:-1].H,v0)))
	return y







