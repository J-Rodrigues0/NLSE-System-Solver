'''

Tests the solver for a single Schrodinger equation of the type H psi = E psi with H = P²/2 + V

'''

from EQUATION import *
from SYSTEM import *
from TERM import *

import scipy as sc
import numpy as np
from scipy import special

def V(x,y,V0,R0):
	V_ = np.zeros((len(x),len(y)))
	index = np.where(x**2 + y**2 > R0**2)
	V_[index] = np.sqrt(x[index]**2 + y[index]**2)
	return V_

def V_harmonic(x,y):
	return x**2/2

def psi_0(X,Y):
	return np.exp(-(((X-0.)**2 + (Y-0.)**2)/0.1**2))*np.exp(1.j*100.*X)

def Dx(sigma,dh):
	return 0.5/dh*(np.roll(sigma,1,axis=0) - np.roll(sigma,-1,axis=0))

def Dy(sigma,dh):
	return 0.5/dh*(np.roll(sigma,1,axis=1) - np.roll(sigma,-1,axis=1))

def D2y(sigma,dh):
	return (np.roll(sigma,1,axis=1) - 2*sigma + np.roll(sigma,-1,axis=1))/dh**2

def sigma_x(x,y):
	sig_x = np.zeros((len(x),len(y)))
	index = np.where(abs(x) >= 0.8)
	sig_x[index] = 10*x[index]**2
	return sig_x

def sigma_y(x,y):
	sig_y = np.zeros((len(x),len(y)))
	index = np.where(abs(y) >= 0.8)
	sig_y[index] = 10*y[index]**2
	return sig_y

Z = 8
N = 2**Z

L = 1.
dh = 2.*L/(N-1)

X, Y = np.mgrid[-L:L+dh:dh,-L:L+dh:dh]


k_line = np.fft.fftfreq(N,dh)
k_x, k_y = np.meshgrid(k_line,k_line)
K2 = k_x**2 + k_y**2

V0 = 1e5
R0 = 0.8

dt = 1e-4
N_stride = 100
N_steps = 100

#Initializing system and equation
system = SYSTEM()
eq1 = EQUATION('eq1',psi_0(X,Y))

system.add_equation(eq1) #Adding equation to the system

#Binding Potential Term
term2 = TERM(V(X,Y,V0,R0),'Position','Binding Potential')
#eq1.add_term(term2)

#Gross-Pitaevskii Term
f = lambda eq: np.abs(eq.solution)**2
kwargs = {'Function': f,'Variables':{'eq':eq1}}

term3 = TERM(f(eq1),'Position','Psi squared',True,**kwargs)
eq1.add_term(term3)

#Kinetic Energy Term
sig_x = sigma_x(X,Y)
sig_y = sigma_y(X,Y)

dx_sig_x = Dx(sig_x,dh)
dy_sig_y = Dy(sig_y,dh)

V_ = V(X,Y,V0,R0)

term4 = TERM(0.5*K2,'Momentum','P_squared')
eq1.add_term(term4)

term5 = TERM(1.j*sig_x,'Position','Sigma x')
eq1.add_term(term5)

phi = lambda eq: sig_x*(V_ + np.abs(eq.solution)**2 + D2y(eq.solution,dh))
kwargs = {'Function':phi,'Variables':{'eq':eq1}}
term6 = TERM(phi(eq1),'Position','Phi',False,True,**kwargs)
eq1.add_term(term6)

system.solve(X,Y,dt,N_stride,N_steps)