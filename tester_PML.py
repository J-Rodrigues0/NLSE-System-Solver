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
	V_ = np.zeros((len(x),len(y))) + V0
	index = np.where(x**2 + y**2 < R0**2)
	V_[index] = 0
	return V_

def V_harmonic(x,y):
	return x**2/2

def psi_0(X,Y):
	return np.exp(-(((X-0.)**2 + (Y-0.)**2)/0.1**2))*np.exp(1.j*0.*X)

def Dx(sigma,dh):
	return 0.5*dh*(np.roll(sigma,1,axis=0) - np.roll(sigma,-1,axis=0))

def Dy(sigma,dh):
	return 0.5*dh*(np.roll(sigma,1,axis=1) - np.roll(sigma,-1,axis=1))

def sigma_x(x,y):
	sig_x = np.zeros((len(x),len(y)))
	index = np.where(x >= 0.9)
	sig_x[index] = -1e10
	return sig_x

def sigma_y(x,y):
	sig_y = np.zeros((len(x),len(y)))
	index = np.where(y >= 0.9)
	sig_y[index] = -1e10
	return sig_y

Z = 6
N = 2**Z

L = 1.
dh = 2.*L/(N-1)

X, Y = np.mgrid[-L:L+dh:dh,-L:L+dh:dh]


k_line = np.fft.fftfreq(N,dh)
k_x, k_y = np.meshgrid(k_line,k_line)
K2 = k_x**2 + k_y**2

V0 = 1e10
R0 = 0.4

dt = 1e-4
N_stride = 100
N_steps = 100

#Initializing system and equation
system = SYSTEM()
eq1 = EQUATION('eq1',psi_0(X,Y))

system.add_equation(eq1) #Adding equation to the system

#Binding Potential Term
term2 = TERM(V(X,Y,V0,R0),'Position','Binding Potential')
eq1.add_term(term2)

#Gross-Pitaevskii Term
f = lambda eq: np.abs(eq.solution)**2
kwargs = {'Function': f,'Variables':{'eq':eq1}}

term3 = TERM(f(eq1),'Position','Psi squared',True,**kwargs)
eq1.add_term(term3)

#Non time-derivative terms (NTD)
sig_x = sigma_x(X,Y)
sig_y = sigma_y(X,Y)

dx_sig_x = Dx(sig_x,dh)
dy_sig_y = Dy(sig_y,dh)

term4_0 = TERM(0.5j*(sig_x + sig_y),'Position','NTD1')
eq1.add_term(term4_0)

term4_1 = TERM(-0.25*sig_y*dx_sig_x*k_x + -0.25*sig_x*dy_sig_y*k_y,'Momentum','NTD2')
eq1.add_term(term4_1)

#Time derivative terms (TD)

phi0 = lambda eq: 0.25j*sig_x*sig_y*eq.solution
kwargs0 = {'Function':phi0,'Variables':{'eq':eq1}}
term5_0 = TERM(phi0(eq1),'Position','Phi0',False,True,**kwargs0)
eq1.add_term(term5_0)

#Terms in V
V_ = V(X,Y,V0,R0)
phi1 = lambda eq: 0.5*(sig_x + sig_y)*V_*eq.solution
kwargs1 = {'Function':phi1,'Variables':{'eq':eq1}}
term5_1 = TERM(phi1(eq1),'Position','Phi1',False,True,**kwargs1)
eq1.add_term(term5_1)

phi2 = lambda eq: 0.25*sig_x*sig_y*V_*eq.solution
kwargs2 = {'auxiliary':True,'Function':phi2,'Variables':{'eq':eq1}}
term5_2 = TERM(phi2(eq1),'Position','Phi2',False,True,**kwargs2)
eq1.add_term(term5_2) #Auxiliary term

phi3 = lambda: term5_2.matrix
kwargs3 = {'Function':phi3,'Variables':{}}
term5_3 = TERM(phi3(),'Position','Phi3',False,True,**kwargs3)
eq1.add_term(term5_3)

#Terms in psi²
phi4 = lambda eq: 0.5*(sig_x + sig_y)*np.abs(eq.solution)**2*eq.solution
kwargs4 = {'Function':phi4,'Variables':{'eq':eq1}}
term5_4 = TERM(phi4(eq1),'Position','Phi1',False,True,**kwargs4)
eq1.add_term(term5_4)

phi5 = lambda eq: 0.25*sig_x*sig_y*np.abs(eq.solution)**2*eq.solution
kwargs5 = {'auxiliary':True,'Function':phi5,'Variables':{'eq':eq1}}
term5_5 = TERM(phi5(eq1),'Position','Phi2',False,True,**kwargs5)
eq1.add_term(term5_5)

phi6 = lambda: term5_5.matrix
kwargs6 = {'Function':phi6,'Variables':{}}
term5_6 = TERM(phi6(),'Position','Phi3',False,True,**kwargs6)
eq1.add_term(term5_6)

#Time derivative term (??)

system.solve(X,Y,dt,N_stride,N_steps)