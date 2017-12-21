'''

Defines the class EQUATION which will be used as an object in the system of NLSE we want to solve.

'''

import numpy as np
import pylab as pl
import matplotlib.pyplot as plt

class EQUATION:
	def __init__(self,ground_state = None,name = None):
		self.n_terms = 0
		self.terms = {}
		self.name = name
		self.L = None #Part to be solved in the momentum representation
		self.N = None #Part to be solved in the position representation
		self.solution = ground_state

	def __str__(self):
		terms_str = '\n'
		for i in self.terms:
			terms_str += str(self.terms[i]) + '\n'
		return '--EQUATION-- \nName: %s \nNº of Terms: %s \nTerms: %s' %(self.name,self.n_terms,terms_str)

	def add_term(self,term):
		self.terms[self.n_terms] = term  
		self.n_terms += 1

	def parts(self):
		for i in self.terms:
			term = self.terms[i]
			if term.representation == 'Momentum':
				if self.L == None:
					self.L = term.matrix
				else:
					self.L += term.matrix
			else:
				if self.N == None:
					self.N = term.matrix
				else:
					self.N += term.matrix

	def solve(self,X,Y,dt,N_stride,N_steps):
		self.parts()

		psi0 = self.solution
		psi = psi0
		psi = np.fft.fft2(psi)
		psi = np.exp(0.5j*self.L*dt)*psi

		t = 0
		for n in range(N_stride):
			print(n)
			#Iteration through full steps
			for m in range(N_steps):
				t += dt
				psi = np.fft.ifft2(psi)
				
				psi = np.exp(1.j*self.N*dt)*psi

				psi = np.fft.fft2(psi)
				psi = np.exp(1.j*self.L*dt)*psi

			#Last full step in N
			psi = np.fft.ifft2(psi)
			psi = np.exp(1.j*self.N*dt)*psi

			#Last half step in L
			psi = np.fft.fft2(psi)
			psi = np.exp(0.5j*self.L*dt)*psi

			psi = np.fft.ifft2(psi)

			prob = np.abs(psi)**2

			#Make the bonecade
			pl.figure()
			pl.contour(X,Y,self.N.real,10)
			pl.contourf(X,Y,prob)
			pl.colorbar()
			plt.savefig('/home/joaorodrigues/Desktop/MCE/Final/figs/%s.png' %n,dpi=200)
			pl.close()

			#Continue the run
			psi = np.fft.fft2(psi)
			psi = np.exp(0.5j*self.L*dt)*psi

			self.solution = psi

