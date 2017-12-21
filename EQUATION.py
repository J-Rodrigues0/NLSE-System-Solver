'''

Defines the class EQUATION which will be used as an object in the system of NLSE we want to solve.

'''

import numpy as np

class EQUATION:
	def __init__(self,name = '',ground_state = None):
		self.n_terms = 0
		self.terms = {}
		self.name = name
		self.L = [] #Part to be solved in the momentum representation
		self.N = [] #Part to be solved in the position representation
		self.solution = ground_state

	def __str__(self):
		terms_str = '\n'
		for i in self.terms:
			terms_str += str(self.terms[i]) + '\n'
		return '--EQUATION-- \nName: %s \nNº of Terms: %s \nTerms: %s' %(self.name,self.n_terms,terms_str)

	def add_term(self,term):
		self.terms[self.n_terms] = term  
		self.n_terms += 1

	def parts(self): #Updates the L and N parts of the equation
		momentum_matrix = []
		position_matrix = []
		for i in self.terms:
			term = self.terms[i]
			if term.time_variant:
				term.matrix = term.function(**term.variables).astype(np.float64)
			if term.representation == 'Momentum':
				if momentum_matrix == []:
					momentum_matrix = term.matrix
				else:
					momentum_matrix += term.matrix
			else:
				if position_matrix == []:
					position_matrix = term.matrix
				else:
					position_matrix += term.matrix
		self.L = momentum_matrix
		self.N = position_matrix

	def step(self,dt):
		psi = self.solution

		#Half step in L
		psi = np.fft.fft2(psi)
		psi = np.exp(0.5j*self.L*dt)*psi

		#Full step in N
		psi = np.fft.ifft2(psi)
		psi = np.exp(1.j*self.N*dt)*psi

		#Half step in L
		psi = np.fft.fft2(psi)
		psi = np.exp(0.5j*self.L*dt)*psi

		psi = np.fft.ifft2(psi)

		self.solution = psi

