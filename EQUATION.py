'''

Defines the class EQUATION which will be used as an object in the system of NLSE we want to solve.

'''

import numpy as np

class EQUATION:
	def __init__(self,name = '',ground_state = None):
		self.n_terms = 0 #Number of terms in the equation
		self.terms = {} #Dictionary of terms in the equation
		self.name = name #Name of the equation (optional)
		self.L = [] #Part to be solved in the momentum representation
		self.N = [] #Part to be solved in the position representation
		self.V = [] #Binding potential of the equation
		self.solution = ground_state #Solution of the equation

	def __str__(self): #String conversion
		terms_str = '\n'
		for i in self.terms:
			terms_str += str(self.terms[i]) + '\n'
		return '--EQUATION-- \nName: %s \nNº of Terms: %s \nTerms: %s' %(self.name,self.n_terms,terms_str)

	def add_term(self,term): #Adds a term to the equation
		self.terms[self.n_terms] = term  
		self.n_terms += 1

	def parts(self): #Updates the L and N parts of the equation
		momentum_matrix = [] #Matrix with the momentum representation parts
		position_matrix = [] #Matrix with the position representation parts
		for i in self.terms: #Iterate through the terms
			term = self.terms[i]
			if term.time_variant: #If term is time variant, recalculate its value
				term.matrix = term.function(**term.variables).astype(np.float64)
			if term.representation == 'Momentum':
				if momentum_matrix == []:
					momentum_matrix = np.copy(term.matrix)
				else:
					momentum_matrix += term.matrix
			else:
				if position_matrix == []:
					position_matrix = np.copy(term.matrix)
				else:
					position_matrix += term.matrix
				if term.name == 'Binding Potential': #If there is a binding potential to represent in the figures
					self.V = term.matrix
		self.L = momentum_matrix
		self.N = position_matrix

	def step(self,dt): #Performs a step of the SSSFM
		psi = self.solution

		#Half step in L
		psi = np.fft.fft2(psi)
		psi = np.exp(0.5j*self.L*dt)*psi

		#Full step in N
		psi = np.fft.ifft2(psi)
		psi = np.exp(1.0j*self.N*dt)*psi

		#Half step in L
		psi = np.fft.fft2(psi)
		psi = np.exp(0.5j*self.L*dt)*psi

		psi = np.fft.ifft2(psi)

		self.solution = psi