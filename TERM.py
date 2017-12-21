'''

Defines the class TERM which is a term in an EQUATION object in the system of NLSE equations we want to solve.

The TERM has an optional name, a defining matrix and the representation in which it should be solved 
(either positions or momentums).

'''

class TERM:
	def __init__(self,matrix,representation,name = '',time_variant = False,**kwargs):
		self.matrix = matrix
		self.representation = representation
		self.name = name
		self.time_variant = time_variant

		if time_variant:			
			self.function = kwargs['Function']
			self.variables = kwargs['Variables']

	def __str__(self):
		return '--TERM-- \nName: %s \nRepresentation: %s' %(self.name,self.representation)
