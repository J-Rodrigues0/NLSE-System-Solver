'''

Defines the class TERM which is a term in an EQUATION object in the system of NLSE equations we want to solve.

The TERM has an optional name, a defining matrix and the representation in which it should be solved 
(either positions or momentums).

'''

class TERM:
	def __init__(self,matrix,representation,name = '',time_variant = False,**kwargs):
		self.matrix = matrix #Matrix representation of the term
		self.representation = representation #Representation of the term. Either 'Momentum' or 'Position'
		self.name = name #Name of the term (optional)
		self.time_variant = time_variant #If the term if time variant

		if time_variant:			
			self.function = kwargs['Function'] #Function to call to update the term if time variant
			self.variables = kwargs['Variables'] #kwargs to pass to self.function to update the term

	def __str__(self): #String conversion
		return '--TERM-- \nName: %s \nRepresentation: %s' %(self.name,self.representation)
