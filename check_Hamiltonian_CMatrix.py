import numpy as np
from reading_CMatrix_dat_files import *

class Graph():  
	def __init__(self, vertices):  
		self.graph = [[0 for column in range(vertices)] 
							for row in range(vertices)]  
		self.V = vertices  
  
	 # Check if this vertex is an adjacent vertex of the previously added vertex and is not included in the path earlier 
	def isSafe(self, v, pos, path):  
		# Check if current vertex and last vertex in path are adjacent  
		if self.graph[ path[pos-1] ][v] == 0:  
			return False
		# Check if current vertex not already in path  
		for vertex in path:  
			if vertex == v:  
				return False
  
		return True
  
	# Define recursive function here  
	def hamCycleUtil(self, path, pos):  
  
		# base case: if all vertices are included in the path  
		if pos == self.V:  
			# Last vertex must be adjacent to the first vertex in path to make a cyle  
			if self.graph[ path[pos-1] ][ path[0] ] >= 1: # Might need to change here so that I can use values other than 1 to denote connectivity 
				# need to put condition greater then equal to one since the last list could be connected to the first more than once
				return True
			else:  
				return False
  
		# Try different vertices as a next candidate in Hamiltonian Cycle. We don't try for 0 as we included 0 as starting point in hamCycle()  
		for v in range(1,self.V):  
			if self.isSafe(v, pos, path) == True:
				path[pos] = v  
				if self.hamCycleUtil(path, pos+1) == True:  
					return True
				# Remove current vertex if it doesn't lead to a solution  
				path[pos] = -1
  
		return False
  
	def hamCycle(self):  
		path = [-1] * self.V  
  		# Can start with zero since undirected graph. Can start anywhere since it is cycle
		path[0] = 0
  
		if self.hamCycleUtil(path,1) == False:
			# print ("Solution does not exist\n")
			soln_writer = open('hamsolno.txt','w')
			soln_writer.writelines(('Solution does not exist\n'))
			soln_writer.writelines(str(np.array(CMatrix))) 
			return False
  
		self.printSolution(path)  
		return True
  
	def printSolution(self, path):  
		# print ("Solution Exists: Following is one Hamiltonian Cycle") 
		soln_writer = open('hamsolyes.txt', 'w') 
		soln_writer.writelines(('Solution Exists: Following is one Hamiltonian Cycle\n'))
		for vertex in path:  
			# print (vertex + 1, end = " ")
			soln_writer.writelines(str(vertex + 1) + " ")
		# print (path[0]+1, "\n")
		soln_writer.writelines('\n')
		# print(CMatrix)
		f = open('CMatrix.dat','wb')
		dill.dump(CMatrix, f)
		f.close()
		soln_writer.writelines(str(np.array(CMatrix)))




# fetch CMatrix in current folder
CMatrix = get_CMatrix()
g1 = Graph(len(CMatrix[0]))
g1.graph = CMatrix
g1.hamCycle()

# check truth value of Hamiltonian cycle
print(g1.hamCycle())