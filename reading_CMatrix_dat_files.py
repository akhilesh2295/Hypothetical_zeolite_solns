# Reading connectivity Matrices into python matrices
import dill

def get_CMatrix():
	f = open('CMatrix.dat', 'rb') 
	CMatrix = dill.load(f)
	# print(CMatrix)
	f.close()
	return CMatrix
# get_CMatrix()
