# from time import perf_counter
import os
# import subprocess
import shutil
import dill
import ast
import itertools
os.chdir(os.getcwd())
from Cif_to_Unitcell import *
from convert_xyz_to_csv import *
# from Make_matrix_general import *
from Opti_Make_matrix_general import *
import sqlite3
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
from math import gcd
import math
from functools import reduce
pd.options.mode.chained_assignment = None  # default='warn'
# t0= perf_counter()
sufficient_algorithmic_solutions = 10

# rect_box_size = (3,3,3)
np.set_printoptions(threshold=sys.maxsize)

solve_algo_approach = False
make_algo_plot = False
algo_plot_mod = False
gams_model_solve = True
gams_model_plot = True

global list_u
global MatConn
global starter_node
global max_Num

def find_gcd(list):
    # x = reduce(gcd, list)
    return reduce(gcd, list)

def plot_point_with_text(i):
    P = df_new.loc[df_new['node_num'] == i]
    # l = find_in_list_u(i)
    # threedee.scatter(P.x, P.y, P.z, c = (float(P.z) - df['z'].min())/(df['z'].max() - df['z'].min()), cmap = 'viridis', s = 70 )
    threedee.scatter(P.x, P.y, P.z, c = 'red', s = 70 )
    threedee.text(float(P.x), float(P.y), float(P.z), int(P.id)+1, size=14,zdir=None)
    # threedee.annotate( l, xyz = (float(P.x), float(P.y), float(P.z)), size=14,zdir=None)

# def find_in_list_u(target):
#     for i,lst in enumerate(list_u):
#         for j,color in enumerate(lst):
#             if color == target:
#                 return i
#                 # print(i)
#     return None

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=256):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

def plot_point(i):
    P = df_new.loc[df_new['node_num'] == i]
    # l = find_in_list_u(i)
    threedee.scatter(P.x, P.y, P.z, s = 70, c='red')
    # return threedee
    # threedee.text(float(P.x), float(P.y), float(P.z), l, size=14,zdir=None)

def plot_line(i,j):
    P1 = df_new.loc[df_new['node_num'] == i]
    P2 = df_new.loc[df_new['node_num'] == j]
    # P1 = SRU_df.loc[i]
    # P2 = SRU_df.loc[j]
    threedee.plot([float(P1.x),float(P2.x)],[float(P1.y),float(P2.y)],[float(P1.z),float(P2.z)], color = '0.8')
    # return threedee

def plot_point_gray(i):
    P = df_new.loc[df_new['node_num'] == i]
    # l = find_in_list_u(i)
    threedee.scatter(P.x, P.y, P.z, c = '0.7', alpha = 0.5, s = 70)
    # threedee.text(float(P.x), float(P.y), float(P.z), l, size=14,zdir=None)


def plot_line_gray(i,j):
    P1 = df_new.loc[df_new['node_num'] == i]
    P2 = df_new.loc[df_new['node_num'] == j]
    threedee.plot([float(P1.x),float(P2.x)],[float(P1.y),float(P2.y)],[float(P1.z),float(P2.z)],color = '0.9', alpha = 0.5)

def makeCMatrix(filename):
	f = open(filename, 'rb')
	global list_u
	list_u = dill.load(f)
	MatConn = dill.load(f)
	starter_node = dill.load(f)
	# max_Num = dill.load(f)
	f.close()
	global df_new
	global df
	df_file = open('df_new.dat','rb')
	df_new = dill.load(df_file)
	df = dill.load(df_file)
	df_file.close()
	node_array = df_new.node_num.to_numpy()
	CMatrix = np.zeros((len(list_u)-1,len(list_u)-1))
	# for i in range(1,len(list_u)):
	# 	for j in range(1,len(list_u)):
	# 		somelists = [list_u[i], list_u[j]]
	# 		for element in itertools.product(*somelists):
	# 			if MatConn[element[0]][element[1]]:
	# 				CMatrix[i-1][j-1] = 1
	# 				break
	for i in range(1,len(list_u)):
		for j in range(1,len(list_u)):
			# somelists = [list_u[i], list_u[j]]
			list1 = list_u[i]
			list2 = list_u[j]
			for elem1 in list1:
				a=0
				for elem2 in list2:
					if MatConn[elem1][elem2]:
						a +=1
				if CMatrix[i-1][j-1] < a:
					CMatrix[i-1][j-1] = a
			# for element in itertools.product(*somelists):
			# 	if MatConn[element[0]][element[1]]:
			# 		CMatrix[i-1][j-1] = 1
			# 		break
	f = open('CMatrix.dat','wb')
	dill.dump(CMatrix, f)
	f.close()
	return CMatrix

def listid(n):
	# global list_u
	for u in list_u:
		if n in u:
			# print(list_u.index(u))
			return(list_u.index(u))

def getB(n):
	return(np.where(MatConn[n,:] == 1)[0].tolist())

class node(object):
	address = 1
	def __init__(self,parent,val):
		self.parent = parent
		self.val = val
		self.child = []
		self.llid = [listid(self.val)]
		self.plist = [self.val]
		self.address = node.address
		# node.address = node.address + 1

	def appendNode(self,val):
		# print(self.llid, listid(val), self.llid.count(listid(val)))
		# Following condition imposed for weighted approach. Use above commented line for only single node from each list
		# if listid(val) not in self.llid:
		if self.llid.count(listid(val)) < weights_for_lists[listid(val)-1]:
			a = node(self.address,val)
			self.child = self.child + [a]
			a.llid = self.llid + [listid(a.val)]
			a.plist = a.plist + self.plist
			# if (len(a.llid) == max_Num and MatConn[a.plist[-1],a.plist[0]]):   #Constraint to get hamiltonian paths in the solution
			if (len(a.llid) == max_Num and CMatrix[a.llid[-1]-1,a.llid[0]-1]):   #Constraint to ensure hamiltonian cycle exists on lists
			# if (len(a.llid) == max_Num):
				# print (a.plist)
				soln_file = open('hamiltonian_soln.txt', 'a+')
				soln_file.writelines(str(a.plist) + '\n')
				global sol_found_till_now
				sol_found_till_now = sol_found_till_now + 1
				# print('soln found')
				if sol_found_till_now >= sufficient_algorithmic_solutions:
					global sufficient_sol_found 
					sufficient_sol_found = True
					# print(sufficient_sol_found)
				soln_file.close()
			return a
			# print(a.plist)
	# def __repr__(self):
	# 	return str(self.val)

	# def __str__(self):
	# 	return str(self.val)

def appendall(basenode,Blist):
	# print('in appendall')
	for i in Blist:
		basenode.appendNode(i)

def exploreTree(baseN,num):
	# sufficient_sol_found = False
	# print(baseN.val, getB(baseN.val), baseN.plist, baseN.llid)
	global sufficient_sol_found
	# print(sufficient_sol_found)
	if num == 0 or sufficient_sol_found:
		return
	else:
		appendall(baseN,getB(baseN.val))
		# print('valid childs',[i.val for i in baseN.child])
		if baseN.child:
			for i in baseN.child:
				if sufficient_sol_found:
					return
				exploreTree(i,num-1)

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


# def load_data(filename):
# 	global list_u
# 	global MatConn
# 	global starter_node
# 	global max_Num
# 	# global dat_name
# 	f = open(filename, 'rb')
# 	list_u = dill.load(f)
# 	MatConn = dill.load(f)
# 	starter_node = dill.load(f)
# 	max_Num = dill.load(f)
# 	f.close()

def convert_to_gams_data(filename, max_Num, node_in_unit_cell):
	f = open(filename, 'rb')
	list_u = dill.load(f)
	MatConn = dill.load(f)
	starter_node = dill.load(f)
	# max_Num = dill.load(f)
	f.close()
	# global df_new
	# global df
	df_file = open('df_new.dat','rb')
	df_new = dill.load(df_file)
	df = dill.load(df_file)
	df_file.close()
	node_array = df_new.node_num.to_numpy()
	CMat = open('CMatrix.dat','rb')
	CMatrix = dill.load(CMat)
	CMat.close()

	gams_data = open('data.txt', 'w') 
	gams_data.writelines(('Set' + '\n')) 
	gams_data.writelines(("n     'node list'         / 1*" + str(node_array[-1]) + "/" + '\n'))
	gams_data.writelines(("u     'number of unique nodes'        / 1*" + str(len(weights_for_lists)) + "/" + '\n'))
	gams_data.writelines(("l     'number of nodes total by accounting weights'        / 1*" + str(max_Num) + "/" + '\n'))
	string = ""
	for counter in range(1,len(list_u)):
		string = string + str(counter) + "." + str(list_u[counter]).replace("[","(").replace("]",")") + ","
	string = string.rstrip(',')
	gams_data.writelines(("list(u,n) 'list of nodes in each unique list'        /" + string + "/" + '\n'))
	string = ""
	for i,l1 in enumerate(MatConn):
		for j,val in enumerate(l1):
			if val == 1:
				string = string + "(" + str(i) + "." + str(j) + "),"  
	string = string.rstrip(',')
	gams_data.writelines(("M(n,n)  'Connectivity Matrix'        /" + string + "/" + '\n'))
	# Writing CMatrix
	string = ""
	for i,r1 in enumerate(CMatrix):
		for j,val in enumerate(r1):
			if val == 1:
				string = string + "(" + str(i+1) + "." + str(j+1) + "),"  
	string = string.rstrip(',')
	gams_data.writelines(("CMatrix(u,u)  'Adjacency Matrix'        /" + string + "/;" + '\n'))
	gams_data.writelines(("Parameter max_Num 'number of nodes to be selected';        max_Num = " + str(max_Num) + ";" + '\n'))
	# Writing weights as parameters
	string = ""
	for i in range(len(weights_for_lists)):
		string = string + str(i+1) + " = " + str(weights_for_lists[i]) + ","  	
	string = string.rstrip(',')
	gams_data.writelines(("Parameter weight(u) 'weights to be given to lists' /" + string + "/;" + '\n'))
	gams_data.close()
	# weight_eqn_data.open('weighted_eqn.txt','w')
	# gams_data.writelines(("equation list_weights(u); list_weights(u) .. sum(n$list(u,n),y(n))=e=weight(u); " + str(max_Num) + ";" + '\n'))
	# weight_eqn_data.close()

	init_guess = open('init_guess.txt', 'w')
	# sol_file = open('hamiltonian_soln.txt', 'r') 
	init_guess.writelines('y.fx(\'' + str(node_in_unit_cell) + '\') = 1;')
	# Lines = sol_file.readlines()
	# Lines = Lines[0]
	# string = ''
	# algorithmic_solution =  ast.literal_eval(Lines) 
	# for i in algorithmic_solution:
	#     string = string + "y.l(\'" + str(i) + "\') = 1;"
	# init_guess.writelines(string)
	# sol_file.close()    
	init_guess.close()

# t1 = perf_counter() - t0


global starter_node
global max_Num
parent_dir=os.getcwd()
os.chdir(parent_dir)

cmap = plt.get_cmap('Reds_r')
new_cmap = truncate_colormap(cmap, 0.2, 0.6)
cmap = plt.get_cmap('Greys_r')
cmap_grey = truncate_colormap(cmap, 0.3, 0.7)

# for filename in ['tmp_17_1_10_W9gfvB.cif']:
# list_of_158_no_err_results = ['CHA.cif', 'AVL.cif', 'EWO.cif', 'ITH.cif', 'CSV.cif', 'SBN.cif', 'CAN.cif', 'CZP.cif', 'SZR.cif', 'ITW.cif', 'OSI.cif', 'THO.cif', 'MTT.cif', 'BRE.cif', 'ATS.cif', 'SFF.cif', 'SOD.cif', 'MEP.cif', 'AWW.cif', 'MRT.cif', 'MEL.cif', 'AWO.cif', 'RWR.cif', 'MTN.cif', 'SOS.cif', 'CAS.cif', 'SSY.cif', 'MTF.cif', 'DFT.cif', 'CDO.cif', 'BPH.cif', 'KFI.cif', 'NSI.cif', 'RTH.cif', 'BOF.cif', 'IFR.cif', 'OFF.cif', 'STI.cif', 'CFI.cif', 'MFS.cif', 'AFR.cif', 'ASV.cif', 'AFN.cif', 'TON.cif', 'JRY.cif', 'NAT.cif', 'AEL.cif', 'OBW.cif', 'PHI.cif', 'LTJ.cif', 'UFI.cif', 'SGT.cif', 'EEI.cif', 'SBT.cif', 'RRO.cif', 'ATN.cif', 'MAZ.cif', 'EPI.cif', 'SVV.cif', 'UOZ.cif', 'WEI.cif', 'OKO.cif', 'BOG.cif', 'NAB.cif', 'OSO.cif', 'BEC.cif', 'PON.cif', 'SFH.cif', 'SOF.cif', 'LIO.cif', 'AFS.cif', 'MER.cif', 'ETV.cif', 'BOZ.cif', 'NES.cif', 'GME.cif', 'PUN.cif', 'SFS.cif', 'ATO.cif', 'SAS.cif', 'SAO.cif', 'NPO.cif', 'JSN.cif', 'ZON.cif', 'AFI.cif', 'TER.cif', 'LOV.cif', 'LTA.cif', 'ITE.cif', 'APC.cif', 'SFE.cif', 'GON.cif', 'PWW.cif', 'AFX.cif', 'MVY.cif', 'AEN.cif', 'LTL.cif', 'RWY.cif', 'UEI.cif', 'ATT.cif', 'EDI.cif', 'NPT.cif', 'JSW.cif', 'MEI.cif', 'NON.cif', 'GIS.cif', 'PCR.cif', 'IFY.cif', 'USI.cif', 'LAU.cif', 'ERI.cif', 'BIK.cif', 'SFN.cif', 'ABW.cif', 'APD.cif', 'GOO.cif', 'SAF.cif', 'IHW.cif', '_CON.cif', 'RTE.cif', 'IFO.cif', 'ACO.cif', 'STF.cif', 'EZT.cif', 'JOZ.cif', 'OWE.cif', 'MOR.cif', 'MON.cif', 'BCT.cif', 'AFO.cif', 'JNT.cif', 'ETR.cif', 'DAC.cif', 'AEI.cif', 'UTL.cif', 'CGS.cif', 'YUG.cif', 'PTY.cif', 'SFO.cif', 'AFV.cif', 'AET.cif', 'ANA.cif', 'UOS.cif', 'AHT.cif', 'MTW.cif', 'IWV.cif', 'JBW.cif', 'IWR.cif', 'VSV.cif', 'ATV.cif', 'FER.cif', 'IRR.cif', 'SAV.cif', 'VET.cif', 'DON.cif', 'PWO.cif', 'AST.cif', 'RHO.cif']
debug_list = ['SOD.cif']

# for filename in debug_list:
for filename in os.listdir(parent_dir):
	if ('.cif' in filename):
		zeolite_name = filename.replace(".cif","")
		print(zeolite_name)
		os.chdir(parent_dir)
		inside_folder = os.path.join(parent_dir, zeolite_name)
		if os.path.exists(inside_folder):
			pass
		else:
			os.mkdir(inside_folder)
		cif_name = zeolite_name + '.cif'
		xyz_name = zeolite_name + '.xyz'
		csv_name = zeolite_name + '.csv'
		dat_name = zeolite_name + '.dat'
		sol_name = zeolite_name + '.sol'

		# t2 = perf_counter() - t0

		shutil.copy(os.path.join(parent_dir, cif_name), os.path.join(inside_folder, cif_name))
		os.chdir(inside_folder)
		os.system('obabel -icif '+str(filename)+' -oxyz -O '+str(xyz_name)+' --fillUC strict')
		# subprocess.run('obabel -icif '+str(filename)+' -oxyz -O '+str(xyz_name)+' --fillUC strict')

		# Important parameters that generate the entire lattice from the CIF file
		La, Lb, Lc, ax, bx, by, cx, cy, cz = read_and_report(cif_name)
		
		shutil.copy(os.path.join(parent_dir, 'SRU_v2.gms'), os.path.join(inside_folder, 'SRU_v2.gms'))

# 		os.replace(parent_dir +'/'+ csv_name, path+'/'+ csv_name)
# 		os.chdir(path)
# 		# t3 = perf_counter() - t0
		# print(zeolite_name, round(ax,5), round(bx,5), round(cx,5), round(by,5), round(cy,5), round(cz,5))
		# break
		global df_new
		global df
		global list_u

		# if False:
		if os.path.exists('df_new.dat'):
			print('loading data from dat file')
			f = open(dat_name, 'rb')
			list_u = dill.load(f)
			MatConn = dill.load(f)
			starter_node = dill.load(f)
			# max_Num = dill.load(f)
			f.close()
			df_file = open('df_new.dat','rb')
			df_new = dill.load(df_file)
			df = dill.load(df_file)
			df_file.close()
		else:
			print('expanding cif to csv lattice')
			xyz_to_csv(zeolite_name)
			expand_to_lattice_size(zeolite_name, round(ax,6), round(bx,6), round(cx,6), round(by,6), round(cy,6), round(cz,6))
			print('computing vector tuples')
			list_u, MatConn, starter_node, df_new, df = Make_matrix_general(zeolite_name)
		# print((zeolite_name, round(ax,6), round(bx,6), round(cx,6), round(by,6), round(cy,6), round(cz,6)))
		# break
		# print(list_u)
		# break

		# This part of code is to get weights which will be useful for larger zeolites 
		# It will tell me how many from each list shall be required to generate the SRU
		# Note that this has no impact on the adjacency matrix and the hamiltonian nature of it
		# single_unit_cell_nodes = df['x'].between(ax+bx*abs(df['y']/by)+cx*abs(df['z']/cz), -0.0001+2*ax+bx*abs(df['y']/by)+cx*abs(df['z']/cz), inclusive=True)
		# single_unit_cell_nodes = df['y'].between(by+cy*abs(df['z']/cz), -0.0001+2*by+cy*abs(df['z']/cz), inclusive=True)
		# single_unit_cell_nodes = df['z'].between(cz, -0.0001+2*cz, inclusive=True)
		

		single_unit_cell_nodes = df['x'].between(round(ax+bx*abs(df['y']/by)+cx*abs(df['z']/cz),4), -0.0001+round(2*ax+bx*abs(df['y']/by)+cx*abs(df['z']/cz),4), inclusive=True) & df['y'].between(round(by+cy*abs(df['z']/cz),4), -0.0001+round(2*by+cy*abs(df['z']/cz),4), inclusive=True) & 	df['z'].between(round(cz,4), -0.0001+round(2*cz,4), inclusive=True)
		# single_unit_cell_nodes = df['x'].between(round(ax+bx*abs(df['y']/by)+cx*abs(df['z']/cz),4), -0.001+round(2*ax+bx*abs(df['y']/by)+cx*abs(df['z']/cz),4), inclusive=True) & df['y'].between(round(by+cy*abs(df['z']/cz),4), -0.001+round(2*by+cy*abs(df['z']/cz),4), inclusive=True) & 	df['z'].between(round(cz,4), -0.001+round(2*cz,4), inclusive=True)
		single_unit_cell_node_num = df.loc[single_unit_cell_nodes].node_num.tolist()
		weights_for_lists = np.zeros(len(list_u)-1).tolist()
		# print(single_unit_cell_node_num)
		# break
		print(len(weights_for_lists))
		for i in single_unit_cell_node_num:
		    weights_for_lists[int(df_new.loc[df_new['node_num'] == i,'id'])-1] += 1
		    # weights_for_lists[find_in_list_u(i)-1] += 1
		    # print(find_in_list_u(i) - df_new.loc[df_new['node_num'] == i,'id'])
		    # print(df_new.loc[df_new['node_num'] == i,'id'])
		# res = list(itertools.compress(range(len(single_unit_cell_nodes )), single_unit_cell_nodes )) 
		# break
		weights_for_lists = [round(x) for x in weights_for_lists]    
		weights_for_lists = [round(x/find_gcd(weights_for_lists)) for x in weights_for_lists]
		# print(weights_for_lists)
		# break
		# total = 0
		# for ele in range(0, len(weights_for_lists)): 
		# 	total = total + weights_for_lists[ele] 
		# set max_Num to total of weights obtained from list for the algorithmic approach calculations
		max_Num = sum(weights_for_lists)
		if any(ele == 0 for ele in weights_for_lists):
			print('ERROR: Weight for element ==0. Decimal precision error, manual intervention needed. Skipping job')
			break
		# print('max num == ', max_Num)
		# break
		# t4 = perf_counter() - t0
		# break
		# if not os.path.exists(path + '/CMatrix.dat'):
		# break
		print('making CMatrix')
		CMatrix = makeCMatrix(dat_name)

		g1 = Graph(len(list_u)-1)
		g1.graph = CMatrix
		g1.hamCycle()

		
		# t5 = perf_counter() - t0
		# print(path)
		if os.path.exists(inside_folder + '/hamsolyes.txt'):
			if solve_algo_approach:
				print('CMatrix is hamiltonian, looking for structural hamiltonian solutions')
				sol_file = open('hamiltonian_soln.txt','w+')
				
				no_soln_found = True
				global sol_found_till_now
				global sufficient_sol_found 
				sufficient_sol_found = False

				iteration_num = 0
				node_array = df_new.node_num.to_numpy()
				while (no_soln_found and (iteration_num < len(node_array))):
					# sol_file = open('hamiltonian_soln.txt', 'r') 
					Lines = sol_file.readlines() 
					# print(list_u)
					sol_found_till_now = 0
					if Lines == []:
						starter_node = node_array[iteration_num]
						# update_starter_node(iteration_num)
						iteration_num = iteration_num + 1
						# print('in for loop')
						tree = node(0,starter_node)
						exploreTree(tree,max_Num)
						# if (iteration_num == len(node_array)):
						# 	print('last_try')
					else:
						no_soln_found = False
					# sol_file.close()
				sol_file.close()

		# t6 = perf_counter() - t0
			# if not os.path.exists(os.path.join(path,'/data.txt')):
			print('generating data for optimization model')
			convert_to_gams_data(dat_name, max_Num, single_unit_cell_node_num[0])

			# global threedee

			# threedee = plt.figure(dpi=300).gca(projection='3d')

			# if not os.path.exists(os.path.join(path,'/gray_bg.dat')):
			# 	# print('gray_bg exists')
			# 	for i in node_array:
			# 		plot_point_gray(i)
			# 		for j in node_array:
			# 			d = distance(df_new.loc[i-1],df_new.loc[j-1])
			# 			if (d < dist_max and d > dist_min):
			# 				plot_line_gray(i,j)					
			# 	# for line in Lines:
			# 	# 	line_num = line_num + 1
			# 	f = open('gray_bg.dat','wb')
			# 	dill.dump(threedee, f)
			# 	f.close()
			# 	plt.close()

			if make_algo_plot:
				sol_file = open('hamiltonian_soln.txt', 'r') 
				Lines = sol_file.readlines()
				Lines = Lines[:4]
				line_num=0
				for line in Lines:
					line_num = line_num + 1
					algorithmic_solution =  ast.literal_eval(line) 
					
					# f = open('gray_bg.dat', 'rb')
					# threedee2 = dill.load(f)
					# f.close() 

					# algorithmic_solution =  ast.literal_eval(line) 
					# # plot_red(threedee2)
					# for i in algorithmic_solution:
					# 	plot_point(i)
					# 	for j in algorithmic_solution:
					# 		d = distance(df_new.loc[i-1],df_new.loc[j-1])
					# 		if (d < dist_max and d > dist_min):
					# 			plot_line(i,j)

					# threedee2.grid(False)
					# # Hide axes ticks
					# threedee2.set_xticks([])
					# threedee2.set_yticks([])
					# threedee2.set_zticks([])
					# # threedee2.auto_scale_xyz 
					# threedee2.axis('off')
					# plt.savefig('algo_plot' + str(line_num) + '.jpg')

				
					threedee2 = plt.figure(dpi=300, figsize=plt.figaspect(1)).gca(projection='3d')
					SRU_df = pd.DataFrame(columns = ['x' , 'y', 'z', 'list_u'])
					for i in algorithmic_solution:
						P = df_new.loc[df_new['node_num'] == i]
						new_row = {'x':float(P.x), 'y':float(P.y), 'z':float(P.z), 'list_u':P.id}
						SRU_df = SRU_df.append(new_row, ignore_index=True)
						threedee2.text(float(P.x), float(P.y), float(P.z), int(P.id)+1, size=9,zdir=None)
					threedee2.scatter(SRU_df['x'], SRU_df['y'], SRU_df['z'], c=SRU_df['z'], cmap = new_cmap, s = 70)
					for i in SRU_df.index:
						for j in SRU_df.index:
							P1 = SRU_df.loc[i]
							P2 = SRU_df.loc[j]
							d = distance(P1,P2)
							if (d < dist_max and d > dist_min):
							    threedee2.plot([float(P1.x),float(P2.x)],[float(P1.y),float(P2.y)],[float(P1.z),float(P2.z)], color = '0.8')
					threedee2.grid(False)
					threedee2.set_xticks([])
					threedee2.set_yticks([])
					threedee2.set_zticks([])
					threedee2.axis('off')
					# threedee2.auto_scale_xyz
					plt.savefig('algo_plot_SRU' + str(line_num) + '.jpg')
					plt.close()


				# Modified plot with only the nearest gray neighbours
				if algo_plot_mod:
					threedee3 = plt.figure(dpi=300, figsize=plt.figaspect(1)).gca(projection='3d')
					# SRU_df = pd.DataFrame(columns = ['x' , 'y', 'z', 'list_u'])
					local_filter = df['x'].between(SRU_df['x'].min()-1.7*dist_max, SRU_df['x'].max()+1.7*dist_max, inclusive=False) & df['y'].between(SRU_df['y'].min()-1.7*dist_max, SRU_df['y'].max()+1.7*dist_max, inclusive=False) & df['z'].between(SRU_df['z'].min()-1.7*dist_max, SRU_df['z'].max()+1.7*dist_max, inclusive=False)

					algo_bool = np.zeros(len(df), dtype=bool)
					for i in algorithmic_solution:
						algo_bool[i-1] = True
					SRU_filter = df.loc[algo_bool]

					localized_df = df.loc[local_filter]
					local_wo_SRU_df = df.loc[local_filter^algo_bool]

					local_wo_SRU_df = local_wo_SRU_df.sort_values(['z'])
					SRU_filter = SRU_filter.sort_values(['z'])
					# localized_df = np.subtract(local_filter, algo_bool, dtype=bool)
					threedee3.scatter(SRU_df['x'], SRU_df['y'], SRU_df['z'], c=SRU_df['z'], cmap = new_cmap, s = 70, zorder=SRU_df['z'])
					threedee3.scatter(local_wo_SRU_df['x'], local_wo_SRU_df['y'], local_wo_SRU_df['z'], c = local_wo_SRU_df['z'], cmap = cmap_grey, s = 70, zorder=local_wo_SRU_df['z'])
					for i in localized_df.index:
						for j in localized_df.index:
							P1 = localized_df.loc[i]
							P2 = localized_df.loc[j]
							d = distance(P1,P2)
							if (d < dist_max and d > dist_min):
							    threedee3.plot([float(P1.x),float(P2.x)],[float(P1.y),float(P2.y)],[float(P1.z),float(P2.z)], color = '0.9', alpha=0.5)
					for i in SRU_filter.index:
						for j in SRU_filter.index:
							P1 = SRU_filter.loc[i]
							P2 = SRU_filter.loc[j]
							d = distance(P1,P2)
							if (d < dist_max and d > dist_min):
							    threedee3.plot([float(P1.x),float(P2.x)],[float(P1.y),float(P2.y)],[float(P1.z),float(P2.z)], color = 'r', alpha=0.5)
					threedee3.grid(False)
					threedee3.set_xticks([])
					threedee3.set_yticks([])
					threedee3.set_zticks([])
					threedee3.axis('off')
					plt.savefig('algo_plot_mod' + str(line_num) + '.jpg')
				plt.close()

				sol_file.close()
			
			# t7 = perf_counter() - t0
			
			if gams_model_solve:
				print('solving gams model')
				os.system('gams SRU_v2.gms')
			
			# # t8 = perf_counter() - t0

			if gams_model_plot:
				database_fix = os.path.join(inside_folder, 'SRU_v2.db')
				dat_fix = sqlite3.connect(database_fix)
				sol_arr = []
				# y = (pd.read_sql_query("SELECT * FROM y", dat_fix)).to_numpy()
				y = (pd.read_sql_query("SELECT * FROM y WHERE level = 1", dat_fix)).to_numpy()

				for line in y:
					if line[1]==1:
						sol_arr.append(int(line[0]))
				# print(sol_arr)
				# pd.read_sql_query("SELECT * FROM Variable", dat_fix)
				print('writing gams solution to file')
				f = open('opti_soln.txt','w')
				f.writelines(str(sol_arr)) 
				f.close()
				sol_file = open('opti_soln.txt', 'r') 
				Lines = sol_file.readlines() 
				line_num=0
				for line in Lines:
					line_num = line_num + 1
					# global threedee
					threedee = plt.figure().gca(projection='3d')
					algorithmic_solution =  ast.literal_eval(line) 
					for i in algorithmic_solution:
					    # print(i)
					    # plot_point(i)
					    plot_point_with_text(i)
					for i in algorithmic_solution:
					    for j in algorithmic_solution:
					        d = distance(df_new.loc[i-1],df_new.loc[j-1])
					        if (d < dist_max and d > dist_min):
					            plot_line(i,j)
					threedee.set_xlabel('x')
					threedee.set_ylabel('y')
					threedee.set_zlabel('z')
					# threedee.auto_scale_xyz 
					threedee.axis('off')
					plt.savefig('opti_plot' + str(line_num) + '.jpg')
					plt.close()
				sol_file.close()

			# t9 = perf_counter() - t0



# print(t1 , 'function definitions')
# print(t2 ,'in for loop, before reading CIF file')
# print(t3 , 'csv file generated and moved, starting Make_matrix_general if data file not present')
# print(t4 , 'Make_matrix_general completed, creating Adjacency Matrix and checking hamiltonian nature')
# print(t5 , 'gams model solved, starting algorithmic search')
# print(t6 , 'starting to write data for gams model')
# print(t7 , 'starting gams model solving')
# print(t8 , 'algorithmic search completed, plotting begins')
# print(t9 , 'plotting done, code ends')

# print(t10 + 'import and function definitions')
# print(t1 + 'import and function definitions')

# print(t1, t2, t3, t4, t5, t6, t7, t8, t9)
