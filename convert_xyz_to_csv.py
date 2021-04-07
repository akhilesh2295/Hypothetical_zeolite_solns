# zeolite_name = 'CHA'
import pandas as pd

def xyz_to_csv(file_name):
	# takes just the name of the zeolite as input and converts xyz file format to csv with only tetrahedral nodes included
	xyz_name = file_name + '.xyz'
	csv_name = file_name + '.csv'

	file_xyz = open(xyz_name, 'r') 
	Lines = file_xyz.readlines() 
	Lines = list(dict.fromkeys(Lines))  #Required to remove duplicates from xyz to csv
	csv_file = open(csv_name, 'w') 
	csv_file.writelines(('x,y,z' + '\n')) 

	for line in Lines:
		if line[0]=='T' or line[0]=='S' :
			# print("Line{}: {}".format(count, line.strip()))
			a = str(line.split()[1:4])
			a = a.strip("'[]")
			a = a.replace("'","")
			csv_file.writelines(a + '\n')
	file_xyz.close()
	csv_file.close() 


def expand_to_lattice_size(file_name, ax, bx, cx, by, cy, cz):
	# takes just the name of the zeolite as input and converts xyz file format to csv with only tetrahedral nodes included
	csv_name = file_name + '.csv'
	df = pd.read_csv(csv_name)
	# x = ax + bx +cx
	# xyz_name = file_name + '.xyz'
	# print(df)
	df_100 = pd.DataFrame(list(zip(df.x+ax, df.y, df.z)), columns=['x','y','z'])
	df_200 = pd.DataFrame(list(zip(df.x+2*ax, df.y, df.z)), columns=['x','y','z'])
	x_axis = [df, df_100, df_200]
	df_x_axis = pd.concat(x_axis, ignore_index=True)


	df_x10 = pd.DataFrame(list(zip(df_x_axis.x+bx, df_x_axis.y+by, df_x_axis.z)), columns=['x','y','z'])
	df_x20 = pd.DataFrame(list(zip(df_x_axis.x+2*bx, df_x_axis.y+2*by, df_x_axis.z)), columns=['x','y','z'])
	xy_axis = [df_x_axis, df_x10, df_x20]
	df_xy_axis = pd.concat(xy_axis, ignore_index=True)

	df_xy1 = pd.DataFrame(list(zip(df_xy_axis.x+cx, df_xy_axis.y+cy, df_xy_axis.z+cz)), columns=['x','y','z'])
	df_xy2 = pd.DataFrame(list(zip(df_xy_axis.x+2*cx, df_xy_axis.y+2*cy, df_xy_axis.z+2*cz)), columns=['x','y','z'])
	xyz_axis = [df_xy_axis, df_xy1, df_xy2]
	df_xyz_axis = pd.concat(xyz_axis, ignore_index=True)
	# print(len(df_xyz_axis))
	df_condensed = df_xyz_axis.drop_duplicates()
	# print(len(df_condensed))
	# print(df_condensed)
	df_condensed.to_csv(csv_name,index=False, float_format='%g')


# import os
# zeolite_name = 'SOD'
# parent_dir=os.getcwd()
# os.chdir(parent_dir)
# path = os.path.join(parent_dir, zeolite_name)
# os.chdir(path)
# xyz_to_csv(zeolite_name)
# expand_to_lattice_size('SOD', 8.965, 0.0, 0.0, 8.965, 0.0, 8.965)