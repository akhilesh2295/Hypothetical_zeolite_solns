# from typing import Set
from numpy import average
import pandas as pd
import dill
import os
os.chdir(os.getcwd())
from Cif_to_Unitcell import *
from convert_xyz_to_csv import *
from Opti_Make_matrix_general import *
import random


from math import gcd

from functools import reduce

def find_gcd(list):
    return reduce(gcd, list)


parent_dir = os.getcwd()

filename = 'SOD.cif'

zeolite_name = filename.replace(".cif","")
print(zeolite_name)
os.chdir(parent_dir)
inside_folder = os.path.join(parent_dir, zeolite_name)

os.chdir(inside_folder)

cif_name = zeolite_name + '.cif'
xyz_name = zeolite_name + '.xyz'
dat_name = zeolite_name + '.dat'


dat_name = 'SOD.dat'
cif_name = 'SOD.cif'



f = open(dat_name, 'rb')
list_u = dill.load(f)
MatConn = dill.load(f)
starter_node = dill.load(f)
f.close()

df_file = open('df_new.dat','rb')
df_new = dill.load(df_file)
df = dill.load(df_file)
df_file.close()



La, Lb, Lc, ax, bx, by, cx, cy, cz = read_and_report(cif_name)



single_unit_cell_nodes = df['x'].between(round(ax+bx*abs(df['y']/by)+cx*abs(df['z']/cz),4), -0.0001+round(2*ax+bx*abs(df['y']/by)+cx*abs(df['z']/cz),4), inclusive='both') & df['y'].between(round(by+cy*abs(df['z']/cz),4), -0.0001+round(2*by+cy*abs(df['z']/cz),4), inclusive='both') & 	df['z'].between(round(cz,4), -0.0001+round(2*cz,4), inclusive='both')
# single_unit_cell_nodes = df['x'].between(round(ax+bx*abs(df['y']/by)+cx*abs(df['z']/cz),4), -0.001+round(2*ax+bx*abs(df['y']/by)+cx*abs(df['z']/cz),4), inclusive=True) & df['y'].between(round(by+cy*abs(df['z']/cz),4), -0.001+round(2*by+cy*abs(df['z']/cz),4), inclusive=True) & 	df['z'].between(round(cz,4), -0.001+round(2*cz,4), inclusive=True)
single_unit_cell_node_num = df.loc[single_unit_cell_nodes].node_num.tolist()
weights_for_lists = np.zeros(len(list_u)-1).tolist()
# print(single_unit_cell_node_num)
# break
# print(len(weights_for_lists))
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


gb = df_new.loc[df_new.is_in_unit_cell==1].groupby('id')


def point_2_centroid_dist(x,y,z, point):
    return (((x-point.x)**2+(y-point.y)**2+(z-point.z)**2)**0.5)

# def swap(node1, node2):
#     if distance()

def objective(node_num_list):
    curr_nodes = df_new.loc[df_new.node_num.isin(node_num_list)]
    x_centroid = average(curr_nodes['x'])
    y_centroid = average(curr_nodes['y'])
    z_centroid = average(curr_nodes['z'])
    obj = 0
    for j in node_num_list:
        obj += point_2_centroid_dist(x_centroid, y_centroid, z_centroid, df_new.loc[j-1])
    return obj


sol_list = {-1}
for id, w in enumerate(weights_for_lists):
# print(random.sample(list(gb.get_group(id)['node_num']), w)[0])
    sol_list.add(random.sample(list(gb.get_group(id)['node_num']), w)[0])
sol_list.remove(-1)

print(sol_list)

ret = objective(sol_list)
print(ret)



import matplotlib.pyplot as plt
import matplotlib.colors as colors

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=256):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

cmap = plt.get_cmap('Reds_r')
new_cmap = truncate_colormap(cmap, 0.2, 0.6)
cmap = plt.get_cmap('Greys_r')
cmap_grey = truncate_colormap(cmap, 0.3, 0.7)

# for i in range(20):

SRU_df = pd.DataFrame(columns = ['x' , 'y', 'z', 'list_u'])
for i in sol_list:
    P = df_new.loc[df_new['node_num'] == i]
    new_row = pd.DataFrame([{'x':float(P.x), 'y':float(P.y), 'z':float(P.z), 'list_u':P.id}])
    # SRU_df = SRU_df.append(new_row, ignore_index=True)
    SRU_df = pd.concat([SRU_df, new_row], ignore_index=True)

if 1<0:
    # sol_file = open('hamiltonian_soln.txt', 'r') 
    # Lines = sol_file.readlines()
    # Lines = Lines[:4]
    # line_num=0
    # for line in Lines:
    #     line_num = line_num + 1
    #     algorithmic_solution =  ast.literal_eval(line) 

    # algorithmic_solution = list(sol_list)
    fig2 = plt.figure(dpi=300, figsize=plt.figaspect(1))
    threedee2 = fig2.add_subplot(projection='3d')
    # SRU_df = pd.DataFrame(columns = ['x' , 'y', 'z', 'list_u'])
    # for i in sol_list:
    #     P = df_new.loc[df_new['node_num'] == i]
    #     new_row = pd.DataFrame([{'x':float(P.x), 'y':float(P.y), 'z':float(P.z), 'list_u':P.id}])
    #     # SRU_df = SRU_df.append(new_row, ignore_index=True)
    #     SRU_df = pd.concat([SRU_df, new_row], ignore_index=True)
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
    # plt.savefig('algo_plot_SRU' + str(line_num) + '.jpg')
    plt.show()
    # plt.close()


    # Modified plot with only the nearest gray neighbours




if True:
    fig3 = plt.figure(dpi=300, figsize=plt.figaspect(1))
    threedee3 = fig3.add_subplot(projection='3d')
    local_filter = df['x'].between(SRU_df['x'].min()-1.7*dist_max, SRU_df['x'].max()+1.7*dist_max, inclusive='both') & df['y'].between(SRU_df['y'].min()-1.7*dist_max, SRU_df['y'].max()+1.7*dist_max, inclusive='both') & df['z'].between(SRU_df['z'].min()-1.7*dist_max, SRU_df['z'].max()+1.7*dist_max, inclusive='both')

    algo_bool = np.zeros(len(df), dtype=bool)
    for i in sol_list:
        algo_bool[i-1] = True
    SRU_filter = df.loc[algo_bool]

    localized_df = df.loc[local_filter]
    local_wo_SRU_df = df.loc[local_filter^algo_bool]

    local_wo_SRU_df = local_wo_SRU_df.sort_values(['z'])
    SRU_filter = SRU_filter.sort_values(['z'])
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
    # plt.savefig('algo_plot_mod' + str(line_num) + '.jpg')
    plt.show()
    plt.close()

# sol_file.close()
