import itertools
import dill
import numpy as np
import pandas as pd
import math

def distance(P1, P2):
    return (((P1.x-P2.x)**2+(P1.y-P2.y)**2+(P1.z-P2.z)**2)**0.5)

def vector(P1,P2,d):
    return np.round((P2.x-P1.x, P2.y-P1.y, P2.z-P1.z)/d,2)+0.0
    # return np.round((math.floor((P2.x-P1.x)*10000/d)/10000, math.floor((P2.y-P1.y)*10000/d)/10000, math.floor((P2.z-P1.z)*10000/d)/10000),3)+0.0

def assign_vectors(i,j,no_of_close_pts):
    d = distance(df.loc[i],df.loc[j])
    if (d < dist_max and d > dist_min):
        no_of_close_pts = no_of_close_pts + 1
        v = vector(df.loc[i],df.loc[j],d)
        arr = np.array([v[0], v[1], v[2]])
        col_name = 'vector'+str(no_of_close_pts)
        df_arr = pd.DataFrame({col_name: [arr]})
        df.loc[[i],[col_name]] = df_arr
        df.at[i,col_name] = arr
    return no_of_close_pts

dist_min = 2.6
dist_max = 3.4

def Make_matrix_general(zeolite_name):

    csv_name = zeolite_name + '.csv'
    dat_name = zeolite_name + '.dat'
    global df
    df = pd.read_csv(csv_name)
    df = df.drop_duplicates() #Line required to drop duplicates generated in lattice
    # print(df)
    df['node_num'] = df.index + 1

    df["vector1"] = np.nan
    df["vector2"] = np.nan
    df["vector3"] = np.nan
    df["vector4"] = np.nan

    
    for i in df.index:
        no_of_close_pts = 0
        filter_xyz = df['x'].between(df.loc[i].x-dist_max, df.loc[i].x+dist_max, inclusive=False) & df['y'].between(df.loc[i].y-dist_max, df.loc[i].y+dist_max, inclusive=False) & df['z'].between(df.loc[i].z-dist_max, df.loc[i].z+dist_max, inclusive=False)
        if len(df.loc[filter_xyz]) >= 4:
            for j in df.loc[filter_xyz].node_num.to_list():
                j=j-1   # since accessing elements in df using index for the rest of this loop, reduce 1 from j
                if no_of_close_pts==4:    #Break loop if all 4 neighbours found
                    break
                no_of_close_pts = assign_vectors(i,j,no_of_close_pts)

    # pd.set_option('display.max_columns', None)
    # pd.set_option('display.max_rows', None)
    df_new = df.dropna()
    df_new.set_index('node_num')
    # global node_array
    node_array = df_new.node_num.to_numpy()
    df_new.loc[:,"all_vectors"] = ""
    for i in df_new.index:
        tuple_vec = []
        tuple_vec.append(df_new.loc[i].vector1.tolist())
        tuple_vec.append(df_new.loc[i].vector2.tolist())
        tuple_vec.append(df_new.loc[i].vector3.tolist())
        tuple_vec.append(df_new.loc[i].vector4.tolist())
        tuple_vec.sort()
        df_new.loc[i,'all_vectors'] = str(tuple_vec)
    df_new = df_new.assign(id=(df_new['all_vectors']).astype('category').cat.codes)
    df_new.drop(columns='all_vectors')

    grouped_lists = df_new.groupby('id').agg(lambda col: col.tolist())['node_num'].reset_index().node_num.tolist()
    MatConn = np.zeros((node_array[-1]+1,node_array[-1]+1))

    # Output of connectivities after null values dropped. This gives the reduced connectivitiy matrix
    for i in node_array:
        filter_xyz = df['x'].between(df.loc[i-1].x-dist_max, df.loc[i-1].x+dist_max, inclusive=False) & df['y'].between(df.loc[i-1].y-dist_max, df.loc[i-1].y+dist_max, inclusive=False) & df['z'].between(df.loc[i-1].z-dist_max, df.loc[i-1].z+dist_max, inclusive=False)
        for j in node_array:
            d = distance(df_new.loc[i-1],df_new.loc[j-1])
            if (d < dist_max and d > dist_min):
                MatConn[i,j]=1

    list_u = [[0]]+grouped_lists

    # max_Num = len(grouped_lists)

    starter_node = node_array[0]

    f = open(dat_name,'wb')
    dill.dump(list_u, f)
    dill.dump(MatConn, f)
    dill.dump(starter_node, f)
    # dill.dump(max_Num, f)
    f.close()

    f = open('df_new.dat','wb')
    dill.dump(df_new, f)
    dill.dump(df, f)
    f.close()

    return list_u, MatConn, starter_node, df_new, df




# import os
# zeolite_name = 'SOD'
# parent_dir=os.getcwd()
# os.chdir(parent_dir)
# path = os.path.join(parent_dir, zeolite_name)
# os.chdir(path)
# Make_matrix_general(zeolite_name)