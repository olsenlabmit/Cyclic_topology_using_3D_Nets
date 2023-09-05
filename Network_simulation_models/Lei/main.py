"""
#######################################
#                                     #
#-------------- Generating polymer network using Lei algorithm --------------#
#-------- Lei and Liu, Journal of Applied Physics 132, no. 13 (2022) --------# 
#------------------------- Author: Devosmita Sen  ---------------------------#
#                                     #
#######################################

"""

## scipy Delaunay with PBC included- in 3D
# Model: "A network mechanics method to study the mechanism of the large-deformation fracture of elastomers." Journal of Applied Physics 132, no. 13 (2022).
##Modification: PBC has been implemented

import numpy as np
import random
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
import networkx as nx
from itertools import combinations
from timeit import default_timer as timer
#import param as p
import subprocess



## *********** This is now defined at the top***********
L=10
Lx = L; Ly = L; Lz = L
n_nodes=736 #736 for 5006 nodes
#*******************************************************

def dist(lnk_1,lnk_2, Lx, Ly, Lz):
    dx=lnk_1[0]-lnk_2[0]
    dy=lnk_1[1]-lnk_2[1]
    dz=lnk_1[2]-lnk_2[2]
    if(abs(int(round(dx/Lx))*Lx>0) or abs(int(round(dy/Ly))*Ly>0)or abs(int(round(dz/Lz))*Lz>0)):
        print('lnk_1',lnk_1)
        print('lnk_2',lnk_2)
        
    dx = dx - int(round(dx/Lx))*Lx
    dy = dy - int(round(dy/Ly))*Ly
    dz = dz - int(round(dz/Lz))*Lz
    d=[dx,dy,dz]
    dist = LA.norm(d)

    return dist

def tile_points(points, N):

    # in 3D- create more points all around the box
    # 27 copies of smallest unit cell
    d=3 # dimension
    num_copies=3**d
    # 27 tiles surrounding individual coordinates
    point_tile = np.zeros((num_copies * N, 3)) 
    shift_arr=[0,-L,L]
    cnt=0
    for s_x in shift_arr: ## s_x=shift in x
        for s_y in shift_arr:
            for s_z in shift_arr:
                point_tile[N*cnt:N*(cnt+1),0]=points[:,0]+s_x
                point_tile[N*cnt:N*(cnt+1),1]=points[:,1]+s_y
                point_tile[N*cnt:N*(cnt+1),2]=points[:,2]+s_z
                cnt=cnt+1

    return point_tile
    

start=timer()
random.seed(a=None,version=2)


# Array sizes and initializations
nodes   = np.zeros((n_nodes,3), dtype = float)


for i in range (0, n_nodes):
    x = random.uniform(0, Lx)
    y = random.uniform(0, Ly)
    z = random.uniform(0, Lz)
    nodes[i,:] = [x, y, z]
    
old= nodes    
print ('Before', np.shape(nodes))
N=n_nodes
nodes=tile_points(nodes, N)
print ('after',np.shape(nodes))

## find the Delaunay tetrahedrons
tri = Delaunay(nodes)

ax = plt.axes(projection ="3d")
 
# Creating plot
##ax.scatter3D(nodes[:,0], nodes[:,1], nodes[:,2], color = "red",s=1.5)
##ax.scatter3D(old[:,0], old[:,1], old[:,2], color = "blue")
plt.title('Showing all points. BLUE are original points and RED are tiled points')
#plt.legend(loc=2)

##plt.show

##plt.triplot(points[:,0], points[:,1], tri.simplices)
##plt.plot(points[:,0], points[:,1], 'o')

indices = tri.simplices
indices_old=indices
indices_new=[]
for i in indices:
    a=sorted([x%N for x in i])
    if(a not in indices_new):
        core=np.arange(0,N) # only the core points
        if((i[0] in core) or (i[1] in core) or (i[2] in core) or (i[3] in core)):  # any one of the vertex of the tetrahedron belongs to the main cube, and so, that connection is a valid connection
            indices_new.append(a) # N=total number of points
indices=indices_new

print(len(indices))

cells = [("quad", indices_old)]

vertices = nodes[indices]

num_tetra=np.shape(vertices)[0] # number of tetrahedrons= number of crosslinkers


print('Number of tetrahedra (crosslinkers)=',num_tetra)
##plt.show()
chains=[]
G=nx.Graph()
Gmult=nx.MultiGraph()

## Find the coplanar tetrahedrons (tetrahedrons having a common edge)
for t1 in range(num_tetra): # this is index of tetrahedron (node)

    t1_f_nodes=(indices[t1])

    t1_f_nodes_3set=list(combinations(t1_f_nodes, 3))

    for t1_f_nodes in t1_f_nodes_3set: # each plane of tetrahedron considered seperately
        for t2 in range(0,num_tetra):
            if(t1!=t2):# and t2 not in nodes_visited):
                t2_f_nodes=(indices[t2])
    ##            stop
                if(len(set(t1_f_nodes)& set(t2_f_nodes))==3): # any 3 elements are common- means that the face is common- and hence coplanar
                    
                    if(G.has_edge(t1,t2)==False): 
                        chains.append([t1,t2])
                        G.add_edge(t1,t2)
                        Gmult.add_edge(t1,t2)

for i in G:

    if(G.degree[i]!=4):
        print('degree=4 NOT SATISFIED!!')
        stop
np.savetxt('network.txt',np.array(list(G.edges())),delimiter=',',header='Node1,Node2')
print('Number of chains=',len(chains))



## find body centres
# this is just to get the spacial positions, nothing to do with topology
atoms=np.zeros((num_tetra,3))# crosslinks are the centroid (body centres) of each tetrahedron
for i in range(0,num_tetra):
    atoms[i]=np.mean(vertices[i,:],axis=0)

end=timer()
print('Time taken for network generation=',end-start,' seconds')

