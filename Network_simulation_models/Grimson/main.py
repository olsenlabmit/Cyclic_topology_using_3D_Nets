#!/usr/local/bin/env python
# -*- coding: utf-8 -*-

"""
#######################################
#                                     #
#---- Generating polymer network using Grimson algorithm ----#
#-------- Grimson, M. J. Mol. Phys. 1991, 74, 1097. ---------# 
#------------------ Author: Devosmita Sen  ------------------#
#                                     #
#######################################

"""

import numpy as np
import networkx as nx
from timeit import default_timer as timer   
import random

random.seed(a=None,version=2)
start=timer()

G=nx.Graph()
Gmult=nx.MultiGraph()
GDi=nx.DiGraph()
class node:
    def __init__(self,node_pos,node_num):
        self.pos=node_pos
        self.num=node_num

a=17 # divisions on each edge- after implementing periodic boundary condition
# the node with coordinnate a+1 is the same as the node with coordinate 0
# periodic bc implemented
tot_nodes=(a+1)**3
print('tot_nodes',tot_nodes)
nodes=[None]*tot_nodes
node_cnt=0
for x in range(0,a+1):
    for y in range(0,a+1):
        for z in range(0,a+1):
            pos=[x,y,z]
            n=node(pos,node_cnt)
            nodes[node_cnt]=n
            node_cnt=node_cnt+1

node_cnt=0
for x in range(0,a+1):
    for y in range(0,a+1):
        for z in range(0,a+1):

            
            posx1=[(x+1)%(a+1),y%(a+1),z%(a+1)]
            posx2=[(x-1)%(a+1),y%(a+1),z%(a+1)]
            posy1=[x%(a+1),(y+1)%(a+1),z%(a+1)]
            posy2=[x%(a+1),(y-1)%(a+1),z%(a+1)]
            posz1=[x%(a+1),y%(a+1),(z+1)%(a+1)]
            posz2=[x%(a+1),y%(a+1),(z-1)%(a+1)]

            for tmp in range(tot_nodes):
                if(tmp!=node_cnt):
                    if(nodes[tmp].pos in [posx1,posx2,posy1,posy2,posz1,posz2]):

                        if(Gmult.has_edge(node_cnt,tmp)==False): 
                            G.add_edge(node_cnt,tmp)
                            Gmult.add_edge(node_cnt,tmp)
                            GDi.add_edge(node_cnt,tmp)
            node_cnt=node_cnt+1

##stop                        
one=0
two=0
three=0
four=0
five=0
six=0
for i in Gmult:
    if((Gmult.degree(i))==1):
        one=one+1
    if((Gmult.degree(i))==2):
        two=two+1
    if((Gmult.degree(i))==3):
        three=three+1
    if((Gmult.degree(i))==4):
        four=four+1
    if((Gmult.degree(i))==5):
        five=five+1
    if((Gmult.degree(i))==6):
        six=six+1
print('Gmult:')
print('Number of nodes:',len(Gmult))
print('Number of edges',len(Gmult.edges))
print('one',one)
print('two',two)
print('three',three)
print('four',four)
print('five',five)
print('six',six)


one=0
two=0
three=0
four=0
five=0
six=0
for i in G:
    if((G.degree(i))==1):
        one=one+1
    if((G.degree(i))==2):
        two=two+1
    if((G.degree(i))==3):
        three=three+1
    if((G.degree(i))==4):
        four=four+1
    if((G.degree(i))==5):
        five=five+1
    if((G.degree(i))==6):
        six=six+1
print('G:')
print('Number of nodes:',len(G))
print('Number of edges',len(G.edges))
print('one',one)
print('two',two)
print('three',three)
print('four',four)
print('five',five)
print('six',six)

end=timer()

print('Total time taken for making initial lattice =',end-start)

        


netgen_complete=False
nodes_to_cleave=[] # nodes for which the number of neighbors is greater than 4
## corresponds to the nodes for which an associated bond can be cleaved



for node in G.nodes:
    if(len(G[node])>4):
        nodes_to_cleave.append(node)
while(netgen_complete==False):

    
    rnd_node=random.choice(nodes_to_cleave)
    rnd_edge_of_rnd_node=random.choice(list(G.edges(rnd_node)))
    G.remove_edge(rnd_edge_of_rnd_node[0],rnd_edge_of_rnd_node[1])
    
    if(len(G[rnd_edge_of_rnd_node[0]])<=4): # remove first node of edge
        nodes_to_cleave.remove(rnd_edge_of_rnd_node[0])
        
    if(rnd_edge_of_rnd_node[1] in nodes_to_cleave):
        if(len(G[rnd_edge_of_rnd_node[1]])<=4): # remove second node of edge
            nodes_to_cleave.remove(rnd_edge_of_rnd_node[1])
    if(nodes_to_cleave==[]):
        netgen_complete=True


print('Final number of chains in network= ',len(G.edges))
np.savetxt('network.txt',np.array(list(G.edges())),delimiter=',',header='Node1,Node2')


zero=0
one=0
two=0
three=0
four=0

for j in range(0,len(G)):
    i=G[j]
    if(len(i)==0):
        zero=zero+1
    elif(len(i)==1):
        one=one+1
    elif(len(i)==2):
        two=two+1
    elif(len(i)==3):
        three=three+1
    elif(len(i)==4):
        four=four+1

print('four fraction', four/len(G))
print('three fraction', three/len(G))
print('two fraction', two/len(G))
print('one fraction', one/len(G))
print('zero fraction', zero/len(G))

print('bond fraction',len(G.edges)/len(Gmult.edges))
np.savetxt('node_connecitivity_and_bond_fraction.txt',np.array([zero/len(G),one/len(G),two/len(G),three/len(G),four/len(G),len(G.edges)/len(Gmult.edges)]))

