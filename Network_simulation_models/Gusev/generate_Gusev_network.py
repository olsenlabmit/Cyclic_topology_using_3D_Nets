#!/usr/local/bin/env python
# -*- coding: utf-8 -*-

"""
#######################################
#                                     #
#-- Generating polymer network using Gusev algorithm --#
#-------- Macromolecules 2019, 52, 9, 3244–3251 -------#
#----------------- Author: Devosmita Sen  -------------#
#                                     #
#######################################

"""
import os.path
import sys
import time
import math
import random
import matplotlib
matplotlib.use('Agg')
import numpy as np
from timeit import default_timer as timer
import os
from matplotlib.pyplot import cm
import multiprocessing as mp

import networkx as nx
from itertools import combinations
from numpy import linalg as LA
from scipy.optimize import fsolve
from matplotlib import pyplot as plt
import param as p
import shutil

class node_part: 
    def __init__(self, M, y): 
        self.M = M 
        self.y = y


def generate_graph(netgen_flag,test, graph_algo, lattice_type): # specify lattice type only if usig lattice, else: None
    directory = './original_files_netgen/'
    orig_dir = os.path.dirname(directory)
    files=os.listdir(orig_dir)

    directory = './'#+'network'+'/'
    if not os.path.isdir(directory):
      os.mkdir(directory)

    for fname in files:
     
    # copying the files to the
    # destination directory
       shutil.copy2(os.path.join(orig_dir,fname), directory)

    # now add path to frac_weak directory
    file_dir = os.path.dirname(directory)
    sys.path.append(file_dir)
    import ioLAMMPS
    import netgen


    G=nx.Graph()
    if(test==False):
        if(graph_algo=='Gusev'):
            
            if(netgen_flag==0):
              print(graph_algo)
              vflag = 0
            ##   N = 12   
              print('--------------------------')   
              print('----Reading Network-------')   
              print('--------------------------')
              
              filename = "network.txt"
              file_path = os.path.join(directory, filename)
              if not os.path.isdir(directory):
                 os.mkdir(directory)  
              [xlo, xhi, ylo, yhi, zlo, zhi, n_atoms, n_bonds, atoms, bonds, 
                          atom_types, bond_types, mass, loop_atoms,num_secondary_loops] = ioLAMMPS.readLAMMPS_into_graph(G,file_path, vflag, frac_weak)

            elif(netgen_flag==1): # Gusev network generation

              func=p.func
              l0   = 1
              prob = 1.0
              n_chains  = p.n_chains
              n_links   = int(2*n_chains/func)
              L=p.L
            
              num_secondary_loops=netgen.generate_network(G,prob, func, parameters,L, l0, n_chains, n_links, frac_weak)
              directory = './'#+str(int(100*frac_weak))+'/' # while generating network- separate folder created so that existing network data doesn't get deleted
              filename = 'network.txt'
              file_path = os.path.join(directory, filename)
              if not os.path.isdir(directory):
                 os.mkdir(directory)  
              
              [xlo, xhi, ylo, yhi, zlo, zhi, n_atoms, n_bonds, atoms, bonds, 
                      atom_types, bond_types, mass, loop_atoms] = ioLAMMPS.readLAMMPS(file_path,0,frac_weak)

            else:
                print('Invalid network generation flag for graph algorithm: ',graph_algo)

        elif(graph_algo=='KMC'):
            
            if(netgen_flag==0):
              print(graph_algo)

              vflag = 0
            
              print('--------------------------')   
              print('----Reading Network-------')   
              print('--------------------------')
              
              filename = "network_cs=1.1408A3000B1500_00.txt"
              file_path = os.path.join(directory, filename)
              if not os.path.isdir(directory):
                 os.mkdir(directory)  
              ioLAMMPS.readLAMMPS_into_graph_from_bond_only(G,file_path)#, vflag, frac_weak)

        elif(graph_algo=='lattice'):
            
            if(netgen_flag==0):
              print(graph_algo)

              vflag = 0
              
              print('--------------------------')   
              print('----Reading Network-------')   
              print('--------------------------')
              
              filename = lattice_type+"_connectivity.txt"
              file_path = os.path.join(directory, filename)
              if not os.path.isdir(directory):
                 os.mkdir(directory)  
            
              ioLAMMPS.readLAMMPS_into_graph_lattice(G,file_path)#, vflag, frac_weak)    

            else:
              print('Invalid network generation flag for graph algorithm: ',graph_algo)

    elif(test==True):
####        TEST GRAPHS
        G = nx.Graph()
        '''
        nx.add_cycle(G,[1,2,5,4])
        nx.add_cycle(G,[1,4,3,2])
        nx.add_cycle(G,[1,4,6,2])
        '''

        G.clear()
        #square grid
        nx.add_cycle(G,[0,1,4,5])
        nx.add_cycle(G,[1,4,3,2])
        nx.add_cycle(G,[4,3,8,7])
        nx.add_cycle(G,[4,7,6,5])


##        nx.add_cycle(G,[1,2,5,4])
##        nx.add_cycle(G,[1,4,3,2])
##        nx.add_cycle(G,[1,4,6,2])
##
##        G.add_edge(1,2)
##        G.add_edge(1,2)
##        print(list(G.neighbors(1)))
##
##        # 4-7 ring- text example
##        nx.add_cycle(G,[1,2,3,4,5,6,7])
##        nx.add_cycle(G,[7,6,9,8])
##        nx.add_cycle(G,[6,5,10,9])
##        nx.add_cycle(G,[14,13,12,11,5,6,7])
##        G.add_edge(5,11)
##        G.add_edge(11,9)
##
##        nx.add_cycle(G,[6,5,10,9,8,7])
##
##        ##   Testing example 1
##        nx.add_cycle(G,[1,2,3,4])
##        nx.add_cycle(G,[1,5,3,6])
##        G.add_edge(1,3)
##
##        ##   Testing example 2
##        nx.add_cycle(G,[1,2,4,3])#[1,2,3,4]or [1,2,3,4]
##        nx.add_cycle(G,[3,5,2,4])
##        G.add_edge(1,3)
##
##        new test case:
##        nx.add_cycle(G,[1,2,3,4,5,6,7,8])
##        nx.add_cycle(G,[6,5,9,7])
##        G.add_edge(1,10)
##        G.add_edge(1,13)
##        G.add_edge(3,11)
##        G.add_edge(11,12)


##        DRAW NETWORK
        nx.draw(G, with_labels=True,node_color="tab:pink")
        plt.savefig('network.png',bbox_inches='tight')



        plt.show()
    
    return G


if __name__=='__main__':
    
    ###############GENERATE NETWORK#######################
    start1=timer()
    random.seed(a=None, version=2)
    print('First random number of this seed: %d'%(random.randint(0, 10000))) 
    # This is just to check whether different jobs have different seeds
    ##global parameters
    parameters=np.zeros([2,6]) # N, b, K, fit_param, E_b,U0
    parameters[0,:]=np.array([p.N_low,p.b_low,p.K_low,p.fit_param_low,p.E_b_low,p.U0_low])
    parameters[1,:]=np.array([p.N_high,p.b_high,p.K_high,p.fit_param_high,p.E_b_high,p.U0_high])

    frac_weak=0.0


    netgen_flag = 1
    test=False
    graph_algo='Gusev'
    lattice_type='None'


    #######################################
    G=generate_graph(netgen_flag,test, graph_algo,lattice_type)# specify lattice type only if usig lattice, else: None
    #######################################
    end1 = timer()

