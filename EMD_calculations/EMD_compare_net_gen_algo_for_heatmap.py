import networkx as nx
import matplotlib.pyplot as plt
import shutil
import sys
import os.path
import pyomo.environ as pyo
import numpy as np
from timeit import default_timer as timer

start=timer()

#EMD linear optimization code adapted from https://jckantor.github.io/ND-Pyomo-Cookbook/notebooks/03.01-Transportation-Networks.html

max_runs=1
run_1_range=[1,2,3]
run_2_range=[1,2,3]
for run1 in run_1_range:
    for run2 in run_2_range:
        
##        KMC=np.genfromtxt('./KMC/30000/dimlessC_10/Run'+str(run2)+'/lumped_symbols.txt',dtype=int)
##        Gusev=np.genfromtxt('./Gusev/10000/dimlessC_10/Run'+str(run2)+'/lumped_symbols.txt',dtype=int)
        Annable=np.genfromtxt('./Annable/10000/dimlessC_10/Run'+str(run1)+'/lumped_symbols_final.txt',dtype=int)
##        Eichinger=np.genfromtxt('./Eichinger/10000/dimlessC_10/Run'+str(run2)+'/lumped_symbols.txt',dtype=int)
##        Grimson=np.genfromtxt('./Grimson/10000/Run'+str(run1)+'/lumped_symbols.txt',dtype=int)
##        Lei=np.genfromtxt('./Lei/10000/Run'+str(run2)+'/lumped_symbols.txt',dtype=int)
        DWS_1=np.genfromtxt('./Dia_WS/30000/dimlessC_10/Run'+str(run2)+'/lumped_symbols35000.txt',dtype=int)
##        DWS_2=np.genfromtxt('./Dia_WS_Nearest_Neighbor_cnt_by_2/10000/dimlessC_10/Run'+str(run1)+'/lumped_symbols11650.txt',dtype=int)
##        DWS_3=np.genfromtxt('./Dia_WS_Nearest_Neighbor_cnt_by_3/10000/dimlessC_10/Run'+str(run1)+'/lumped_symbols11650.txt',dtype=int)
##        DWS_4=np.genfromtxt('./Dia_WS_Nearest_Neighbor_cnt_by_4/10000/dimlessC_10/Run'+str(run1)+'/lumped_symbols11650.txt',dtype=int)


        net_algo_list=[Annable,DWS_1]#[KMC,Gusev,Annable,Eichinger, Grimson, Lei, DWS_1, DWS_2, DWS_3, DWS_4]
        net_algo_str=['Annable','DWS_1']#['KMC','Gusev','Annable','Eichinger', 'Grimson', 'Lei', 'DWS_1', 'DWS_2', 'DWS_3', 'DWS_4']
        num_algo=len(net_algo_list)
        EMD_mat=np.zeros((num_algo,num_algo))# similarity EMD scores stored in array, symmetric matrix, diagonals 0


        for net1 in range(0,num_algo):
            # network 1
            lumped_symbols_with_freq=net_algo_list[net1] # network_1
            network_symbols_list=lumped_symbols_with_freq[:,0:-1]
            network_symbols_weight_list=lumped_symbols_with_freq[:,6]
            
            for net2 in range(net1+1,num_algo):
                # network_2
                lumped_symbols_with_freq=net_algo_list[net2]  # network_2
                network2_symbols_list=lumped_symbols_with_freq[:,0:-1]
                network2_symbols_weight_list=lumped_symbols_with_freq[:,6]


                lattice_symbols_list=network2_symbols_list
                lattice_symbols_weight_list=network2_symbols_weight_list


                # obtain the length of lists
                network_symbols_list_length = len(network_symbols_list)

                lattice_symbols_list_length = len(lattice_symbols_list)
                ##stop
                # check whether the query is identical to the target
                if np.array_equal((network_symbols_list),(lattice_symbols_list)):
                    network_symbols_array = np.array(network_symbols_list)
                    network_symbols_level_array = np.array(network_symbols_weight_list)
                ##        inds_query = query_smiles_array.argsort()
                    sorted_network_symbols_array = network_symbols_array#[inds_query]
                    sorted_network_symbols_level_array = network_symbols_level_array#[inds_query]

                    lattice_symbols_array = np.array(lattice_symbols_list)
                    lattice_symbols_level_array = np.array(lattice_symbols_weight_list)
                ##        inds_target = target_smiles_array.argsort()
                    sorted_lattice_symbols_array = lattice_symbols_array#[inds_target]
                    sorted_lattice_symbols_level_array = lattice_symbols_level_array#[inds_target]

                    if np.array_equal(
                        network_symbols_array, sorted_lattice_symbols_array
                    ) and np.array_equal(
                        network_symbols_level_array, sorted_lattice_symbols_level_array
                    ):
                ##        return 0.0

                        print('both are same')
                        stop





                #### WEIGHTS ARE NORMALIZED, DISTANCES ARE NOT !!! - WILL DO LATER IF NEEDED

                ## define the required three sets
                Network = {}  # query
                Lattice = {}  # target
                T = {}  # transport flow # Distance matrix

                # define the weight for SMILES in query, w_i

                network_weight_sum = sum(network_symbols_weight_list)
                for i in range(0, network_symbols_list_length):
                    Network[str(network_symbols_list[i])] = network_symbols_weight_list[i]/network_weight_sum #/ query_weight_sum # currently not normalizing , will see later

                # define weight for SMILES in target, w_j
                lattice_weight_sum = sum(lattice_symbols_weight_list)
                for j in range(0, lattice_symbols_list_length):
                    Lattice[str(lattice_symbols_list[j])] = lattice_symbols_weight_list[j]/lattice_weight_sum #target_smiles_weight_list[j] / target_weight_sum


                ## define the distance matrix
                distance_matrix=np.zeros((network_symbols_list_length,lattice_symbols_list_length))

                ## Loop edit distance calculation    
                for i in range(0, network_symbols_list_length):
        ##            if(i%100==0):
        ##                print('i=',i)
                    Vi=network_symbols_list[i]
                    len_symbol=len(Vi)
                    for j in range(0, lattice_symbols_list_length):
                ##        print('j=',j)
                        Vj=lattice_symbols_list[j]
                        dist=0
                        for tmp in range(len_symbol): #Vj should be same length as Vi=6
                            dist=dist+abs(Vi[tmp]-Vj[tmp])
                        distance_matrix[i,j]=dist
                        T[str(network_symbols_list[i]),str(lattice_symbols_list[j])] =  dist

                        
                # linear optimization with pyomo

                # step 0: Create an instance of the model
                model = pyo.ConcreteModel()
                model.dual = pyo.Suffix(direction=pyo.Suffix.IMPORT)

                # step 1: define index sets
                CUS = list(Network.keys())
                SRC = list(Lattice.keys())


                # step 2: define the decision
                model.x = pyo.Var([str(x) for x in network_symbols_list], [str(x) for x in lattice_symbols_list], domain=pyo.NonNegativeReals)

                # step 3: define objective
                model.Cost = pyo.Objective(
                    expr=sum([T[c, s] * model.x[c, s] for c in CUS for s in SRC]), sense=pyo.minimize)

                # step 4: define constraints
                model.src = pyo.ConstraintList()
                for s in SRC:
                    model.src.add(sum([model.x[c, s] for c in CUS]) == Lattice[s])

                model.dmd = pyo.ConstraintList()

                for c in CUS:
                    model.dmd.add(sum([model.x[c, s] for s in SRC]) == Network[c])


                # step 5: solve the model
                results = pyo.SolverFactory("cbc").solve(model)

                if "ok" == str(results.Solver.status):
                    emd=model.Cost()
                    EMD_mat[net1,net2]=emd
                    EMD_mat[net2,net1]=emd
        ##            print('lattice_symbols_list= ',lattice_symbols_list, ', emd= ',emd)
                ##    return emd

                else:
                    print("No Valid Solution Found")
                ##    return False

                    
        np.savetxt('EMD_mat_'+net_algo_str[0]+'_'+str(run1)+'_'+net_algo_str[1]+'_'+str(run2)+'.txt',EMD_mat)
        ##np.savetxt('EMD_mat'+str(net1)+'_'+str(net2)+'.txt',EMD_mat)
        end=timer()
        print('Time taken= ', end-start)
        print(EMD_mat)
        print(net_algo_str)

