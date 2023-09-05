#!/usr/local/bin/env python
# -*- coding: utf-8 -*-

"""
#######################################
#                                     #
#-- ---------------Generating polymer network using Annable algorithm ------------------------#
#--- T. Annable et al./Colloids Surfaces A." Physicochem. Eng. Aspects 112 (1996) 97-116 -----#
#------------------------------  Author: Devosmita Sen  --------------------------------------#
#                                     #
#######################################

"""

import time
import math
import random
import netgen 
import ioLAMMPS
import matplotlib
import numpy as np
from exchange import Optimizer
from numpy import linalg as LA
from scipy.optimize import fsolve
from matplotlib import pyplot as plt
import param as p
from timeit import default_timer as timer

#-------------------------------------#
#       Simulation Parameters         #
#-------------------------------------#
start=timer()
#N  = 12

K  = p.K #1.0
r0 = 0.0
##U0  = 1
tau = 1
del_t = p.del_t#0.005

wrt_step = p.wrt_step#100
func = p.func#4
N    = p.N #12
Nb = N
rho  = 3
l0   = 1
prob = 1.0
n_chains  = p.n_chains#1000
n_links   = int(2*n_chains/func)
L = p.L #32.5984
U0=p.U0_low#56.7578
netgen_flag = 1
swell = 0
vflag = 0


##random.seed(a=500)
random.seed(a=None, version=2)
print('First random number of this seed: %d'%(random.randint(0, 10000))) 
# This is just to check whether different jobs have different seeds

bond_type=1 # 1=no exchange, 2=exchange
if(netgen_flag==0):

   
##   N = 12   
   print('--------------------------')   
   print('----Reading Network-------')   
   print('--------------------------')   
   [xlo, xhi, ylo, yhi, zlo, zhi, n_atoms, n_bonds, atoms, bonds, 
           atom_types, bond_types, mass, loop_atoms] = ioLAMMPS.readLAMMPS("network.txt", bond_type, vflag)
   
   print('xlo, xhi',xlo, xhi) 
   print('ylo, yhi',ylo, yhi) 
   print('zlo, zhi',zlo, zhi) 
   print('n_atoms', n_atoms) 
   print('n_bonds', n_bonds) 
   print('atom_types = ', atom_types) 
   print('bond_types = ', bond_types) 
   print('mass = ', mass) 
   print('primary loops = ', len(loop_atoms)) 
   print('--------------------------')   

elif(netgen_flag==1):

   
   
   netgen.generate_network(prob, func, N,bond_type, L, l0, n_chains, n_links)

   [xlo, xhi, ylo, yhi, zlo, zhi, n_atoms, n_bonds, atoms, bonds, 
           atom_types, bond_types, mass, loop_atoms] = ioLAMMPS.readLAMMPS("network.txt", bond_type, vflag)

else:
   print('Invalid network generation flag')


fstr=open('stress','w')
fstr.write('#Lx, Ly, Lz, lambda, FE, deltaFE, st[0], st[1], st[2], st[3], st[4], st[5]\n') 

flen=open('strand_lengths','w')
flen.write('#lambda, ave(R), max(R)\n') 

fkmc=open('KMC_stats','w')
fkmc.write('#lambda, init bonds, final bonds\n') 

#-------------------------------------#
#       Bond Exchanges                #
#-------------------------------------#
num_reacted=np.zeros(len(atoms)) # number of reacted junction sites on each atom (crosslinker)

for i in loop_atoms: ## add the primary loop chains(bonds) into the bonds list
##   stop
   bonds=np.vstack((bonds,np.array([1,1,int(i),int(i)])))

for bond_cnt in range(len(bonds)):  
   link_1=bonds[bond_cnt,2]
   link_2=bonds[bond_cnt,3]
   num_reacted[link_1-1]=num_reacted[link_1-1]+1
   num_reacted[link_2-1]=num_reacted[link_2-1]+1


num_reacted=num_reacted/2.0 ## because there has been double counting in earlier loop   
mymin = Optimizer(atoms, bonds, xlo, xhi, ylo, yhi, zlo, zhi, K, r0, N,num_reacted, 'Mao')

filename = 'network_initial_for_count.txt' #%d.txt' %(i+1)
ioLAMMPS.writeLAMMPSafterexchange(filename, mymin.xlo, mymin.xhi, mymin.ylo, mymin.yhi, mymin.zlo, mymin.zhi,
                                       mymin.atoms, mymin.bonds, atom_types, 2, mass, loop_atoms)

print(np.shape(mymin.bonds))
num_exchange_steps=p.num_exchange_steps#n_chains*100
bond_index_list=[]
avg_r_arr=[]
avg_r2_arr=[]
loop_frac_arr=[]
secondary_loop_frac_arr=[]
accept_arr_last_few_iter=[]
accept_frac_arr=[]
accept_rate_arr=[]
n_accept=0 # number of accepted moves currently
MC_step_arr=[]

filename = 'restart_network_initial.txt' #%d.txt' %(i+1)
ioLAMMPS.writeLAMMPSafterexchange(filename, mymin.xlo, mymin.xhi, mymin.ylo, mymin.yhi, mymin.zlo, mymin.zhi,
                                       mymin.atoms, mymin.bonds, atom_types, 2, mass, loop_atoms)

num_primary_loops=0
n_bonds = len(mymin.bonds[:,0])

for i in range(0,n_bonds):
 if(mymin.bonds[i,2]==mymin.bonds[i,3]):
    num_primary_loops=num_primary_loops+1

i=-1
while(i<num_exchange_steps):
   
   bond_index,num_reacted,accept,no_Ctmp=mymin.KMCbondexchange(U0, tau, del_t, 0, i+1,func)

   if(no_Ctmp==True): # no available sites for exchange

      continue
      
   else:
      i=i+1
      bond_index_list.append(bond_index)
      accept_arr_last_few_iter.append(int(accept))
      if((i+1)%wrt_step==0): 
         filename = 'restart_network_%d.txt' %(i+1)
   ##      print(mymin.bonds[bond_index,1])
         num_primary_loops, num_secondary_loops=ioLAMMPS.writeLAMMPS(filename, mymin.xlo, mymin.xhi, mymin.ylo, mymin.yhi, mymin.zlo, mymin.zhi,
                                              mymin.atoms, mymin.bonds, atom_types, 2, mass, loop_atoms)
  
         
         print('iteration_step',i+1)
         for bond_index in bond_index_list:
            mymin.bonds[bond_index,1]=1
         bond_index_list=[]
         dist = mymin.bondlengths()
         r=dist[:,3]
         avg_r=np.mean(r) # avg end to end distance
         avg_r_arr.append(avg_r)
         avg_r2=(np.square(r)).mean()
         avg_r2_arr.append(avg_r2)
         loop_frac_arr.append(num_primary_loops/p.n_chains)#len(mymin.bonds))
         secondary_loop_frac_arr.append(num_secondary_loops/p.n_chains)#len(mymin.bonds))
         n_accept_this_set=np.sum(accept_arr_last_few_iter)
         n_accept=n_accept+n_accept_this_set
         accept_frac_arr.append(n_accept_this_set/len(accept_arr_last_few_iter))
         accept_rate_arr.append(n_accept/(i+1))
         accept_arr_last_few_iter=[]
         MC_step_arr.append(i+1)

   
plt.figure()
plt.plot(MC_step_arr,avg_r_arr,'bo-',label='mean_r')      
plt.legend()
plt.savefig('mean_r')

plt.figure()
plt.plot(MC_step_arr,avg_r2_arr,'bo-',label='mean_r2')#,'bo-')
plt.legend()
plt.savefig('mean_r2')

plt.figure()
plt.plot(MC_step_arr,loop_frac_arr,'bo-',label='primary_loop_frac')#,'bo-')
plt.legend()
plt.savefig('primary_loop_frac')


plt.figure()
plt.plot(MC_step_arr,secondary_loop_frac_arr,'bo-',label='secondary_loop_frac')#,'bo-')
plt.legend()
plt.savefig('secondary_loop_frac')

plt.figure()
plt.plot(MC_step_arr,accept_rate_arr,'bo-',label='Acceptance rate')#,'bo-')
plt.legend()
plt.savefig('acceptance_rate')


plt.figure()
plt.plot(MC_step_arr,accept_frac_arr,'bo-',label='Fraction accepted over last '+str(wrt_step)+' steps')#,'bo-')
plt.legend()
plt.savefig('frac_accept')

np.savetxt('all_data',np.transpose(np.array([MC_step_arr,avg_r_arr,avg_r2_arr,loop_frac_arr,secondary_loop_frac_arr,accept_rate_arr,accept_frac_arr])),header='MC_step_arr,avg_r_arr,avg_r2_arr,loop_frac_arr,secondary_loop_frac_arr,accept_rate_arr,accept_frac_arr')


filename = 'restart_network_final.txt' #%d.txt' %(i+1)
ioLAMMPS.writeLAMMPSafterexchange(filename, mymin.xlo, mymin.xhi, mymin.ylo, mymin.yhi, mymin.zlo, mymin.zhi,
                                       mymin.atoms, mymin.bonds, atom_types, 2, mass, loop_atoms)


fstr.close()
flen.close()
fkmc.close()

end=timer()

time=end-start
print('time for ',str(len(mymin.bonds)),'chains is = ',time,' seconds')
plt.show()
