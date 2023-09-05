import math
import time
import random
import numpy as np
import scipy.optimize as opt
from numpy import linalg as LA
from scipy.optimize import fsolve
import param as p


class Optimizer(object):

    def __init__(self, atoms, bonds, xlo, xhi, ylo, yhi, zlo, zhi, K, r0, N, num_reacted, ftype):

        self.atoms = atoms
        self.bonds = bonds
        self.num_reacted=num_reacted # number of reacted groups on each atom (crosslinker)
        self.xlo = xlo
        self.xhi = xhi
        self.ylo = ylo
        self.yhi = yhi
        self.zlo = zlo
        self.zhi = zhi
        self.K = K
        self.r0 = r0          
        self.N = N          
        self.ftype = ftype

    def bondlengths(self):
     
        atoms = self.atoms
        bonds = self.bonds
        Lx = self.xhi - self.xlo
        Ly = self.yhi - self.ylo
        Lz = self.zhi - self.zlo
        n_atoms = len(self.atoms[:,0])
        n_bonds = len(self.bonds[:,0])

        dist = np.zeros((n_bonds,4), dtype=float)

        for i in range (0, n_bonds):
          
              lnk_1 = bonds[i,2]-1
              lnk_2 = bonds[i,3]-1
              delr = atoms[lnk_1,:] - atoms[lnk_2,:]
              
              delr[0] = delr[0] - int(round(delr[0]/Lx))*Lx
              delr[1] = delr[1] - int(round(delr[1]/Ly))*Ly
              delr[2] = delr[2] - int(round(delr[2]/Lz))*Lz
                   
              dist[i,0:3] = delr
              dist[i,3] = LA.norm(delr)
    
        return dist
##
    def bondlength_lnk1_lnk2(self,lnk_1,lnk_2):
     
        atoms = self.atoms
        bonds = self.bonds
        Lx = self.xhi - self.xlo
        Ly = self.yhi - self.ylo
        Lz = self.zhi - self.zlo
        n_atoms = len(self.atoms[:,0])
        n_bonds = len(self.bonds[:,0])

        dist = np.zeros(4, dtype=float)

##        for i in range (0, n_bonds):
          
##              lnk_1 = bonds[i,2]-1
##              lnk_2 = bonds[i,3]-1
        delr = atoms[lnk_1-1,:] - atoms[lnk_2-1,:]
          
        delr[0] = delr[0] - int(round(delr[0]/Lx))*Lx
        delr[1] = delr[1] - int(round(delr[1]/Ly))*Ly
        delr[2] = delr[2] - int(round(delr[2]/Lz))*Lz
               
        dist[0:3] = delr
        dist[3] = LA.norm(delr)
    
        return dist
##
##    
    def invlangevin(self, x):
        return x*(2.99942 - 2.57332*x + 0.654805*x**2)/(1-0.894936*x - 0.105064*x**2)

    def kuhn_stretch(self, lam, E_b):
       
        def func(x, lam, E_b):
            y = lam/x
            beta = self.invlangevin(y)
            return E_b*np.log(x) - lam*beta/x
   
        if lam == 0:
           return 1
        else:
           lam_b = opt.root_scalar(func,args=(lam, E_b),bracket=[lam,lam+1],x0=lam+0.05)
           return lam_b.root
##
    def get_bondforce(self, r):

        K  = self.K
        r0 = self.r0
        Nb = self.N # b = 1 (lenght scale of the system)
        E_b = 1200
 
        x = (r-r0)/Nb
        if(x<0.90):
           lam_b = 1.0
           fbkT  = self.invlangevin(x) ## fb/kBT
           fbond = -K*fbkT/r
        elif(x<1.4):
           lam_b = self.kuhn_stretch(x, E_b)
           fbkT  = self.invlangevin(x/lam_b)/lam_b ## fb/kBT
           fbond = -K*fbkT/r
        else:
           lam_b = x + 0.05
           fbkT  = 325 + 400*(x-1.4)    ## fb/kBT        
           fbond = -K*fbkT/r
 
        return fbond, lam_b  



    def KMCbondexchange(self, U0, tau, delta_t, pflag, index,func):
        num_primary_loops=0
        n_bonds = len(self.bonds[:,0])

        for i in range(0,n_bonds):
          if(self.bonds[i,2]==self.bonds[i,3]):
             num_primary_loops=num_primary_loops+1
        

        
        dist = self.bondlengths()
        bond_index    = random.randint(0, len(self.bonds)-1)

        old_r=dist[bond_index,3]
        r=old_r
        
        if(r > 0):
              [fbond, lam_b] = self.get_bondforce(r)
        else: fbond = 0.0

        fit_param = 1
        fbkT = -fbond*r/self.K
##            U0=56.7578
        old_rate = math.exp(- fbkT*fit_param) #rate
##        print('old_fbkT',fbkT)


        
        bond_end=random.randint(0, 1)
##        print((1-bond_end)+2)
        C1=self.bonds[bond_index,bond_end+2] # crosslinker_1- the crosslinker with which the chain end is associated
        C2=self.bonds[bond_index,(1-bond_end)+2] # crosslinker 2- chain is currently associated, but will be broken
        potential_Ctmp=[] # should be partially reacted, and new distance should be <= contour length
        # the condition for partial reaction is because for exchanging, we will need another crosslinker with which to exchange, which has to be partially reacted
##        partially_reacted=[]
        for i in range(1,len(self.num_reacted)+1): # i is crosslinker_number, crosslinker_index=crosslinker_number-1
            if self.num_reacted[i-1]<func:
                
                dist=self.bondlength_lnk1_lnk2(C1,i) # i is potential Ctmp
                r=dist[3]
                if(r<=p.N*p.b):
                    potential_Ctmp.append(i)

        no_Ctmp=True
        if(len(potential_Ctmp)>0):
##            stop
            no_Ctmp=False
            C_index=random.randint(0, len(potential_Ctmp)-1)
            Ctmp=potential_Ctmp[C_index]
       
            dist=self.bondlength_lnk1_lnk2(C1,Ctmp)
            new_r=dist[3]
            r=new_r
    
            if(r > 0):
                  [fbond, lam_b] = self.get_bondforce(r)
            else: fbond = 0.0

            fit_param = 1
            fbkT = -fbond*r/self.K
    
            new_rate = math.exp(-fbkT*fit_param)  #rate 
   
            ratio=new_rate/old_rate ## r2/r1
            
            accept=False
            if(ratio>=1):
                accept=True
            else:
                rnd_num  = random.uniform(0,1)
                if(ratio>=rnd_num):
                    accept=True
                else:
                    accept=False # reject move
 
            if(accept==True):
   
                if(self.bonds[bond_index,(1-bond_end)+2]!=C2):
                    stop
                self.bonds[bond_index,(1-bond_end)+2]=Ctmp # exchange done!
                self.bonds[bond_index,1]=2
    
                self.num_reacted[Ctmp-1]=self.num_reacted[Ctmp-1]+1
                self.num_reacted[C2-1]=self.num_reacted[C2-1]-1
                
            return bond_index, self.num_reacted, accept, no_Ctmp
        else:
            return False, False, False, no_Ctmp
        
 
