"""
#######################################
#                                     #
#------------------Generating polymer network using Eichinger algorithm ------------------------#
#Leung, Yu‐Kwan, and B. E. Eichinger. "Computer simulation of end‐linked elastomers. I. Trifunctional networks cured in the bulk." The Journal of chemical physics 80.8 (1984): 3877-3884.
#------  Author: Devosmita Sen  --------#
#                                     #
#######################################

The growth of molecules proceeds by joining available neighboring reactive groups. Formation of an A-B bond takes
place in the order of increasing rAB- the distance between functional groups A and B.

Discrete molecules are treated as graphs where a vertex is defined as a condensed point of a graph
and the degree of a vertex equals the number of primary chains attached to it.

Reaction between the chosen A and B points occurs provided that the following conditions are satisfied:
(1) the A group is unreacted
(2) the degree of the cross linker to which the B group belongs is less than the maximum allowable number f
(3) the corresponding rAB is less than the set distance parameter d.

Thus, the extent of reaction of A functional groups PAis
controlled by the value d.



"""
## A2Bf type polymer network
import numpy as np
import random
import networkx as nx
from numpy import linalg as LA
from timeit import default_timer as timer
import param as p

start=timer()
random.seed(a=None, version=2)
##random.seed(a=100)
n= p.N# number of skeletal bonds in a polymer (number of Kuhn segments) ##CHECK THIS!!!- OKAY
########M0= 37*1e-3 #kg/mol # g/mol# average MW of one bond unit
########rho=0.97*1e-3/1e-6 #kg/m3 #g/cm3# density of network
########Na=6.022e23# Avogadro number
Np=p.n_chains# number of A2 prepolymers
# cubical box length L
########L=(n*M0*Np/(rho*Na))**(1/3.0)
########L=L*1e10 #Amstrong
L=p.L
#The distribution of end-to-end distances of the primary chains assumes a three-dimensional Gaussian function with
#a one-dimensional variance sigma2
########l= 1.64e-10 #m# A (Amstrong)# length of one backbone bond
b=p.b #l
########Cn=6.3# characteristic ratio
########sigma2=(Cn*n*l**2)/3
sigma2=(n*b**2/3) ###CHECK THIS!!- CORRECT!
sigma=np.sqrt(sigma2)
########sigma=sigma*1e10 #Amstrong
f=p.func# functionality
r= 1# stoichiomatric ratio
Nc=int(2*r*Np/f)# number of crosslinkers

d=n*b #m # contour length


def last_value_list(t):
    return t[-1]


def dist(I_x,I_y,I_z,J_x,J_y,J_z, Lx, Ly, Lz):
    dx=I_x-J_x
    dy=I_y-J_y
    dz=I_z-J_z
    dx = dx - int(round(dx/Lx))*Lx
    dy = dy - int(round(dy/Ly))*Ly
    dz = dz - int(round(dz/Lz))*Lz
    d=[dx,dy,dz]
    dist = LA.norm(d)

    return dist


#Index A functional groups by I.
#The first ends of the primary chains are.represented by 1= 2N - 1 and the second ends by I = 2N for 1<.N<.Np .
coord_I_x = np.zeros(2*Np) #coordinates of odd I (one end of chain)
coord_I_y = np.zeros(2*Np)
coord_I_z = np.zeros(2*Np)

for tmp in range(1,2*Np,2): # odd
    coord_I_x[tmp]=random.uniform(0,L)
for tmp in range(1,2*Np,2):
    coord_I_y[tmp]=random.uniform(0,L)
for tmp in range(1,2*Np,2):
    coord_I_z[tmp]=random.uniform(0,L)

##coord_even_I_x = np.zeros(Np) #coordinates of even I (other end of chain)
##coord_even_I_y = np.zeros(Np)
##coord_even_I_z = np.zeros(Np)

for tmp in range(0,2*Np,2):
    coord_I_x[tmp]=random.gauss(0,sigma)+coord_I_x[tmp+1]
for tmp in range(0,2*Np,2):
    coord_I_y[tmp]=random.gauss(0,sigma)+coord_I_y[tmp+1]
for tmp in range(0,2*Np,2):
    coord_I_z[tmp]=random.gauss(0,sigma)+coord_I_z[tmp+1]


#Index f-functional cross-linking units by J: 1<.J<.Nc
coord_J_x=np.zeros(Nc)
coord_J_y=np.zeros(Nc)
coord_J_z=np.zeros(Nc)

for tmp in range(0,Nc):
    coord_J_x[tmp]=random.uniform(0,L)
for tmp in range(0,Nc):
    coord_J_y[tmp]=random.uniform(0,L)
for tmp in range(0,Nc):
    coord_J_z[tmp]=random.uniform(0,L)


#STEP (6)
n_seg=7 # number of regions each side of the box is divided into

bins=np.linspace(0,L,n_seg)

# assign regions

region_I_x=np.zeros(2*Np) # region on x axis to which end_1 of chains belong to
region_I_y=np.zeros(2*Np) # similarly y
region_I_z=np.zeros(2*Np) # similarly z


##region_x_I_odd=[]
##region_y_I_odd=[]
##region_z_I_odd=[]
##
##region_x_I_even=[]
##region_y_I_even=[]
##region_z_I_even=[]
##
##region_x_J=[]
##region_y_J=[]
##region_z_J=[]



for tmp in range(1,2*Np,2):  ## odd
    # tmp denotes chain end/crosslinker index
    # assign region for end_1 of chains
    for k in range(0,len(bins)): # index over elements of bins
        if(k==len(bins)-1):
            if(coord_I_x[tmp]>=bins[k]):# and coord_I_x[tmp]<bins[k+1]):
                region_I_x[tmp]=k
            if(coord_I_y[tmp]>=bins[k]):# and coord_I_x[tmp]<bins[k+1]):
                region_I_y[tmp]=k
            if(coord_I_z[tmp]>=bins[k]):# and coord_I_x[tmp]<bins[k+1]):
                region_I_z[tmp]=k
        else:        
            if(coord_I_x[tmp]>=bins[k] and coord_I_x[tmp]<bins[k+1]):
                region_I_x[tmp]=k
            if(coord_I_y[tmp]>=bins[k] and coord_I_y[tmp]<bins[k+1]):
                region_I_y[tmp]=k
            if(coord_I_z[tmp]>=bins[k] and coord_I_z[tmp]<bins[k+1]):
                region_I_z[tmp]=k
    
    
##region_I_even_x=np.zeros(Np)# region on x axis to which end_2 of chains belong to
##region_I_even_y=np.zeros(Np)
##region_I_even_z=np.zeros(Np)

for tmp in range(0,2*Np,2): # even
    # tmp denotes chain end/crosslinker index
    # assign region for end_2 of chains
    for k in range(0,len(bins)): # index over elements of bins
        if(k==len(bins)-1):
            if(coord_I_x[tmp]>=bins[k]):# and coord_I_x[tmp]<bins[k+1]):
                region_I_x[tmp]=k
            if(coord_I_y[tmp]>=bins[k]):# and coord_I_x[tmp]<bins[k+1]):
                region_I_y[tmp]=k
            if(coord_I_z[tmp]>=bins[k]):# and coord_I_x[tmp]<bins[k+1]):
                region_I_z[tmp]=k
        else:
            if(coord_I_x[tmp]>=bins[k] and coord_I_x[tmp]<bins[k+1]):
                region_I_x[tmp]=k
            if(coord_I_y[tmp]>=bins[k] and coord_I_y[tmp]<bins[k+1]):
                region_I_y[tmp]=k
            if(coord_I_z[tmp]>=bins[k] and coord_I_z[tmp]<bins[k+1]):
                region_I_z[tmp]=k

    
region_J_x=np.zeros(Nc) # region on x axis to which crosslinker belong to
region_J_y=np.zeros(Nc)
region_J_z=np.zeros(Nc)

for tmp in range(0,Nc):
    # assign region for crosslinkers
    for k in range(0,len(bins)): # index over elements of bins
        if(k==len(bins)-1):
            if(coord_J_x[tmp]>=bins[k]):# and coord_I_x[tmp]<bins[k+1]):
                region_J_x[tmp]=k
            if(coord_J_y[tmp]>=bins[k]):# and coord_I_x[tmp]<bins[k+1]):
                region_J_y[tmp]=k
            if(coord_J_z[tmp]>=bins[k]):# and coord_I_x[tmp]<bins[k+1]):
                region_J_z[tmp]=k
        else:
            if(coord_J_x[tmp]>=bins[k] and coord_J_x[tmp]<bins[k+1]):
                region_J_x[tmp]=k
            if(coord_J_y[tmp]>=bins[k] and coord_J_y[tmp]<bins[k+1]):
                region_J_y[tmp]=k
            if(coord_J_z[tmp]>=bins[k] and coord_J_z[tmp]<bins[k+1]):
                region_J_z[tmp]=k


PAR=np.zeros(2*Np,dtype='int32')#parity of the chain end, which indicates whether it has reacted or not, is designated as PAR
# PAR[I]=0 equivalent to conn_I[I][1-I%2]=-1
DEG=np.zeros(Nc,dtype='int32') # degree of visited ends
#DEG[J] equivalent to len(conn_J[J])
IND=np.zeros(Nc,dtype='int32') # connecting index = indices of chain ends to which the connections are made through
# not really required here     #the shared cross linker


# Connectivity data:
conn_J=[] #np.ones((Nc,f)) # connections of J: conn_J[j]= contains list of I indices (chain end indices to which it is connected)
for j in range(0,Nc):
    conn_J.append([])

##conn_J=conn_J*-1 # if -1: then not connected 
conn_I=np.ones((2*Np)) # connections of I: conn_I[i,0]= J_index to which end_1 of chain i is connected
                                         # conn_I[i,1]= J_index to which end_2 of chain i is connected
conn_I=conn_I*-1



G1=nx.Graph()


##STEP (14) -Sweep over regions
for region_x in range(0,n_seg):
    print(region_x)
    for region_y in range(0,n_seg):
        for region_z in range(0,n_seg):
            J_list=[] # list of Js belonging to current region
            I_list=[] # list of I's belonging to current region and neighboring regions
            for Jtmp in range(0,Nc):
                if(region_J_x[Jtmp]==region_x and region_J_y[Jtmp]==region_y and region_J_z[Jtmp]==region_z): ## this line could be replaced by coordinates and the previous few steps about assigning regions to each element can be omitted- if running into emmeory or time issues
                    J_list.append(Jtmp)
            for Itmp in range(0,2*Np): # odd I
                if(region_I_x[Itmp]in [region_x,region_x-1,region_x+1] and region_I_y[Itmp] in [region_y,region_y-1,region_y+1] and region_I_z[Itmp] in [region_z,region_z-1,region_z+1]):
                    # currently defining neighboring regions as the entire cubic shell around the main region- will change later if required
                    I_list.append(Itmp)
                    
            dist_set=[]#{'I':[],'J':[]}
            
            r_IJ=np.zeros((len(I_list),len(J_list)))
            for I_cnt in range(len(I_list)):
                for J_cnt in range(len(J_list)):
                    I=I_list[I_cnt]
                    J=J_list[J_cnt]
                    r_IJ[I_cnt][J_cnt]=dist(coord_I_x[I],coord_I_y[I],coord_I_z[I],coord_J_x[J],coord_J_y[J],coord_J_z[J],L,L,L)
                    if(r_IJ[I_cnt][J_cnt]<=d):
                        dist_set.append([I,J,r_IJ[I_cnt][J_cnt]])
                        
            dist_set=sorted(dist_set,key=last_value_list)
##            print(dist_set)
            
            #STEP (13)
            for idx in range(len(dist_set)):
                #STEP (11)
                
                I=dist_set[idx][0] # chain end not chain itself!
                
                J=dist_set[idx][1]
##                if(J==28):
##                    print('J=28')
##                    print('I',I)
##                    print('DEG[28]',DEG[28])
##                    print('PAR[I]',PAR[I])
##                    stop
##                print(J)
                if(PAR[I]==0 and DEG[J]<f):
##                    if(J==28):
##                        print('J=28 inside')
                    
##                    stop
                    PAR[I]=1 #should this be PAR[I] or PAR[J}???
                    DEG[J]=DEG[J]+1
                    conn_J[J].append(I)
##                    I_ch=int(I/2) # chain number, not chain end
##                    conn_I[I_ch][1-I%2]=J  #if odd: 1-I%2=0, if even: 1-I%2=1
                    conn_I[I]=J
##                    stop
########                    I1=I
########                    if(I1%2==0): # even 
########                        I2=I1+1
########                    else:
########                        I2=I1-1
########                    J2=conn_I[I2][1-I2%2]
########                    J1=J
########                    if(J2!=-1):
########                        if(J1!=J2):
########                            G1.add_edge(J1,J2)
########                        else:
########                            G1.add_node(J1)
##                            loop_atoms.append(J1)
##        Gmult.add_edge(J1,J2)
##                    if(I%2==1):# odd
##                        conn_I[I][0]=J
##                    else:# even
##                        conn_I[I][1]=J
                #STEP (12)
                if(IND[J]==0):
                    IND[J]=I
                else:
                    IND[J]=min(I,IND[J])
                

# this might not be required for my analysis
# STEP (15)
NEW=np.zeros(2*Np,dtype='int32')
for I in range(0,2*Np):
    for J in range(0,Nc):
        if(IND[J]==I): # is this equivalent to saying that Ith end is attached to the Jth crosslinker???
            NEW[I]=IND[J]
    if(PAR[I]==0):
        NEW[I]=I

#STEP (16)
END=np.zeros((2,Np),dtype='int32') # end number of chain
for N in range(0,Np):
    END[0,N]=min(NEW[2*N-1],2*N-1)
    END[1,N]=min(NEW[2*N],2*N)
    
# converting to graph object
##G=nx.Graph()
Gmult=nx.MultiGraph()

loop_atoms=[]
bonds_arr=[]

# look at each chain- and then see which junctions each of the chain ends are connected to
# form a edge between the two junctions J1 and J2
# this method will be able to avoid double counting of secondary loops

for I_ch in range(0,Np):
    I1=I_ch*2+1# odd
    I2=I_ch*2# even
    J1=conn_I[I1]
    J2=conn_I[I2]
    if((J1 !=-1) and (J2 !=-1)):
        bonds_arr.append([J1,J2])
        Gmult.add_edge(J1,J2)
    
##J_arr=list(range((0,Nc))) ## array from which J will be chosen for graph object
##visited_I=[]
##for J1 in range(0,Nc):
##    for I1 in conn_J[J1]: # I1 is chain end index
##        #1-I1%2: if I1%2=0, then even (1) and if I1%2=1, then odd (0) # end 1
##        # do for end 2: if I1%2=0, then end2 is odd(0) and if I1%2=1, then end2 is even(1)
##        # this means: odd-even index of end2 =I1%2
####        if(I1%2==0): # even 
####            I2=I1+1
####        else:
####            I2=I1-1
####        I_ch=int(I1/2)# chain index
##        if(I1%2==0):
##            I2=I1+1
##        else:
##            I2=I1-1
##        J2=conn_I[I2]
####        if(J2!=-1):
##        if(J2 !=-1):
##            if(J1!=J2):
##                G.add_edge(J1,J2)
##                
##            else:
##                G.add_node(J1)
##                loop_atoms.append(J1)
##            Gmult.add_edge(J1,J2)
##            bonds_arr.append([J1,J2])
####        else:
####            stop

    
bonds_arr=np.array(bonds_arr)
np.savetxt('network.txt',bonds_arr,delimiter=',',header='Node1,Node2')
##print('DEG',DEG)
##print(G.degree(np.arange(1,len(G))))
end=timer()
print('len(np.where(DEG==0)[0])',len(np.where(DEG==0)[0]))
print('len(Gmult)',len(Gmult))
print('time',end-start)

num_reacted=0
for i in range(len(conn_I)):
    if(conn_I[i]!=-1):
        num_reacted=num_reacted+1
print(num_reacted/(Np*2))
