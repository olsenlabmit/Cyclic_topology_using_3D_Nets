U0_low=10.0
U0_high=40.0

N_low=12.0
N_high=12.0
b_low=1.0
b_high=1.0

N=N_low
b=b_low


cR3=10 # dimless conc
n_chains=10000

conc=cR3/(N*b**2)**1.5 #(chains/nm3)
L=(n_chains/conc)**(1/3)

print('L',L)
print('C_mM', C_mM)

