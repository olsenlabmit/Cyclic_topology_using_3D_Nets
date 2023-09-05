U0=56.7578
U0_low=U0
U0_high=U0
N_low=12.0
N_high=12.0
b_low=1.0
b_high=1.0

N=N_low
b=b_low

##frac_weak=0.0

#L=10 
##conc=2.0 #(chains/nm3)
cR3=10 # dimless conc
n_chains=10000
func=4

conc=cR3/(N*b**2)**1.5 #(chains/nm3)
L=(n_chains/conc)**(1/3)

