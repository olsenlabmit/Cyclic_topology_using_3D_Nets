U0_low=10.0
U0_high=40.0

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

conc=cR3/(N*b**2)**1.5 #(chains/nm3)
L=(n_chains/conc)**(1/3)







K_low=1.0
K_high=1.0

fit_param_low=1.0
fit_param_high=1.0

E_b_low=1200.0
E_b_high=1200.0

func=4

tol=0.01
max_itr = 100000
write_itr = 10000
wrt_step = 5
