U0=56.7578
U0_low=U0
U0_high=U0
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
C_mM=conc/0.6022 

K_low=1.0
K_high=1.0
K=K_low
fit_param_low=1.0
fit_param_high=1.0

E_b_low=1200.0
E_b_high=1200.0

func=4
del_t=0.04
num_exchange_steps=n_chains*50
wrt_step = int(num_exchange_steps/50)
