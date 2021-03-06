import pandas as pd
import numpy as np
import radvel

# radvel fit -s 120066.py -d results/master
# radvel mcmc -s 120066.py -d results/master --maxGR 1.001 --minsteps 1000 --nsteps 100000 --minpercent 100 --thin 25

"""
"keywords"
"""
linearP = False  # fit using uniform priors
informative_per_prior = False # include A Vandenburg's informative period prior

vary_dvdt = False # include a trend
"""
"""

starname = 'HD120066'
nplanets = 2   
instnames = ['k', 'j', 'a', 'm']    
ntels = len(instnames)     
fitting_basis = 'logper tc secosw sesinw logk'
if linearP:
    fitting_basis = 'per tc secosw sesinw k'
bjd0 = 2440000.

# stellar mass & error
stellar = dict(mstar=1.16, mstar_err=.12)

# load in data
data_cps = pd.read_csv('data/120066.txt')
data_mcd = pd.read_csv('data/HD120066_McD.ALL',
	names=['time','mnvel','errvel', 'SVAL','sval_err'], header=None,sep='\s+'
)
data_mcd['tel']='m'
data_mcd['time'] -= 40000.
data = pd.concat([data_cps, data_mcd], ignore_index=True)

baseline = np.max(data.time.values) - np.min(data.time.values)
time_base = np.median(data.time.values) 

def initialize_params():
    params = radvel.Parameters(1,basis='per tp e w k')
    params['per1'] = radvel.Parameter(value=25962.)
    params['tp1'] = radvel.Parameter(value=18134.)
    params['e1'] = radvel.Parameter(value=0.84)
    params['w1'] = radvel.Parameter(value=-0.26)
    params['k1'] = radvel.Parameter(value=38.)
    params['dvdt'] = radvel.Parameter(value=0, vary=vary_dvdt)
    params['curv'] = radvel.Parameter(value=0, vary=False)
    
    # Convert input orbital parameters into the fitting basis
    params = params.basis.to_any_basis(params,fitting_basis)
    
    return params

# initialize the orbit parameters and the orbit model
params = initialize_params()
params['gamma_j'] = radvel.Parameter(value=-44.57)
params['jit_j'] = radvel.Parameter(value=2.56)
params['gamma_k'] = radvel.Parameter(value=-45.33)
params['jit_k'] = radvel.Parameter(value=3.36)
params['gamma_a'] = radvel.Parameter(value=-38.76)
params['jit_a'] = radvel.Parameter(value=4.03)
params['gamma_m'] = radvel.Parameter(value=-5.07)
params['jit_m'] = radvel.Parameter(value=6.13)

priors = [
    radvel.prior.EccentricityPrior( 1 ), # Keeps eccentricity < 1
    radvel.prior.HardBounds('jit_k', 0.0, 10.0),
    radvel.prior.HardBounds('jit_j', 0.0, 10.0),
    radvel.prior.HardBounds('jit_a', 0.0, 10.0),
    radvel.prior.HardBounds('jit_m', 0.0, 10.0),
] 


if informative_per_prior:
   priors.append(radvel.prior.InformativeBaselinePrior('logper1', baseline))
