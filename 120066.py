import pandas as pd
import numpy as np
import radvel

# radvel fit -s 120066.py -d 2019-2-19/master
# radvel mcmc -s 120066.py -d 2019-2-19/master --maxGR 1.001 --minsteps 1000 --nsteps 10000 --minpercent 100

"""
"keywords"
"""
linearP = False  # fit in linear P
informative_per_prior = False # include A Vandenburg's informative period prior
informative_per_prior_with_dur = False # same as above, but assume event duration is 2015-mid-2018 (3.5yrs)
fit_recentpoints_only = False # do a fit to only data taken after 2017
vary_dvdt = False # include a trend
longP = True # enforce a 200 yr period to show that likelihood is lower
"""
"""

starname = 'HD120066'
nplanets = 2   
instnames = ['k', 'j', 'a', 'm']    
ntels = len(instnames)     
fitting_basis = 'logper tc secosw sesinw logk'
if linearP:
    fitting_basis = 'per tc secosw sesinw logk'
bjd0 = 2440000.

# stellar mass & error
stellar = dict(mstar=1.16, mstar_err=.12)

# load in data
data_cps = pd.read_csv('~/Dropbox/planet-pi/data/120066.txt')
data_mcd = pd.read_csv('~/Dropbox/planet-pi/data/HD120066_McD.ALL',
	names=['time','mnvel','errvel', 'SVAL','sval_err'], header=None,sep='\s+'
)
data_mcd['tel']='m'
data_mcd['time'] -= 40000.
data = pd.concat([data_cps, data_mcd], ignore_index=True)

if fit_recentpoints_only:
    data = data[data.time >= 17754.]
    instnames = ['j', 'a', 'm'] 

baseline = np.max(data.time.values) - np.min(data.time.values)
time_base = np.median(data.time.values) 

def initialize_params():
    params = radvel.Parameters(1,basis='per tp e w k')
    if longP:
        params['per1'] = radvel.Parameter(value=200.*365.)
    else:
        params['per1'] = radvel.Parameter(value=21860.)
    params['tp1'] = radvel.Parameter(value=18134.)
    params['e1'] = radvel.Parameter(value=0.84)
    params['w1'] = radvel.Parameter(value=-0.26)
    params['k1'] = radvel.Parameter(value=38.)
    params['dvdt'] = radvel.Parameter(value=0, vary=vary_dvdt)
    params['curv'] = radvel.Parameter(value=0, vary=False)
    
    # Convert input orbital parameters into the fitting basis
    params = params.basis.to_any_basis(params,fitting_basis)

    if longP:
        params['logper1'].vary = False
    
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
elif informative_per_prior_with_dur:
    duration = 3.5*365.25
    priors.append(radvel.prior.InformativeBaselinePrior('logper1', baseline, duration))