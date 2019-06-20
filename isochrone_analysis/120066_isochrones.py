import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm as nm

import corner
from isochrones import StarModel
from isochrones.mist import MIST_Isochrone

# params from Gaia archive, accessed Dec 2018
Teff = (4444.2, 160) 
logg = (3.0, 0.5)
feh = (0.0, 0.25)
parallax = (31.92,0.05)

# age from BJ's Specmatch-Gaia results
age = [9.9,0.1] # log10(age in yr)

mist = MIST_Isochrone()
model  = StarModel(mist, Teff=Teff, logg=logg, feh=feh, parallax=parallax, age=age, use_emcee=True)
model.fit(niter=1000)

# add in prior on age
age_post = model.samples.age_0

age_prior_probs = nm(age[0],age[1]).pdf(age_post)
compare = np.random.uniform(size=len(age_post))

new_samples = model.samples.loc[age_prior_probs > compare]
new_samples = new_samples[['mass_0_0','radius_0_0','feh_0','age_0','distance_0']]

# make figure
corner.corner(
	new_samples,
	labels=['Mass [Msol]','Radius [Rsol]','Fe/H','log10(age [Gyr])','distance [pc]'],
	quantiles=[0.16, 0.5, 0.84],
    show_titles=True
)
plt.savefig('isochrones_corner.png', dpi=250)

# save results
model.save_hdf('120066.hdf5', overwrite=True)


