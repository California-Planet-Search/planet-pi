import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
import os

from galpy.orbit import Orbit
from galpy.potential import MWPotential2014

# compute this many orbits
num_samples = 50

# params from Gaia archive, accessed Dec 2018
par_med = 31.75676467890995
ra_med = 206.73579924515266
dec_med = 6.349899363461747
pmra_med = -510.4469988734486
pmdec_med = -110.22479327853173
rv_med = -30.41781739548299

par = np.random.normal(par_med,0.03904992068059236, size=num_samples) # [mas]
ra = np.random.normal(ra_med,0.034181020450823614, size=num_samples) # [deg]
dec = np.random.normal(dec_med,0.029049724029680198, size=num_samples) # [deg]
pmra = np.random.normal(pmra_med,0.07061308080526092, size=num_samples) # [mas/yr]
pmdec = np.random.normal(pmdec_med,0.06383967745385508, size=num_samples) # [mas/yr]
rv = np.random.normal(rv_med,0.19533398113754827, size=num_samples) # km/s

ts = np.linspace(0.,2.,2500.)*u.Gyr

o1 = Orbit(vxvv=[ra_med,dec_med,1/par_med,pmra_med,pmdec_med,rv_med], radec=True)
o1.integrate(ts, MWPotential2014)
o1.plot(color='red', linewidth=1.0)

for i in np.arange(num_samples):

	o = Orbit(vxvv=[ra[i],dec[i],1/par[i],pmra[i],pmdec[i],rv[i]], radec=True)
	o.integrate(ts, MWPotential2014)

	o.plot(color='grey', overplot=True, alpha=0.1, linewidth=0.3)

plt.plot([o1.R()],[o1.z()],'b*', zorder=2, markersize=20)

plt.savefig(
	'{}/Dropbox/Apps/Overleaf/120066/plots/galactic_orbit.png'.format(os.path.expanduser('~')),
	dpi=250
)
