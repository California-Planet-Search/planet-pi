import numpy as np
from astropy.time import Time
import matplotlib.pyplot as plt

import radvel
from orbitize.kepler import _calc_ecc_anom

def calculate_true_anomaly(epoch, sma, ecc, tau, mtot):
	""" Calculate true anomaly

	Args:
		epoch (float): epoch at which to calculate contrast [mjd]
		sma (float): semimajor axis [au]
		ecc (float): eccentricity
		tau (float): epoch of periastron, in mjd, divided by period
		mtot (float): total mass [solar masses]

	Returns: true anomaly [float]

	History:
		2018: shamelessly stolen from orbitize.kepler.calc_orbit() by Sarah Blunt
	"""

	# compute period (from Kepler's third law) and mean motion
	period = np.sqrt((sma**3)/mtot)*365.25
	mean_motion = 2*np.pi/(period) # in rad/day

	# compute mean anomaly
	manom = (mean_motion*epoch - 2*np.pi*tau) % (2.0*np.pi)

	# compute eccentric anomalies
	eanom = _calc_ecc_anom(manom, ecc)

	# compute the true anomalies
	tanom = 2.*np.arctan(np.sqrt( (1.0 + ecc)/(1.0 - ecc))*np.tan(0.5*eanom))

	return tanom


def calculate_lambertian_contrast(epoch, sma, ecc, aop, inc, tau, mtot, proj_sep, radius=11.2):
	""" Calculate flux ratio for an object using the Lambertian disk approximation

	Args:
		epoch (astropy.time.Time): epoch at which to calculate contrast
		sma (float): semimajor axis [au]
		ecc (float): eccentricity
		aop (float): argument of periastron [radians]
		inc (float): inclination [radians]
		tau (float): epoch of periastron, in mjd, divided by period
		mtot (float): total mass [solar masses]
		proj_sep (float): projected separation between planet and star [au]
		radius (float [optional]): radius of planet, in Earth radii [1 R_J = 11.2 R_E]

	Returns: 
		tuple containing:
			flux ratio (float)
			phase angle (float) [radians]

	History:
		< 2018: written by Eric Nielsen
		2018: ported to Python by Sarah Blunt
	"""

	epoch = epoch.jd - 2400000.5

	albedo = 0.5 # Bob Brown's value for a Jupiter in the visible.

	true_anomaly = calculate_true_anomaly(epoch, sma, ecc, tau, mtot) # [radians]

	orbit_radius = sma*(1. - ecc**2) / (1. + ecc * np.cos(true_anomaly))

	z = -np.cos(aop) * np.sin(inc) * np.sin(true_anomaly) - np.cos(true_anomaly) * np.sin(inc) * np.sin(aop)
	z = z * orbit_radius # line-of-sight distance between star and planet [au]

	sep = np.sqrt(proj_sep**2 + z**2)  # total (not projected) separation between star and planet [au]

	out_beta = np.arctan2(-proj_sep,z) + np.pi  # angle between your eyeball, the planet, and the host star

	phase_function = (np.sin(out_beta) + (np.pi - out_beta) * np.cos(out_beta)) / np.pi # Lambert phase function

	fluxratio = ((6.371e8 * radius)**2) * albedo / (1.496e13 * sep)**2. * phase_function

	return fluxratio, out_beta

if __name__=='__main__':
	# 51 Eri b params
	epochs = Time(np.linspace(2020,2120,500),format='decimalyear')
	epp = Time(2016.0,format='decimalyear').jd - 2400000.5
	period = 60.68*365.25
	plx=33.57
	sma, ecc, inc, argp, tau, mtot, sep = (18., 0.28, np.radians(126.7), np.radians(4.2), epp/period, 1.75, 400./plx)
	fr, _ = calculate_lambertian_contrast(epochs, sma, ecc, argp, inc, tau, mtot, sep, radius=1.)
	plt.plot(fr)
	plt.show()

