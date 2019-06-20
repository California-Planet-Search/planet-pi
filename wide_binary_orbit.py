import numpy as np
import os.path

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from astropy import units as u
from astropy.constants import G

mas2deg = 2.7777777777778e-7

"""
Reference: Sect 3.4 of Jack Wisdom's Lecture Notes for MIT Class 12.603, Spring 2019
"""

num_samples = int(1e7)

# params from Gaia archive, accessed Dec 2018
par_1 = np.random.normal(31.75676467890995,0.03904992068059236, size=num_samples)*u.mas
par_2 = np.random.normal(31.916619970689208,0.04588075116786399, size=num_samples)*u.mas
ra1 = np.random.normal(206.73579924515266,0.034181020450823614*mas2deg, size=num_samples)*u.deg.to('', equivalencies=u.dimensionless_angles())
ra2 = np.random.normal(206.86779239249455,0.041854907527526505*mas2deg, size=num_samples)*u.deg.to('', equivalencies=u.dimensionless_angles())
dec1 = np.random.normal(6.349899363461747,0.029049724029680198*mas2deg, size=num_samples)*u.deg.to('', equivalencies=u.dimensionless_angles())
dec2 = np.random.normal(6.315179418424761,0.027597433308644042*mas2deg, size=num_samples)*u.deg.to('', equivalencies=u.dimensionless_angles())
pmra1 = np.random.normal(-510.4469988734486,0.07061308080526092, size=num_samples)*u.mas.to('', equivalencies=u.dimensionless_angles())/u.yr
pmra2 = np.random.normal(-509.44024258805166,0.0833659380642554, size=num_samples)*u.mas.to('', equivalencies=u.dimensionless_angles())/u.yr
pmdec1 = np.random.normal(-110.22479327853173,0.06383967745385508, size=num_samples)*u.mas.to('', equivalencies=u.dimensionless_angles())/u.yr
pmdec2 = np.random.normal(-111.02152731289831,0.06053443317129055, size=num_samples)*u.mas.to('', equivalencies=u.dimensionless_angles())/u.yr
rv1 = np.random.normal(-30.41781739548299,0.19533398113754827, size=num_samples)*u.kilometer/u.second
rv2 = np.random.normal(-30.668597020721364,0.1500618955235784, size=num_samples)*u.kilometer/u.second

# params from isochrones & specmatch-Gaia analysis
m1 = np.random.normal(1.07,0.04, size=num_samples)*u.M_sun
m2 = np.random.normal(0.67,0.05, size=num_samples)*u.M_sun

dist_1 = par_1.to(u.parsec, equivalencies=u.parallax())
dist_2 = par_2.to(u.parsec, equivalencies=u.parallax())

du_z = rv1 - rv2
dz = dist_1 - dist_2

dx = (ra1 - ra2)*dist_1
dy = (dec1 - dec2)*dist_1

du_x = (pmra1 - pmra2)*dist_1
du_y = (pmdec1 - pmdec2)*dist_1

# compute mu & reduced mass
mu = G*m1*m2
m = 1./(1./m1 + 1./m2)

# compute specific angular momentum (|r x v|)
h = np.sqrt((dy*du_z - dz*du_y)**2 + (dz*du_x - dx*du_z)**2 + (dx*du_y - dy*du_x)**2)

# compute radius
r = np.sqrt(dx**2 + dy**2 + dz**2) # [cm]

# compute total velocity
v = np.sqrt(du_x**2 + du_y**2 + du_z**2)

# METHOD 1:
# =========

# compute semi-major axis
k = mu/m
a = 1./(2./r - v**2/k)
a_au = a.to(u.au)

# compute eccentricity
p = (h**2)/k
rdot = (dx*du_x + dy*du_y + dz*du_z)/r
ecostrueanom = p/r - 1
esintrueanom = h*rdot/k
e = np.sqrt(ecostrueanom**2 + esintrueanom**2)

# compute period
n = np.sqrt(mu/(m*(a**3)))
per = (2.*np.pi/n).to(u.megayear)
# =========



# METHOD 2:
# =========

# mu = G*(m1+m2)

# # compute specific energy
# E = 0.5*v**2 - mu/r

# # compute semi-major axis
# a = -0.5*mu/E
# a_au = a.to(u.au)

# # compute eccentricity
# e = np.sqrt(1 - (h**2)/(a*mu))

# # compute period
# per = np.sqrt(4.*(np.pi**2)*(a**3)/mu).to(u.megayear)
# =========


print('Percentage of orbits that are bound: {}%'.format(100.*len(e[e<1])/num_samples))

# compute inclination
i = np.arccos(((dx*du_y - dy*du_x)/h).decompose()).to(u.degree)

fig = plt.figure(figsize=(4,11))
fig.subplots_adjust(hspace=.35)
gs = gridspec.GridSpec(4,1, figure=fig)
plt.rcParams.update({'font.size': 12})
bins = 200
color = 'grey'

ax0 = plt.subplot(gs[0])
max_a = 205
ax0.hist(a_au.value[e<1]/1e3, bins=bins, range=(0,max_a), color=color)
ax0.set_xlim(0,max_a)
ax0.set_yticks([])
ax0.set_xlabel('a [x10$^{3}$ au]')

ax1 = plt.subplot(gs[1])
ax1.hist(e.value[e<1], bins=bins, color=color)
ax1.set_yticks([])
ax1.set_xlabel('e')

ax2 = plt.subplot(gs[2])
ax2.hist(i.value[e<1], bins=bins, color=color)
ax2.set_yticks([])
ax2.set_xlim(50,150)
ax2.set_xlabel('inc [$^{\circ}$]')

ax3 = plt.subplot(gs[3])
n, bins, _ = ax3.hist(a_au.value[e<1]*(1-e.value[e<1])/1e3, bins=bins, color=color)
ax3.set_yticks([])
ax3.set_xlabel('a(1-e) [x10$^{3}$ au]')
max_minsep = 50
ax3.set_xlim(0,max_minsep)

print('Peak of the a(1-e) histogram occurs at {} au.'.format(1e3*bins[np.argmax(n)]))
print('Median of the a(1-e) histogram occurs at {} au.'.format(
	np.median(a_au.value[e<1]*(1-e.value[e<1]))
))

fig = fig.tight_layout()

plt.savefig(
	'{}/Dropbox/Apps/Overleaf/120066/plots/wide_companion.png'.format(os.path.expanduser('~')),
	dpi=250
)


