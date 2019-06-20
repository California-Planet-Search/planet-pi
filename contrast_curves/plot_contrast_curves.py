import numpy as np 
import os.path
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 12})


def cont2mag(cont):
	mag = np.log10(cont)*(100**(1/5))
	return mag

hst_conts = np.transpose(np.loadtxt('HST_contrasts.tdt.txt', skiprows=1, usecols=(1,6,7)))
rad_h = hst_conts[0]
w6 = hst_conts[1]
w1 = hst_conts[2]

naco_conts = np.transpose(np.loadtxt('NaCo_contrasts.dat'))
rad_n = naco_conts[0]
naco_cont = naco_conts[1]/5.

"""
Note: I divide NaCo contrasts by 5 because their data reduction pipeline
assumes that a point source needs to be 5-simga above the noise to be detectable.
For this paper, we just plot 1-sigma detection limits.
"""

fig, ax = plt.subplots(figsize=(5,4))
plt.scatter(rad_h, w6, label='HST (W6)', s=25, facecolors='none', edgecolors='r', marker='o')
plt.scatter(rad_h, w1, label='HST (W1)', s=25, facecolors='none', edgecolors='k', marker='s')
plt.scatter(rad_n, naco_cont, label='NaCo', s=25, facecolors='none', edgecolors='b', marker='v')

ax.set_yscale('log')
ax.set_ylim(3e-9,3e-5)

ax2 = ax.twinx()
mn, mx = ax.get_ylim()
ax2.set_ylim(cont2mag(mn),cont2mag(mx))
ax2.set_ylabel('$\\Delta$ Magnitude (mag)')

ax.legend()
ax.set_ylabel('Contrast Limit')
ax.set_xlabel('Angular Separation [arcsec]')

fig.tight_layout()
plt.savefig('{}/Dropbox/Apps/Overleaf/120066/plots/contrast_curves.png'.format(os.path.expanduser('~')), dpi=250)