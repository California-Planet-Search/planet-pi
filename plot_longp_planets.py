import pandas as pd
import matplotlib.pylab as plt
from numpy import nan
import os.path

"""
Written by Lauren Weiss, 2018. Edited a bit by Sarah.
"""


def reformat_author_string(string_list): # useful for making citations
    authorstrings=[]
    for entry in string_list:
        authorsyear = entry.split('=')[1].replace('href','').replace('_ET_AL_','').split('_')
        while len(authorsyear) > 3:
            authorsyear.remove('')
        #    print authorsyear
        author = authorsyear[0].title() # first author, not all caps
        year = authorsyear[-1] # year
        authorstring="\citep{%s%s}"%(author,year)
        # print authorstring
        authorstrings.append(authorstring) # add this string
    return authorstrings

nea = pd.read_csv('data/planets_2023.csv',skiprows=94)
# print(nea.keys())
discmeths = pd.unique(nea.discoverymethod)
# hd12_dict = {"pl_hostname":["HR 5183 b"],
#                      "pl_orbper":73*365.25,
#                      "pl_orbper_unc":14*365.25,
#                      "pl_bmassj":3.23, # Mjup
#                      "pl_bmassj_unc":0.17,
#                      "pl_orbsmax":18, # AU
#                      "pl_orbsmax_unc":2,
#                      "pl_orbeccen": 0.83,
#                      "pl_orbeccen_unc":0.02,
#                      "K":"\semiamp", # m/s
#                      "pl_discmethod":"Radial Velocity",
#                      "disc_authorstring":'This paper',
#                      "def_authorstring":'This paper'}
# hd12 = pd.DataFrame(hd12_dict)
# hd12 = pd.DataFrame(hd12_dict)

#### Update other planets

# HR 8799, Wang et al. 2018
# nea.loc[nea.pl_name=='HR 8799 b','pl_orbsmax'] = 70.8
# nea.loc[nea.pl_name=='HR 8799 b','pl_orbeccen'] = 0.018
# nea.loc[nea.pl_name=='HR 8799 b','pl_bmassj'] = 5.8
# nea.loc[nea.pl_name=='HR 8799 c','pl_orbsmax'] = 43.1
# nea.loc[nea.pl_name=='HR 8799 c','pl_orbeccen'] = 0.022
# nea.loc[nea.pl_name=='HR 8799 c','pl_bmassj'] = 7.2

# nea.loc[nea.pl_name=='HR 8799 d','pl_orbsmax'] = 26.2
# nea.loc[nea.pl_name=='HR 8799 d','pl_orbeccen'] = 0.129
# nea.loc[nea.pl_name=='HR 8799 d','pl_bmassj'] = 7.2

# nea.loc[nea.pl_name=='HR 8799 e','pl_orbsmax'] = 16.2
# nea.loc[nea.pl_name=='HR 8799 e','pl_orbeccen'] = 0.118
# nea.loc[nea.pl_name=='HR 8799 e','pl_bmassj'] = 7.2


# recompute orbital periods with Kepler's law
# HR8799_mstar = 1.52
# nea.loc[nea.pl_name=='HR 8799 b','pl_orbper'] = pow(pow(nea.loc[nea.pl_name=='HR 8799 b','pl_orbsmax'],3.) / HR8799_mstar, 1/2.) * 365.25
# nea.loc[nea.pl_name=='HR 8799 c','pl_orbper'] = pow(pow(nea.loc[nea.pl_name=='HR 8799 c','pl_orbsmax'],3.) / HR8799_mstar, 1/2.) * 365.25
# nea.loc[nea.pl_name=='HR 8799 d','pl_orbper'] = pow(pow(nea.loc[nea.pl_name=='HR 8799 d','pl_orbsmax'],3.) / HR8799_mstar, 1/2.) * 365.25
# nea.loc[nea.pl_name=='HR 8799 e','pl_orbper'] = pow(pow(nea.loc[nea.pl_name=='HR 8799 e','pl_orbsmax'],3.) / HR8799_mstar, 1/2.) * 365.25

# update citation!
# nea.loc[nea.pl_hostname=='HR 8799','def_authorstring'] = "\citep{Wang2018}"

#### Make a table with the longest-period planets (> 20 years)
# longp = nea[nea.pl_orbper >= 365.25*20.] # 20 years and longer

#### Make a table with planets with semi-major axes between 10 and 100
# longp = nea[(nea.pl_orbsmax >= 10) & (nea.pl_orbsmax <= 100)] # 20 years and longer

# ### redo the author strings from nea format to latex format for long period planets
# longp['disc_authorstring'] = reformat_author_string(longp.pl_disc_reflink)
# longp['def_authorstring'] = reformat_author_string(longp.pl_def_reflink)
# longp.loc[longp.pl_hostname=='HR 8799','def_authorstring'] = "\citep{Wang2018}"

###


### Merge HD 120066 into the table
# longp = pd.concat([longp,hd12])

### Add some notes
# longp['notes'] = ['']*len(longp)
# longp.loc[longp.disc_authorstring=="\citep{Qian2010}",'notes'] = 'Binary star'

# ### Output table
# columns_of_interest=['pl_name','pl_orbper','pl_orbsmax','pl_bmassj','pl_orbeccen','pl_discmethod','disc_authorstring','def_authorstring','notes']
# print longp.loc[:,columns_of_interest].head()
# longp.loc[:,columns_of_interest].to_csv('long-period-planets.csv',index=False)


#### Make some nice plots! 

import itertools
marker = itertools.cycle(['o', '^','d', 'P', 's', 'v','p','H',"<",">"]) 
#color = itertools.cycle

### Plot 1: msini vs. zemimajor axis ###
import numpy as np
counter = 0
for meth in discmeths:
    if meth in ['Radial Velocity', 'Imaging', 'Transit']:
        planets = nea[nea.discoverymethod==meth]
        # planets = planets[planets.st_age < 0.5]
        planets = planets.dropna(subset=['pl_orbsmax','pl_bmassj'])

        # remove planets with large period uncertainty
        # if meth == 'Transit' or meth == 'Radial Velocity':
        #     planets = planets[planets['pl_orbper']/planets['pl_orbpererr1'] > 1.5]
        # planets = planets[np.abs(planets['st_age']/planets['st_ageerr1']) > 1.5]
        # print meth, len(planets)
        # print(planets['pl_bmassj'])
        color = 'C{}'.format(counter)
        if len(planets) > 1:
            counter += 1
            if counter == 3:
                counter += 1
            plt.scatter(#planets.pl_orbsmax,
                    planets.pl_orbsmax,
                    planets.pl_bmassj,
    #                 xerr=planets.pl_orbsmaxerr1,
    #                 yerr=planets.pl_bmassjerr1,fmt='o',
                    label = meth,
                    marker = marker.next(),
                    color=color, s=20, alpha=0.5)

plt.xlim([1e-2,1e4])
plt.ylim([2e-3,20])
# plt.scatter(hd12.pl_orbsmax,hd12.pl_bmassj, marker='*',color='C3',label='HR 5183 b',s=250)
plt.yscale('log')
plt.xscale('log')
plt.xlabel('Semi-major Axis [au]')
# plt.xlabel('Orbital Period [days]')
plt.ylabel('Mass or Msini [M$_{\\mathrm{{J}}}$]')
plt.legend()
plt.savefig('planets.png', dpi=250)
# plt.savefig('{}/Dropbox/Apps/Overleaf/120066/plots/nea_planets.pdf'.format(os.path.expanduser('~')),format='pdf',dpi=250)
#plt.close()
plt.clf()

### Plot 2: eccentricity vs. orbital period ###
# marker = itertools.cycle(['o', '^','d', 'P', '*', 'v','p','H',"<",">"]) 

# for meth in discmeths:
#     planets = nea[nea.discoverymethod==meth]
#     planets = planets.dropna(subset=['pl_orbper','pl_orbeccen'])
#     print meth, len(planets)
#     if len(planets) > 3:
#         plt.scatter(#planets.pl_orbsmax,
#                  planets.pl_orbper,
#                  planets.pl_orbeccen,
# #                 xerr=planets.pl_orbsmaxerr1,
# #                 yerr=planets.pl_bmassjerr1,fmt='o',
#                  label = meth,
#                  marker = marker.next())
#                  # color=color)
# plt.xlim([1,1e7])
# plt.ylim([0,1])
# plt.errorbar(hd12.pl_orbper,hd12.pl_orbeccen,hd12.pl_orbeccen_unc,hd12.pl_orbper_unc, marker='s',color='k',label='HR 5183 b')#,s=50)
# plt.xscale('log')
# #plt.xlabel('Semi-major Axis [AU]')
# plt.xlabel('Orbital Period [days]')
# plt.ylabel('Eccentricity')
# plt.legend(loc='lower right')
# # plt.savefig('nea_planets_ecc.pdf',format='pdf')
# #plt.show()
# #plt.close()
# plt.clf()
