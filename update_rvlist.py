import cpsutils.io
import pandas as pd
import os

# sync velocities from cadence to local computer
_ = os.system('cps sync -t vel')

# sync S values from cadence to local computer
_ = os.system('cps sync -t sval')

# read in velocities
vst = cpsutils.io.loadcps('120066', getS=True)
vst_a = cpsutils.io.loadcps(
	'120066', getS=True, hires_rk=False, hires_rj=False, apf=True
)

# set 'tel' for each velocity: 'a' for apf, 'k' for pre-2004 HIRES, 'j' for post-2004 HIRES
hires_tel = [obnm[1] for obnm in vst.obnm]
vst['tel'] = hires_tel
vst_a['tel'] = 'a'

time_offset_from_jd = 2440000.

# concatenate HIRES & APF velocities
vst = pd.concat([vst,vst_a])
vst['time'] = vst.jd - time_offset_from_jd
vst = vst.sort_values(by='time')

# send velocities to csv file
vst.loc[:,['time','mnvel','errvel','tel','SVAL']].to_csv('data/120066.txt',index=False)