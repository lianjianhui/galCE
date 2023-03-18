'''
.. moduleauthor:: Jianhui Lian (ljh520hw@gmail.com)

General purpose:
................

The script to set up the parameters and run galce.py
'''

from pylab import *
from astropy.io import ascii
from astropy.io import fits
import galce 
import acc_setup as setup
import time

t0 = time.time()

# Gas accretion history
twidth = 0.02
ageUni = 13.7
gasAcc = setup.gasAcc(twidth,ageUni)

acc_initial = 0.8
acc_secular = 0.
acc_second = 0.0316
acc_post = 0.0
sfe_initial = 1 #3
sfe_secular = 0.631
sfe_second = 0.79
sfe_postburst = 0.25
time0 = 0.4
time1 = 7
dt = 1
tau_sfe = 0.5
tau_initial = 0.1
tau_second = 100

gasAcc.mw_multi_burst(acc_initial, tau_initial, acc_secular, acc_second, tau_second, acc_post, sfe_initial, tau_sfe, sfe_secular, sfe_second, sfe_postburst, time0, time1, dt)
#print (gasAcc.sfeH)

#Run the chemical evolution code
mgas0 = 0
imf = 'kr'
outflow = 1
fsn1a = 0.012
nks = 1.5
mgenhance = 0
yields = 'l18'
mass_1a_up = 4.5
dtd_1a_power = -1
dtd_1a_min = 0.15
masscal = 0
outputFile = '/home/jianhui/projects/galce/tracks/multi-burst/'+'%.4f' % acc_initial+'-%.4f' % acc_secular+'-%.4f' % acc_second+'-%.4f' % acc_post+'-%.2f' % sfe_initial+'-%.3f' % sfe_secular+'-%.2f' % sfe_second+'-%.2f' % sfe_postburst+'-%.1f' % time0+'-%.1f' % time1+'-%.1f' % dt+'-%.1f' % tau_sfe+'-%.1f' % tau_second+'-%.1f' % outflow+'.fits'
cemodel = galce.galCE(outputFile, gasAcc, mgas0, imf, outflow, twidth, nks, yields, dtd_1a_power, dtd_1a_min)
cemodel.run()
print (outputFile)

print ('time used: ',time.time()-t0)
