'''
.. moduleauthor:: Jianhui Lian

General purpose:
................

The script to set up the parameters of the GCE model tuned for the Milky Way's bulge based on Lian et al. 2018c (arXiv.2007.12179).

'''

from pylab import *
from astropy.io import ascii
from astropy.io import fits
import galce 
import acc_setup as setup
import time
import subprocess

t0 = time.time()

# Set gas accretion history
twidth = 0.02
ageUni = 13.7
gasAcc = setup.gasAcc(twidth,ageUni)

acc_initial = 1
acc_secular = 0.
acc_second = 0.
acc_post = 0.0
sfe_initial = 0.6 #3
sfe_secular = 0.1
sfe_second = 1
sfe_postburst = 1
time0 = 1.5
time1 = 14
dt = 1
tau_sfe = 0.5
tau_initial = 100
tau_second = 100

gasAcc.mw_multi_burst(acc_initial, tau_initial, acc_secular, acc_second, tau_second, acc_post, sfe_initial, tau_sfe, sfe_secular, sfe_second, sfe_postburst, time0, time1, dt)

#Run the chemical evolution code
mgas0 = 0
imf = 'kr'
outflow = 0
nks = 1.5
mgenhance = 0
yields = 'l18'
dtd_1a_power = -1
dtd_1a_min = 0.15

filename = 'model1.fits'
outputFile = '/home/jianhui/projects/galce/tracks/'+filename
cemodel = galce.galCE(outputFile, gasAcc, mgas0, imf, outflow, twidth, nks, yields, dtd_1a_power, dtd_1a_min)
cemodel.run()
print (outputFile)

res = subprocess.check_output(['python3','mock-single.py',filename])
print ('time used: ',time.time()-t0)
