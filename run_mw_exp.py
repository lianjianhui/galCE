'''
.. moduleauthor:: Jianhui Lian (732592593@qq.com)

General purpose:
................

The script to set up the parameters of the general GCE model and run galce.py
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

acc_initial = 0.1
acc_secular = 0.
acc_second = 0.1
acc_post = 0.0
sfe_initial = 1.8 #3
sfe_secular = 0.3
sfe_second = 0.8
sfe_postburst = 0.2
time0 = 1.5
time1 = 7
dt = 1
tau_sfe = 0.5
tau_initial = 7
tau_second = 100

gasAcc.exp_single(acc_initial, tau_initial, sfe_initial)

#Run the chemical evolution code
mgas0 = 0
imf = 'kr'
outflow = 0
nks = 1.5
mgenhance = 0
yields = 'l18'
dtd_1a_power = -1
dtd_1a_min = 0.15
dtd_nsm_min = 0.01
dtd_nsm_power = -1
mrsn_f = 0.004

filename = 'model3.fits'
outputFile = '/home/jianhui/projects/galce/tracks/'+filename
cemodel = galce.galCE(outputFile, gasAcc, mgas0, imf, outflow, twidth, nks, yields, dtd_1a_power, dtd_1a_min,dtd_nsm_min,dtd_nsm_power,mrsn_f)
cemodel.run()
print (outputFile)

res = subprocess.check_output(['python3','mock-single.py',filename])
print ('time used: ',time.time()-t0)
