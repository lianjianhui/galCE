'''
.. moduleauthor:: Jianhui Lian

General purpose:
................
The script to generate the gas accretion and star formation history. 
'''

from astropy.io import ascii
from pylab import *

class gasAcc:
	'''
	class object of gas accretion and star formation history. 
	Evolution history of gas outflow will be included in the near future
	'''

	def __init__(self, twidth=0.1, ageUni=13.7):
		self.twidth = twidth
		self.ageUni = ageUni

	def mw_multi_burst(self, acc_initial, tau_initial, acc_secular, acc_second, tau_second, acc_post, sfe_initial, tau_sfe, sfe_secular, sfe_second, sfe_postburst, time0=2, time1=7, dt=1):
		'''
		This method uses the gas accretion and star formation history framework that consists of four phases taken from Lian et al. 2020, MNRAS, 497, 2371	
		Parameters:
			tau_sfe: the timescale of SFE declining during the first quenching episode. 
			time0: transition time when initial gas accretion ended.
			time1: transition time when second gas accretion started
			time2: transition time when second gas accretion ended.
		'''
		t = np.arange(0,self.ageUni,self.twidth)
		self.accH = np.zeros(len(t)) #Gas accretion hitory
		self.sfeH = np.zeros(len(t)) #Evolutionary history of star formation efficiency

		#Initial gas accretion and star burst
		id0 =  t<=time0
		self.accH[id0] = acc_initial*exp(-1*(t[id0])/tau_initial)
		self.sfeH[id0] = sfe_initial

		#Quenching + secular evolution phase		
		id1 = (t>time0)&(t<=time1)
		self.accH[id1] = acc_secular#*exp(-1*(t[id1]-time0)/5)
		if sfe_secular<sfe_initial:
			sfeH_id1 = sfe_initial*exp(-1*(t[id1]-time0)/tau_sfe) #Assuming SFE decrease exponentially until reaching the threshold set for the secular phase. 
			idsfe1 = np.where(sfeH_id1<sfe_secular)
			sfeH_id1[idsfe1] = sfe_secular
		else:
			sfeH_id1 = sfe_secular
		self.sfeH[id1] = sfeH_id1

		#Second accretion and burst 
		id2 = (t>time1)&(t<=time1+dt)
		self.accH[id2] = acc_second#*exp(-1*(t[id2]-time1)/tau_second)
		self.sfeH[id2] = sfe_second

		#Post starburst after the second accretion
		id3 = (t>time1+dt)
		self.accH[id3] = acc_post
		self.sfeH[id3] = sfe_postburst
		
