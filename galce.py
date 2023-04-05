'''
.. moduleauthor:: Jianhui Lian (732592593@qq.com)

General purpose:
................

The class galCE is dedicated to run the chemical evolution calculation. It takes the gas accretion and star formation history parameters from 
acc_setup.py and other parameters from run_mw.py.

'''
from astropy.io import ascii
from pylab import *
from kimf import kimf
from astropy.table import Table
#from estimations_3d import estimation
from acc_setup import *
from astropy.io import fits

class galCE: 
 '''
  Notes
  -----

     .. note::
  *This is how it proceeds :*
  # reads the gas accretion and star formation history from galce_setup.py
  # read in the yields table and organize the yields as the same function of time with the model grid. 
  # integrate the metals released by each ejecta on the above time grid. 
  # iterate the chemical evolution, star formation, and gas accretion calculation
  # finally writes the output file.  
        # evolution time in Gyr
 '''

 def __init__(self, outputFile, gasAcc, mgas0=0, imf='kr', outflow=0, twidth=0.1, nks=1.5, yields_ccsn='l18', dtd_1a_power=-1, dtd_1a_min=0.15, 
              dtd_nsm_min=0.15, dtd_nsm_power=-1, mrsn_f=0):
  self.gasAcc = gasAcc
  self.outputFile = outputFile
  self.imf = imf
  self.fsn1a = 0.012
  self.twidth = twidth        
  self.age_uni = 13.7
  self.t = np.arange(0,self.age_uni,twidth)
  self.nks = nks
  self.mrange = np.arange(1200)*0.1+0.1 # to get more accurate estimate of total number of stars 
  self.mass_lft = np.append(np.append(np.arange(160)*0.025+1,np.arange(30)*0.5+5),np.arange(21)*5+20)
  self.area = 1.e6 #per square kpc^2
  self.mgenhance = 0
  zgrid_yields = np.array([1.4e-5,2.e-5,5.e-5,1.e-4,3.e-4,1.e-3,2.e-3,3.e-3,6.e-3,8.e-3,1.e-2,1.4e-2,2.e-2])
  self.zgrid = log10(zgrid_yields/0.014)
  self.yields_ccsn = yields_ccsn
  self.mgas0 = mgas0
  self.masscal = 0 #Boolean to calculate real-time stellar mass and integrated [M/H] and [alpha/Fe] or not 
  self.dtd_1a_min = dtd_1a_min
  self.dtd_1a_max = 21 
  self.dtd_1a_power = -1
  self.dtd_nsm_min = dtd_nsm_min
  self.dtd_nsm_max = 1.e6
  self.nsm_per_m = 2.e-5
  self.nsm_m_ej = 2.5e-2
  self.dtd_nsm_power = dtd_nsm_power
  self.mrsn_f = mrsn_f
  self.solar = ascii.read('yields/solar_G07.txt',data_end=84)
  self.outflow = outflow

  #SN-Ia and NSM DTD normalisation
  temp_sampling = np.arange(self.dtd_1a_min,self.dtd_1a_max,0.01)
  temp_sampling = 10**np.arange(log10(self.dtd_1a_min),log10(self.dtd_1a_max),0.01)
  dt_sampling = temp_sampling[1:]-temp_sampling[:-1]
  self.norm_1a_dtd = 1./np.sum(temp_sampling[:-1]**self.dtd_1a_power*dt_sampling)
  temp_sampling = 10**np.arange(log10(self.dtd_nsm_min),log10(self.dtd_nsm_max),0.01)
  dt_sampling = temp_sampling[1:]-temp_sampling[:-1]
  self.norm_nsm_dtd = 1./np.sum(temp_sampling[:-1]**self.dtd_nsm_power*dt_sampling)

 def get_yields(self):
  self.agb = fits.getdata('yields/AGB-Cristallo15-cube-galce.fits')

  if self.yields_ccsn=='l18':
   self.snii = fits.getdata('yields/CCSN-LC18-cube-R-ave-galce.fits')

  self.sn1a = ascii.read('yields/yields_sn1a_i99_list.txt')

  self.nsm = ascii.read('yields/yields_nsm_arnould07.txt')

  self.mrsn = ascii.read('yields/yields_mrsn_nishimura15.txt')

 def streamline_yields(self,dmass,t,m_h):
  if self.imf=='kr':
   slope1 = 1.3
   slope2 = 2.3
   nk = kimf(self.mrange,1,slope1,slope1,slope2)
   kk = sum(nk)
   numtot = kk*dmass
   num = kimf(self.mass_lft,dmass,slope1,slope1,slope2)
  id316 = (self.mrange>=3)&(self.mrange<=16)
  f316 = sum(nk[id316])

  #lifetime and delay time in time grid 
  lft = fits.getdata('yields/lifetime-cube-230226.fits')/1.e9
  idz_grid = np.argmin(abs(self.zgrid-m_h))
  lft_z = lft[:,idz_grid]
  self.lft_z = lft_z
  lft_tid = (lft_z[lft_z<self.age_uni-t]/self.twidth).astype(int) #can be improved to consider non even time interval 
  dtd_1a_id = np.arange(round((self.age_uni-t)/self.twidth))
  dtd_1a_sp = dtd_1a_id*self.twidth+self.twidth/2
  dtd_1a_valid = np.ones(len(dtd_1a_sp))
  dtd_1a_valid[dtd_1a_sp<=self.dtd_1a_min] = 0
  dtd_nsm_id = np.arange(round((self.age_uni-t)/self.twidth))

  #sample NSM DTD in a dense grid
  dtd_nsm_sp = np.arange(0,13.7-t,0.001)
  dtd_nsm_valid = np.ones(len(dtd_nsm_sp))
  dtd_nsm_valid[dtd_nsm_sp<=self.dtd_nsm_min] = 0

  t_id = round(t/self.twidth)

  weight_sn1a = numtot*self.fsn1a*f316*self.norm_1a_dtd*dtd_1a_sp**self.dtd_1a_power*self.twidth*dtd_1a_valid
  weight_nsm_dense = dmass*self.nsm_per_m*self.norm_nsm_dtd*dtd_nsm_sp**self.dtd_nsm_power*0.001*dtd_nsm_valid
  weight_nsm_dense[np.isnan(weight_nsm_dense)] = 0
  weight_nsm = np.zeros(round(self.age_uni/self.twidth)-t_id)

  id_mrsn = (self.mass_lft>=13)&(self.mass_lft<=25)
  num_mrsn = num
  num_mrsn[~id_mrsn] = 0
  weight_mrsn = num_mrsn*self.mrsn_f

  t_multi = int(self.twidth/0.001)
  for m in range(round(self.age_uni/self.twidth)-t_id):
      weight_nsm[m] = np.sum(weight_nsm_dense[t_multi*m:t_multi*(m+1)])*self.nsm_m_ej
  
  agbej = np.zeros((round(self.age_uni/self.twidth),84))
  sniiej = np.zeros((round(self.age_uni/self.twidth),84))
  sn1aej = np.zeros((round(self.age_uni/self.twidth),84))
  nsmej = np.zeros((round(self.age_uni/self.twidth),84))
  mrsnej = np.zeros((round(self.age_uni/self.twidth),84))

  for el in range(84):
   weights_agb = self.agb[el,:,idz_grid]*num
   agbej[t_id:t_id+np.max(lft_tid)+1,el] = np.bincount(lft_tid, weights=weights_agb[lft_z<self.age_uni-t])  

   weights_snii = self.snii[el,:,idz_grid]*num*(1-self.mrsn_f)
   sniiej[t_id:t_id+np.max(lft_tid)+1,el] = np.bincount(lft_tid, weights=weights_snii[lft_z<self.age_uni-t])

   sn1aej[t_id:,el] = np.bincount(dtd_1a_id, weights=self.sn1a['yields'][el]*weight_sn1a)

   nsmej[t_id:,el] = np.bincount(dtd_nsm_id, weights=self.nsm['yields'][el]*weight_nsm)
   
   mrsnej[t_id:t_id+np.max(lft_tid)+1,el] = np.bincount(lft_tid, weights=self.mrsn['yields'][el]*weight_mrsn[lft_z<self.age_uni-t])

  self.agbej = self.agbej + agbej
  self.sniiej = self.sniiej + sniiej
  self.sn1aej = self.sn1aej + sn1aej
  self.nsmej = self.nsmej + nsmej
  self.mrsnej = self.mrsnej + mrsnej

  #calculate the rate of events
  agb_mass = np.zeros(len(self.mass_lft))
  agb_mass[self.agb[0,:,0]>0] = 1
  agbrate = np.zeros(round(self.age_uni/self.twidth))
  weights_agbrate = agb_mass*num
  agbrate[t_id:t_id+np.max(lft_tid)+1] = np.bincount(lft_tid, weights=weights_agbrate[lft_z<self.age_uni-t])

  snii_mass = np.zeros(len(self.mass_lft))
  snii_mass[self.snii[0,:,0]>0] = 1
  sniirate = np.zeros(round(self.age_uni/self.twidth))
  weights_sniirate = snii_mass*num
  sniirate[t_id:t_id+np.max(lft_tid)+1] = np.bincount(lft_tid, weights=weights_sniirate[lft_z<self.age_uni-t])

  sn1arate = np.zeros(round(self.age_uni/self.twidth))
  sn1arate[t_id:] = np.bincount(dtd_1a_id, weights=weight_sn1a)

  nsmrate = np.zeros(round(self.age_uni/self.twidth))
  nsmrate[t_id:] = np.bincount(dtd_nsm_id, weights=weight_nsm/self.nsm_m_ej)

  mrsnrate = np.zeros(round(self.age_uni/self.twidth))
  mrsnrate[t_id:t_id+np.max(lft_tid)+1] = np.bincount(lft_tid, weights=weight_mrsn[lft_z<self.age_uni-t])

  self.agbrate = self.agbrate + agbrate
  self.sniirate = self.sniirate + sniirate
  self.sn1arate = self.sn1arate + sn1arate
  self.nsmrate = self.nsmrate + nsmrate
  self.mrsnrate = self.mrsnrate + mrsnrate

 def run(self):
  t0 = time.time()
  self.mgas = self.mgas0
  self.fe_h = -10
  self.mass_el = 0
  self.x_h = np.zeros(84)-10

  atomic_mass = ascii.read('yields/elements_mass.csv')
  
  #Mass-to-light ratio table from Claudia Maraston
  '''mlr = ascii.read('/Users/jianhui/project/manga/sfh/data/ML-ratio.txt')
  age = mlr['col2'][600:900]
  zh = mlr['col1'][600:900]
  ml = mlr['col4'][600:900]
  e_ml = estimation(age,zh,ml)'''

  output = np.zeros((len(self.t)-1,),dtype=[('t',float),('gasmass',float),('stellarmass',float),('sfr',float),('fe_h',float),('agb_rate',float),('snii_rate',float),('sn1a_rate',float),('nsm_rate',float),('mrsn_rate',float)]) 
  output_abun = np.zeros((len(self.t)-1,84))

  #set initial abundance
  self.x_h = np.zeros(84)-10
  self.get_yields()

  self.agbej = np.zeros((round(self.age_uni/self.twidth),84))
  self.sniiej = np.zeros((round(self.age_uni/self.twidth),84))
  self.sn1aej = np.zeros((round(self.age_uni/self.twidth),84))
  self.nsmej = np.zeros((round(self.age_uni/self.twidth),84))
  self.mrsnej = np.zeros((round(self.age_uni/self.twidth),84))
  self.agbrate = np.zeros(round(self.age_uni/self.twidth))
  self.sniirate = np.zeros(round(self.age_uni/self.twidth))
  self.sn1arate = np.zeros(round(self.age_uni/self.twidth))
  self.nsmrate = np.zeros(round(self.age_uni/self.twidth))
  self.mrsnrate = np.zeros(round(self.age_uni/self.twidth))

  for i in range(len(self.t)-1):
 
   #calculate star formation rate
   sur_gas = self.mgas/self.area
   sfr = self.gasAcc.sfeH[i]*2.5*1.e-4*0.75**self.nks*sur_gas**self.nks*self.area*1.e-6 #KS law, Kennicutt et al. 1998
   if sfr*(self.t[i+1]-self.t[i])*1.e9>self.mgas:
    sfr = self.mgas/((self.t[i+1]-self.t[i])*1.e9)
   dmass_star = sfr*(self.t[i+1]-self.t[i])*1.e9

   #calculate real-time stellar mass, light- ans mass-weighted stellar metallicity 
   #always calculate the final stellar mass
   avgz = np.zeros((i+1,4))
   if (self.masscal==0)&(i>=len(self.t)-2):
    for j in range(i):
     idy = self.lft_z <= ((self.t[i]-self.t[j])*1.e9) #identify the mass of a star with evolution time shorter than the age of Universe at t[i] 
     numz = kimf(self.mass_lft,10**output['sfr'][j]*(self.t[j+1]-self.t[j])*1.e9,1.3,1.3,2.3) #Kroupa IMF
     myoung = sum(numz[idy]*self.mass_lft[idy])
     avgz[j,0] = 10**output['sfr'][j]*(self.t[j+1]-self.t[j])*1.e9-myoung
    avgz[0,1] = 0
    avgz[0,2] = 0

   self.streamline_yields(dmass_star,self.t[i],self.fe_h)
   output['t'][i] = self.t[i]
   output['gasmass'][i] = round(log10(self.mgas),4)
   idinf = ~np.isinf(avgz[:,0])
   output['stellarmass'][i] = round(log10(np.nansum(avgz[idinf,0])),4)
   output['sfr'][i] = round(log10(sfr),4)
   output['agb_rate'][i] = self.agbrate[i]
   output['snii_rate'][i] = self.sniirate[i]
   output['sn1a_rate'][i] = self.sn1arate[i]
   output['nsm_rate'][i] = self.nsmrate[i]
   output['mrsn_rate'][i] = self.mrsnrate[i]

   #if self.masscal:
   # for j in range(i):
   #  idy = self.agb['lt'] <= log10((self.t[i]-self.t[j])*1.e9) #identify the mass of a star with evolution time shorter than the age of Universe at t[i] 
   #  numz = kimf(self.agb['Mi'],10**output['sfr'][j]*(self.t[j+1]-self.t[j])*1.e9,1.3,1.3,2.3) #Kroupa IMF
   #  myoung = sum(numz[idy]*self.agb['Mi'][idy])
   #  avgz[j,0] = 10**output['sfr'][j]*(self.t[j+1]-self.t[j])*1.e9-myoung
   #  avgz[j,1] = output['m/h'][j]
   #  avgz[j,2] = output['a/fe'][j]

     #Calculate luminosity of each ssp to estimate light-weighted properties
   '''z_table = output['m/h'][j]
     if z_table<1.e-10:
       z_table = 1.e-10
     ml_j = e_ml.estimate(self.t[i]-self.t[j],log10(z_table))
     avgz[j,3] = 10**output['sfr'][j]*(self.t[i+1]-self.t[i])*1.e9/ml_j'''
   # avgz[0,1] = 0
   # avgz[0,2] = 0


   #mass of each element returned to the ISM 
   rtn = self.agbej[i,:]+self.sniiej[i,:]+self.sn1aej[i,:]+self.nsmej[i,:]+self.mrsnej[i,:]
   #rtn = self.sniiej[i,:]+self.sn1aej[i,:]
   #print (self.t[i],self.sniiej[i,62],self.sn1aej[i,62],self.agbej[i,62],self.nsmej[i,62],self.fe_h)

   #current abundance in the ISM, atomic number ratio!
   if self.mgas>0:
    self.x_h = self.mass_el/self.mgas/0.75/atomic_mass['AtomicMass'][:84]
   
   dgas = self.gasAcc.accH[i]*(self.t[i+1]-self.t[i])*1.e9
   dmass_gas = (-1)*sfr*(self.t[i+1]-self.t[i])*1.e9+dgas+np.nansum(rtn)
   self.mgas = self.mgas+dmass_gas

   #mass and abundance of each element after enrichment
   self.mass_el = self.mass_el-self.mass_el/self.mgas*dmass_star+rtn

   #update the metal content after outflow
   self.mass_el = self.mass_el-self.outflow*sfr*(self.t[i+1]-self.t[i])*1.e9*self.mass_el/self.mgas
   output_abun[i,:] = log10(self.x_h)-self.solar['col2'][:85]+12

   #update the gas content after outflow
   self.mgas = self.mgas-self.outflow*sfr*(self.t[i+1]-self.t[i])*1.e9

   self.fe_h = output_abun[i,25]
   output['fe_h'][i] = self.fe_h

   #if self.masscal:
   # output['m/h_mw'][i] = round(np.nansum(avgz[:,0]*(avgz[:,1]))/np.sum(avgz[:,0]),4)
   # output['m/h_lw'][i] = round(np.nansum(avgz[:,3]*(avgz[:,1]))/np.sum(avgz[:,3]),4)
   # output['a/fe_mw'][i] = round(np.nansum(avgz[:,0]*(avgz[:,2]))/np.sum(avgz[:,0]),4)
   # output['a/fe_lw'][i] = round(np.nansum(avgz[:,3]*(avgz[:,2]))/np.sum(avgz[:,3]),4)

  image_hdu = fits.PrimaryHDU()
  hdu1 = fits.BinTableHDU(data=output)
  hdu2 = fits.ImageHDU(data=output_abun)
  hdu_list = fits.HDUList([image_hdu,hdu1,hdu2])
  hdu_list.writeto(self.outputFile,overwrite=True)
