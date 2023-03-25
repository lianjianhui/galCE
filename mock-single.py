'''
.. moduleauthor:: Jianhui Lian

General purpose:
................

The script to generate a mock catalog of stars for a given GCE model.
'''

from astropy.io import ascii
from pylab import *
from kimf import kimf
import sys
from astropy.io import fits

#read in GCE model track
acch = 'multi-burst'
filename = sys.argv[1]
file = fits.open('tracks/'+acch+'/'+filename)
model = file[1].data
abun = file[2].data
feh = abun[:,25]
mgh = abun[:,11]
mgfe = mgh-feh
dt = model['t'][1:]-model['t'][:-1]

mass_lft = np.append(np.append(np.arange(180)*0.025+0.5,np.arange(30)*0.5+5),np.arange(21)*5+20)

#mass range of stars to be simulated
#caveat: this range might affect the sampling of stars 
mvalid = mass_lft[(mass_lft>=0.5)&(mass_lft<40)] 
lft = fits.getdata('yields/lifetime-cube-230313-lowm.fits')/1.e9
zgrid_yields = np.array([1.4e-5,2.e-5,5.e-5,1.e-4,3.e-4,1.e-3,2.e-3,3.e-3,6.e-3,8.e-3,1.e-2,1.4e-2,2.e-2])
zgrid = log10(zgrid_yields/0.014)

#probability distribution as a function of mass and age
pdf_ori = np.zeros((len(mvalid),len(model)))
for i in range(len(model)-1):
 numz = kimf(mvalid,10**model['sfr'][i]*dt[i]*1.e9,1.3,1.3,2.3)
 idz_lft = np.argmin(abs(feh[i]-zgrid))
 mlife = lft[:,idz_lft][(mass_lft>=0.5)&(mass_lft<40)]
 idy = (mlife >= (13.7-model['t'][i]))
 if (len(mlife[idy])>0):
  pdf_ori[idy,i] = numz[idy]
pdf_ori = pdf_ori/np.sum(pdf_ori)

#probability distribution as a function of evolution time (i.e. age) after considering stellar evolution
pdf_age = np.sum(pdf_ori,axis=0)
pdf_age_norm = pdf_age/np.sum(pdf_age)


mocksz = 5000 #number of stars in the mock catalog

#mock catalog
mock = np.zeros((mocksz,),dtype=[('age',float),('mass',float),('fe/h',float),('mg/fe',float)])

#sampling age distribution according to the probability distribution
randices_age = np.random.choice(np.arange(len(model['t'])),mocksz,replace=True,p=pdf_age_norm)

#uncertainties applied to the mock properties
sigma_age = 0.2
sigma_feh = 0.02
sigma_mgfe = 0.03

#age of mock stars
ages = 13.7-model['t'][randices_age]
mock['age'] = ages*10**np.random.normal(0,sigma_age,len(mock))

#resampling stars with age older than the universe to ensure all stars younger than the Universe
idiv = np.where(mock['age']>13.7)
for i in range(len(idiv[0])):
 oldp = exp(-0.5*(log10(13.7-model['t'])-log10(ages[idiv[0]][i]))**2/sigma_age**2)
 oldp = oldp/sum(oldp)
 newage = model['t'][np.random.choice(np.arange(len(model['t'])),1,replace=True,p=oldp)]
 mock['age'][idiv[0][i]] = 13.7-newage[0]

#element abunances of mock stars
model_feh = feh
model_mgfe = mgfe
mock['fe/h'][:] = feh[randices_age]+np.random.normal(0,sigma_feh,len(mock))
mock['mg/fe'][:] = mgfe[randices_age]+np.random.normal(0,sigma_mgfe,len(mock))

#mass of mock stars
birth_star = 13.7-mock['age']
for i in range(len(model)-1):
 idt = ((birth_star>=model['t'][i])&(birth_star<model['t'][i+1]))
 if len(birth_star[idt])>0:
  pdf_mass = pdf_ori[:,i]
  numz = kimf(mvalid,10**model['sfr'][i]*dt[i]*1.e9,1.3,1.3,2.3)
  idy = (mlife >= log10((13.7-model['t'][i])*1.e9))
  if (len(mlife[idy])>0):
   pdf_mass[idy] = numz[idy]
  if np.sum(pdf_mass)>0:
   pdf_mass = pdf_mass/np.sum(pdf_mass)
   mass = np.random.choice(mvalid,len(birth_star[idt]),replace=True,p=pdf_mass)
   mock['mass'][idt] = mass

#save the mock catalog
modelname = filename.split('fits')[0]
ascii.write(mock,'mocks/'+acch+'/'+modelname+'txt',overwrite=1)

