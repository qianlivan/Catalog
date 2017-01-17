import pyfits
from pylab import *

import math
import os,sys
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import mpl_toolkits
from matplotlib.patches import Ellipse


hdulist = pyfits.open('t13_new.fits')
image = hdulist[0].data
nx = hdulist[0].header['naxis1']
ny = hdulist[0].header['naxis2']
nz = hdulist[0].header['naxis3']
crvalx = hdulist[0].header['crval1'] 
cdeltax = hdulist[0].header['cdelt1'] 
crpixx = hdulist[0].header['crpix1'] 
crvaly = hdulist[0].header['crval2'] 
cdeltay = hdulist[0].header['cdelt2'] 
crpixy = hdulist[0].header['crpix2'] 
crvalz = hdulist[0].header['crval3'] 
cdeltaz = hdulist[0].header['cdelt3'] 
crpixz = hdulist[0].header['crpix3'] 

x = np.arange(-crpixx*cdeltax+crvalx,(nx-1-crpixx)*cdeltax+crvalx,cdeltax)
y = np.arange(-crpixy*cdeltay+crvaly,(ny-1-crpixy)*cdeltay+crvaly,cdeltay)


vfield=np.zeros([ny,nx])
vfieldtemp=np.load('sigma1_0_499.npy')
vfield=vfield+vfieldtemp
vfieldtemp=np.load('sigma1_500_999.npy')
vfield=vfield+vfieldtemp
vfieldtemp=np.load('sigma1_1000_1499.npy')
vfield=vfield+vfieldtemp
vfieldtemp=np.load('sigma1_1500_2068.npy')
vfield=vfield+vfieldtemp
print size(vfieldtemp[0,:]),size(vfieldtemp[:,0])


vfield[vfield<0]=0
#vfield[vfield>12.0]=0

os.system('rm -f sigma1.fits')
hduout=pyfits.PrimaryHDU(vfield)
hdulistout=pyfits.HDUList([hduout])
hdulistout.writeto('sigma1.fits')

ax = plt.subplot(111)
#im = plt.imshow(vfield, cmap=cm.gist_heat
#im = plt.imshow(vfield, cmap=cm.rainbow
im = plt.imshow(vfield, cmap=cm.spectral
                ,origin='lower', aspect='equal'
                ,interpolation='none')
xlabel('RA')
ylabel('Dec')
plt.colorbar(im,orientation='vertical')
savefig('sigmafield1.eps')
plt.show()
