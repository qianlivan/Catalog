import pyfits
from pylab import *

import math
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
z = np.arange(-crpixz*cdeltaz+crvalz,(nz-1-crpixz)*cdeltaz+crvalz,cdeltaz)
z=np.arange(nz)*1.0

#np.savetxt('freq.txt',z/1000.0,fmt='%s')
np.savetxt('freq.txt',z,fmt='%s')

#nx=10
#ny=10

for ix in range(nx):
    for iy in range(ny):
        spec = image[:,iy,ix]
        filename=str(ix)+'_'+str(iy)+'.txt'
        np.savetxt(filename,spec,fmt='%s')

