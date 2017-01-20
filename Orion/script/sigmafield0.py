import pyfits
from pylab import *
import ephem

import math
import os,sys
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import mpl_toolkits
from matplotlib.patches import Ellipse

hdulist = pyfits.open('nh3_22.fits')
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

hdulist1 = pyfits.open('Tk.fits')
image1 = hdulist1[0].data
nx1 = hdulist1[0].header['naxis1']
ny1 = hdulist1[0].header['naxis2']
nz1 = hdulist1[0].header['naxis3']
crvalx1 = hdulist1[0].header['crval1'] 
cdeltax1 = hdulist1[0].header['cdelt1'] 
crpixx1 = hdulist1[0].header['crpix1'] 
crvaly1 = hdulist1[0].header['crval2'] 
cdeltay1 = hdulist1[0].header['cdelt2'] 
crpixy1 = hdulist1[0].header['crpix2'] 
crvalz1 = hdulist1[0].header['crval3'] 
cdeltaz1 = hdulist1[0].header['cdelt3'] 
crpixz1 = hdulist1[0].header['crpix3'] 

x1 = np.arange(-crpixx1*cdeltax1+crvalx1,(nx1-1-crpixx1)*cdeltax1+crvalx1,cdeltax1)
y1 = np.arange(-crpixy1*cdeltay1+crvaly1,(ny1-1-crpixy1)*cdeltay1+crvaly1,cdeltay1)


#print x
#print y
#print x1
#print y1

vfield=np.load('sigma0.npy')

for i in range(nx):
   for j in range(ny):
       xtemp=x[i]-crvalx
       ytemp=y[j]-crvaly
       #print xtemp,ytemp
       temprho=np.sqrt((xtemp*np.pi/180.0)**2+(ytemp*np.pi/180.0)**2)
       tempc=np.arcsin(temprho)
       dec=np.arcsin(np.cos(tempc)*np.sin(crvaly*np.pi/180.0)+ytemp*np.pi/180.0*np.sin(tempc)*np.cos(crvaly*np.pi/180.0)/temprho)
       ra=crvalx*np.pi/180.0+np.arctan((xtemp*np.pi/180.0*np.sin(tempc))/(temprho*np.cos(crvaly*np.pi/180.0)*np.cos(tempc)-ytemp*np.pi/180.0*np.sin(crvaly*np.pi/180.0)*np.sin(tempc)))
       #star = ephem.FixedBody()
       #star._ra=ra
       #star._dec=dec
       #star._epoch=ephem.J2000
       #star.compute(epoch=ephem.B1950)
       #star._epoch=ephem.B1950
       #star.compute(epoch=ephem.J2000)
       #ra1=float(star.ra)
       #dec1=float(star.dec)
       new=ephem.Equatorial(ra,dec,epoch=ephem.J2000)
       old=ephem.Equatorial(new,epoch=ephem.B1950)
       ra1=float(old.ra)
       dec1=float(old.dec)
       print ra1,dec1
       xtemp1=np.cos(dec1)*np.sin(ra1-crvalx1*np.pi/180.0)
       ytemp1=np.cos(crvaly1*np.pi/180.0)*np.sin(dec1)-np.sin(crvaly1*np.pi/180.0)*np.cos(dec1)*np.cos(ra1-crvalx1*np.pi/180.0)
       xtemp1=xtemp1/np.pi*180.0+crvalx1
       ytemp1=ytemp1/np.pi*180.0+crvaly1
       itemp=int((xtemp1-crvalx1)/cdeltax1+crpixx1)
       jtemp=int((ytemp1-crvaly1)/cdeltay1+crpixy1)
       #print ra,dec,ra1,dec1
       #print ra/np.pi*180/15,dec/np.pi*180,ra1/np.pi*180/15,dec1/np.pi*180
       print ra/np.pi*180,dec/np.pi*180,ra1/np.pi*180,dec1/np.pi*180
       print itemp,jtemp

ax = plt.subplot(111)
im = plt.imshow(vfield, cmap=cm.spectral
                ,origin='lower', aspect='equal'
                ,interpolation='none'
                ,extent=[max(x),min(x),min(y),max(y)])
xlabel('RA')
ylabel('Dec')
plt.colorbar(im,orientation='vertical')
savefig('sigmafield0.eps')
plt.show()
