#!/usr/bin/env python

import pyfits
from pylab import *

import math
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import mpl_toolkits
from matplotlib.patches import Ellipse
#from matplotlib.patches import Wedge

#hdulist = pyfits.open('/mnt/data1/lqian/clumps/t13_new.fits')
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
#X, Y = np.meshgrid(x, y)

# image[z,x,y]!!!! The first dimension is velocity channel!
Z = image[0,:,:]
Z = Z*0.0
for i in range(0,nz):
  Z = Z+image[i,:,:]

# set negative value to zero
Z = (Z+abs(Z))/2.0

ax = plt.subplot(111)
im = plt.imshow(Z, cmap=cm.gist_yarg,
                origin='lower', aspect='equal',
                interpolation='none',
                extent=[max(x),min(x),min(y),max(y)])
locs, labels = plt.xticks()
xtv = locs
for i in range(0,len(locs)):
  xtv[i] = int((locs[i]-int(locs[i]/15.0)*15.0)/15.0*60.0)
#print xtv
x_ticks=[str(xtv[0]),str(xtv[1]),str(xtv[2]),str(xtv[3]),
          str(xtv[4]),'04h '+str(xtv[5])+' m' ,'04h '+str(xtv[6])+' m']
plt.gca().set_xticklabels(x_ticks)

xlabel('RA')
ylabel('Dec')

'''
ra,dec,major,minor,angle = np.loadtxt('/mnt/data1/lqian/clumps/newcut_modified_2012March16/txt/catalog.txt',unpack=True,usecols=[16,17,18,19,20])
xr = ra
yd = dec
angle=90.0-angle
for i in range(len(ra)):
    yd[i] = dec[i]
    xr[i] = (ra[i]-crvalx)*math.cos(dec[i]*3.1415926/180.0)+crvalx

ells = [Ellipse(xy=[xr[i],yd[i]],width=major[i]*2,height=minor[i]*2,angle=angle[i])
        for i in xrange(len(ra))]

for e in ells:
    ax.add_artist(e)
    e.set_alpha(1.0)
    e.set_clip_box(ax.bbox)
#    e.set_facecolor([0.0,1.0,0.0])
    e.set_edgecolor([0.0,1.0,0.0])
    e.set_fill(False)
    e.set_linewidth(0.5)

#plt.colorbar(im,orientation='horizontal')
plt.colorbar(im,orientation='vertical')

'''
plt.show()

