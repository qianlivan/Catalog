#!/usr/bin/env python

#import msvcrt
import pyfits
from pylab import *

import math
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import mpl_toolkits
from matplotlib.patches import Ellipse
from lmfit.models import GaussianModel,ExponentialModel
import sys
import matplotlib.pyplot as plt



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

z = np.arange(-crpixz*cdeltaz+crvalz,(nz-crpixz)*cdeltaz+crvalz,cdeltaz)

x = z/1000.0

#print size(x)

velocity0=np.zeros([ny,nx])
velocity1=np.zeros([ny,nx])
velocity2=np.zeros([ny,nx])

amplitude0=np.zeros([ny,nx])
amplitude1=np.zeros([ny,nx])
amplitude2=np.zeros([ny,nx])

sigma0=np.zeros([ny,nx])
sigma1=np.zeros([ny,nx])
sigma2=np.zeros([ny,nx])


outfile=file('Taurus_velocity.txt','a')

for ix in range(nx):
#for ix in range(10):
    print ix
    g0_c=0.0
    g0_cmin=-3.0
    g0_cmax=13.0
    g0_sig=0.5
    g0_sigmin=0.1
    g0_sigmax=1.0
    g0_amp=1.0
    g0_ampmin=0.1
    g1_c=1.0
    g1_cmin=-3.0
    g1_cmax=13.0
    g1_sig=0.5
    g1_sigmin=0.1
    g1_sigmax=1.0
    g1_amp=1.0
    g1_ampmin=0.1
    g2_c=2.0
    g2_cmin=-3.0
    g2_cmax=13.0
    g2_sig=0.5
    g2_sigmin=0.1
    g2_sigmax=1.0
    g2_amp=1.0
    g2_ampmin=0.1
    for iy in range(ny):
    #for iy in range(10):
	print ix,iy
        y = image[:,iy,ix]
	#print size(y)
	noise=np.std(y[0:10])
	peak=np.max(y)
	index_y=y>(peak-0.1)
	tempy=y[index_y]
	tempx=x[index_y]
        g0_amp=tempy[0]
        g0_c=tempx[0]
        g1_amp=tempy[0]+np.random.rand()
        g1_c=tempx[0]+np.random.rand()
        g2_amp=tempy[0]+np.random.rand()
        g2_c=tempx[0]+np.random.rand()
	if(peak<5*noise):
	    velocity0[iy,ix]=-5.0
	    velocity1[iy,ix]=-5.0
	    velocity2[iy,ix]=-5.0
            amplitude0[iy,ix]=0.0
            amplitude1[iy,ix]=0.0
            amplitude2[iy,ix]=0.0
            sigma0[iy,ix]=0.0
            sigma1[iy,ix]=0.0
            sigma2[iy,ix]=0.0
        else:
            gauss0 = GaussianModel(prefix='g0_')
            pars = gauss0.make_params()
            pars['g0_center'].set(g0_c, min=g0_cmin, max=g0_cmax)
            pars['g0_sigma'].set(g0_sig, min=g0_sigmin, max=g0_sigmax)
            pars['g0_amplitude'].set(g0_amp, min=g0_ampmin)
            
            gauss1 = GaussianModel(prefix='g1_')
            pars.update(gauss1.make_params())
            pars['g1_center'].set(g1_c, min=g1_cmin, max=g1_cmax)
            pars['g1_sigma'].set(g1_sig, min=g1_sigmin, max=g1_sigmax)
            pars['g1_amplitude'].set(g1_amp, min=g1_ampmin)
            
            gauss2  = GaussianModel(prefix='g2_')
            pars.update(gauss2.make_params())
            pars['g2_center'].set(g2_c, min=g2_cmin, max=g2_cmax)
            pars['g2_sigma'].set(g2_sig, min=g2_sigmin, max=g2_sigmax)
            pars['g2_amplitude'].set(g2_amp, min=g2_ampmin)
            
            mod = gauss0+gauss1+gauss2
            init = mod.eval(pars,x=x)
            #plt.plot(x,init,'k--')
  
            out = mod.fit(y,pars,x=x)
            #help(out)
	    #print out.params
            #sys.exit()
	    temp00=out.params['g0_amplitude'].value
	    temp01=out.params['g0_center'].value
	    temp02=out.params['g0_sigma'].value
	    temp10=out.params['g1_amplitude'].value
	    temp11=out.params['g1_center'].value
	    temp12=out.params['g1_sigma'].value
	    temp20=out.params['g2_amplitude'].value
	    temp21=out.params['g2_center'].value
	    temp22=out.params['g2_sigma'].value
	    #temp=[[temp01,temp00],[temp11,temp10],[temp21,temp20]]
	    temp=[[temp00,temp01,temp02],[temp10,temp11,temp12],[temp20,temp21,temp22]]
	    temp.sort(reverse=True)
	    print temp
	    print peak
	    g0_amp=temp[0][0]
	    g0_c=temp[0][1]
	    g0_sig=temp[0][2]
	    g1_amp=temp[1][0]
	    g1_c=temp[1][1]
	    g1_sig=temp[1][2]
	    g2_amp=temp[2][0]
	    g2_c=temp[2][1]
	    g2_sig=temp[2][2]
	    if (g0_amp>0.15 and g0_sig<0.95):
	        velocity0[iy,ix]=g0_c
                amplitude0[iy,ix]=g0_amp
                sigma0[iy,ix]=g0_sig
	    else:
	        velocity0[iy,ix]=-5.0
                amplitude0[iy,ix]=0.0
                sigma0[iy,ix]=0.0
	    if (g1_amp>0.15 and g1_sig<0.95):
	        velocity1[iy,ix]=g1_c
                amplitude1[iy,ix]=g1_amp
                sigma1[iy,ix]=g1_sig
	    else:
	        velocity1[iy,ix]=-5.0
                amplitude1[iy,ix]=0.0
                sigma1[iy,ix]=0.0
	    if (g2_amp>0.15 and g2_sig<0.95):
	        velocity2[iy,ix]=g2_c
                amplitude2[iy,ix]=g2_amp
                sigma2[iy,ix]=g2_sig
	    else:
	        velocity2[iy,ix]=-5.0
                amplitude2[iy,ix]=0.0
                sigma2[iy,ix]=0.0
	    
            #plt.plot(x,y)
            #plt.plot(x,out.best_fit,'r-')
            #plt.show()
	    #next=msvcrt.getch()
	    #print next
	    #if (next=='x'):
	    #    sys.exit()
        outstr=str(ix)+'    '+str(iy)+'    '+str(velocity0[iy,ix])+'    '+str(velocity1[iy,ix])+'    '+str(velocity2[iy,ix])+'    '+str(sigma0[iy,ix])+'    '+str(sigma1[iy,ix])+'    '+str(sigma2[iy,ix])+'    '+str(amplitude0[iy,ix])+'    '+str(amplitude1[iy,ix])+'    '+str(amplitude2[iy,ix])+'\n'
        outfile.write(outstr)

outfile.close()

np.save("velocity0.npy",velocity0)
np.save("velocity1.npy",velocity1)
np.save("velocity2.npy",velocity2)
np.save("sigma0.npy",sigma0)
np.save("sigma1.npy",sigma1)
np.save("sigma2.npy",sigma2)
np.save("amplitude0.npy",amplitude0)
np.save("amplitude1.npy",amplitude1)
np.save("amplitude2.npy",amplitude2)
