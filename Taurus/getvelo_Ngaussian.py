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

ngaussian=6


hdulist = pyfits.open('../t13_new.fits')
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

createVar=locals()

for i in range(ngaussian):
    createVar['velocity'+str(i)] = np.zeros([ny,nx])    
    createVar['amplitude'+str(i)] = np.zeros([ny,nx])    
    createVar['sigma'+str(i)] = np.zeros([ny,nx])    

residual=np.zeros([nz,ny,nx])


nx1str=sys.argv[1]
nx2str=sys.argv[2]
nx1=int(nx1str)
nx2=int(nx2str)
nx2=min([nx,nx2])


for ix in range(nx1,nx2+1):
    print ix
    for i in range(ngaussian):
        createVar['g'+str(i)+'_c'] = 0.0
        createVar['g'+str(i)+'_cmin'] = -3.0
        createVar['g'+str(i)+'_cmax'] = 13.0
        createVar['g'+str(i)+'_sig'] = 0.5
        createVar['g'+str(i)+'_sigmin'] = 0.1
        createVar['g'+str(i)+'_sigmax'] = 1.0
        createVar['g'+str(i)+'_amp'] = 1.0
        createVar['g'+str(i)+'_ampmin'] = 0.1
    for iy in range(ny):
	print ix,iy
        y = image[:,iy,ix]
	noise=np.std(y[0:10])
	peak=np.max(y)
	index_y=y>(peak-0.1)
	tempy=y[index_y]
	tempx=x[index_y]

        for i in range(ngaussian):
            locals()['g'+str(i)+'_amp'] = tempy[0]+np.random.rand() 
            locals()['g'+str(i)+'_c'] = tempx[0]+np.random.rand()
	if(peak<5*noise):
            for i in range(ngaussian):
                locals()['velocity'+str(i)][iy,ix] = -5.0
                locals()['amplitude'+str(i)] = 0.0
                locals()['sigma'+str(i)] = 0.0
        else:
            gauss0 = GaussianModel(prefix='g0_')
            i=0
            tempg_c=locals()['g'+str(i)+'_c']
            tempg_cmin=locals()['g'+str(i)+'_cmin']
            tempg_cmax=locals()['g'+str(i)+'_cmax']
            tempg_sig=locals()['g'+str(i)+'_sig']
            tempg_sigmin=locals()['g'+str(i)+'_sigmin']
            tempg_sigmax=locals()['g'+str(i)+'_sigmax']
            tempg_amp=locals()['g'+str(i)+'_amp']
            tempg_ampmin=locals()['g'+str(i)+'_ampmin']
            pars = gauss0.make_params()
            pars['g0_center'].set(tempg_c, min=tempg_cmin, max=tempg_cmax)
            pars['g0_sigma'].set(tempg_sig, min=tempg_sigmin, max=tempg_sigmax)
            pars['g0_amplitude'].set(tempg_amp, min=tempg_ampmin)
            mod=gauss0
            for i in range(1,ngaussian):
                prestr='g'+str(i)+'_'
                centerstr='g'+str(i)+'_center'
                sigmastr='g'+str(i)+'_sigma'
                amplitudestr='g'+str(i)+'_amplitude'
                
                tempg_c=locals()['g'+str(i)+'_c']
                tempg_cmin=locals()['g'+str(i)+'_cmin']
                tempg_cmax=locals()['g'+str(i)+'_cmax']
                tempg_sig=locals()['g'+str(i)+'_sig']
                tempg_sigmin=locals()['g'+str(i)+'_sigmin']
                tempg_sigmax=locals()['g'+str(i)+'_sigmax']
                tempg_amp=locals()['g'+str(i)+'_amp']
                tempg_ampmin=locals()['g'+str(i)+'_ampmin']

                createVar['gauss'+str(i)] = GaussianModel(prefix=prestr)
                pars.update(locals()['gauss'+str(i)].make_params())
                pars['g'+str(i)+'_center'].set(tempg_c, min=tempg_cmin, max=tempg_cmax)
                pars['g'+str(i)+'_sigma'].set(tempg_sig, min=tempg_sigmin, max=tempg_sigmax)
                pars['g'+str(i)+'_amplitude'].set(tempg_amp, min=tempg_ampmin)
                mod = mod+locals()['gauss'+str(i)]

            
            init = mod.eval(pars,x=x)
  
            out = mod.fit(y,pars,x=x)
            residual[:,iy,ix]=out.residual
	    temp00=out.params['g0_amplitude'].value
	    temp01=out.params['g0_center'].value
	    temp02=out.params['g0_sigma'].value
	    temp10=out.params['g1_amplitude'].value
	    temp11=out.params['g1_center'].value
	    temp12=out.params['g1_sigma'].value
	    temp20=out.params['g2_amplitude'].value
	    temp21=out.params['g2_center'].value
	    temp22=out.params['g2_sigma'].value
	    temp30=out.params['g3_amplitude'].value
	    temp31=out.params['g3_center'].value
	    temp32=out.params['g3_sigma'].value
	    temp40=out.params['g4_amplitude'].value
	    temp41=out.params['g4_center'].value
	    temp42=out.params['g4_sigma'].value
	    temp50=out.params['g5_amplitude'].value
	    temp51=out.params['g5_center'].value
	    temp52=out.params['g5_sigma'].value
	    temp=[[temp00,temp01,temp02],[temp10,temp11,temp12],[temp20,temp21,temp22],[temp30,temp31,temp32],[temp40,temp41,temp42],[temp50,temp51,temp52]]
	    temp.sort(reverse=True)
	    print temp
	    print peak
            for i in range(ngaussian):
                locals()['g'+str(i)+'_amp']=temp[i][0]
                locals()['g'+str(i)+'_c']=temp[i][1]
                locals()['g'+str(i)+'_sig']=temp[i][2]
	        if (locals()['g'+str(i)+'_amp']>0.15 and locals()['g'+str(i)+'_sig']<0.95):
	            locals()['velocity'+str(i)][iy,ix]=locals()['g'+str(i)+'_c']
                    locals()['amplitude'+str(i)][iy,ix]=locals()['g'+str(i)+'_amp']
                    locals()['sigma'+str(i)][iy,ix]=locals()['g'+str(i)+'_sig']
	        else:
	            locals()['velocity'+str(i)][iy,ix]=-5.0
                    locals()['amplitude'+str(i)][iy,ix]=0.0
                    locals()['sigma'+str(i)][iy,ix]=0.0

filename='velocity0_'+nx1str+'_'+nx2str+'.npy'
np.save(filename,velocity0)
filename='velocity1_'+nx1str+'_'+nx2str+'.npy'
np.save(filename,velocity1)
filename='velocity2_'+nx1str+'_'+nx2str+'.npy'
np.save(filename,velocity2)
filename='velocity3_'+nx1str+'_'+nx2str+'.npy'
np.save(filename,velocity3)
filename='velocity4_'+nx1str+'_'+nx2str+'.npy'
np.save(filename,velocity4)
filename='velocity5_'+nx1str+'_'+nx2str+'.npy'
np.save(filename,velocity5)
filename='sigma0_'+nx1str+'_'+nx2str+'.npy'
np.save(filename,sigma0)
filename='sigma1_'+nx1str+'_'+nx2str+'.npy'
np.save(filename,sigma1)
filename='sigma2_'+nx1str+'_'+nx2str+'.npy'
np.save(filename,sigma2)
filename='sigma3_'+nx1str+'_'+nx2str+'.npy'
np.save(filename,sigma3)
filename='sigma4_'+nx1str+'_'+nx2str+'.npy'
np.save(filename,sigma4)
filename='sigma5_'+nx1str+'_'+nx2str+'.npy'
np.save(filename,sigma5)
filename='amplitude0_'+nx1str+'_'+nx2str+'.npy'
np.save(filename,amplitude0)
filename='amplitude1_'+nx1str+'_'+nx2str+'.npy'
np.save(filename,amplitude1)
filename='amplitude2_'+nx1str+'_'+nx2str+'.npy'
np.save(filename,amplitude2)
filename='amplitude3_'+nx1str+'_'+nx2str+'.npy'
np.save(filename,amplitude3)
filename='amplitude4_'+nx1str+'_'+nx2str+'.npy'
np.save(filename,amplitude4)
filename='amplitude5_'+nx1str+'_'+nx2str+'.npy'
np.save(filename,amplitude5)

filename='residual_'+nx1str+'_'+nx2str+'.npy'
np.save(filename,residual)
