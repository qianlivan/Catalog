import pyfits
from pylab import *

import math
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import mpl_toolkits
from matplotlib.patches import Ellipse

#vfield=np.load('velocity0_0_499.npy')
vfield=np.load('velocity0_500_999.npy')

vfield[vfield<0]=0
vfield[vfield>12.0]=0

ax = plt.subplot(111)
im = plt.imshow(vfield, cmap=cm.gist_heat
                ,origin='lower', aspect='equal'
                ,interpolation='none')
xlabel('RA')
ylabel('Dec')
plt.colorbar(im,orientation='vertical')
savefig('velocityfield.eps')
plt.show()
