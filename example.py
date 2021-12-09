import meaShear
import numpy as np
import astropy.io.fits as pyfits
import numpy.lib.recfunctions as rfn

# Read PSF image
psfData=pyfits.getdata('data/psf_test_1.fits')
# Read GAL image
print('Reading the simulated image ditorted by g_1= 0.02, g2=0.00')
galDatAll=pyfits.getdata('data/gal_test_1.fits')
imgList=[galDatAll[i*64:(i+1)*64,0:64] for i in range(4)]

ellRes=[]
for i in range(4):
    ellRes.append(meaShear.measureEllip(imgList[i],psfData))

ellRes =   rfn.stack_arrays(ellRes,usemask=False)

g1_est=np.average(ellRes['fpfs_e1'])/np.average(ellRes['fpfs_RE'])
g2_est=np.average(ellRes['fpfs_e2'])/np.average(ellRes['fpfs_RE'])
print('estimated shear is: g1= %.5f, g2= %.5f' %(g1_est,g2_est))
