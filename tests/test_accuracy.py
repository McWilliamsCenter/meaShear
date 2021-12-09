import meaShear
import numpy as np
import astropy.io.fits as pyfits
import numpy.lib.recfunctions as rfn

def run_noiseless_gal(itest):
    # Read PSF image
    psfData=pyfits.getdata('../data/psf_test_%d.fits' %(itest))
    # Read GAL image
    galDatAll=pyfits.getdata('../data/gal_test_%d.fits' %(itest))
    imgList=[galDatAll[i*64:(i+1)*64,0:64] for i in range(4)]


    ellRes=[]
    for i in range(4):
        ellRes.append(meaShear.measureEllip(imgList[i],psfData))

    ellRes =   rfn.stack_arrays(ellRes,usemask=False)

    g1_est=np.average(ellRes['fpfs_e1'])/np.average(ellRes['fpfs_RE'])
    g2_est=np.average(ellRes['fpfs_e2'])/np.average(ellRes['fpfs_RE'])

    np.testing.assert_almost_equal(g1_est, 0.02, 5)
    # 1e-5 accuracy
    np.testing.assert_almost_equal(g2_est, 0.00, 5)
    # 1e-5 accuracy
    return

def test_noiseless_gals():
    run_noiseless_gal(itest=1)
    run_noiseless_gal(itest=2)
    return

