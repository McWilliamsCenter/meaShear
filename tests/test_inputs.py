import meaShear
import numpy as np
import astropy.io.fits as pyfits
import numpy.lib.recfunctions as rfn
from nose.tools import assert_raises

def test_input_types():
    # Read PSF image
    psfData=pyfits.getdata('../data/psf_test_1.fits')
    # Read GAL image
    galData=pyfits.getdata('../data/gal_test_1.fits')[0:64,0:64]
    assert_raises(TypeError,meaShear.measureEllip,'hello',psfData,0.7,100)
    assert_raises(TypeError,meaShear.measureEllip,galData,'world',0.7,100)
    assert_raises(TypeError,meaShear.measureEllip,galData,psfData,True,100)
    return

