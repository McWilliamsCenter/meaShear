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
    assert_raises(TypeError,meaShear.measureEllip,'hello',psfData,scale_par=0.7,weight_par=100)
    assert_raises(TypeError,meaShear.measureEllip,galData,'world',0.7,100)
    #TODO: Add tests for ensuring weight_par is real value
    return

def test_input_values():
    # Read PSF image
    psfData=pyfits.getdata('../data/psf_test_1.fits')
    # Read GAL image
    galData=pyfits.getdata('../data/gal_test_1.fits')[0:64,0:64]
    assert_raises(ValueError,meaShear.measureEllip,galData,psfData,-0.7,100)
    assert_raises(ValueError,meaShear.measureEllip,galData,psfData,7,100)
    #TODO: Add tests for preventing nan, inf values in galData and psfData
    return

