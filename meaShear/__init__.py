from .base import *
from . import imgutil

__version__="0.0.1"

def measureEllip(galData,psfData,scale_par=0.85,weight_par=100.):
    """
    Measure the FPFS ellipticity and its shear response
    Parameters:
        galData:    2D array of galaxy image    [ndarray]
        psfData:    2D array of PSF image       [ndarray]
        scale_par:  The typical measurement -- 1/scale_par*seeing
        weight_par: The weighting parameter changes the relative weights between galaxies

    Returns :
        galaxy ellipticity (e) and shear response (R)
    """
    # type control
    if type(galData) is not np.ndarray:
        raise TypeError('galData must be a numpy ndarray')
    if type(psfData) is not np.ndarray:
        raise TypeError('psfData must be a numpy ndarray')
    if type(scale_par) not in [int,float]:
        raise TypeError('scale_par must be a real nmber >0 and <1')

    # value control
    if scale_par>=1 or scale_par<=0:
        raise ValueError('scale_par must be a real nmber >0 and <1')

    moments =   measureMoments(galData,psfData,scale_par)
    ellRes  =   fpfsM2E(moments,weight_par)
    return ellRes
