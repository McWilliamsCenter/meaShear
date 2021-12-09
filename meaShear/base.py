# Simplified version of FPFS shear estimator
# Copyright 20210905 Xiangchong Li.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#
# python lib

from . import imgutil
import numpy as np
import numpy.lib.recfunctions as rfn

def deconvolvePow(data,psfPow,indK,order=1.):
    """
    Deconvolve the galaxy power with the PSF power
    Parameters:
        data :      galaxy power, centerred at middle
        psfPow :    PSF power, centerred at middle
        indK:       the k-ranges to use
        order:      set to 1 for deconvolving PSF; 2 for
                    estimating covariance of shapelet mdoes
    Returns :
        Deconvolved galaxy power (truncated at rlim)

    """

    out =   np.zeros(data.shape,dtype=np.float64)
    out[indK]=data[indK]/psfPow[indK]**order
    return out

def get_kscale(powImg,sigma):
    """
    Get rlim, the pixles outside rlim are sufficiently supressed by the shaplet
    Gaussian kernel, so that we truncate at this radius in Fourier space and
    set the outside pixels to 0
    Parameters:
        powImg :    power function in Fourier space
        sigma  :    scale radius of the Gaussian kernel of shapelet basis vectors
    Returns:
        indK:       the k-ranges to use, outside this ranges are set to 0
    """
    ngrid   =   powImg.shape[0]
    thres   =   1.e-3
    for dist in range(ngrid//5,ngrid//2-1):
        ave =  abs(np.exp(-dist**2./2./sigma**2.)/powImg[ngrid//2+dist,ngrid//2])
        ave +=  abs(np.exp(-dist**2./2./sigma**2.)/powImg[ngrid//2,ngrid//2+dist])
        ave =   ave/2.
        if ave<=thres:
            rlim=   dist
            break
    indX=np.arange(ngrid//2-rlim,ngrid//2+rlim+1)
    indY=indX[:,None]
    ind2D=np.ix_(indX,indX)
    return ind2D

def shapelets_transform(data,chi):
    """
    Project image onto shapelet basis vectors.
    This is equivalent to fitting with a uniform covariance. By doing so, we
    loss some precision but the shapelet basis vectors are kept orthogonal and
    we avoid entangling other shaelet modes in the fitting. This avoids model
    bias.
    process.
    Parameters:
        data:   image to transfer

    Returns:
        projection in shapelet space

    """

    # Only uses M00, M22 (real and img) and M40
    _indC  =   np.array([0,12,20])
    # Moments
    M       =   np.sum(data[None,:,:]*chi[_indC,:,:],axis=(1,2))
    types   =   [('fpfs_M00','>f8'),\
                ('fpfs_M22c','>f8'),('fpfs_M22s','>f8'),\
                ('fpfs_M40','>f8')\
                ]
    out     =   np.array((M.real[0],\
                M.real[1],M.imag[1],\
                M.real[2]),dtype=types)
    return out

def fpfsM2E(moments,const=1.):
    """
    Estimate FPFS ellipticities from fpfs moments

    Parameters:
        moments:    input FPFS moments     [float array]
        const:      the weighting Constant [float]

    Returns:
        an array of (FPFS ellipticities, FPFS ellipticity response, FPFS flux
        ratio, and FPFS selection response)

    """
    #Get weight
    weight  =   moments['fpfs_M00']+const
    #Ellipticity
    e1      =   moments['fpfs_M22c']/weight
    e2      =   moments['fpfs_M22s']/weight
    e1sq    =   e1*e1
    e2sq    =   e2*e2
    #FPFS flux ratio
    s0      =   moments['fpfs_M00']/weight
    s4      =   moments['fpfs_M40']/weight
    #FPFS sel Respose (part1)
    e1sqS0  =   e1sq*s0
    e2sqS0  =   e2sq*s0

    eSq     =   e1sq+e2sq
    eSqS0   =   e1sqS0+e2sqS0
    #Response factor
    RE      =   1./np.sqrt(2.)*(s0-s4+e1sq+e2sq)
    types   =   [('fpfs_e1','>f8'),('fpfs_e2','>f8'),('fpfs_RE','>f8'),\
                ('fpfs_s0','>f8'), ('fpfs_eSquare','>f8'), ('fpfs_RS','>f8')]
    ellDat  =   np.array(np.zeros(moments.size),dtype=types)
    ellDat['fpfs_e1']   =   e1
    ellDat['fpfs_e2']   =   e2
    # In Li et. al (2018) there is a minus sign difference
    ellDat['fpfs_RE']   =   -1.*RE
    ellDat['fpfs_s0']   =   s0
    ellDat['fpfs_eSquare']  =   eSq
    ellDat['fpfs_RS']   =   (eSq-eSqS0)/np.sqrt(2.)
    return ellDat

def measureMoments(galData,psfData,beta=0.85):
    """
    Measure the shapelets modes (moments) from galaxy's Fourier power after
    deconvolving PSF.
    Parameters:
        galData:    2D array of galaxy image
        psfData:    2D array of PSF image
        beta   :    scale parameter

    Returns :
        moments after PSF deconvolution
    """
    ngrid   =   psfData.shape[0]
    psfPow  =   imgutil.getFouPow(psfData)
    # PSF scale
    sigmaPsf    =   imgutil.getRnaive(psfPow)
    # shapelet scale
    sigma   =   max(min(sigmaPsf*beta,4.),1.)
    # get the k-scale in Fourier space to truncate
    ind2D   =   get_kscale(psfPow,sigma)
    galPow  =   imgutil.getFouPow(galData)

    decPow  =   deconvolvePow(galPow,psfPow,ind2D,order=1.)
    # Preparing shapelets (reshaped)
    chi     =   imgutil.shapelets2D(ngrid,4,sigma).reshape((25,ngrid,ngrid))
    mm      =   shapelets_transform(decPow,chi)
    return mm

