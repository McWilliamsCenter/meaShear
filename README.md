# Unit Testing

In case you have not prepared your code to implement unit test, please use the
following example code for shear estimation.

## Download and Install

```shell
git clone https://github.com/mr-superonion/meaShear.git
cd meaShear
pip install .
```

## About the code

This code implements a simple estimator of the weak gravitational lensing shear -- i.e., the distortion of light from distant source galaxies as they pass through the large-scale structure of the Universe on their way to us.  The deflection of light rays by the matter distribution causes changes in galaxy brightness, size, and shape.  Coherent changes in shape can be used to measure those deflections and learn about the distribution of dark matter in the Universe.

The code in this repository implements one step in that process: measurement of galaxy shapes of images, after accounting for blurring of the images due to the light passing through the atmosphere and telescope optics.

### shear estimation
<p align="center">
<img src="fig/shear_distort.png" alt="shear" width="800">
</p>

```python

import meaShear
import numpy as np
import astropy.io.fits as pyfits
import numpy.lib.recfunctions as rfn
# Read PSF image
psfData=pyfits.getdata('data/psf_test.fits')
# Read GAL image
print('Reading the simulated image distorted by g_1= 0.02, g2=0.00')
galDatAll=pyfits.getdata('data/gal_test.fits')
imgList=[galDatAll[i*64:(i+1)*64,0:64] for i in range(4)]

ellRes=[]
for i in range(4):
    # input image and PSF to estimate e and R
    eR=meaShear.measureEllip(imgList[i],psfData)
    ellRes.append()

ellRes =   rfn.stack_arrays(ellRes,usemask=False)

# average over response to estimate shear
g1_est=np.average(ellRes['fpfs_e1'])/np.average(ellRes['fpfs_RE'])
g2_est=np.average(ellRes['fpfs_e2'])/np.average(ellRes['fpfs_RE'])
print('estimated shear is: g1= %.5f, g2= %.5f' %(g1_est,g2_est))
```

### ring test
Ring thest is widely applied to check the accuracy of a shear estimator.


## Test Development List

+ Code performance (accuracy)
    -   Use ring test to test the accuracy of the estimator on noiseless
        galaxies
+ Inputs Control (interact correctly with users)
    -   Test to make sure the code requires the inputs galaxies and PSF images
        are ndarrays
    -   Test to make sure that if the input image has values that are
        problematic, like NaN or Inf, the code does something reasonable
    -   Test to make sure the tuning parameters in reasonable ranges
        (0<scale_par<1), (weight_par>0)

### Example
You can find a testing example [here](./tests/test_accuracy.py), and you can run the test by

```shell
cd tests
nosetests -v --with-coverage --cover-package=meaShear
```
