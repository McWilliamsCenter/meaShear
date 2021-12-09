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
### shear estimation
<p align="center">
<img src="fig/shear_distort.png" alt="shear" width="800">
</p>

### ring test
+ Simulate a group galaxies with intrinsic ellipticity average to zero and
distort them by known shear to test the shear estimation.

```python

import meaShear
import numpy as np
import astropy.io.fits as pyfits
import numpy.lib.recfunctions as rfn
# Read PSF image
psfData=pyfits.getdata('data/psf_test.fits')
# Read GAL image
print('Reading the simulated image ditorted by g_1= 0.02, g2=0.00')
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


## Test Development List

+ Code performance (accuracy)
    -   Use ring test to test the accuracy of the estimator on noiseless
        galaxies
+ Inputs Control (interact correctly with users)
    -   type control:  Inputs galaxies and PSF images are ndarrays
    -   type control:  The tuning parameters are real
    -   value control: If the input image has values that are
        problematic, like NaN or Inf, the code does something reasonable
    -   value control: the tuning parameters in reasonable ranges
        (0<scale_par<1), (weight_par>0)

### Example
You can find a testing example [here](./tests/test_accuracy.py), and you can run the test by

```shell
cd tests
nosetests -w tests -v --with-coverage --cover-package=meaShear
```

### Useful functions


```python
np.testing.assert_almost_equal(out_come, expectation, 5) # 1e-5 accuracy
assert_raises(TypeError,your_function,input_1,input_2,input_3,input_4)
assert_raises(ValueError,your_function,input_1,input_2,input_3,input_4)
```
