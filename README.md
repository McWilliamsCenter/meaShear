# Unit Testing

If you have not prepared your code to implement unit test, please use the
following example.

## Download and Install

```shell
git clone https://github.com/mr-superonion/meaShear.git
cd meaShear
pip install .
```

## Example
You can find a testing example [here](./tests/test_accuracy.py), and you can run the test by

```shell
cd tests
nosetests -v --with-coverage --cover-package=meaShear
```

## Test Development List
+ Inputs Control (interact correctly with users)
    -   Test to make sure the code requires the inputs galaxies and PSF images
        are ndarrays
    -   Test to make sure that if the input image has values that are
        problematic, like NaN or Inf, the code does something reasonable
    -   Test to make sure the tuning parameters in reasonable ranges
        ($0<$scale_par$<1$), (weight_par$>0$)
+ Code performance (accuracy)
    -   Use ring test to test the accuracy of the estimator on noiseless
        galaxies
