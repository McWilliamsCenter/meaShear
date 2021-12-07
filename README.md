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

## Things to do
As you can see, the code suffers many kinds of errors
+ TypeError
+ Parameter ranges
+ ...
