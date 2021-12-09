from setuptools import setup, find_packages
import numpy

setup(
    name='meaShear',
    version='0.0.1',
    description='shear estimation',
    author='Xiangchong Li et al.',
    author_email='mr.superonion@hotmail.com',
    python_requires='>=3.6',
    install_requires=[
        'numpy',
        'astropy',
    ],
    include_dirs=numpy.get_include(),
    packages=find_packages(),
    include_package_data=True,
    zip_safe=False,
)
