# reduceccd
Basic Image CCD reduction using CCDPROC package

*reduceccd* is a small pipeline that performs basic automatic image CCD reduction on a set of files: 

+ Create master bias
+ Create master flat for each filter
+ Performs bias and flat fielding correction of the science images
+ Cosmic ray removal
+ Sky subtraction
+ Aligning and combination of images of the same filter

# Instalation instructions

The following scripts require these packages to be installed in order to run:

+ numpy
+ scipy
+ astropy
+ [ccdproc](https://github.com/astropy/ccdproc)
+ [photutils](https://github.com/astropy/photutils)

Download the package and inside the directory type: 
```python
python setup.py install
```

## Instructions

To perform basic data reductions, follow these steps:
