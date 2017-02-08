# reduceccd
Basic Image CCD reduction using CCDPROC package

*reduceccd* performs a basic automatic image CCD reduction on a set of files: master bias and master flat (for each filter) creation, bias and flat fielding correction of science images, cosmic ray removal, sky subtraction, aligning and combination of images of the same filter.

# Instalation instructions

The following scripts require these packages to be installed in order to run:

+ numpy
+ scipy
+ astropy
+ [ccdproc](https://github.com/astropy/ccdproc)
+ [photutils](https://github.com/astropy/photutils)

## Instructions

To perform basic data reductions, follow these steps:
