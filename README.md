# reduceccd
Basic Image CCD reduction using CCDPROC package

*reduceccd* is a small pipeline that performs basic automatic image CCD reduction on a set of files: 

+ Create master bias
+ Create master flat for each filter
+ Performs bias and flat fielding correction of the science images
+ Cosmic ray removal
+ Sky subtraction
+ Align and combination of images of the same filter

# Installation instructions

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

If the images contain the basic information in the header, a standar reduction script would be something like:

```python
import reduceccd as rc
import warnings
import logging

logging.disable(logging.INFO)
warnings.simplefilter('ignore')

# Your raw data
path = '/User/me/route/to/my/raw/data/'
# Output directory
dout = 'redudir/'

# If you want to reduce all filters, just write "None" 
# If you want just to reduce a few filters, provide a list or string (one filter)
filters = None # ['Ha6563', 'Ha6607', 'V', 'R', 'B']
# Trim section (None otherwise)
fits_section ='[:,:995]'

# CCD info
gain = 1.5 # +- 0.02 e / ADU
readnoise = 8.23 # +- 0.12 e rms

# If you want to reduce all objects, just write "None" 
# If you want just to reduce a particular object(s), provide a list or string (one object)
objects = None # 'ngc7465'

rc.reduceNight(path, filters, fits_section=fits_section, dout=dout, gain=gain, create_flat=True,
        create_bias=True, sky=True, correct_images=True, objects=objects, align=True)
```
