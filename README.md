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

To install the package type in a terminal:

```
pip install https://github.com/rgbIAA/reduceccd/archive/master.zip
```

Or download the package and inside the directory type: 

```python
python setup.py install
```

## Instructions

If the images contain the basic information in the header, a standard reduction script would be something like:

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

# If you want to reduce all filters, leave "None"
# If you just want to reduce a few filters, provide a list or a string
filters = None # ['Ha6563', 'Ha6607', 'V', 'R', 'B']
# Trim section (None otherwise)
fits_section ='[:,:995]'

# CCD info
gain = 1.5 # +- 0.02 e / ADU
readnoise = 8.23 # +- 0.12 e rms

# If you want to reduce all objects, just write "None" 
# If you just want to reduce a particular set of objects, provide a list or string (one object)
objects = None # 'ngc7465'

rc.reduceNight(path, filters, fits_section=fits_section, dout=dout, gain=gain, create_flat=True,
        create_bias=True, sky=True, correct_images=True, objects=objects, align=True)
```

By default, `reduceNight` looks into the header of the images to sort the images by type in the following way:

+ Bias: keyword `imagetyp` = 'BIAS'
+ Flat: keyword `imagetyp` = 'FLAT'
+ Object: keyword `imagetyp` = 'LIGHT'

You can modify the way `reduceNight` looks for the different types of images. Let's say one night you have the keyword `imagetyp` = 'FLAT' in the header of some flats, and in others you have the string "flat" in the name of the file, but NOT in the header. You can give an OR search command in the form of a list of dictionaries:

```python
dfilter_flat = [{'imagetyp':'FLAT'}, {'find':'flat'}]
rc.reduceNight(path, filters, fits_section=fits_section, dout=dout, gain=gain, dfilter_flat=dfilter_flat)
```
`reduceNight` will take as flats files that *either* have the keyword `imagetyp` = 'BIAS' in the header OR files that contain the string "flat" in their names ('find' is the dictionary keyword to search for a string in a file name).

You can also provide AND instances for searching files just adding more keywords in the SAME dictionary:

```python
dfilter_flat = {'imagetyp':'FLAT', 'find':'flat'}
rc.reduceNight(path, filters, fits_section=fits_section, dout=dout, gain=gain, dfilter_flat=dfilter_flat)
```
In this case, `reduceNight` will consider flats files that have the keyword `imagetyp` = 'BIAS' in the header AND files that contain the string "flat" in their names.

There is also a `dfilter_bias` and `dfilter_images` for bias and science images, respectively. 
