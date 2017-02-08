############################################################################
#                  Image reduction with CCDPROC package                    #
#                                                                          #
#                             VERSION: 0.0.3                               #
#                                                                          #
#                          Ruben Garcia-Benito                             #
#                                RGB@IAA                                   #
#                                                                          #
#                       Last Change: 2017/01/31                            #
############################################################################

from photutils import Background2D, SigmaClip, MedianBackground, MeanBackground
from photutils import SExtractorBackground, BiweightLocationBackground
from photutils import ModeEstimatorBackground, MMMBackground
from skimage.feature import register_translation
from astropy.stats import sigma_clipped_stats 
from scipy.ndimage import interpolation
from ccdproc import ImageFileCollection
from photutils import make_source_mask
import astropy.io.fits as pyfits
from astropy.table import Table
from astropy import units as u
from ccdproc import CCDData
from copy import deepcopy
#from pyraf import iraf 
import numpy as np
import ccdproc
import sys
import os
import re

################################ VERSION ###################################
VERSION = '0.0.3'                                                          #
############################################################################


# ------------------------------ FUNCTIONS -----------------------------
def join_path(name, directory=None):
    if directory is not None:
        name = os.path.join(directory, os.path.basename(name))
    return name

def fits2CCDData(lfits, key_unit='BUNIT', key_file='FILENAME', unit=None, single=False):
    lccd = []
    if not isinstance(lfits, (tuple, list)):
        lfits = [lfits]
    for fits in lfits:
        fits_unit = unit
        if os.path.exists(fits):
           hdr = pyfits.getheader(fits)
        else:
           print ('>>> WARNING: File "%s" NOT found' % os.path.basename(fits))
           continue
        if key_unit is not None and key_unit in hdr:
            try:
                fits_unit = eval('u.%s' % hdr[key_unit])
            except:
                pass
        if fits_unit is None:
            if key_unit is not None:
                sys.exit('>>> Units NOT found in header ("%s") of image "%s". Specify one in "unit" variable' % (key_unit, os.path.basename(fits)))
            else:
                print ('>>> WARNING: "key_unit" not specified')
        ccd = CCDData.read(fits, unit=fits_unit)
        if key_file is not None and not key_file in ccd.header:
            ccd.header[key_file] = os.path.basename(fits)
        lccd.append(ccd)
    if len(lccd) == 0:
        print ('>>> WARNING: NO files found!')
        return
    if single and len(lccd) == 1:
        lccd = lccd[0]
    return lccd

def getObjects(image_file_collection):
    table = image_file_collection.summary
    pattern = re.compile('flat|bias|bd\+|feige|land')
    objects = [item.split('-')[0].split('.')[-1] for item in table['file'] if pattern.search(item.lower()) is None]
    return np.unique(objects).tolist()

def getFilters(image_file_collection, key_filter='filter'):
    table = image_file_collection.summary
    filters = np.unique(table[key_filter].data)
    if isinstance(filters, np.ma.MaskedArray):
       filters = np.sort(filters[~filters.mask])
    return filters.tolist()

def getFilename(ccdfile, key='FILENAME'):
    namefile = None
    if isinstance(ccdfile, str):
        namefile = ccdfile
    elif isinstance(ccdfile, CCDData):
        if key in ccdfile.header:
            namefile = ccdfile.header[key]
    return namefile

def addKeyHdr(hdr, key, value):
    if value is not None and key is not None:
        hdr[key] = value
    return hdr

def ImageFileCollectionFilter(image_file_collection, ldfilter={}, dirname=None, abspath=True, 
	key_find='find', invert_find=False, return_mask=False, key='file', 
	dkeys=None, copy=True):
    dirname = image_file_collection.location if dirname is None else dirname
    if not abspath:
        dirname = None
    if not isinstance(ldfilter, (tuple, list)):
        ldfilter = [ldfilter]
    if dkeys is not None:
        ldfilter = addKeysListDict(ldfilter, dkeys, copy=True)
    mask  = None
    lmask = []
    for dfilter in ldfilter:
        if dfilter is None or len(dfilter) == 0:
            continue
        lfmask = []
        if key_find in dfilter:
            key_string = dfilter.pop(key_find)
            lfiles_find_mask = np.char.find(image_file_collection.summary[key].data, key_string) > -1
            if invert_find:
                lfiles_find_mask = np.invert(lfiles_find_mask)
            lfmask.append(lfiles_find_mask)
        lfiles = image_file_collection.files_filtered(**dfilter)
        lfiles_mask = np.in1d(image_file_collection.summary['file'], lfiles)
        lfmask.append(lfiles_mask)
        fmask = np.logical_and.reduce(lfmask)
        lmask.append(fmask)
    if len(lmask) > 0:
        mask = np.logical_or.reduce(lmask)
        list_files = image_file_collection.summary['file'][mask]
    else:
        list_files = image_file_collection.summary['file']
    list_files = [join_path(image, dirname) for image in list_files]
    if return_mask:
        return list_files, mask
    else:
        return list_files

def getListFiles(list_files, dfilter=None, mask=None, key_find='find', invert_find=False):
    if isinstance(list_files, ImageFileCollection):
        if mask is None:
                list_files = ImageFileCollectionFilter(list_files, dfilter, key_find=key_find, invert_find=invert_find)
        else:
                list_files = [join_path(image, list_files.location) for image in list_files.summary['file'][mask]]
    return list_files

def getTimeTable(list_files, key_time='EXPTIME', key_file='file', ext=0, abspath=False, mask=-1, clean=True, sort=True, **kwargs):
    dt = {key_file: [], key_time: []}
    for fits in list_files:
        fits = fits if abspath else os.path.basename(fits)
        dt[key_file].append(fits)
        hdr = pyfits.getheader(fits, ext=ext, **kwargs)
        if key_time in hdr:
            dt[key_time].append(hdr[key_time])
        else:
            dt[key_time].append(mask)
    table = Table(dt)
    if clean:
        id_time = table[key_time] == mask
        table = table[np.invert(id_time)]
    if sort:
        table.sort(key_time)
    return table

def addKeysListDict(ldict, dkeys, copy=True, force=False):
    if ldict is not None and dkeys is not None:
        if not isinstance(ldict, (list, tuple)):
             ldict = [ldict]
        if copy:
            ldict = deepcopy(ldict)
        for i in range(len(ldict)):
             if ldict[i] is not None and isinstance(ldict[i], dict):
                 for key in dkeys:
                     if key is not None:
                         if not key in ldict[i] or (key in ldict[i] and force):
                             ldict[i][key] = dkeys[key]
    return ldict

def ammendHeader(header):
    if '' in header:
        del header['']
    return header

def cleanCosmic(ccd, mbox=15, rbox=15, gbox=11, sigclip=5, cleantype="medmask", cosmic_method='lacosmic'):
    ctype = cosmic_method.lower().strip()
    ctypes = ['lacosmic', 'median']
    if not ctype in ctypes:
        print ('>>> Cosmic ray type "%s" NOT available [%s]' % (ctype, ' | '.join(ctypes)))
        return
    if ctype == 'lacosmic':
        ccd = ccdproc.cosmicray_lacosmic(ccd, sigclip=sigclip, cleantype=cleantype)
    elif ctype == 'median':
        ccd = ccdproc.cosmicray_median(ccd, mbox=mbox, rbox=rbox, gbox=gbox)
    if isinstance(ccd, CCDData):
        ccd.header['COSMIC'] = ctype.upper()
    return ccd

def create_master_bias(list_files, fitsfile=None, fits_section=None, gain=None, method='median', 
	dfilter={'imagetyp':'bias'}, mask=None, key_find='find', invert_find=False, sjoin=','):
    if gain is not None and not isinstance(gain, u.Quantity):
        gain = gain * u.electron / u.adu
    lbias = []
    list_files = getListFiles(list_files, dfilter, mask, key_find=key_find, invert_find=invert_find)
    for filename in list_files:
        ccd = CCDData.read(filename, unit= u.adu)
        trimmed = True if fits_section is not None else False
        ccd = ccdproc.trim_image(ccd, fits_section=fits_section, add_keyword={'trimmed': trimmed})
        if gain is not None:
            ccd = ccdproc.gain_correct(ccd, gain)
        lbias.append(ccd)
    combine = ccdproc.combine(lbias, method=method)
    if gain is not None and not 'GAIN' in combine.header:
        combine.header.set('GAIN', gain.value, gain.unit)
    combine.header['CGAIN'] = True if gain is not None else False
    combine.header['IMAGETYP'] = 'BIAS'
    combine.header['CCDVER'] = VERSION
    if sjoin is not None:
        combine.header['LBIAS'] = sjoin.join([os.path.basename(fits) for fits in list_files])
    combine.header['NBIAS'] = len(list_files)
    if fitsfile is not None:
        combine.header['FILENAME'] = os.path.basename(fitsfile)
        combine.write(fitsfile, clobber=True)
    return combine

def create_master_flat(list_files, flat_filter=None, fitsfile=None, bias=None, fits_section=None, gain=None, 
	method='median', key_filter='filter', dfilter={'imagetyp':'FLAT'}, mask=None, key_find='find', 
	invert_find=False, sjoin=','):
    if gain is not None and not isinstance(gain, u.Quantity):
        gain = gain * u.electron / u.adu
    lflat = []
    if dfilter is not None and key_filter is not None and flat_filter is not None:
        dfilter = addKeysListDict(dfilter, {key_filter: flat_filter})
    list_files = getListFiles(list_files, dfilter, mask, key_find=key_find, invert_find=invert_find)
    for filename in list_files:
        ccd = CCDData.read(filename, unit= u.adu)
        trimmed = True if fits_section is not None else False
        ccd = ccdproc.trim_image(ccd, fits_section=fits_section, add_keyword={'trimmed': trimmed})
        if gain is not None:
            ccd = ccdproc.gain_correct(ccd, gain)
        if bias is not None:
            if isinstance(bias, str):
                bias = fits2CCDData(bias, single=True)
            ccd = ccdproc.subtract_bias(ccd, bias)
        lflat.append(ccd)
    combine = ccdproc.combine(lflat, method=method)
    if gain is not None and not 'GAIN' in combine.header:
        combine.header.set('GAIN', gain.value, gain.unit)
    combine.header['CGAIN'] = True if gain is not None else False
    combine.header['IMAGETYP'] = 'FLAT'
    combine.header['CCDVER'] = VERSION
    addKeyHdr(combine.header, 'MBIAS', getFilename(bias))
    if sjoin is not None:
        combine.header['LFLAT'] = sjoin.join([os.path.basename(fits) for fits in list_files])
    combine.header['NFLAT'] = len(list_files)
    if fitsfile is not None:
        combine.header['FILENAME'] = os.path.basename(fitsfile)
        combine.write(fitsfile, clobber=True)
    return combine

def create_master_flat_from_dict(list_files, dflat, verbose=True, **kwargs):
    dflat_ccd = {}
    for flat_filter in dflat:
        master_flat = dflat[flat_filter]
        if verbose:
            print ('>>> Creating FLAT: %s' % os.path.basename(master_flat))
        dflat_ccd[flat_filter] = create_master_flat(list_files, flat_filter, master_flat, **kwargs)
    return dflat_ccd

def ccdproc_images_filter(list_files, image_filter=None, master_flat=None, master_bias=None, fits_section=None, gain=None, readnoise=None, 
	error=False, sky=True, dout=None, cosmic=False, mbox=15, rbox=15, gbox=11, cleantype="medmask", cosmic_method='lacosmic', 
	sigclip=5, key_filter='filter', dfilter={'imagetyp':'LIGHT'}, mask=None, key_find='find', invert_find=False, **kwargs):
    if gain is not None and not isinstance(gain, u.Quantity):
        gain = gain * u.electron / u.adu
    if readnoise is not None and not isinstance(readnoise, u.Quantity):
        readnoise = readnoise * u.electron 
    if dfilter is not None and key_filter is not None and image_filter is not None:
        dfilter = addKeysListDict(dfilter, {key_filter: image_filter})
    list_files = getListFiles(list_files, dfilter, mask, key_find=key_find, invert_find=invert_find)
    dccd = {}
    for filename in list_files:
        ccd = CCDData.read(filename, unit= u.adu)
        nccd = ccdproc.ccd_process(ccd, trim=fits_section, gain=gain, master_bias=master_bias, master_flat=master_flat, readnoise=readnoise, error=error)
        for key in ccd.header:
            if not key in nccd.header:
                nccd.header[key] = ccd.header[key]
        # Better get rid of the cosmic rays BEFORE subtracting the global sky background
        if cosmic:
            nccd = cleanCosmic(nccd, mbox=mbox, rbox=rbox, gbox=gbox, sigclip=sigclip, cleantype=cleantype, cosmic_method=cosmic_method)
        if sky:
            nccd = subtract_sky_ccd(nccd, **kwargs)
        addKeyHdr(nccd.header, 'MBIAS', getFilename(master_bias))
        addKeyHdr(nccd.header, 'MFLAT', getFilename(master_flat))
        filename = 'c%s' % os.path.basename(filename)
        dccd[filename] = nccd
        filename = join_path(filename, dout)
        nccd.header['FILENAME'] = os.path.basename(filename)
        nccd.header['CCDVER'] = VERSION
        nccd.header = ammendHeader(nccd.header)
        nccd.write(filename, clobber=True)
    return dccd

def ccdproc_images(list_files, dmaster_flat, master_bias=None, fits_section=None, gain=None, sky=True, dout=None, verbose=True, **kwargs):
    dccd_images = {}
    if gain is not None and not isinstance(gain, u.Quantity):
        gain = gain * u.electron / u.adu
    for filt in dmaster_flat:
        if verbose:
            print ('>>> Reducing Filter: %s' % filt)
        dccd = ccdproc_images_filter(list_files, filt, master_flat=dmaster_flat[filt], master_bias=master_bias, fits_section=fits_section, gain=gain, sky=sky, dout=dout, **kwargs)
        dccd_images[filt] = dccd
    return dccd_images

def subtract_sky_ccd_percentile(ccd, skysub='SKYSUB', skyval='SKYVAL', lower=1., upper=95., check_skysub=True, max_val=10000):
    if check_skysub and skysub is not None and skysub in ccd.header and ccd.header[skysub]:
        return ccd
    if not check_skysub and skysub is not None and skysub in ccd.header and ccd.header[skysub]:
        print ('WARNING: Image already sky subtracted! Subtracting AGAIN the sky! [set check_skysub = True for checking]')
    data = ccd.data[np.isfinite(ccd.data) & (ccd.data > 0.)].flatten()
    percentile = np.percentile(data, [lower, upper])
    sky = np.median(percentile)
    if max_val is not None and sky >= max_val:
        if 'FILENAME' in ccd.header:
            print ('WARNING: File "%s" has high sky level (%s >= %s). Correction NOT applied' % (ccd.header['FILENAME'], sky, max_val))
        else:
            print ('WARNING: File has high sky level (%s >= %s). Correction NOT applied' % (sky, max_val))
        return ccd
    ccd.data -= sky
    ccd.header[skysub] = True
    ccd.header[skyval] = sky
    ccd.header['SKYTYPE'] = '1D'
    return ccd

def subtract_sky_ccd_mask(ccd, skysub='SKYSUB', skyval='SKYVAL', snr=3, npixels=5, dilate_size=11, sigma=3, check_skysub=True, max_val=10000):
    if check_skysub and skysub is not None and skysub in ccd.header and ccd.header[skysub]:
        return ccd
    if not check_skysub and skysub is not None and skysub in ccd.header and ccd.header[skysub]:
        print ('WARNING: Image already sky subtracted! Subtracting AGAIN the sky! [set check_skysub = True for checking]')
    mask = make_source_mask(ccd.data, snr=sigma, npixels=npixels, dilate_size=dilate_size)
    mean, median, std = sigma_clipped_stats(ccd.data, sigma=sigma, mask=mask)
    sky = median
    if max_val is not None and sky >= max_val:
        if 'FILENAME' in ccd.header:
            print ('WARNING: File "%s" has high sky level (%s >= %s). Correction NOT applied' % (ccd.header['FILENAME'], sky, max_val))
        else:
            print ('WARNING: File has high sky level (%s >= %s). Correction NOT applied' % (sky, max_val))
        return ccd
    ccd.data -= sky
    ccd.header[skysub] = True
    ccd.header[skyval] = sky
    ccd.header['SKYTYPE'] = '1D'
    return ccd

def subtract_sky_ccd_2d(ccd, skysub='SKYSUB', skyval='SKYVAL', box_size=(50,50), filter_size=(3,3), method='sextractor', 
	check_skysub=True, max_val=10000, fitsky=None, sigma=3., iters=10, edge_method='crop'):
    # Although edge_method='pad' is recommended, it does give errors "total size of new array must be unchanged" for v0.3
    if check_skysub and skysub is not None and skysub in ccd.header and ccd.header[skysub]:
        return ccd
    if not check_skysub and skysub is not None and skysub in ccd.header and ccd.header[skysub]:
        print ('WARNING: Image already sky subtracted! Subtracting AGAIN the sky! [set check_skysub = True for checking]')

    dsky = {'median': MedianBackground(), 'mean': MeanBackground(), 'sextractor': SExtractorBackground(), 
	'biweight': BiweightLocationBackground(), 'mode': ModeEstimatorBackground(), 'mm': MMMBackground()}

    if not method in dsky:
        print ('WARNING: Background type "%s" NOT available!! Selecting "median" from [%s]' % (bkg_type, ' | '.join(dksy.keys())))
        bkg_estimator = dksy['median']
    else:
        bkg_estimator = dsky[method]

    sigma_clip = SigmaClip(sigma=sigma, iters=iters) if sigma is not None and iters is not None else None

    bkg = Background2D(ccd.data, box_size, filter_size=filter_size, sigma_clip=sigma_clip, bkg_estimator=bkg_estimator, edge_method=edge_method)

    if max_val is not None and bkg.background_median >= max_val:
        if 'FILENAME' in ccd.header:
            print ('WARNING: File "%s" has high sky level (%s >= %s). Correction NOT applied' % (ccd.header['FILENAME'], bkg.background_median, max_val))
        else:
            print ('WARNING: File has high sky level (%s >= %s). Correction NOT applied' % (bkg.background_median, max_val))
        return ccd
    ccd.data -= bkg.background
    ccd.header[skysub] = True
    ccd.header[skyval] = bkg.background_median
    ccd.header['SKYTYPE'] = '2D'
    ccd.header['BKGTYPE'] = method
    ccd.bkg = bkg
    if fitsky is not None:
        pyfits.writeto(fitsky, bkg.background, clobber=True)
    return ccd

def subtract_sky_ccd(ccd, method='2d', **kwargs):
    method = method.strip().lower()
    methods = ['percentile', 'mask', '2d']
    if method == 'percentile':
        return subtract_sky_ccd_percentile(ccd, **kwargs)
    elif method == 'mask':
        return subtract_sky_ccd_mask(ccd, **kwargs)
    elif method =='2d':
        return subtract_sky_ccd_2d(ccd, **kwargs)
    else:
        print ('>>> Method "%s" NOT found! Selecting "2d" from [%s]' % (method, ' | '.join(methods)))
        return subtract_sky_ccd_2d(ccd, **kwargs)

def subtract_sky_ccds(lccd, method='2d', **kwargs):
    if not isinstance(lccd, (tuple, list)):
        lcdd = [lccd]
    nlccd = []
    for ccd in lccd:
        nlccd.append(subtract_sky_ccd(ccd, method, **kwargs))
    return nlccd

def align_combine_images(list_files, fitsfile, ref_image_fits=None, precision=100, dout=None, hexp='EXPTIME', force=False, 
	minmax_clip=False, minmax_clip_min=0.0, sigma_clip=True, date=None, clip_extrema=False, func=np.ma.median, sigclip=5, 
	cosmic=False, mbox=15, rbox=15, gbox=11, cleantype="medmask", cosmic_method='lacosmic', sky=False, dict_sky={}, 
	dict_combine={}, suffix=None, hfilter='FILTER', key_file='file', hobj='OBJECT', ext=0, method='average', 
	align=True):
    if ref_image_fits is None:
        tobj = getTimeTable(list_files, key_time=hexp, key_file=key_file, ext=ext, abspath=True, mask=-1, sort=True, clean=True)
        ref_image_fits, list_files = tobj[key_file][-1], tobj[key_file][:-1].data.tolist()
    ref_image = CCDData.read(ref_image_fits)
    ref_image = ccdproc.cosmicray_lacosmic(ref_image, sigclip=sigclip)
    lccd = [ref_image]
    lexp = [ref_image.header[hexp]]
    all_images = [os.path.basename(ref_image_fits)]
    for img in list_files:
        image, hdr = pyfits.getdata(img, header=True)
        image_suffix = suffix if suffix is not None else img.split('.')[-1]
        if 'REFIMA' in hdr and 'IMAGES' in hdr and not force:
            continue
        if cosmic:
            image, _ = cleanCosmic(nccd, mbox=mbox, rbox=rbox, gbox=gbox, sigclip=sigclip, cleantype=cleantype, cosmic_method=cosmic_method)
        if align:
            shift, error, diffphase = register_translation(image, ref_image.data, precision)
            offset_image = interpolation.shift(image, (-shift[0], -shift[1]))
            img_shifted = '%s_shift.%s' % (img.split('.fit')[0], image_suffix)
            all_images.append(os.path.basename(img_shifted))
            img_shifted = join_path(img_shifted, dout)
            pyfits.writeto(img_shifted, offset_image, clobber=True)
            lccd.append(CCDData(offset_image, unit=ref_image.unit))
        else:
            all_images.append(os.path.basename(img))
            lccd.append(CCDData(image, unit=ref_image.unit))
        lexp.append(hdr[hexp])
        #iraf.unlearn("imshift")
        #iraf.imshift.interp_type = "spline3"
        #iraf.imshift(img, img, shift[1], shift[0])
    # Combine images
    lexp = np.array(lexp)
    scale_func = None # 1. / lexp
    combine = ccdproc.combine(lccd, method=method, scale=scale_func, minmax_clip=minmax_clip, func=func,
		minmax_clip_min=minmax_clip_min, sigma_clip=sigma_clip, clip_extrema=clip_extrema, **dict_combine)
    if sky:
        combine = subtract_sky_ccd(combine, **dict_sky)
    combine.header['IMAGES'] = str(' | '.join(all_images))
    combine.header['REFIMA'] = os.path.basename(ref_image_fits)
    combine.header['IMGSEXP'] = ' | '.join(map(str,lexp[1:].tolist() + [lexp[0]]))
    dir_out = os.path.dirname(ref_image_fits) if dout is None else dout
    fitsfile = join_path(fitsfile, dir_out)
    combine.header['FILENAME'] = os.path.basename(fitsfile)
    combine.header['CCDVER'] = VERSION
    combine.write(fitsfile, clobber=True)

def align_combine(image_file_collection, filters=None, objects=None, dout=None, suffix=None, date=None, key_find='file', 
	key_filter='filter', invert_find=False, dfilter={'imagetyp': 'LIGHT'}, verbose=True, force=False, find_obj=True,
	find_filter=True, **kwargs):
    if filters is None:
        filters = getFilters(image_file_collection)
    if objects is None:
        objects = getObjects(image_file_collection)
    if not isinstance(filters, (tuple, list)):
        filters = [filters]
    if not isinstance(objects, (tuple, list)):
        objects = [objects]
    table = image_file_collection.summary
    ddata = image_file_collection.location
    for obj in objects:
        dkeys = {key_find: obj} if find_obj else None
        mask_object = ImageFileCollectionFilter(image_file_collection, dfilter, dkeys=dkeys, invert_find=invert_find, return_mask=True)[-1]
        if not np.any(mask_object):
             print ('WARNING: Object "%s" does NOT exists!!' % obj)
             continue
        for image_filter in filters:
            dkeys = {} 
            if find_obj:
                dkeys[key_find] = obj
            if find_filter:
                dkeys[key_filter] = image_filter
            mask_table = ImageFileCollectionFilter(image_file_collection, dfilter, dkeys=dkeys, invert_find=invert_find, return_mask=True)[-1]
            if 'refima' in table.colnames and 'images' in table.colnames and not force:  
                mask = np.logical_and(table['refima'].mask, table['images'].mask)           
                mask_table = np.logical_and(mask_table, mask)
            tobj = table[mask_table]
            tobj.sort('exposure')
            if tobj['file'].data.size == 0:
                print ('WARNING: Object "%s" does NOT contain filter "%s"' % (obj, image_filter))
                continue
            if verbose:
                print ('>>> OBJECT: %s (Filter %s)' % (obj, image_filter))
            ref_image_fits = join_path(tobj['file'][-1], image_file_collection.location)
            list_files = [join_path(fits, image_file_collection.location) for fits in tobj['file'][:-1].data]
            image_suffix = suffix if suffix is not None else ref_image_fits.split('.')[-1]
            fitsfile = '%s_%s.%s' % (obj, image_filter, image_suffix) if date is None or len(date) == 0 else '%s_%s_%s.%s' % (obj, image_filter, date, image_suffix)
            align_combine_images(list_files, fitsfile, ref_image_fits=ref_image_fits, dout=dout, force=force, **kwargs)

def reduceNight(path, filters=None, fits_section=None, date=None, dout=None, create_bias=True, create_flat=True, 
	correct_images=True, gain=None, readnoise=None, sky=True, lower=1.0, upper=95., combine=True, align=True,
	objects=None, cosmic=True, mbox=15, rbox=15, gbox=11, cleantype="medmask", sky_after=True, dict_sky={}, 
	dict_combine={}, method='median', mask_bias=None, dfilter_bias={'imagetyp':'bias'}, lfits_bias=None, 
	invert_find_bias=False, dfilter_flat={'imagetyp':'FLAT'}, mask_flat=None, cosmic_method='lacosmic', 
	invert_find_flat=False, dfilter_images={'imagetyp':'LIGHT'}, mask_images=None, key_find='find', 
	key_filter='filter', find_obj=True, invert_find_images=False, find_filter=True, 
	invert_find_align=False, suffix=None, dict_align_combine={}, verbose=True):

    if date is None:
        date = path.split('/')[-2] if path.endswith('/') else path.split('/')[-1]
        if not date.isdigit():
            date = None
    if date is None:
        print ('>>> WARNING: You better provide a date!')

    if dout is not None and not os.path.exists(dout):
        dout = None
        print ('>>> WARNING: output directory does NOT exists! Saving to local working directory!')

    # Find all images
    ic_all = ImageFileCollection(path, keywords=None)

    # Find suffix
    if suffix is None:
        suffix = ic_all.summary['file'][0].split('.')[-1]
   
    
    if filters is None:
        filters = getFilters(ic_all)
        print ('>>> Filters found: %s' % ' | '.join(filters))

    if objects is None:
        objects = getObjects(ic_all)
        print ('>>> Objects found: %s' % ' | '.join(objects))

    # Sky option
    if sky:
        sky_before = not sky_after
    else:
        sky_before = False
        sky_after  = False

    # Create name master bias
    master_bias = 'master_bias_%s.%s' % (date, suffix) if date is not None else 'master_bias.%s' % suffix
    master_bias = join_path(master_bias, dout)

    # Create name master Flats
    dflat = {}
    for filt in filters:
        master_flat = 'master_flat_%s_%s.%s' % (filt, date, suffix) if date is not None else 'master_flat_%s.%s' % (filt, suffix)
        dflat[filt] = join_path(master_flat, dout)

    #--------------- Master bias file ---------------------
    ccd_master_bias = None
    if create_bias and master_bias is not None:
        lfits_bias = ic_all if lfits_bias is None else lfits_bias
        if verbose:
            print ('>>> Creating BIAS: %s' % os.path.basename(master_bias))
        ccd_master_bias = create_master_bias(lfits_bias, master_bias, fits_section=fits_section, gain=gain, method=method, dfilter=dfilter_bias, mask=mask_bias, key_find=key_find, invert_find=invert_find_bias)
    if not create_bias and master_bias is not None:
        ccd_master_bias = fits2CCDData(master_bias, single=True)

    # -------- Create Master File for each filter ---------
    dccd_master_flat = None
    if create_flat:
        dccd_master_flat = create_master_flat_from_dict(ic_all, dflat, bias=ccd_master_bias, fits_section=fits_section, gain=gain, method=method, dfilter=dfilter_flat, mask=mask_flat, key_find=key_find, invert_find=invert_find_flat, verbose=verbose)
    else:
        dccd_master_flat = {}
        for key in dflat:
            dccd_master_flat[key] = fits2CCDData(dflat[key], single=True)

    # ----------- Correct all science and standard stards images -------------
    if correct_images:
        ccdproc_images(ic_all, dccd_master_flat, master_bias=ccd_master_bias, fits_section=fits_section, gain=gain, dout=dout, sky=sky_before, 
		cosmic=cosmic, mbox=mbox, rbox=rbox, gbox=gbox, cleantype=cleantype, cosmic_method=cosmic_method, key_filter=key_filter,
        	dfilter=dfilter_images, mask=mask_images, key_find=key_find, invert_find=invert_find_images, verbose=verbose, **dict_sky)

    # ----------- Align and combine -------------------
    if combine:
        align_file_collection = image_file_collection if dout is None else ImageFileCollection(dout, keywords=None)
        align_combine(align_file_collection, filters, objects, dout=dout, sky=sky_after, dict_sky=dict_sky, dict_combine=dict_combine, key_find=key_find, invert_find=invert_find_align, suffix=suffix, dfilter=dfilter_images, find_obj=find_obj, find_filter=find_filter, align=align, verbose=verbose, **dict_align_combine)
# ----------------------------------------------------------------------
