# This note is for JACKS data reduction


## 0. Download the uncalibrated data
One suggested structure of the directory for your working space.
- _your working space_
  - _raw_ (put raw.ms here)
  - _calibrated_ (put pip_script.py here)
    - _cont_ (split the continuum spw here and do the imaging of continuum)
    - _water_ (split the water maser spw here and generate the selfcal table)
    - _nh3_ (put the mask file here)
      - _11_
      - _22_
      - _33_
      - _44_
      - _55_
The default running environment is casa 6.5.4 unless it is noted with **python**. 

#### NOTE! The IDs of spw are different between 25A-303 and 22B-193
For 25A-303:
   - 7 -- NH3 1,1
   - 6 -- NH3 2,2
   - 5  -- NH3 3,3
   - 4  -- NH3 4,4
   - 3  -- NH3 5,5
   - 2  -- H64a
   - 21 -- H2O
   - 20 -- H66a
   - 19 -- NH3 5,4 5( 4)0a- 5( 4)0s
   - 18 -- NH3 4,3 4( 3)0a- 4(-3)0s
   - 17 -- NH3 3,2 3( 2)0a- 3( 2)0s
   - 16 -- NH3 7,6 7( 6)0a- 7(-6)0s
   - 15 -- NH3 2,1 2( 1)0a- 2( 1)0s
   - 8\~14,22\~28 -- continuum

For 22B-193:
   - 11 -- NH3 1,1
   - 10 -- NH3 2,2
   - 9  -- NH3 3,3
   - 8  -- NH3 4,4
   - 7  -- NH3 5,5
   - 6  -- H64a
   - 23 -- H2O
   - 22 -- H66a
   - 21 -- NH3 5,4 5( 4)0a- 5( 4)0s
   - 20 -- NH3 4,3 4( 3)0a- 4(-3)0s
   - 19 -- NH3 3,2 3( 2)0a- 3( 2)0s
   - 18 -- NH3 7,6 7( 6)0a- 7(-6)0s
   - 17 -- NH3 2,1 2( 1)0a- 2( 1)0s
   - 2\~5,12\~16,24\~30 -- continuum #NOTE! There are 2 spws missing in some observations. Please check the listobs.txt carefully!


## 1. Rerun the pipeline and check the weblog

Run the script in `calibrated`. Remember to revise **myvis** to your ms file.


```python
#example for the pipeline
context = h_init()
context.set_state('ProjectSummary', 'proposal_code', 'VLA/22B-193')
context.set_state('ProjectSummary', 'proposal_title', 'JACKS')
context.set_state('ProjectSummary', 'piname', 'EACMills')
context.set_state('ProjectSummary', 'observatory', 'Karl G. Jansky Very Large Array')
context.set_state('ProjectSummary', 'telescope', 'EVLA')
context.set_state('ProjectStructure', 'ppr_file', 'PPR.xml')
context.set_state('ProjectStructure', 'recipe_name', 'hifv_cal')
# Define the raw visibility data, which should be retrieved (and unpacked) from the NRAO archive
myvis = '../raw/25A-303.sb48641652.eb48925296.60819.25704964121.ms'
try:
    hifv_importdata(vis=myvis, session=['session_1'])
    #hifv_hanning(pipelinemode="automatic") # Skip Hanning smoothing
    hifv_flagdata(hm_tbuff='1.5int', fracspw=0.01, intents='*POINTING*,*FOCUS*,*ATMOSPHERE*,*SIDEBAND_RATIO*, *UNKNOWN*, *SYSTEM_CONFIGURATION*, *UNSPECIFIED#UNSPECIFIED*')
    hifv_vlasetjy(pipelinemode="automatic")
    hifv_priorcals(pipelinemode="automatic")
    hifv_syspower(pipelinemode="automatic")
    hifv_testBPdcals(pipelinemode="automatic")
    #hifv_checkflag(checkflagmode='bpd-vla') # Skip RFI flagging for the bandpass calibrator
    hifv_semiFinalBPdcals(pipelinemode="automatic")
    #hifv_checkflag(checkflagmode='allcals-vla') # Skip RFI flagging for the phase calibrator
    hifv_solint(pipelinemode="automatic")
    hifv_fluxboot(pipelinemode="automatic")
    hifv_finalcals(pipelinemode="automatic")
    hifv_applycals(pipelinemode="automatic")
    #hifv_checkflag(checkflagmode='target-vla') # Skip RFI flagging for the targets
    #hifv_statwt(pipelinemode="automatic") # Skip adjusting statistical weights
    hifv_plotsummary(pipelinemode="automatic")
    hif_makeimlist(intent='PHASE,BANDPASS', specmode='cont')
    hif_makeimages(hm_masking='centralregion')
finally:
    h_save()

```

Then check the weblog and flag the abnormal scans using **plotms** and **flagdata**.

## 2. Split the water maser spw and generate the selfcal table

Split the water maser spw in `water`. 


```python
#generate the listobs
myviscal='../25A-303.sb48641652.eb48925296.60819.25704964121.ms'
listobs(myviscal, listfile=myviscal+'.listobs.txt', overwrite=True)

#check the fields
tb.open(myviscal+'/FIELD')
fieldlist = tb.getcol('NAME')[2:] # Skip the first two sources, which are calibrators
mosaicname = fieldlist[-1][:-1]

#split the spw of water maser
split(vis=myviscal,outputvis=mosaicname+'_H2O.ms',spw='21',datacolumn='corrected',keepflags=False,field='2~15')
```

Now you have `*_H2O.ms`. You can check the emission of water maser by **plotms** or making a dirty map (preferred).


```python
#make a dirty cube and find the channels and fields with emission
import numpy as np
import analysisUtils as au

inputvis = mosaicname + '_H2O.ms'

uvcontsub_results  =  uvcontsub (vis = inputvis, outputvis = inputvis+'.contsub', fitorder=1, fitspec = "0:201~800;1050~1350;1450~1850") #fitspec is the line-free channel checked from plotms

#split the H2O spw and do uvcontsub
inputvis=inputvis+'.contsub'
imname=mosaicname+'_H2O.dirty'

#set the mapcenter for the mosaic
au.plotmosaic(inputvis, figfile=True, showplot=False)
ptgcenters = au.getRADecForFields(inputvis)
racenter = ( np.max(ptgcenters[0]) + np.min(ptgcenters[0]) )*0.5 # radian
if racenter < 0:
    racenter = 2*np.pi + racenter
deccenter = ( np.max(ptgcenters[1]) + np.min(ptgcenters[1]) )*0.5 # radian
mapcenter = 'J2000 '+str(racenter)+'rad '+str(deccenter)+'rad'

#run the tclean task. let parallel=True if you use mpicasa
default(tclean)
tclean(vis = inputvis,
  imagename = imname,
  specmode = 'cube',
  niter = 0,
  threshold = '12mJy',
  deconvolver = 'hogbom',
  gridder = 'mosaic',
  restfreq = '22.23508GHz',
  interactive = False,
  imsize = [3600,2500],
  cell = '0.2arcsec',
  weighting = 'briggs',
  robust = 0.5,
  outframe = 'lsrk',
  pblimit = 0.1,
  restoringbeam = 'common',
  phasecenter = mapcenter,
  perchanweightdensity = False,
#  parallel = True,
  mosweight = True)
```

You may find some water masers in some channels and fields. Split them and generate the selfcal table.


```python
#split the peak channel of water maser in one field. Here is an example for maser in field 10 and channel 432.
split(vis=myviscal, outputvis=mosaicname+'_H2O_F10_c432.ms', field='10', spw='21:432', datacolumn='corrected')

#do selfcal for the splited ms files (2 phase cal and 2 amp cal, or more if needed)
inputvis=mosaicname+'_H2O_F10_c432.ms'
listobs(vis=inputvis, listfile=inputvis+'.listobs.txt')
default(tclean)
tclean(vis = inputvis,
  imagename = mosaicname+'_H2O_F10_c432_initial',
  specmode = 'cube',
  niter = 1000000,
  threshold = '10mJy',
  deconvolver = 'hogbom',
  gridder = 'mosaic',
  restfreq = '22.23508GHz',
  interactive = True,
  imsize = 1024,
  cell = '0.2arcsec',
  weighting = 'briggs',
  robust = 0.5,
  outframe = 'lsrk',
  pblimit = 0.1,
  restoringbeam = 'common',
  perchanweightdensity = False,
  mosweight = True,
  usepointing = False,
  savemodel = 'modelcolumn')
gaincal(vis=inputvis, caltable='CMZ_28_H2O_F10_c432.pha0.cal',solint='inf',calmode='p',refant='ea09',gaintype='T', minsnr=3.0)
applycal(vis=inputvis,gaintable=['CMZ_28_H2O_F10_c432.pha0.cal'],calwt=False,applymode='calonly')

default(tclean)
tclean(vis = inputvis,
  imagename = mosaicname+'_H2O_F10_c432_pha0',
  specmode = 'cube',
  niter = 1000000,
  threshold = '10mJy',
  deconvolver = 'hogbom',
  gridder = 'mosaic',
  restfreq = '22.23508GHz',
  interactive = True,
  imsize = 1024,
  cell = '0.2arcsec',
  weighting = 'briggs',
  robust = 0.5,
  outframe = 'lsrk',
  pblimit = 0.1,
  restoringbeam = 'common',
  perchanweightdensity = False,
  mosweight = True,
  usepointing = False,
  savemodel = 'modelcolumn')
gaincal(vis=inputvis, caltable='CMZ_28_H2O_F10_c432.pha1.cal',solint='15s',calmode='p',refant='ea09',gaintype='T', minsnr=3.0)
applycal(vis=inputvis,gaintable=['CMZ_28_H2O_F10_c432.pha1.cal'],calwt=False,applymode='calonly')

default(tclean)
tclean(vis = inputvis,
  imagename = mosaicname+'_H2O_F10_c432_pha1',
  specmode = 'cube',
  niter = 1000000,
  threshold = '10mJy',
  deconvolver = 'hogbom',
  gridder = 'mosaic',
  restfreq = '22.23508GHz',
  interactive = True,
  imsize = 1024,
  cell = '0.2arcsec',
  weighting = 'briggs',
  robust = 0.5,
  outframe = 'lsrk',
  pblimit = 0.1,
  restoringbeam = 'common',
  perchanweightdensity = False,
  mosweight = True,
  usepointing = False,
  savemodel = 'modelcolumn')
gaincal(vis=inputvis,caltable='CMZ_28_H2O_F10_c432.amp.cal',solint='inf',calmode='ap',refant='ea09',gaintype='T',solnorm=True,gaintable='CMZ_28_H2O_F10_c432.pha1.cal')
applycal(vis=inputvis,gaintable=['CMZ_28_H2O_F10_c432.pha1.cal','CMZ_28_H2O_F10_c432.amp.cal'],calwt=False,applymode='calonly')

default(tclean)
tclean(vis = inputvis,
  imagename = mosaicname+'_H2O_F10_c432_amp',
  specmode = 'cube',
  niter = 1000000,
  threshold = '10mJy',
  deconvolver = 'hogbom',
  gridder = 'mosaic',
  restfreq = '22.23508GHz',
  interactive = True,
  imsize = 1024,
  cell = '0.2arcsec',
  weighting = 'briggs',
  robust = 0.5,
  outframe = 'lsrk',
  pblimit = 0.1,
  restoringbeam = 'common',
  perchanweightdensity = False,
  mosweight = True,
  usepointing = False,
  savemodel = 'modelcolumn')
gaincal(vis=inputvis,caltable='CMZ_28_H2O_F10_c432.amp2.cal',solint='60s',calmode='ap',refant='ea09',gaintype='T',solnorm=True,gaintable='CMZ_28_H2O_F10_c432.pha1.cal')
applycal(vis=inputvis,gaintable=['CMZ_28_H2O_F10_c432.pha1.cal','CMZ_28_H2O_F10_c432.amp2.cal'],calwt=False,applymode='calonly')

default(tclean)
tclean(vis = inputvis,
  imagename = mosaicname'_H2O_F10_c432_amp2',
  specmode = 'cube',
  niter = 1000000,
  threshold = '10mJy',
  deconvolver = 'hogbom',
  gridder = 'mosaic',
  restfreq = '22.23508GHz',
  interactive = True,
  imsize = 1024,
  cell = '0.2arcsec',
  weighting = 'briggs',
  robust = 0.5,
  outframe = 'lsrk',
  pblimit = 0.1,
  restoringbeam = 'common',
  perchanweightdensity = False,
  mosweight = True,
  usepointing = False,
  savemodel = 'modelcolumn')
```

Now you have the selfcal table: `*_H2O_F*_c*.pha1.cal` and `*_H2O_F*_c*.amp2.cal`. You can repeat the selfcal process if there are more masers in other fields.

## 3. Imaging the continuum

Split the continuum spws and run the following script in `cont`


```python
#split the cont spws
myviscal='../25A-303.sb48641652.eb48925296.60819.25704964121.ms'
tb.open(myviscal+'/FIELD')
fieldlist = tb.getcol('NAME')[2:] # Skip the first two sources, which are calibrators
mosaicname = fieldlist[-1][:-1]
split(vis=myviscal, outputvis=mosaicname+'_CONT.ms', spw='8~14,22~28', datacolumn='corrected', keepflags=False, field='2~15')
inputvis = mosaicname+'_CONT.ms'
listobs(vis=inputvis,listfile=inputvis+'.listobs.txt', overwrite=True)

#set the mapcenter for mosaic
au.plotmosaic(inputvis, figfile=True, showplot=False)
ptgcenters = au.getRADecForFields(inputvis)
racenter = ( np.max(ptgcenters[0]) + np.min(ptgcenters[0]) )*0.5 # radian
if racenter < 0:
    racenter = 2*np.pi + racenter
deccenter = ( np.max(ptgcenters[1]) + np.min(ptgcenters[1]) )*0.5 # radian

#Then, apply selfcal table and re-weight different spws. Use the "Name" of fields (can be find in the listobs.txt) to specify the field which to be applycaled. Note: for observations in 22B-193, spwmap = [[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]]

applycal(vis=inputvis, gaintable=['../water/CMZ_28_H2O_F6_c926.pha1.cal','../water/CMZ_28_H2O_F6_c926.amp2.cal'], spwmap=[[0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0]], field='CMZ_28e', calwt=False, applymode='calonly', flagbackup=False)
applycal(vis=inputvis, gaintable=['../water/CMZ_28_H2O_F14_c926.pha1.cal','../water/CMZ_28_H2O_F14_c926.amp2.cal'], spwmap=[[0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0]], field='CMZ_28m', calwt=False, applymode='calonly', flagbackup=False)
applycal(vis=inputvis, gaintable=['../water/CMZ_28_H2O_F3_c944.pha1.cal','../water/CMZ_28_H2O_F3_c944.amp2.cal'], spwmap=[[0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0]], field='CMZ_28b', calwt=False, applymode='calonly', flagbackup=False)
applycal(vis=inputvis, gaintable=['../water/CMZ_28_H2O_F12_c944.pha1.cal','../water/CMZ_28_H2O_F12_c944.amp2.cal'], spwmap=[[0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0]], field='CMZ_28k', calwt=False, applymode='calonly', flagbackup=False)
applycal(vis=inputvis, gaintable=['../water/CMZ_28_H2O_F11_c1469.pha1.cal','../water/CMZ_28_H2O_F11_c1469.amp2.cal'], spwmap=[[0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0]], field='CMZ_28j', calwt=False, applymode='calonly', flagbackup=False)
applycal(vis=inputvis, gaintable=['../water/CMZ_28_H2O_F5_c432.pha1.cal','../water/CMZ_28_H2O_F5_c432.amp2.cal'], spwmap=[[0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0]], field='CMZ_28d', calwt=False, applymode='calonly', flagbackup=False)

statwt(vis=inputvis, minsamp=2, datacolumn='data')

#Run tclean to produce the dirty map and the final image.
default(tclean)
imname = mosaicname+'_dirty.CONT'
mapcenter = 'J2000 '+str(racenter)+'rad '+str(deccenter)+'rad'
tclean(vis = inputvis,
  imagename = imname,
  specmode = 'mfs',
  deconvolver = 'mtmfs',
  nterms = 2,
  scales = [0,5,15,50,150],
  niter = 0,
  threshold = '0.07mJy', # measured rms is 0.023 mJy, here we use a threshold of 3sigma
  gridder = 'mosaic',
  interactive = False,
  imsize = [3600,2500],
  cell = '0.2arcsec',
  weighting = 'briggs',
  robust = 0.5,
  pblimit = 0.2,
  mosweight = True,
  usepointing = False,
  phasecenter = mapcenter,
  usemask = 'auto-multithresh', # See https://casaguides.nrao.edu/index.php/Automasking_Guide
  sidelobethreshold = 2.0,
  noisethreshold = 4.25,
  minbeamfrac = 0.3,
  lownoisethreshold = 1.5,
  negativethreshold = 0.,
  nsigma = 3.,
  cyclefactor = 3.0 # To suppress divergence
  )
 
default(tclean)
imname = mosaicname+'.CONT'
tclean(vis = inputvis,
  imagename = imname,
  specmode = 'mfs',
  deconvolver = 'mtmfs',
  nterms = 2,
  scales = [0,5,15,50,150],
  niter = 100000,
  threshold = '0.07mJy', # measured rms is 0.023 mJy, here we use a threshold of 3sigma
  gridder = 'mosaic',
  interactive = False,
  imsize = [3600,2500],
  cell = '0.2arcsec',
  weighting = 'briggs',
  robust = 0.5,
  pblimit = 0.2,
  mosweight = True,
  usepointing = False,
  phasecenter = mapcenter,
  usemask = 'auto-multithresh', # See https://casaguides.nrao.edu/index.php/Automasking_Guide
  sidelobethreshold = 2.0,
  noisethreshold = 4.25,
  minbeamfrac = 0.3,
  lownoisethreshold = 1.5,
  negativethreshold = 0.,
  nsigma = 3.,
  cyclefactor = 3.0 # To suppress divergence
  )

#export to fits file
exportfits(imagename=imname+'.image.tt0', fitsimage=imname+'.image.tt0.fits', dropdeg=True)
impbcor(imagename=imname+'.image.tt0',pbimage=imname+'.pb.tt0',outfile=imname+'.image.pbcor.tt0')
exportfits(imagename=imname+'.image.pbcor.tt0', fitsimage=imname+'.image.pbcor.tt0.fits', dropdeg=True)
exportfits(imagename=imname+'.pb.tt0', fitsimage=imname+'.pb.tt0.fits', dropdeg=True)
exportfits(imagename=imname+'.alpha', fitsimage=imname+'.alpha.fits', dropdeg=True)
exportfits(imagename=imname+'.alpha.error', fitsimage=imname+'.alpha.error.fits', dropdeg=True)
```

## 4. Imaging the spectral cubes (Taking NH3 (1,1) as an example)

To make a NH3 cube, you prefer to use a mask file from SWAG data. We provide a large mask file, "NH3_3-3_large.fits", on our website. You can cut it into smaller size that fit to your mosaic field via the following codes in **python**.


```python
from astropy.io import fits
from astropy.wcs import WCS
from reproject import reproject_interp
import numpy as np

# --- Configuration ---
input_fits = 'NH3_3-3_large.fits'
output_fits = 'NH3_3-3_small_cut.fits'

# Your J2000 mapcenter (can be find in previous mapcenter)
ra_deg = 266.350807 
dec_deg = -29.416454 

# Your image size
a = 3600
b = 2500

# 1. Load the Galactic Mask
hdul = fits.open(input_fits)
data_gal = hdul[0].data
header_gal = hdul[0].header
wcs_gal = WCS(header_gal)

# 2. Create the Target J2000 Header (The "Desired" Grid)
# We build the header that defines the axb box aligned with RA/Dec
target_header = header_gal.copy()
target_header['WCSAXES'] = 2
target_header['NAXIS1'] = a
target_header['NAXIS2'] = b
target_header['CTYPE1'] = 'RA---SIN'
target_header['CTYPE2'] = 'DEC--SIN'
target_header['CRVAL1'] = ra_deg
target_header['CRVAL2'] = dec_deg
target_header['CRPIX1'] = a / 2.0
target_header['CRPIX2'] = b / 2.0
# Ensure the cell size matches (0.2 arcsec)
target_header['CDELT1'] = -5.55555555555555e-05
target_header['CDELT2'] = 5.55555555555555e-05
target_header['CUNIT1'] = 'deg'
target_header['CUNIT2'] = 'deg'

# Remove Galactic-specific keys that cause rotation conflicts
for key in ['PV2_1', 'PV2_2', 'LONPOLE', 'LATPOLE']:
    if key in target_header: del target_header[key]

# 3. Perform the Reprojection (The Rotation)
# This rotates the Galactic pixels into the J2000 grid
print("Reprojecting and rotating Galactic mask to J2000...")
new_data, footprint = reproject_interp((data_gal, wcs_gal), target_header)

# 4. Cleanup and Save
# Ensure the mask remains binary (0 or 1) after interpolation
new_data = np.where(new_data > 0.5, 1.0, 0.0).astype(np.float32)

new_hdu = fits.PrimaryHDU(data=new_data, header=target_header)
new_hdu.writeto(output_fits, overwrite=True)
print(f"Successfully saved rotated J2000 mask: {output_fits}")
```

Then you need to transform the fits file in casa format:


```python
importfits(fitsimage='NH3_3-3_small_cut.fits', imagename='*_NH3.mask', overwrite=True)
```

Now, you can start the imaging process of NH3 (1,1) cube.


```python
import numpy as np
import analysisUtils as au

#Split the NH3 (1,1) spw
myviscal='../25A-303.sb48641652.eb48925296.60819.25704964121.ms'
split(vis=myviscal, outputvis=mosaicname+'_NH3_11.ms',field='2~15',spw='7',datacolumn='corrected')

#check the listobs
inputvis = mosaicname+'NH3_11.ms'
listobs(vis = inputvis, listfile = inputvis+'.listobs')

#apply the selfcal table
applycal(vis=inputvis, gaintable=['../../water/CMZ_28_H2O_F6_c926.pha1.cal','../../water/CMZ_28_H2O_F6_c926.amp2.cal'], spwmap=[[0],[0]], field='CMZ_28e', calwt=False, applymode='calonly', flagbackup=False)
applycal(vis=inputvis, gaintable=['../../water/CMZ_28_H2O_F14_c926.pha1.cal','../../water/CMZ_28_H2O_F14_c926.amp2.cal'], spwmap=[[0],[0]], field='CMZ_28m', calwt=False, applymode='calonly', flagbackup=False)
applycal(vis=inputvis, gaintable=['../../water/CMZ_28_H2O_F3_c944.pha1.cal','../../water/CMZ_28_H2O_F3_c944.amp2.cal'], spwmap=[[0],[0]], field='CMZ_28b', calwt=False, applymode='calonly', flagbackup=False)
applycal(vis=inputvis, gaintable=['../../water/CMZ_28_H2O_F12_c944.pha1.cal','../../water/CMZ_28_H2O_F12_c944.amp2.cal'], spwmap=[[0],[0]], field='CMZ_28k', calwt=False, applymode='calonly', flagbackup=False)
applycal(vis=inputvis, gaintable=['../../water/CMZ_28_H2O_F11_c1469.pha1.cal','../../water/CMZ_28_H2O_F11_c1469.amp2.cal'], spwmap=[[0],[0]], field='CMZ_28j', calwt=False, applymode='calonly', flagbackup=False)
applycal(vis=inputvis, gaintable=['../../water/CMZ_28_H2O_F5_c432.pha1.cal','../../water/CMZ_28_H2O_F5_c432.amp2.cal'], spwmap=[[0],[0]], field='CMZ_28d', calwt=False, applymode='calonly', flagbackup=False)

#rebin the channels
mstransform(vis = inputvis, outputvis = inputvis+'.rebin', datacolumn = 'all', chanaverage = True, chanbin = 2)

#substract the continuum
inputvis = inputvis+'.rebin'
uvcontsub_result = uvcontsub(vis = inputvis, outputvis = inputvis+'.contsub', fitspec = '0:150~650;1700~1900', fitorder = 1)
inputvis = inputvis+'.contsub'

#set the mapcenter
au.plotmosaic(inputvis, figfile=True, showplot=False)
ptgcenters = au.getRADecForFields(inputvis)
racenter = ( np.max(ptgcenters[0]) + np.min(ptgcenters[0]) )*0.5 # radian
if racenter < 0:
    racenter = 2*np.pi + racenter
deccenter = ( np.max(ptgcenters[1]) + np.min(ptgcenters[1]) )*0.5 # radian
mapcenter = 'J2000 '+str(racenter)+'rad '+str(deccenter)+'rad'

#do the tclean. let parallel = True if you use mpicasa.
default(tclean)
tclean(vis = inputvis,
  imagename = imname,
  specmode = 'cube',
  nchan = 900,
  width = '0.4km/s',
  start = '-180km/s',
  niter = 100000000,
  phasecenter = mapcenter,
  threshold = '12mJy',
  deconvolver = 'multiscale',
  scales = [0,5,15,50,150],
  gridder = 'mosaic',
  restfreq = '23.6944955GHz',
  interactive = False,
  imsize = [3600,2500],
  cell = '0.2arcsec',
  weighting = 'briggs',
  robust = 0.5,
  outframe = 'lsrk',
  pblimit = 0.1,
  restoringbeam = 'common',
  perchanweightdensity = True,
#  parallel = True,
  usepointing = False,
  mosweight = True,
  
  usemask='user',
  mask='../CMZ28_NH3.mask'
  )

#export the fits file
exportfits(imagename=imname+'.image', fitsimage=imname+'.image.fits', dropdeg=True)
impbcor(imagename=imname+'.image',pbimage=imname+'.pb',outfile=imname+'.image.pbcor')
exportfits(imagename=imname+'.image.pbcor', fitsimage=imname+'.image.pbcor.fits', dropdeg=True)
exportfits(imagename=imname+'.pb', fitsimage=imname+'.pb.fits', dropdeg=True)
```


