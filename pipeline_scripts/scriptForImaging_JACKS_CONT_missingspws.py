
import numpy as np
import os
import analysisUtils as au

##### Notes ###################################################
# In principle all JACKS SBs use the same setup, so the SPW indices should be the same for all SBs
# SPWs:
#   11 -- NH3 1,1
#   10 -- NH3 2,2
#   9  -- NH3 3,3
#   8  -- NH3 4,4
#   7  -- NH3 5,5
#   6  -- H64a
#   23 -- H2O
#   22 -- H66a
#   21 -- NH3 5,4 5( 4)0a- 5( 4)0s
#   20 -- NH3 4,3 4( 3)0a- 4(-3)0s
#   19 -- NH3 3,2 3( 2)0a- 3( 2)0s
#   18 -- NH3 7,6 7( 6)0a- 7(-6)0s
#   17 -- NH3 2,1 2( 1)0a- 2( 1)0s
#   2~5,12~16,24~30 -- continuum
#
###############################################################

#### Image cm continuum ################################
myviscal = '../22B-193_sb42958191_1_1.59903.74775888889.ms'
listobs(myviscal, listfile=myviscal+'.listobs.txt', overwrite=True)
# Two spws are missing!!!

fieldlist = vishead(myviscal, mode='get',hdkey='field') 
mosaicname = fieldlist[0][-1][:-1]

split(vis=myviscal, outputvis=mosaicname+'_CONT.ms', spw='2~5,12~14,22~28', datacolumn='corrected', keepflags=False, field='2~15')
inputvis = mosaicname+'_CONT.ms'
listobs(vis=inputvis,listfile=inputvis+'.listobs.txt', overwrite=True)

# Inspect the mosaic and find the center
au.plotmosaic(inputvis, figfile=True, showplot=False)
ptgcenters = au.getRADecForFields(inputvis)
racenter = ( np.max(ptgcenters[0]) + np.min(ptgcenters[0]) )*0.5 # radian
if racenter < 0:
    racenter = 2*np.pi + racenter
deccenter = ( np.max(ptgcenters[1]) + np.min(ptgcenters[1]) )*0.5 # radian

# Check the spws
xaxis = 'channel'
yaxis = 'amp'
avgtime = '1e8'
spw = '0~13'
antenna = 'ea01'
#plotms(vis=inputvis, xaxis=xaxis, yaxis=yaxis, avgtime=avgtime, spw=spw, antenna=antenna, avgscan=True, iteraxis='spw',coloraxis='corr')

# Flag edge channels of the basebands, and channels around strong lines
mode = 'manual'
spw = '1:30~35,2:6~11,3:0~10,6:0~10;33~38,7:58~63'
# Three strong lines:
## NH3 (1,1) -- 22.694 GHz
## NH3 (3,3) -- 23.870 GHz
## H2O  --  22.223 GHz
## Also noisy edges in spws 3, 6, 7
flagbackup = False
flagdata(vis=inputvis, mode=mode, spw=spw, flagbackup=flagbackup)

## Re-evaluate the weights
## Inputvis is split out from the corrected datacolumn and only has a 'data' column
statwt(vis=inputvis, minsamp=2, datacolumn='data')

## Self calibration using H2O masers?
#applycal(vis='mosaicname_CONT.ms', spwmap=[[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]], gaintable=['../mosaicname_H2O/selfcal/mosaicname.pha1.cal','../mosaicname_H2O/selfcal/mosaicname.amp.cal'], calwt=False, applymode='calonly', flagbackup=False)

##################################
default(tclean)
imname = mosaicname+'.CONT'
mapcenter = 'J2000 '+str(racenter)+'rad '+str(deccenter)+'rad'
tclean(vis = inputvis,
  imagename = imname,
  specmode = 'mfs',
  deconvolver = 'mtmfs',
  nterms = 2,
  scales = [0,5,15,50],
  niter = 100000000000,
  threshold = '0.1mJy', # measured rms is 0.03 mJy, here we use a threshold of 3sigma
  gridder = 'mosaic',
  interactive = False,
  imsize = 4096,
  cell = '0.2arcsec',
  weighting = 'briggs',
  robust = 0.5,
  pblimit = 0.2,
  mosweight = True,
  usepointing = False,
  phasecenter = mapcenter,
  usemask = 'auto-multithresh', # See https://casaguides.nrao.edu/index.php/Automasking_Guide
  sidelobethreshold = 1.25, # Use the 7m parameters as recommended
  noisethreshold = 5.0,
  minbeamfrac = 0.1,
  lownoisethreshold = 2.0,
  negativethreshold = 0.,
  nsigma = 3.,
  pbmask = 0.25, # make it bigger than pblimit (=0.2 by default)
  cyclefactor = 3.5 # To suppress divergence
  )

exportfits(imagename=imname+'.image.tt0', fitsimage=imname+'.image.tt0.fits', dropdeg=True)
impbcor(imagename=imname+'.image.tt0',pbimage=imname+'.pb.tt0',outfile=imname+'.image.pbcor.tt0')
exportfits(imagename=imname+'.image.pbcor.tt0', fitsimage=imname+'.image.pbcor.tt0.fits', dropdeg=True)
exportfits(imagename=imname+'.pb.tt0', fitsimage=imname+'.pb.tt0.fits', dropdeg=True)
exportfits(imagename=imname+'.alpha', fitsimage=imname+'.alpha.fits', dropdeg=True)
exportfits(imagename=imname+'.alpha.error', fitsimage=imname+'.alpha.error.fits', dropdeg=True)
