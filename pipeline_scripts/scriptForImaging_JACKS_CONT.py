
import numpy as np
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
myviscal = '../22B-193.sb43003623.eb43145710.59933.73626988426.ms'
listobs(myviscal, listfile=myviscal+'.listobs.txt', overwrite=True)

tb.open(myviscal+'/FIELD')
fieldlist = tb.getcol('NAME')[2:] # Skip the first two sources, which are calibrators
mosaicname = fieldlist[-1][:-1]

split(vis=myviscal, outputvis=mosaicname+'_CONT.ms', spw='2~5,12~16,24~30', datacolumn='corrected', keepflags=False, field='2~15')
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
#xaxis = 'channel'
#yaxis = 'amp'
#avgtime = '1e8'
#spw = '0~15'
#antenna = 'ea01'
#plotms(vis=inputvis, xaxis=xaxis, yaxis=yaxis, avgtime=avgtime, spw=spw, antenna=antenna, avgscan=True, iteraxis='spw',coloraxis='corr')

## Flag problematic channels
#mode = 'manual'
#spw = '3:0~10,8:0~10;32~38,1:30~35,2:6~11;20~25,6:20~25,7:46~50'
#flagbackup=False
#flagdata(vis=inputvis, mode=mode, spw=spw, flagbackup=flagbackup)
#mode = 'manual'
#spw = '0~15:0~1;62~63,4;9:58~63'
#flagbackup=False
#flagdata(vis=inputvis, mode=mode, spw=spw, flagbackup=flagbackup)

#mode = 'manual'
#spw = '6:0~25,7:0~63'
#antenna = 'ea11&ea22;ea05&ea17'
#flagbackup=False
#flagdata(vis=inputvis, mode=mode, spw=spw, antenna=antenna, flagbackup=flagbackup)
## RFI in spw 6/7/8?
## We should still run hifv_checkflag to remove RFI?

## Re-evaluate the weights
## Inputvis is split out from the corrected datacolumn and only has a 'data' column
statwt(vis=inputvis, minsamp=2, datacolumn='data')

## Self calibration using H2O masers, if possible
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
  #scales = [0,5,15,50,150],
  scales = [0,5,15,50],
  niter = 2000,
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
  sidelobethreshold = 2.0,
  noisethreshold = 4.25,
  minbeamfrac = 0.3,
  lownoisethreshold = 1.5,
  negativethreshold = 0.,
  nsigma = 3.,
  pbmask = 0.25, # make it bigger than pblimit (=0.2 by default)
  cyclefactor = 3.0 # To suppress divergence
  )

exportfits(imagename=imname+'.image.tt0', fitsimage=imname+'.image.tt0.fits', dropdeg=True)
impbcor(imagename=imname+'.image.tt0',pbimage=imname+'.pb.tt0',outfile=imname+'.image.pbcor.tt0')
exportfits(imagename=imname+'.image.pbcor.tt0', fitsimage=imname+'.image.pbcor.tt0.fits', dropdeg=True)
exportfits(imagename=imname+'.pb.tt0', fitsimage=imname+'.pb.tt0.fits', dropdeg=True)
exportfits(imagename=imname+'.alpha', fitsimage=imname+'.alpha.fits', dropdeg=True)
exportfits(imagename=imname+'.alpha.error', fitsimage=imname+'.alpha.error.fits', dropdeg=True)
