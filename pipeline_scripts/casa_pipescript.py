context = h_init()
context.set_state('ProjectSummary', 'proposal_code', 'VLA/null')
context.set_state('ProjectSummary', 'proposal_title', 'unknown')
context.set_state('ProjectSummary', 'piname', 'unknown')
context.set_state('ProjectSummary', 'observatory', 'Karl G. Jansky Very Large Array')
context.set_state('ProjectSummary', 'telescope', 'EVLA')
context.set_state('ProjectStructure', 'ppr_file', 'PPR.xml')
context.set_state('ProjectStructure', 'recipe_name', 'hifv_cal')
# Define the raw visibility data, which should be retrieved (and unpacked) from the NRAO archive
myvis = '22B-193.sb43003623.eb43145710.59933.73626988426'
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
