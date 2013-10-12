pro NIRSPECFIT_INITIALIZE

common NIRSPECFIT, pfolder, calfolder, mfolder, outputfolder, $
 model_names, units, $
 nirspec_lref, nirspec_lrng, nirspec_res0, nirspec_npix, nirspec_pdisp, $
 samp_res, vsini_dv, vsini_min, vsini_max, $
 cvel

; PHYSICAL CONSTANTS
cvel = 3.0e5				; c in km/s

; PROGRAM INFORMAITON
pfolder = '/Users/adam/idl/nirspecfit/'
mfolder = '/Users/adam/models/'
calfolder = pfolder+'/calibrations/'
model_names = ['BT-Settled-08','Saumon-13']
units = 'erg/s/cm2/micron'			; units of models (at stellar surface)

; NIRSPEC PARAMETERS
nirspec_lref = 76.56			; blaze x order
nirspec_lrng = [0.95,2.45]		; wavelength range of interest
nirspec_res0 = 10300.		; R-theta in arcseconds
;nirspec_res0 = 10000.		; R-theta in arcseconds (from Fitzgerald)
nirspec_pdisp = 0.144		; dispersion in spectral dimension in arcseconds 
nirspec_npix = 1024.		; width of detector

; FITTING PARAMETERS
samp_res =100000.			; sampling resolution = 3 km/s
vsini_dv = 1.0				; minimum, maximum and sampling for v sin i
vsini_min = vsini_dv				
vsini_max = 150.

return
end
