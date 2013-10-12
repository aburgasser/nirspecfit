; -----------------------------------------------------------------------------------------------------------------
;
; NIRSPECFIT_CALIBRATIONS
; FUNCTION TO SELECTED, UNZIP, AND READ IN CALIBRATION FILES
;
; -----------------------------------------------------------------------------------------------------------------

FUNCTION NIRSPECFIT_CALIBRATIONS, order, mref, settled=settled, saumon=saumon

common NIRSPECFIT, pfolder, calfolder, mfolder, outputfolder, $
 model_names, units, $
 nirspec_lref, nirspec_lrng, nirspec_res0, nirspec_npix, nirspec_pdisp, $
 samp_res, vsini_dv, vsini_min, vsini_max, $
 cvel

order = order(0) ; one at a time

; select model 
cfile = calfolder+'models_'
case 1 of
 keyword_set(settled): mref = 0
 keyword_set(saumon): mref = 1
 n_elements(mref) eq 0: mref = 0
 mref gt 1: message, 'Only Settled (/settled) Saumon (/saumon) models are available'
 else: mref=mref
endcase

cfile = cfile+model_names(mref)+'_order'
print, 'Fitting to '+model_names(mref)

; select order
cfile = cfile+strtrim(string(order),2)+'.dat'

; unzip calibration file
chk = file_search(cfile+'.gz') 
if (chk ne '') then spawn, 'gunzip '+cfile+'.gz'

; read in calibration file
chk = file_search(cfile) 
if (chk eq '') then message, 'Cannot find calibration file '+cfile

restore, file=cfile
;spawn, 'gzip '+cfile

return, str

end
