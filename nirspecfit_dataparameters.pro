;
; NIRSPECFIT_DATAPARAMETERS
;
; FUNCTION TO READ IN ALL DATA PARAMETERS

PRO nirspecfit_dataparameters, pstr, dstr, reset=reset

common NIRSPECFIT, pfolder, calfolder, mfolder, outputfolder, $
 model_names, units, $
 nirspec_lref, nirspec_lrng, nirspec_res0, nirspec_npix, nirspec_pdisp, $
 samp_res, vsini_dv, vsini_min, vsini_max, $
 cvel


on_error, 0
tb = '	'
monstr = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']


if (n_elements(savefile) eq 0) then savefile = 'data_parameters.dat'
if (file_search(pstr.datafile) eq '') then message, 'Cannot find data input file '+pstr.datafile

if (file_search(outputfolder+savefile) eq '' or keyword_set(reset)) then begin


; RESTORE STRUCTURE FROM NIRSPEC REDUCTION
if (strpos(pstr.datafile,'.dat') gt 0 and file_search(pstr.datafile) ne '') then begin
 restore, file=pstr.datafile
 dstr = str			; rename structure
; fix as nirspec_reduction catches up
 if (max(strpos(tag_names(str),'WAVE_PARAM')) eq -1) then $
   dstr = create_struct(dstr,'wave_param',nirspec_lref/dstr.order*[1.,nirspec_pdisp/nirspec_res0, 0.5*(nirspec_pdisp/nirspec_res0)^2])

; OTHERWISE READ IN DATA
endif else begin

; POPULATE INFO FROM FIT PARAMETER FILE
 dstr = create_struct($
  'name',pstr.name,$
  'coordinate',pstr.coordinate,$
  'designation',pstr.designation,$
  'date',pstr.date,$
  'order',pstr.order,$
  'jd',0L,$
  'airmass',1.0,$
  'snr',20.)

; compute coordinates if needed
if (dstr.coordinate eq '' and dstr.designation ne '') then $
 dstr.coordinate = designation2coord(dstr.designation)
ra = 0.
dec = 0.
if (dstr.coordinate ne '') then begin
 tmp = minsec(dstr.coordinate)
 ra = tmp(0)
 dec = tmp(1)
endif

; READ IN SPECTRAL DATA
; source
; original REDSPEC output format
if (pstr.data_type eq 0) then begin
  readcol, pstr.datafile, pixel, wave, flux, format='i,f,f', /silent, comment='#'
; must fudge noise array for this format
  noise = flux/dstr.snr+1./dstr.snr
; pixel, wave, flux, flux_unc format
endif else $
 readcol, pstr.datafile, pixel, wave, flux, noise, format='i,f,f,f', /silent, comment='#'

scl = median(flux)
flux = flux/scl
noise=noise/scl

; calibrator
if (file_search(calfile) eq '') then message, 'Cannot find data file '+calfile

; original REDSPEC output format
if (pstr.data_type eq 0) then begin
  readcol, pstr.calfile, pixel, wave_cal, flux_cal, format='i,f,f', /silent, comment='#'
; must fudge noise array for this format
  noise_cal = flux_cal/100.+1./100.
endif else $
 readcol, calfile, pstr.pixel, wave, flux_cal, noise_cal, format='i,f,f,f', /silent, comment='#'

scl = median(flux_cal)
flux_cal = flux_cal/scl
noise_cal=noise_cal/scl

w = where(pixel ge pstr.prng(0) and pixel le pstr.prng(1))

dstr = create_struct(dstr,$
 'ra',ra,$
 'dec',dec,$
 'pixel',pixel(w),$
 'wave',wave(w),$
 'spectrum',flux(w),$
 'spectrum_unc',noise(w),$
 'cal_spectrum',flux_cal(w),$
 'cal_spectrum_unc',noise_cal(w))

endelse

; CONVERT OBSDATE TO JDATE IF NEEDED
if (dstr.jd eq 0 and dstr.date ne '') then begin
 if (strlen(dstr.date) eq 6) then $
 dstr.jd = julday(fix(strmid(dstr.date,2,2)),fix(strmid(dstr.date,4,2)),fix('20'+strmid(dstr.date,0,2)))
 if (strlen(str.obsdate) eq 8) then $
 dstr.jd = julday(fix(strmid(dstr.date,4,2)),fix(strmid(dstr.date,6,2)),fix(strmid(dstr.date,0,4)))
endif

; COMPUTE BARYCENTRIC VELOCITY FOR CORRECTION IF NEEDED
if (max(strpos(tag_names(dstr),'VBARY')) eq -1 and dstr.jd ne 0) then begin
  rvconv = [cos(dstr.dec/!radeg)*cos(dstr.ra/!radeg),cos(dstr.dec/!radeg)*sin(dstr.ra/!radeg),sin(dstr.dec/!radeg)]
 baryvel, dstr.jd, 2000,vh,vb
 dstr = create_struct(dstr,$
  'vbary', total(vb*rvconv),$
  'vbary_flag', 1)
endif

; GENERATE MASK AND CALCULATE DOF
mask = dstr.spectrum*0.+1.
w = where(dstr.spectrum le 0,cnt)
if (cnt gt 0) then mask(w) = 0.
cal_mask = dstr.cal_spectrum*0.+1.
w = where(dstr.cal_spectrum le 0,cnt)
if (cnt gt 0) then cal_mask(w) = 0.
dstr = create_struct(dstr,$
 'mask',mask,$
 'cal_mask',cal_mask,$
 'dof',total(mask)-n_elements(pstr.param)-1.,$
 'cal_dof',total(cal_mask)-n_elements(pstr.param)-1.)

; push things back to fit parameter array
pstr.name = dstr.name
if (max(strpos(tag_names(dstr),'DESIGNATION')) gt -1) then pstr.designation = dstr.designation
if (max(strpos(tag_names(dstr),'COORDINATE')) gt -1) then pstr.coordinate = dstr.coordinate
pstr.date = dstr.date
pstr.jdate = dstr.jd
pstr.obsdate = '20'+strmid(dstr.date,0,2)+' '+monstr(fix(strmid(dstr.date,2,2))-1)+' '+strmid(dstr.date,4,2)



save, pstr, dstr, file=outputfolder+savefile
endif else restore, file=outputfolder+savefile

FINISH: return

end
