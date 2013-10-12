; ITEMS TO DO
; improve initial wavelength fit using nirspec's known 
;	wavelength calibration - cross correlate with transmission function
; reduce sizes of stored structures
; correct noise spectra for data (earlier in the extraction process)
; better understanding on how to extract vsin i
; do fit as an ameoba routine?
; eliminate ringing (frequency mask)
; implement noise sampling (NO?)
; THERE IS A STRANGE LOADING ERROR

; MAIN PROGRAM

pro nirspecfit, inputfile, calfile=calfile, parreset=parreset, calreset=calreset, sampnoise=sampnoise, nofit=nofit, standardreset=standardreset, settled=settled, saumon=saumon, model_only=model_only, calonly=calonly, trng=trng, grng=grng, outfolder=outfolder, debug=debug

on_error, 0

common NIRSPECFIT, pfolder, calfolder, mfolder, outputfolder, $
 model_names, units, $
 nirspec_lref, nirspec_lrng, nirspec_res0, nirspec_npix, nirspec_pdisp, $
 samp_res, vsini_dv, vsini_min, vsini_max, $
 cvel

; constants
rvmax = 200.	; maximum RV allowed for intiial cross correlation
nnoise = 50	; number of noise loops to run
tb='	'
f=''

if (n_elements(outfolder) eq 0) then outputfolder ='fit/' else outputfolder=outfolder
if (file_search(outputfolder) eq '') then spawn, 'mkdir '+outputfolder

; initialize parameters
nirspecfit_initialize

; read in user's parameters 
if (n_elements(inputfile) eq 0) then inputfile = 'input.txt'
nirspecfit_fitparameters, inputfile, pstr, trng=trng,grng=grng, /reset
nirspecfit_dataparameters, pstr, dstr, /reset

; add/modify additional parameters
if (keyword_set(debug)) then pstr.debug = 1
pstr = create_struct(pstr,$
 'outputfolder',outputfolder)

; initialize calibrations
cstr = nirspecfit_calibrations(dstr.order,settled=settled,saumon=saumon)
wt = where(cstr.teff ge min(pstr.teffrng) and cstr.teff le max(pstr.teffrng),cntt)
wg = where(cstr.logg ge min(pstr.loggrng) and cstr.logg le max(pstr.loggrng),cntg)
wz = where(cstr.z ge min(pstr.zrng) and cstr.z le max(pstr.zrng),cntz)
wcld = where(cstr.cloud ge min(pstr.cldrng) and cstr.cloud le max(pstr.cldrng),cntcld)

if (keyword_set(nofit)) then goto, PARAMS


; FIRST FIT A0 STAR TO GET WAVELENGTH CALIBRATION

fitname = 'calibrator'
if (keyword_set(standardreset) or file_search(outputfolder+fitname+'.dat') eq '') then begin

 param = pstr.param
 str = create_struct(pstr,$
  'pixel',dstr.pixel,$
  'wave',dstr.wave,$
  'flux',dstr.cal_spectrum,$
  'noise',dstr.cal_spectrum_unc,$
  'mask',dstr.cal_mask,$
  'wave_model',cstr.wave_model,$
  'trans',cstr.trans,$
  'model',cstr.trans*0.+1.,$
  'fitname',fitname,$
  'model_vsini',cstr.vsini,$
  'dof',dstr.dof,$
  'vbary',dstr.vbary,$
  'vbary_flag',dstr.vbary_flag)
 str.nchain = 1000L			; short chain for calibrator

 print, ''
FITCAL: print, 'Fitting calibrator spectrum for '+dstr.name+' for wavelength scale'
 print, ''
    
 nirspecfit_fitcal, str, debug=pstr.debug, /report
 param = str.param
; read, f, prompt='Happy? (Y or n)'
; if (strlowcase(f) eq 'n') then goto, FITCAL 
endif else begin
 restore, file=outputfolder+fitname+'.dat'
 param = str.param
endelse


; plot out spectrum with just correction from A0V star
 plotstr = create_struct(pstr,$
  'pixel',dstr.pixel,$
  'wave',dstr.wave,$
  'flux',dstr.spectrum,$
  'noise',dstr.spectrum_unc,$
  'mask',dstr.mask,$
  'wave_model',cstr.wave_model,$
  'trans',cstr.trans,$
  'model',cstr.trans*0.+1.,$
  'fitname',fitname,$
  'model_vsini',cstr.vsini,$
  'dof',dstr.dof)
plotstr.param = param

nirspecfit_exportspec, plotstr, file=pstr.outputfolder+'nirspec_'+pstr.name+'_'+pstr.date+'_spectrum.txt', /plot

if (keyword_set(calonly)) then goto, FINISH



; NOW FIT SPECTRUM TO ALL MODELS USING THE LEARNED SET

param0 = param

for i=0,cntt-1 do begin
 for j=0,cntg-1 do begin
  for k=0,cntz-1 do begin
   for l=0,cntcld-1 do begin
  
  fitname = dstr.name+'_t'+strtrim(string(fix(cstr.teff(wt(i)))),2)+'_g'+strmid(strtrim(string(cstr.logg(wg(j))),2),0,3)+'_z'+strmid(strtrim(string(cstr.z(wz(k))),2),0,3)+'_cld'+strmid(strtrim(string(cstr.cloud(wcld(l))),2),0,3)
  fitstr = create_struct(pstr,$
   'pixel',dstr.pixel,$
   'flux',dstr.spectrum,$
   'noise',dstr.spectrum_unc,$
   'mask',dstr.mask,$
   'wave_model',cstr.wave_model,$
   'trans',cstr.trans,$
   'model',reform(cstr.models(*,wt(i),wg(j),wz(k),wcld(l),*)),$
   'model_vsini',cstr.vsini,$
   'dof',dstr.dof,$
   'vbary',dstr.vbary,$
   'vbary_flag',dstr.vbary_flag,$
   'fitname',fitname)

; INITIALIZE THE RV WITH A CROSS-CORRELATION
 if (i eq 0 and j eq 0 and k eq 0) then begin
  nirspecfit_model, fitstr, param, mdl, trs, cont, px
;  w = where(px ge min(fitstr.pixel) and px le max(fitstr.pixel))
;  fit = poly_fit(px(w),fitstr.wave_model(w),fitstr.wave_order,/double)
  conv = cvel*fitstr.wave_param(1)/fitstr.wave_param(0)
  pixrng = min([rvmax/conv,n_elements(mdl)/2.1])
  xcorl, fitstr.flux, mdl, ceil(pixrng), sft, /fine, plot=pstr.debug
  param0(n_elements(param)-2) = sft*conv
  param = param0
 endif
  
; adopt source parameters from previous fit
  fitstr.param = param0
;  fitstr.param(n_elements(param)-3:n_elements(param)-1) = param(n_elements(param)-3:n_elements(param)-1)

; modify variations to calibrator parameters
 fitstr.pvar(1:fitstr.wave_order) = 0.			; keep all but wavelength offset fixed
 fitstr.pvar(0) = 5.
 fitstr.pvar(n_elements(param)-4) = 0.		; keep profile broadening fixed
; fitstr.pvar(n_elements(param)-3) = 0.01
; if (i eq 0 and j eq 0 and k eq 0) then fitstr.pvar(n_elements(param)-2) = 5.

; for kkk=0,fitstr.wave_order do fitstr.pvar(kkk) = 3./(median(fitstr.wave_model)^kkk)
 
  print, ''
  print, 'Fitting '+fitstr.fitname
  print, ''
  nirspecfit_fitspec, fitstr, debug=pstr.debug, /report

  param = fitstr.param_best
  
  
  
; now repeat by sampling noise
  if (keyword_set(sampnoise)) then begin
   fitname = dstr.name+'_t'+strtrim(string(fix(cstr.teff(wt(i)))),2)+'_g'+strmid(strtrim(string(cstr.logg(wg(j))),2),0,3)+'_z'+strmid(strtrim(string(cstr.z(wz(k))),2),0,3)+'_cld'+strmid(strtrim(string(cstr.cloud(wcld(l))),2),0,3)+'_noise'
   if (max(strpos(tag_names(fitstr),'PARAM_NOISE')) lt 0) then fitstr = create_struct(pstr,$
    'param_noise',param*0.,$
    'param_noise_unc',param*0.)
   paramarr = fltarr(n_elements(param),nnoise)
   chiarr = fltarr(nnoise)
   flux0 = fitstr.flux
   print, ''
   print, 'Sampling noise for '+fitstr.fitname
   print, ''
   
   for mmm=0,nnoise-1 do begin
    nfitstr = fitstr
    nfitstr.flux = flux0+randomn(sd,n_elements(fitstr.flux))*fitstr.noise
    nfitstr.nchain = 300		; LIMIT THE NUMBER OF FITS
    nirspecfit_fitspec, nfitstr, debug=dstr.debug, /REPORT
    paramarr(*,mmm) = nfitstr.param_best
    chiarr(mmm) = nfitstr.minchi
   endfor
   
;   prob = 2.(1.-f_pdf(chiarr/min(chiarr),str.dof,str.dof))
   for mmm=0,n_elements(param)-1 do begin
;    m1 = total(prob*paramarr(mmm,*))/total(prob)
;    m2 = sqrt(abs(total(prob*(m1^2-paramarr(mmm,*)^2))/total(prob)))
;    pstr.param(mmm) = m1
;    pstr.param_opt(mmm) = m1
;    pstr.param_opt_unc(mmm) = m2
    fitstr.param_noise(mmm) = mean(paramarr(mmm,*))
    fitstr.param_noise_unc(mmm) = stddev(paramarr(mmm,*))
   endfor
;   str = create_struct(str,$
;    'param_noise', fitstr.param_noise,$
;    'param_noise_unc', fitstr.param_noise_unc)
    
   save, fitstr, file=str.outputfolder+fitstr.fitname+'.dat'
   
  endif

   endfor
  endfor
 endfor
endfor



; COMPILE PARAMETERS
PARAMS: fitname = pstr.name+'_t'+strtrim(string(fix(cstr.teff(wt(0)))),2)+'_g'+strmid(strtrim(string(cstr.logg(wg(0))),2),0,3)+'_z'+strmid(strtrim(string(cstr.z(wz(0))),2),0,3)+'_cld'+strmid(strtrim(string(cstr.cloud(wcld(0))),2),0,3)
if (file_search(pstr.outputfolder+fitname+'.dat') eq '') then $
 print, 'Cannot perform mean parameter determination; cannot find file '+pstr.outputfolder+fitname+'.dat' else begin

for i=0,cntt-1 do begin
 for j=0,cntg-1 do begin
  for k=0,cntz-1 do begin
   for l=0,cntcld-1 do begin

   file = pstr.outputfolder+pstr.name+'_t'+strtrim(string(fix(cstr.teff(wt(i)))),2)+'_g'+strmid(strtrim(string(cstr.logg(wg(j))),2),0,3)+'_z'+strmid(strtrim(string(cstr.z(wz(k))),2),0,3)+'_cld'+strmid(strtrim(string(cstr.cloud(wcld(l))),2),0,3)+'.dat'
   restore, file=file
   p = dblarr(n_elements(str.param)+5,n_elements(str.param_mcmc(0,*)))
   p(0:n_elements(str.param),*) = str.param_mcmc
   p(n_elements(str.param)+1,*) = cstr.teff(wt(i))
   p(n_elements(str.param)+2,*) = cstr.logg(wg(j))
   p(n_elements(str.param)+3,*) = cstr.z(wz(k))
   p(n_elements(str.param)+4,*) = cstr.cloud(wcld(l))
   if (i+j+k+l eq 0) then params = p else params = [[params],[p]]
   mn = min(p(n_elements(str.param),*),loc)
   bp = reform(p(*,loc))
   if (i+j+k+l eq 0) then bparams = bp else bparams = [[bparams],[bp]]

   endfor
  endfor
 endfor
endfor

 
; chi-square values
chisq = reform(params(n_elements(str.param),*))/str.dof
bchisq = reform(bparams(n_elements(str.param),*))/str.dof
mn = min(bchisq,loc)
minchi = bchisq(loc)

; populate structure
str = create_struct(str,$
 'param_mcmc_all', params,$
 'chisq_mcmc_all', chisq,$
 'best_params',bparams,$
 'best_chisq',bchisq)
str.minchi = minchi
str.param = reform(str.best_params(0:n_elements(str.param)-1,loc))
l1 = where(cstr.teff eq str.best_params(n_elements(str.param)+1,loc))
l2 = where(cstr.logg eq str.best_params(n_elements(str.param)+2,loc))
l3 = where(cstr.z eq str.best_params(n_elements(str.param)+3,loc))
l4 = where(cstr.cloud eq str.best_params(n_elements(str.param)+4,loc))
mn = min(abs(str.best_params(n_elements(str.param)-1,loc)-str.model_vsini),l5)
str.model = reform(cstr.models(*,l1(0),l2(0),l3(0),l4(0),*))


; SUMMARIZE RESULTS

nirspecfit_summarize, str

endelse



FINISH: return
end



