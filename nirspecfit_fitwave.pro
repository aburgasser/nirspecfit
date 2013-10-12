; NIRSPECFIT
; ROUTINE TO DO THE WAVELENGTH CALIBRATION ON CAL STAR

PRO NIRSPECFIT_FITWAVE, str, debug=debug, initialize=initialize


common NIRSPECFIT, pfolder, calfolder, mfolder, outputfolder, $
 model_names, units, $
 nirspec_lref, nirspec_lrng, nirspec_res0, nirspec_npix, nirspec_pdisp, $
 samp_res, vsini_dv, vsini_min, vsini_max, $
 cvel

param  = str.param
xcorrad0 = 500			; cross-correlation radius for first pixel centering
fitrng1 = 200			; fit range for slope fitting
xcorrad1 = 70			; cross-correlation radius for slope fitting
xcorrad2 = 50			; cross-correlation radius for correcting wavelength calibration
varrad = 3				; radius for computing variance spectrum
smbox = 50			; smoothing kernel
nloop = 20			; loops for chi-fit
f=''

; WORKING SECTION OF SPECTRUM
wsample = where(str.pixel ge str.prng(0) and str.pixel le str.prng(1),cntsample)

; FIRST FIT 2ND ORDER POLYNOMIAL TO DATA, SUBTRACT AND NORMALIZE
fit = poly_fit(str.pixel,smooth(str.flux,smbox),2)
flx = str.flux/poly(str.pixel,fit)

; CREATE A WEIGHTING SPECTRUM BASED ON VARIATION IN DATA
varwt = flx*0.
x = findgen(varrad*2+1)-varrad
for i=varrad,n_elements(varwt)-1-varrad do varwt(i) = stddev(flx(x+i))
varwt(1000:n_elements(varwt)-1) = 0

nirspecfit_model, str, param, mdl, trs, cont, px, flag=flag

if (keyword_set(initialize)) then begin

; FIRST USE GUESS FROM NIRSPEC INSTRUMENT PARAMETERS
 wave_param = nirspec_lref/str.order*[1.,nirspec_pdisp/nirspec_res0, 0.5*(nirspec_pdisp/nirspec_res0)^2]
 wv = poly(str.pixel-0.5*nirspec_npix,str.wave_param) 
 wpar = poly_fit(wv-str.lam0,str.pixel,2)
 param(0:str.wave_order) = [reform(wpar)]
 nirspecfit_model, str, param, mdl, trs, cont, px, flag=flag

 if (keyword_set(debug)) then begin
  window, 0
  !p.color=0
  !p.background=255
  !p.multi=0
  plot, str.pixel, flx/median(flx), xra=[0,max(str.pixel)], /xsty
  oplot, str.pixel, mdl/median(mdl), color=2
  oplot, [0,0]+str.prng(0), linestyle=1
  oplot, [0,0]+str.prng(1), linestyle=1
  print, param(0:str.wave_order)
  wdelete, 0
 endif

; CROSS-CORRELATE WHOLE SPECTRUM TO FIND ZERO POINT
XC_FULL: nirspecfit_model, str, param, mdl, trs, cont, px, flag=flag
 mdl0=mdl
 xcorl, flx, mdl, xcorrad0, sft
; sft = findgen(2*xcorrad0+1)-xcorrad0
; chis = sft*0.
; for jjj=0,2*xcorrad0 do begin
;  p = param
;  p(0) = p(0)+sft(jjj)
;  nirspecfit_model, str, p, mdl, trs, cont, px, flag=flag
;  chis(jjj) = total(varwt*(flx-median(flx))*(mdl-median(mdl)))/(median(flx)*median(mdl))
; endfor

; print, param(0:str.wave_order)
; mx = max(chis,loc)
 param(0) = param(0)+sft
 nirspecfit_model, str, param, mdl, trs, cont, px, flag=flag

if (keyword_set(debug)) then begin
 window, 0
 !p.color=0
 !p.background=255
 !p.multi=[0,1,2]
 xcorl, flx, mdl, xcorrad0, sft, plot=debug
 plot, str.pixel, flx/median(flx), xra=[0,max(str.pixel)], /xsty
 oplot, str.pixel, mdl-0.5, color=2
 oplot, [0,0]+str.prng(0), [0,10], linestyle=1
 oplot, [0,0]+str.prng(1), [0,10], linestyle=1
 !p.multi=0
 print, param(0:str.wave_order)
; read, f , prompt='Hit return to continue '
 wdelete, 0
endif

endif

; PERFORM CROSS CORRELATIONS ON STRONG FEATURES IN DESIRED RANGE 
; select strong peaks that don't significantly overlap
XC_SHARP: for jjj=0,nloop-1 do begin
 mdl0 = mdl
 s = sort(0.-varwt(wsample))
 refpx = str.pixel(wsample(s(0:50)))
 msk = refpx*0
 while i lt n_elements(refpx)-1 do begin
  w = where(abs(refpx-refpx(i)) lt 0.5*xcorrad2,cnt)
  if (cnt gt 0) then pop, refpx, w(1:cnt-1), /index
  i = i+1
 end
 refpx = refpx(0:20<(n_elements(refpx)-1))

 sftpx = refpx*0
 if (debug) then !p.multi=[0,5,5]
 for i=0,n_elements(refpx)-1 do begin
  rng = [(refpx(i)-xcorrad2)>0,(refpx(i)+xcorrad2)<(n_elements(str.pixel)-1)]
  xcorl, flx(rng(0):rng(1)), mdl(rng(0):rng(1)), xcorrad2/2., sft, plot=debug
  sftpx(i) = sft
 endfor
if (debug) then !p.multi=0

; fit curve to residual offset
fit = poly_fit(interpol(str.wave_model-str.lam0,px,refpx),sftpx,str.wave_order)
param(0:2) = param(0:2)+fit
nirspecfit_model, str, param, mdl, trs, cont, px, flag=flag

if (keyword_set(debug)) then begin
 window, 0
 !p.color=0
 !p.background=255
 !p.multi=[0,2,1]
 plot, refpx, sftpx, psym=4, /ysty, /xsty
 plot, str.pixel, flx/median(flx), xra=str.prng, /xsty, yra=[-0.5,1.2], /ysty
 oplot, str.pixel, mdl-0.25, color=4
 oplot, str.pixel, flx/median(flx)-mdl, color=2
 print, param(0:str.wave_order)
; read, f , prompt='Hit return to continue '
 wait, 1
 wdelete, 0
endif

endfor


FINISH: str.param = param

; CONVERT BACK TO WAVELENGTH TO ASSESS RANGE
w = where(px gt 0 and px le nirspec_npix-1)
fit = poly_fit(px,str.wave_model,str.wave_order)
str.wave_param = fit

return
end
