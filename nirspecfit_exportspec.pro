; NIRSPECFIT
; EXPORT CALIBRATED SPECTRUM

PRO nirspecfit_exportspec, str, file=file, calibrator=calibrator, plot=plot

common NIRSPECFIT, outfolder

tb='	'
if (n_elements(file) eq 0) then file = str.outfolder+'spectrum.txt'

; generate model
nirspecfit_model, str, str.param, mdl, trs, cont, px, flag=flag

; determine wavelength calibration
wv = interpol(str.wave_model,px,str.pixel)
wsample = where(str.pixel ge str.prng(0) and str.pixel le str.prng(1),cntsample)
 
; print out spectrum
openw, unit, file, /get_lun
printf, unit, '# Exported NIRSPECFIT extract file for '+str.name+' order '+strtrim(string(str.order),2)+' obtained on '+str.obsdate+ '(Julian Date '+strtrim(string(str.jdate),2)+')'
if (keyword_set(calibrator)) then printf, unit, '# Calibration source'
printf, unit, '# Data column contains telluric absorption and is not flat fielded'
printf, unit, '# Pixel	Wavelength (micron)	Data	Uncertainty	Best Model	Transmission	Continuum	Corrected Data'

for i=0,cntsample-1 do printf, unit, $
  strtrim(string(str.pixel(wsample(i))),2)+tb+$
  strtrim(string(wv(wsample(i))),2)+tb+$
  strtrim(string(str.flux(wsample(i))),2)+tb+$
  strtrim(string(str.noise(wsample(i))),2)+tb+$
  strtrim(string(mdl(wsample(i))),2)+tb+$
  strtrim(string(trs(wsample(i))),2)+tb+$
  strtrim(string(cont(wsample(i))),2)+tb+$
  strtrim(string(str.flux(wsample(i))/cont(wsample(i))/trs(wsample(i))),2)
close, unit
free_lun, unit

; plot spectrum
if (keyword_set(plot)) then nirspecfit_plotexport, file, strrep(file,'.txt','.eps'), name=str.name, /pixel
  
return
end
