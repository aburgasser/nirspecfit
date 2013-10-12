; NIRSPECFIT_MODEL
; ROUTINE TO CREATE MODEL

PRO nirspecfit_model, str, param, mdl, trs, cont, px, flag=flag, nospin=nospin

common NIRSPECFIT, pfolder, calfolder, mfolder, outputfolder, $
 model_names, units, $
 nirspec_lref, nirspec_lrng, nirspec_res0, nirspec_npix, nirspec_pdisp, $
 samp_res, vsini_dv, vsini_min, vsini_max, $
 cvel

ilam = [0,str.wave_order]		
iflx = ilam(1)+1+[0,str.cont_order]
ibroad = iflx(1)+1
itrans = ibroad+1
ivel = itrans+1
ivrot = ivel+1
flag = 0

; reference frame
; wavelength -> pixel
px = poly(str.wave_model-str.lam0,param(ilam(0):ilam(1)))
w = where(px ge min(str.pixel)-100. and px le max(str.pixel)+100.,cnt)
if (cnt eq 0) then begin
 flag =1
 print, 'warning - wavelength scale did not converge'
 return
endif

; transmission
trs = (interpol(nirspecfit_broaden(str.wave_model(w),str.trans(w),param(ibroad),/gaussian),px(w),str.pixel,/quad))^(param(itrans))

; continuum
cont = poly(str.pixel,param(iflx(0):iflx(1)))

; model of source
; wavelength -> pixel
px_s = poly((str.wave_model*(1.+param(ivel)/cvel))-str.lam0,param(ilam(0):ilam(1)))
w = where(px_s ge min(str.pixel)-100. and px_s le max(str.pixel)+100.,cnt)
if (cnt eq 0) then begin
 flag =1
 return
endif

; choose model if necessary
mn = min(abs(param(ivrot)-str.model_vsini),vloc)
param(ivrot) = str.model_vsini(vloc)
mdl0 = reform(str.model(*,(vloc<(n_elements(str.model(0,*))-1))))

; CHANGE THIS TO INTEGRATE OVER PIXEL
; USE ????
mdl = interpol(nirspecfit_broaden(str.wave_model(w),mdl0(w),param(ibroad),/gaussian),px_s(w),str.pixel,/quad)*cont*trs

return
end
