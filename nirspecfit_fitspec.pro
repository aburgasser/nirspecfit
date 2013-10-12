; NIRSPECFIT
; FIT SPECTRA

PRO NIRSPECFIT_FITSPEC, str, debug=debug, report=report, calibrator=calibrator

common NIRSPECFIT, pfolder, calfolder, mfolder, outputfolder, $
 model_names, units, $
 nirspec_lref, nirspec_lrng, nirspec_res0, nirspec_npix, nirspec_pdisp, $
 samp_res, vsini_dv, vsini_min, vsini_max, $
 cvel


tek_color
tb='	'
chicut = 0.95
rvmax = 200.			; maximum absolute RV allowed
vsinimin = 0.1			; minimum vsini allowed
vsinimax = 150.		; maximum vsini allowed
sgmax = 15.			; maximum instrumental broadening allowed
sgmin = 0.1			; minimum instrumental broadening allowed
tamin = 0.2			; minimum transmission exponent

; MCMC for wavelength parameters of model
wvar = where(str.pvar ne 0.,cntvar)
;print, str.param
param = str.param
store = dblarr(n_elements(param)+1,str.nchain*cntvar)

; where we want to measure the wavelength solution
wsample = where(str.pixel ge str.prng(0) and str.pixel le str.prng(1),cntsample)

; initialize model
ilam = [0,str.wave_order]		
iflx = ilam(1)+1+[0,str.cont_order]
ibroad = iflx(1)+1
itrans = ibroad+1
ivel = itrans+1
ivrot = ivel+1

nirspecfit_model, str, param, mdl, trs, cont, px, flag=flag
if (flag eq 1) then message, 'Warning: could not initialize with input wavelength solution'

; INITIALIZE CONTINUUM FIT
w = where(mdl ne 0. and str.mask ne 0)
param(iflx(0):iflx(1)) = poly_fit(str.pixel(w),smooth(str.flux(w)/(mdl(w)/cont(w)),11),str.cont_order,/double)

; initial chi2
var = str.mask*(str.flux-mdl)^2/str.noise^2
var = var(wsample)
s = sort(abs(var))
chi0 = total(var(s(0:((n_elements(var)*str.chicut)<(n_elements(var)-1)))))
;chi0 = total(var)

; start MCMC loop
param0 = param

for i=0L,str.nchain-1L do begin
 for j=0,cntvar-1 do begin
SETPAR: param(wvar(j)) = param(wvar(j))+randomn(sd,/double)*str.pvar(wvar(j))

; force sg if needed
  param(ibroad) = (param(ibroad)>sgmin)<sgmax
; force ta if needed
  param(itrans) = (param(itrans))>tamin
; force rv if needed
  param(ivel) = ((param(ivel))>(-1.*rvmax))<rvmax
; force vsini if needed
;  param(ivrot) = ((param(ivrot))>vsinimin)<vsinimax
; interpolate vsini onto model grid
  mn = min(abs(param(ivrot)-str.model_vsini),vloc)
  param(ivrot) = str.model_vsini(vloc)

; create model
 nirspecfit_model, str, param, mdl, trs, cont, px
 
; compute variance
 diff = str.mask*(str.flux-mdl)
 var = diff^2/str.noise^2
 var = var(wsample)
 s = sort(var)
 chi = total(var(s(0:((n_elements(var)*str.chicut)<(n_elements(var)-1)))))
; chi = total(var)

; check for improvement and iterate parameters
 if (randomu(sd) le exp(-0.5*(chi-chi0))) then begin
   param0 = param
   chi0 = chi
 endif else $
   param = param0

   store(*,i*cntvar+j) = [param0,chi0]

  if (keyword_set(debug) and (chi0 eq chi)) then begin
   print, i, param0, chi0
   window
   !p.multi=[0,1,2]
   !p.color=0   
   nirspecfit_model, str, param0, mdl, trs, cont, px
   plot, str.pixel, str.mask*str.flux, xra=str.prng, /ysty, /xsty
   oplot, str.pixel, mdl, color=3
   oplot, str.pixel, mdl/trs, color=2
   plot, str.pixel, str.mask*diff/str.noise, yra=[-1.,1]*10., xra=str.prng, /xsty, /ysty
   !p.multi=0
   wait, 0.1
  endif

 endfor
endfor


; determine best fit parameters to structure
chisq = reform(store(n_elements(param),*))
best = str.param
optimal = str.param
optimal_e = str.param
for i=0,n_elements(str.param)-1 do begin
 nirspecfit_bestparam, reform(store(i,*)), chisq, str.dof, bst, opt, opte, prob
 best(i) = bst
 optimal(i) = opt
 optimal_e(i) = opte
endfor
minchi = min(chisq)
print, minchi, best

; repopulate derived parameters
str.param = best
str.wave_param = reform(best(ilam(0):ilam(1)))
str.cont_param = reform(best(iflx(0):iflx(1)))
str.sg = reform(best(ibroad))
str.ta = reform(best(itrans))
str.rv = reform(best(ivel))
str.vsini = reform(best(ivrot))

str.param_best = best
str.param_opt = optimal
str.param_opt_unc = optimal_e
if (max(strpos(tag_names(str),'PARAM_MCMC')) eq -1.) then begin
 str = create_struct(str,$
  'param_mcmc',store,$
  'chisq_mcmc', chisq/str.dof) 
 endif
str.param_mcmc = store
str.chisq_mcmc = chisq/str.dof
str.minchi = minchi
 
save, str, file=outputfolder+str.fitname+'.dat'


; REPORTING
if keyword_set(report) then nirspecfit_report, str, /eps
nirspecfit_exportspec, str, file=str.outputfolder+str.fitname+'_spectrum.txt', calibrator=calibrator, /plot
    
return
end
