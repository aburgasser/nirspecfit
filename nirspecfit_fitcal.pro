; NIRSPECFIT
; ROUTINE TO DO THE WAVELENGTH AND TRANSMISSION CALIBRATION ON CAL STAR

PRO NIRSPECFIT_FITCAL, str, debug=debug, report=report


common NIRSPECFIT, pfolder, calfolder, mfolder, outputfolder, $
 model_names, units, $
 nirspec_lref, nirspec_lrng, nirspec_res0, nirspec_npix, nirspec_pdisp, $
 samp_res, vsini_dv, vsini_min, vsini_max, $
 cvel


on_error, 0
nsamp = 1000.
pixborder = 30		; border of pixels for wave->pixel conversion
samppix0 = 400	; radius of pixels to shift for initial test
samppix = 40		; radius of pixels to shift for later tests
delpix = 10		; pixel shift range; must be at most 1/2 samppix
nloop = 20		; number of  big loops to make 
nloopw = 20		; number of  loops to make for wavelength shifts
shortchain = 1000	; MCMC chain length 
pxtresh = 1.5		; treshhold for pixel correction scale
smfact = 50		; smoothing bin size
sftfitorder = 7		; polynomial order for fitting shift array
f=''
mxshift = 20		; maximum shift allowed
tek_color

;debug=1

param = str.param
ilam = [0,str.wave_order]		
iflx = ilam(1)+1+[0,str.cont_order]
ibroad = iflx(1)+1
itrans = ibroad+1
ivel = ibroad+2
wsample = where(str.pixel ge str.prng(0) and str.pixel le str.prng(1),cntsample)

; FIT WAVELENGTH SOLUTION FOR STANDARD
nirspecfit_fitwave, str, debug=0, /initialize
param = str.param
nirspecfit_model, str, param, mdl, trs, cont, px, flag=flag

; FIT REST OF PARAMETERS
POSTWAVE: parr = dblarr(n_elements(param)+1,nloop)
parr(n_elements(param),*) = 1.d9
loop = 0

; INITIALIZE CONTINUUM
SETCONT: w = where(mdl ne 0. and str.mask ne 0)
param(iflx(0):iflx(1)) = poly_fit(str.pixel(w),smooth(str.flux(w)/mdl(w)*cont(w),smfact),str.cont_order,/double)
nirspecfit_model, str, param, mdl, trs, cont, px, flag=flag


; INTIALIZE SCALE FACTOR FOR TRANSMISSION
SETSCALE: trsn = trs^(1./param(itrans))
fact = param(itrans)+findgen(nsamp+1)/nsamp-0.2
w = where(fact gt 0.1,cnt)
fact = fact(w)
chis = fact*0.+1.e9
for i=0L,cnt-1L do chis(i) = total(str.mask*(str.flux-mdl*trsn^(fact(i)))^2/str.noise^2)
mn = min(chis,loc)
minchi = chis(loc)
param(itrans) = fact(loc)
nirspecfit_model, str, param, mdl, trs, cont, px, flag=flag

if (keyword_set(debug)) then begin
 window, 0
 !p.color=0
 !p.background=255
 plot, str.pixel, str.flux, /xsty, /ysty
 oplot, str.pixel, mdl-0.3*median(str.flux), color=2
 oplot, [0,0]+str.prng(0), [0,2.*max(str.flux)], linestyle=1
 oplot, [0,0]+str.prng(1), [0,2.*max(str.flux)], linestyle=1
endif


; GO THROUGH MCMC FOR ALL PARAMETERS


str.pvar(0:str.wave_order) = 0.		; don't change wavelength scale
;for kkk=0,str.lamc_order do str.pvar(kkk) = 3./(median(str.wave_model)^kkk)
;str.pvar(str.lamc_order+1:str.lamc_order+str.flxc_order+1) = 0.		; don't change flux scale
str.pvar(n_elements(param)-2:n_elements(param)-1) = 0.	; don't change RV, vsini
str.param = param
;str.nchain = shortchain

param0 = param

; run through MCMC routine
nirspecfit_fitspec, str, debug=0, /calibrator

param = str.param
nirspecfit_model, str, param, mdl, trs, cont, px, flag=flag
chi2 = total(str.mask*(str.flux-mdl)^2/str.noise^2)
if (chi2 gt minchi) then param = param0 else minchi=chi2

if (keyword_set(debug)) then begin
 !p.multi=0
 plot, str.pixel, str.flux, yra=[-0.5,1.5], /ysty, /xsty, title=strtrim(string(stddev(str.mask*(str.flux-mdl))),2)
 oplot, str.pixel, mdl, color=3
 oplot, str.pixel, str.mask*(str.flux-mdl), color=150
; read, f, prompt='press return'
; wait, 1
endif

mdl0 = mdl
px0 = px
loop=loop+1

plot, str.flux
oplot, mdl, color=2


;if (loop lt nloop) then goto, SETWAVE else goto, OUTLOOP


; CHECK ON FREQUENCY STRUCTURE - INSERT WHEN THIS WORKS
diff = str.mask*(str.flux-mdl)
xdim=n_elements(diff)
nfil=xdim/2+1
freq=findgen(nfil)
xdimfft=250

diffs = median(diff,30)
difffreq = fft(diff-diffs)
yp=abs(difffreq[0:nfil-1])^2
yp=yp/max(yp)
plot,freq,yp,xrange=[0,xdimfft],xstyle=1,ystyle=1, $
    xtitle='Frequency ('+strtrim(string(xdim),2)+' pix)!E-1', $
    ytitle='Power Spectrum',xmargin=[10,3],charsize=1.5
w = where(abs(yp-median(yp)) gt stddev(yp),cnt)
freq2=findgen(nfil/2+1)/(nfil/float(xdim))
fil = (freq2 gt max(freq(w))) or (freq2 lt min(freq(w)))
;fil = (freq gt max(freq(w))) or (freq lt min(freq(w)))
fil=float(fil)
fil=[fil,reverse(fil[1:*])]
fil=float(fft(fil,/inverse))
fil=fil/nfil
fil=shift(fil,nfil/2)
fil=fil*hanning(nfil)

df1 = convol(diff-diffs,fil,/edge_wrap)+diffs
plot, diff
oplot, df1, color=2


OUTLOOP: str.param = param
save, str, file=str.outputfolder+str.fitname+'.dat'

if (keyword_set(debug)) then begin
 !p.multi=0
 nirspecfit_model, str, param, mdl, trs, cont, px, flag=flag
 plot, str.pixel, str.flux, yra=[-0.5,1.5], /ysty, /xsty, title=strtrim(string(mn/(total(str.mask)-n_elements(param))),2)
 oplot, str.pixel, mdl, color=3
 oplot, str.pixel, str.mask*(str.flux-mdl), color=150
; read, f, prompt='Press return to continue'
endif

if keyword_set(report) then nirspecfit_report, str, /eps
nirspecfit_exportspec, str, file=str.outputfolder+str.fitname+'_spectrum.txt', /calibrator, /plot

return
end
