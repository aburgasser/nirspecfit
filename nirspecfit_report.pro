; NIRSPECFIT
; REPORTING ROUTINE FOR INTERMEDIATE DATA PRODUCTS

PRO NIRSPECFIT_REPORT, str, eps=eps

tek_color
!p.font=0
!p.multi=0
!p.thick=4
!x.thick=3
!y.thick=3
tb='	'
trix = [-1,0,1,-1]
triy = [-1,1,-1,-1]
crx = rndoff(10.*cos(findgen(41)*!pi/20.))
cry = rndoff(10.*sin(findgen(41)*!pi/20.))
sqx = [-1,-1,1,1,-1]
sqy = [-1,1,1,-1,-1]

wsample = where(str.pixel ge str.prng(0) and str.pixel le str.prng(1))

if (keyword_set(eps)) then begin
    set_plot, 'ps'
    device, /encapsulated, ysize=24, xsize=28, filename=str.outputfolder+str.fitname+'_fit.eps', /portrait, bits_per_pixel=8, /color
endif else plotps, str.outputfolder+str.fitname+'_fit.ps', /open, /color


; PLOT SPECTRUM FIT WITH MODEL AND TRANSMISSION
 nirspecfit_model, str, str.param, mdl, trs, cont, px, flag=flag
 diff = str.mask*(str.flux-mdl)
 chi = total(diff^2/str.noise^2)/str.dof

 yra = [0.-2.*(max([str.noise(wsample),0.-diff(wsample)])),max([str.flux(wsample),mdl(wsample)]+3.*str.noise(wsample))]
 !p.multi=0
plot, str.pixel(wsample), str.mask(wsample)*str.flux(wsample), yra=yra, /ysty, /xsty, xtitle='Pixels', ytitle='!3Normalized Flux', charsize=2, ymargin=[8,6], xmargin=[8,6]
 oplot, str.pixel(wsample), mdl(wsample), color=3
 oplot, str.pixel(wsample), str.noise(wsample), color=0
 oplot, str.pixel(wsample), 0.-str.noise(wsample), color=0
 oplot, str.pixel(wsample), mdl(wsample)/trs(wsample), color=2
 for j=ceil(min(str.prng)/150.)*150.,floor(max(str.prng)/150.)*150.,150 do $
  xyouts, j, yra(1)+0.02*(yra(1)-yra(0)), strtrim(string(float(interpol(str.wave_model,px,j,/quad))),2), align=0.5, charsize=1.5
 xyouts, mean(str.pixel(wsample)), yra(1)+0.1*(yra(1)-yra(0)), 'Wavelength (!9m!3m)', charsize=1.5, align=0.5
 oplot, str.pixel(wsample), diff(wsample), color=150
  xyouts, max(str.pixel(wsample))-0.05*(max(str.pixel(wsample))-min(str.pixel(wsample))), max(diff(wsample)), 'O-C', align=1, charsize=1.5
  oplot, str.pixel, str.pixel*0., linestyle=1
 for j=0,n_elements(str.param)-1 do xyouts, max(str.pixel(wsample))+0.02*(max(str.pixel(wsample))-min(str.pixel(wsample))), yra(1)-(0.8*(yra(1)-yra(0))/n_elements(str.param))*(n_elements(str.param)-j), strtrim(string(str.param(j)),2), align=0, charsize=0.8

 xyouts, max(str.pixel(wsample))-0.05*(max(str.pixel(wsample))-min(str.pixel(wsample))), yra(1)-0.1*(yra(1)-yra(0)), '!9c!3!U2!N = '+strtrim(string(float(chi)),2), align=1, charsize=1.5

if (max(strpos(tag_names(str),'PARAM_MCMC')) ge 0) then begin

if (keyword_set(eps)) then begin
   device, /close
    device, /encapsulated, ysize=24, xsize=28, filename=str.outputfolder+str.fitname+'_pdist.eps', /portrait, bits_per_pixel=8, /color
endif 
 
 
 
; PLOT DISTRIBUTIONS OF PARAMETERS
  !p.multi=[0,ceil(n_elements(str.param)/fix(sqrt(n_elements(str.param)))),ceil(sqrt(n_elements(str.param)))]
  chi2 = reform(str.param_mcmc(n_elements(str.param),*))
  prob = exp(-0.5*(chi2-min(chi2)))
  wp = where(prob ge 1.e-2)
  for i=0,n_elements(str.param)-1 do begin
   bin = ((max(str.param_mcmc(i,wp))-min(str.param_mcmc(i,wp)))/(1.*str.nbin))>1.e-3
   rng = [-1.,1.]*0.5*str.nbin*bin+min(str.param_mcmc(i,wp))
   x = dindgen(str.nbin)*bin+rng(0)
   y = x*0.
   for j=0,str.nbin-1 do begin
    wb = where(str.param_mcmc(i,wp) ge x(j)-bin/2. and str.param_mcmc(i,wp) lt x(j)+bin/2. and finite(str.param_mcmc(i,wp)),cnt)
;    if (cnt gt 0) then y(j) = total(prob(wp(wb)))
    if (cnt gt 0) then y(j) = total(wp(wb))
   endfor
   y = y/max(y)
   w = where(y ge 0.5)
   fit = gaussfit(x,y,a,nterms=4,est=[max(y,loc),x(loc),0.5*(max(x(w))-min(x(w))),0.01])
   plot, x, y, psym=10, yra=[0,1.2], /ysty, xra = rng+[-1.,1.]*2*bin, /xsty, charsize=2, xtitle=str.pname(i), ytitle='Normalized Distribution'
   oplot, x, fit, color=2
   xyouts, rng(1)-bin, 1.1, a(1), align=1
   xyouts, rng(1)-bin, 1.0, a(2), align=1
  endfor
  !p.multi=0

 endif

if (keyword_set(eps)) then  device, /close else plotps, /close



; REPORT PARAMETERS
openw, unit, str.outputfolder+str.fitname+'_fit.txt', /get_lun
printf, unit, '# Best-fit NIRSPECFIT parameters for '+str.name
printf, unit, '# Observation date: '+str.obsdate+' (Julian Date '+strtrim(string(str.jdate),2)+')'
if (str.vbary_flag eq 1) then printf, unit, '# The following barycentric correction should be added to RV: '+strtrim(string(str.vbary),2)
printf, unit, '# Analysis completed on '+systime()
printf, unit, '# PARAMETER	BEST VALUE	OPTIMAL VALUE	OPTIMAL UNCERTAINTY	NOISE SAMPLING VALUE	NOISE SAMPLING FIT UNCERTAINTY
for i=0,n_elements(str.param)-1 do begin
 f = str.pname(i)
 f=f+tb+strtrim(string(str.param(i)),2)+tb+strtrim(string(str.param_best(i)),2)+tb+strtrim(string(str.param_opt(i)),2)+tb+strtrim(string(str.param_opt_unc(i)),2)
 if (max(strpos(tag_names(str),'PARAM_NOISE')) ge 0) then f=f+tb+strtrim(string(str.param_noise(i)),2)+tb+strtrim(string(str.param_noise_unc(i)),2)
 printf, unit, f
endfor

; compute chi2
printf, unit, 'chi2'+tb+strtrim(string(str.minchi),2)
printf, unit, 'dof'+tb+strtrim(string(str.dof),2)

close, unit
free_lun, unit



return
end
