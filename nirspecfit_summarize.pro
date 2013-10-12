; NIRSPECFIT
; REPORTING ROUTINE FOR FINAL RESULTS

PRO NIRSPECFIT_SUMMARIZE, str

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


; PLOT CHI2 AS A FUNCTION OF TEFF WITH LOG G  SEPARATED OUT
rv = str.param_mcmc_all(n_elements(str.param)-2,*)
vsini = str.param_mcmc_all(n_elements(str.param)-1,*)
teff = str.param_mcmc_all(n_elements(str.param)+1,*)
logg = str.param_mcmc_all(n_elements(str.param)+2,*)
chisq = str.chisq_mcmc_all
bchisq = str.best_chisq
brv = str.best_params(n_elements(str.param)-2,*)
bvsini = str.best_params(n_elements(str.param)-1,*)
bteff = str.best_params(n_elements(str.param)+1,*)
blogg = str.best_params(n_elements(str.param)+2,*)

set_plot, 'ps'
device, /encapsulated, ysize=24, xsize=28, filename=str.outputfolder+'results_chisq.eps', /portrait, bits_per_pixel=8, /color

loadct, 3, /silent
s = sort(blogg)
u = uniq(blogg(s))
bloggref = blogg(s(u))
lincol = 200./n_elements(bloggref)
xrng = [min(teff),max(teff)]+[-1.,1.]*50.
yrng = [min(bchisq)*0.8,max(bchisq)*1.2]
!p.multi=0

plot, [0],[10], xra=xrng, yra=yrng, xtitle='!3T!Deff!N (K)', ytitle='!9c!3!U2!Dr!N', charsize=2, xsty=1, ysty=1, xmargin=[8,8], yticks=5
;usersym, crx, cry, /fill
;oplot, teff, chisq, psym=8, symsize=0.02
for j=0,n_elements(bloggref)-1 do begin
 w = where(blogg eq bloggref(j))
 s = sort(bteff(w))
 oplot, bteff(w(s)), bchisq(w(s)), color=lincol*j
 xyouts, xrng(1)-0.15*(xrng(1)-xrng(0)), yrng(1)-(0.1+j*0.07)*(yrng(1)-yrng(0)), 'log g = '+strmid(strtrim(string(bloggref(j)),2),0,3), color=lincol*j, charsize=1.2, align=0.
 oplot, xrng(1)-[0.18,0.25]*(xrng(1)-xrng(0)), [0,0]+yrng(1)-(0.11+j*0.07)*(yrng(1)-yrng(0)), color=lincol*j
endfor
usersym, crx, cry, /fill
for i=0,n_elements(bteff)-1 do oplot, [bteff(i)],[bchisq(i)],psym=8,symsize=bvsini(i)*0.2/max(bvsini)
vy = (findgen(5)+1)*max(bvsini)/5.
for i=0,4 do begin
 oplot, [xrng(1)-0.18*(xrng(1)-xrng(0))],[yrng(1)-(0.1+n_elements(bloggref)*0.07+i*0.05)*(yrng(1)-yrng(0))], psym=8,symsize=vy(i)*0.2/max(bvsini)
 xyouts, xrng(1)-0.15*(xrng(1)-xrng(0)), yrng(1)-(0.11+n_elements(bloggref)*0.07+i*0.05)*(yrng(1)-yrng(0)), strmid(strtrim(string(rndoff(vy(i))),2),0,3)+' km/s', charsize=1.2, align=0
endfor

oplot, [0,1.e4],[1.,1.]*3., linestyle=1
oplot, [0,1.e4],[1.,1.]*5., linestyle=1
oplot, [0,1.e4],[1.,1.]*10., linestyle=1

device, /close
set_plot, 'x'


; PLOT BEST FIT MODEL
set_plot, 'ps'
device, /encapsulated, ysize=24, xsize=28, filename=str.outputfolder+'results_bestspec.eps', /portrait, bits_per_pixel=8, /color
loadct, 0, /silent
tek_color

mn = min(bchisq,loc)
param = str.best_params(0:n_elements(str.param)-1,loc)
wsample = where(str.pixel ge str.prng(0) and str.pixel le str.prng(1),cntsample)
nirspecfit_model, str, param, mdl, trs, cont, px, flag=flag
diff = str.mask*(str.flux-mdl)
chi = total(diff^2/str.noise^2)/str.dof
print, bchisq(loc), chi

yra = [0.-2.*(max([str.noise(wsample),0.-diff(wsample)])),max([str.flux(wsample),mdl(wsample)]+3.*str.noise(wsample))*1.2]
!p.multi=0

plot, str.pixel(wsample), str.mask(wsample)*str.flux(wsample), yra=yra, /ysty, /xsty, xtitle='Pixels', ytitle='!3Normalized Flux', charsize=2, ymargin=[8,6], xmargin=[8,6]
 oplot, str.pixel(wsample), mdl(wsample), color=3
 oplot, str.pixel(wsample), str.noise(wsample), color=0
 oplot, str.pixel(wsample), 0.-str.noise(wsample), color=0
 oplot, str.pixel(wsample), mdl(wsample)/trs(wsample), color=2
for j=ceil(min(str.prng)/100.)*100.,floor(max(str.prng)/100.)*100.,100 do $
  xyouts, j, yra(1)+0.02*(yra(1)-yra(0)), strtrim(string(float(interpol(str.wave_model,px,j,/quad))),2), align=0.5, charsize=2
xyouts, mean(str.pixel(wsample)), yra(1)+0.1*(yra(1)-yra(0)), 'Wavelength (!9m!3m)', charsize=2, align=0.5
oplot, str.pixel(wsample), diff(wsample), color=150
xyouts, max(str.pixel(wsample))-0.05*(max(str.pixel(wsample))-min(str.pixel(wsample))), max(diff(wsample)), 'O-C', align=1, charsize=1.5
oplot, str.pixel, str.pixel*0., linestyle=1

 xyouts, max(str.pixel(wsample))-0.05*(max(str.pixel(wsample))-min(str.pixel(wsample))), yra(1)-0.1*(yra(1)-yra(0)), str.name, align=1, charsize=1.5
 xyouts, max(str.pixel(wsample))-0.05*(max(str.pixel(wsample))-min(str.pixel(wsample))), yra(1)-0.15*(yra(1)-yra(0)), $
   strtrim(string(fix(bteff(loc))),2)+'/'+$
   strmid(strtrim(string(blogg(loc)),2),0,3)+'/'+$
   strmid(strtrim(string(brv(loc)+str.vbary),2),0,5)+'/'+$
   strtrim(string(fix(bvsini(loc))),2)+'.0', $
   align=1, charsize=1.5, color=2
 xyouts, max(str.pixel(wsample))-0.05*(max(str.pixel(wsample))-min(str.pixel(wsample))), yra(1)-0.2*(yra(1)-yra(0)), '!9c!3!U2!Dr!N = '+strtrim(string(float(bchisq(loc))),2), align=1, charsize=1.5

device, /close
set_plot, 'x'


; REPORT BEST MODEL FIT  
rpnames = ['Teff (K)','log g (cgs)','[M/H]','Heliocentric RV (km/s)','v sin i (km/s)']
wpnames = ['T$_{eff}$ (K)','log~g (cm~s$^{-2}$)','[M/H]','Heliocentric~RV (km~s$^{-1}$)\tablenotemark{a}','v$\sin{i}$ (km~s$^{-1}$)']
rorder = n_elements(str.param)+[1,2,3,-2,-1]

best_parameters = str.param(rorder)
optimal_parameters = str.param(rorder)
optimal_uncertainties = str.param(rorder)
for i=0,n_elements(rorder)-1 do begin
 nirspecfit_bestparam, reform(str.param_mcmc_all(rorder(i),*)), str.chisq_mcmc_all, str.dof, best, opt, opt_e, prob
 best_parameters(i) = best
 optimal_parameters(i) = opt
 optimal_uncertainties(i) = opt_e
endfor

w = where(strpos(rpnames,'RV') gt -1)
best_parameters(w(0)) = best_parameters(w(0))+str.vbary
optimal_parameters(w(0)) = optimal_parameters(w(0))+str.vbary

print, ''
print, ''
print, 'Best fit model: '
for i=0,n_elements(rorder)-1 do print, rpnames(i)+' = '+strtrim(string(best_parameters(i)),2)
print, 'Reduced chi2 = '+strtrim(string(min(str.chisq_mcmc_all)),2)
print, ''
print, ''
print, 'Optimal parameters: '
for i=0,n_elements(rorder)-1 do print, rpnames(i)+' = '+strtrim(string(optimal_parameters(i)),2)+'+/-'+strtrim(string(optimal_uncertainties(i)),2)
print, ''
print, ''

; PRINT OFF TABLE OF RESULTS
openw, unit, str.outputfolder+'results_table.tex', /get_lun
printf, unit, '\begin{deluxetable}{lcc}
printf, unit, '\tabletypesize{\small}
printf, unit, '\tablecaption{NIRSPEC Order '+strtrim(string(fix(str.order)),2)+' Fit Parameters for '+str.name+' \label{tab:nirspecfit}}
printf, unit, '\tablewidth{0pt}
printf, unit, '\tablehead{
printf, unit, '\colhead{Parameter} &
printf, unit, '\colhead{Best Fit Model} & 
printf, unit, '\colhead{Weighted Means} \\
printf, unit, '}
printf, unit, '\startdata
for i=0,n_elements(rorder)-1 do printf, unit, $
 '{'+wpnames(i)+'} & '+$
 strtrim(string(float(best_parameters(i,0))),2)+' & '+$
 strtrim(string(float(optimal_parameters(i))),2)+'$\pm$'+strtrim(string(float(optimal_uncertainties(i))),2)+' \\'
printf, unit, '$\chi^2_r$ & '+strtrim(string(min(float(str.chisq_mcmc_all))),2)+' & \nodata \\'
printf, unit, '\enddata
printf, unit, '\tablenotetext{a}{Includes barycentric motion correction of '+strtrim(string(str.vbary),2)+'~km~s$^{-1}$' 
printf, unit, '\end{deluxetable}
close, unit
free_lun, unit


; ANALYZE PARAMETER DISTRIBUTIONS

for iii=0,3 do begin
 case 1 of
  iii eq 0: begin
   fname = 'teff'
   pname = 'T!Deff!N (K)'
   p = reform(str.param_mcmc_all(n_elements(str.param)+1,*))
   bin = 100.
  end
  iii eq 1: begin
   fname = 'logg'
   pname = 'log g (cm/s!u2!n)'
   p = reform(str.param_mcmc_all(n_elements(str.param)+2,*))
   bin = 0.5
  end
  iii eq 2: begin
   fname = 'rv'
   pname = 'Heliocentric RV (km/s)'
   p = reform(str.param_mcmc_all(n_elements(str.param)-2,*))+str.vbary
   bin = 0.5
  end
 else: begin
   fname = 'vsini'
   pname = 'V sin i (km/s)'
   p = reform(str.param_mcmc_all(n_elements(str.param)-1,*))
   bin = 2
  end
 endcase
 
 rng = [min(p),max(p)]+[-1.,1.]*bin
 x = findgen((max(rng)-min(rng))/bin+1)*bin+min(rng)
 y = x*0.
 for j=0,n_elements(x)-1 do begin
  wb = where(p ge x(j)-bin/2. and p lt x(j)+bin/2.,cnt)
  if (cnt gt 0) then y(j) = total(prob(wb))
 endfor
 y = y/max(y)
 nirspecfit_bestparam, p, str.chisq_mcmc_all, str.dof, bp, op, op_e
 
;   w = where(y ge 0.5)
;   fit = gaussfit(x,y,a,nterms=4,est=[max(y,loc),x(loc),0.5*(max(x(w))-min(x(w))),0.01])

 set_plot, 'ps'
 device, /encapsulated, ysize=24, xsize=24, filename=str.outputfolder+'results_distribution_'+fname+'.eps', /portrait, bits_per_pixel=8, /color
 loadct, 0, /silent
 tek_color
 yrng = [0.1,1.3]
  
 plot, x, y, psym=10, yra=yrng, /ysty, xra = rng, /xsty, charsize=2, xtitle=pname, ytitle='Normalized Distribution'
 oplot, [0,0]+bp, yrng, linestyle=2
 oplot, [0,0]+op, yrng, linestyle=1
 xyouts, rng(0)+0.05*(rng(1)-rng(0)), 1.22, 'Best Fit Value: '+strtrim(string(float(bp)),2), align=0, charsize=1.5
 xyouts, rng(0)+0.05*(rng(1)-rng(0)), 1.15, 'Weighted Value: '+strtrim(string(float(op)),2)+'+/-'+strtrim(string(float(op_e)),2), align=0, charsize=1.5

;   oplot, x, fit, color=2
;   xyouts, rng(1)-bin, 1.1, a(1), align=1
;   xyouts, rng(1)-bin, 1.0, a(2), align=1


 device, /close
 set_plot, 'x'

endfor

 
FINISH: return
end
