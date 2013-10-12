; NIRSPECFIT
; PLOT AN EXPORTED FILE

PRO nirspecfit_plotexport, specfile, pfile, pixel=pixel, name=name

; error checking
if (file_search(specfile) eq '') then begin
 print, 'Cannot locate '+specfile+' for plotting'
 goto, FINISH
endif

; read in data
readcol, specfile, px, wave, raw, noise, mdl, trans, cont, corrected, format='i,f,f,f,f,f,f,f', /silent, comment='#'

; set pfile if needed
if (n_elements(pfile) eq 0) then pfile = strrep(specfile,'.txt','.eps')
if (n_elements(name) eq 0) then name = ''

; plot
set_plot, 'ps'
device, /encapsulated, ysize=24, xsize=28, filename=pfile, /portrait, bits_per_pixel=8, /color
loadct, 0, /silent
tek_color

;yra = [min([raw-noise])<0,max([raw,2.*noise])*1.2+0.5]
yra = [-0.05,2.2]
scl = median(raw)/median(corrected)
!p.multi=0

if (keyword_set(pixel)) then begin
 x = px
 xtit = 'Pixels'
 ym = [6,6]
endif else begin
 x = wave
 xtit = 'Wavelength (!9m!3m)'
 ym = [6,4]
endelse

y1 = corrected/median(corrected)
y2 = raw/median(raw)
y3 = noise/median(raw)
plot, x, y1, yra=yra, /ysty, /xsty, xtitle=xtit, ytitle='!3Normalized Flux', charsize=2, ymargin=ym; , title=name
 oplot, x, y2, color=150
 oplot, x, y3, color=0
 oplot, x, trans+1., color=0, linestyle=1
 oplot, x, x*0., linestyle=1
 w = where(x ge max(x)-0.25*(max(x)-min(x)))
 xyouts, max(x)-0.1*(max(x)-min(x)), 2.1, ' Transmission + 1 ', align=0.5, charsize=1.2, color=0
 xyouts, max(x)-0.1*(max(x)-min(x)), max([y1(w),y2(w)])+0.1, ' Raw ', align=0, charsize=1.2, color=150
 xyouts, max(x)-0.1*(max(x)-min(x)), max([y1(w),y2(w)])+0.1, ' Corrected ', align=1, charsize=1.2, color=0
 w = where(x le min(x)+0.25*(max(x)-min(x)))
 xyouts, min(x)+0.1*(max(x)-min(x)), max(y3(w))+0.1, ' Noise ', align=0.5, charsize=1.2, color=0

if (keyword_set(pixel)) then begin
 for j=ceil(min(px)/150.)*150.,floor(max(px)/150.)*150.,150 do $
  xyouts, j, yra(1)+0.02*(yra(1)-yra(0)), strtrim(string(float(interpol(wave,px,j,/quad))),2), align=0.5, charsize=2
 xyouts, mean(px), yra(1)+0.1*(yra(1)-yra(0)), 'Wavelength (!9m!3m)', charsize=2, align=0.5
endif

device, /close
set_plot, 'x'

FINISH: return
end
