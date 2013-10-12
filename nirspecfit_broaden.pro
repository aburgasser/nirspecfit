; NIRSPECFIT
; FUNCTION TO BROADEN SPECTRUM FOR ROTATION OR GAUSSIAN

FUNCTION nirspecfit_broaden, wave, flux, vbroad, rotate=rotate, gaussian=gaussian

cvel = 3.e5

vres = cvel*median(abs((wave-shift(wave,1))/wave))

case 1 of 
 keyword_set(rotate): begin
  kern = lsf_rotate(vres,vbroad)
 end
 keyword_set(gaussian): begin
  x = findgen(ceil(20.*vbroad/vres)+1)
  x = (x/max(x)-0.5)*20.
  kern = exp(-0.5*x^2)
 end
 else: begin
  x = findgen(ceil(20.*vbroad/vres)+1)*10.
  kern = exp(-0.5*x^2)
 end
endcase

kern = kern/total(kern)

return, convol(flux,kern)
end
