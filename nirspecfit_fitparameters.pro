;
; NIRSPECFIT_FITPARAMETERS
;
; FUNCTION TO READ IN ALL PROGRAM PARAMETERS AND INPUT DRIVING FILE

PRO nirspecfit_fitparameters, file, pstr, trng=trng, grng=grng, cldrng=cldrng, savefile=savefile, reset=reset

common NIRSPECFIT, pfolder, calfolder, mfolder, outputfolder, $
 model_names, units, $
 nirspec_lref, nirspec_lrng, nirspec_res0, nirspec_npix, nirspec_pdisp, $
 samp_res, vsini_dv, vsini_min, vsini_max, $
 cvel


on_error, 0
tb = '	'
if (n_elements(savefile) eq 0) then savefile = 'fit_parameters.dat'
if (file_search(file) eq '') then message, 'Cannot find data input file '+file

; some presets
if (n_elements(trng) eq 0) then trng = [1500,2500]
if (n_elements(grng) eq 0) then grng = [4.5,5.5]

arr = readarr(file,/comment)
inputs = strarr(n_elements(arr),8)
for i=0,n_elements(arr)-1 do begin
 tmp = str_sep(arr(i),tb)
 for j=0,n_elements(tmp)-1 do inputs(i,j) = tmp(j)
endfor

if (file_search(outputfolder+savefile) eq '' or keyword_set(reset)) then begin
; SET UP STRUCTURE
dstr = create_struct($
; command parameters
 'name','',$
 'designation','',$
 'coordinate','',$
 'date','',$
 'jdate',0L,$
 'obsdate','',$
 'datafile','',$
 'calfile','',$
 'data_type',1, $
 'order',0, $
 'debug',0,$
 'reuse_cal',1,$
 'nchain', 2000L,$
 'chicut', 0.95,$
 'nbin', 30,$
 'teffrng', trng,$
 'loggrng', grng,$
 'cldrng', [2.0,2.0],$
 'zrng', [0.0,0.0],$
 'rv', 0.,$
 'rv_var', 5.d0,$
 'vsini', 5.,$
 'vsini_var', 5.d0,$
 'sg', 5.,$
 'sg_var', 1.,$
 'ta', 1.0,$
 'ta_var', 0.25,$
 'cont_param', [1.,0.,0.],$
 'cont_var', [0.3,3.e-4,3.e-7],$
 'lam0', 0., $
 'prng',[5,600])
 
 
; COMMANDING
w = where(strupcase(inputs(*,0)) eq 'DEBUG',cnt)
if (cnt gt 0) then dstr.debug = fix(inputs(w(0),1))

w = where(strupcase(inputs(*,0)) eq 'REUSE_CAL',cnt)
if (cnt gt 0) then dstr.reuse_cal = fix(inputs(w(0),1))


; DATA PARAMETERS
w = where(strupcase(inputs(*,0)) eq 'DATA_TYPE',cnt)
if (cnt gt 0) then dstr.data_type = fix(inputs(w(0),1))

w = where(strupcase(inputs(*,0)) eq 'DATAFILE',cnt)
if (cnt gt 0) then dstr.datafile=inputs(w(0),1)

w = where(strupcase(inputs(*,0)) eq 'CALFILE',cnt)
if (cnt gt 0) then dstr.calfile=inputs(w(0),1)

w = where(strupcase(inputs(*,0)) eq 'ORDER',cnt)
if (cnt gt 0) then dstr.order = fix(inputs(w(0),1))

dstr.lam0 = nirspec_lref/dstr.order

; FITTING PARAMETERS
w = where(strupcase(inputs(*,0)) eq 'TEFFRNG',cnt)
if (cnt gt 0) then begin
 dstr.teffrng(*) = fix(inputs(w(0),1))
 if (inputs(w(0),2) ne '') then dstr.teffrng(1) = fix(inputs(w(0),2))
endif

w = where(strupcase(inputs(*,0)) eq 'LOGGRNG',cnt)
if (cnt gt 0) then begin
 dstr.loggrng(*) = float(inputs(w(0),1))
 if (inputs(w(0),2) ne '') then dstr.loggrng(1) = float(inputs(w(0),2))
endif

w = where(strupcase(inputs(*,0)) eq 'ZRNG',cnt)
if (cnt gt 0) then begin
 dstr.zrng(*) = float(inputs(w(0),1))
 if (inputs(w(0),2) ne '') then dstr.zrng(1) = float(inputs(w(0),2))
endif


; INITIALIZATION PARAMETERS
w = where(strupcase(inputs(*,0)) eq 'RV_INIT',cnt)
if (cnt gt 0) then dstr.rv = float(inputs(w(0),1))

w = where(strupcase(inputs(*,0)) eq 'VSINI_INIT',cnt)
if (cnt gt 0) then dstr.vsini = float(inputs(w(0),1))

w = where(strupcase(inputs(*,0)) eq 'BROAD_INIT',cnt)
if (cnt gt 0) then dstr.sg = float(inputs(w(0),1))

w = where(strupcase(inputs(*,0)) eq 'ALPHA_INIT',cnt)
if (cnt gt 0) then dstr.ta = float(inputs(w(0),1))

w = where(strupcase(inputs(*,0)) eq 'PIXEL_RANGE',cnt)
if (cnt gt 0) then begin
 if (inputs(w(0),1) ne '' and inputs(w(0),2) ne '') then dstr.prng = long(inputs(w(0),1:2))
endif


; MCMC FITTING PARAMETERS
w = where(strupcase(inputs(*,0)) eq 'NCHAIN',cnt)
if (cnt gt 0) then dstr.nchain = long(inputs(w(0),1))

w = where(strupcase(inputs(*,0)) eq 'CHICUT',cnt)
if (cnt gt 0) then dstr.chicut = float(inputs(w(0),1))

w = where(strupcase(inputs(*,0)) eq 'RV_VAR',cnt)
if (cnt gt 0) then dstr.rv_var = double(inputs(w(0),1))

w = where(strupcase(inputs(*,0)) eq 'VSINI_VAR',cnt)
if (cnt gt 0) then dstr.vsini_var = double(inputs(w(0),1))

w = where(strupcase(inputs(*,0)) eq 'BROAD_VAR',cnt)
if (cnt gt 0) then dstr.sg_var = double(inputs(w(0),1))

w = where(strupcase(inputs(*,0)) eq 'ALPHA_VAR',cnt)
if (cnt gt 0) then dstr.ta_var = double(inputs(w(0),1))

w = where(strupcase(inputs(*,0)) eq 'WAVE_PARAM',cnt)
if (cnt gt 0) then dstr = create_struct(dstr,'wave_param',double(inputs(w(0),1:3))) else $
 dstr = create_struct(dstr,'wave_param',dstr.lam0*[1.,nirspec_pdisp/nirspec_res0, 0.5*(nirspec_pdisp/nirspec_res0)^2])

w = where(strupcase(inputs(*,0)) eq 'WAVE_VAR',cnt)
if (cnt gt 0) then dstr = create_struct(dstr,'wave_var',double(inputs(w(0),1:3))) else $
 dstr = create_struct(dstr,'wave_var',dstr.wave_param/10.)


; PARAMETER FITTING PARAMETERS
w = where(strupcase(inputs(*,0)) eq 'NBIN',cnt)
if (cnt gt 0) then dstr.nbin = long(inputs(w(0),1))

dstr = create_struct(dstr,$
 'wave_order',n_elements(dstr.wave_param)-1,$
 'cont_order',n_elements(dstr.cont_param)-1)

; GENERATE FINAL SET OF PARAMETER ARRAYS
dstr = create_struct(dstr,$
 'param',[dstr.wave_param,dstr.cont_param,dstr.sg,dstr.ta,dstr.rv,dstr.vsini],$
 'pvar',[dstr.wave_var,dstr.cont_var,dstr.sg_var,dstr.ta_var, dstr.rv_var,dstr.vsini_var],$
 'pname',['L'+strtrim(sindgen(dstr.wave_order+1),2),'F'+strtrim(sindgen(dstr.cont_order+1),2),'SG','!9a!3','RV','V sin i'])

dstr = create_struct(dstr,$
 'param_best', dstr.param,$
 'param_opt', dstr.param,$
 'param_opt_unc', dstr.param,$
 'minchi', 1.e9)

pstr = dstr		; rename structure

save, pstr, file=outputfolder+savefile
endif else restore, file=outputfolder+savefile

FINISH: return

end
