; NIRSPECFIT_GENERATEMODELSET
; RUN SPARINGLY TO GENERATE MODEL SETS

PRO nirspecfit_generatemodels, orders=orders

common NIRSPECFIT, pfolder, calfolder, mfolder, outputfolder, $
 model_names, units, $
 nirspec_lref, nirspec_lrng, nirspec_res0, nirspec_npix, nirspec_pdisp, $
 samp_res, vsini_dv, vsini_min, vsini_max, $
 cvel

nirspecfit_initialize

tek_color

; solar atlas files
tfold = calfolder+'solatlas/'			; solar atlas from Livinsgston & Wallace 1991

; flags for what models to create
flag_settled08 = 0
flag_saumon13 = 1

; vsini and resolution parameters
model_vsini = findgen(ceil((vsini_max-vsini_min)/vsini_dv)+1)*vsini_dv+vsini_min


; ORDERS AND WAVELENGTH RANGES TO GENERATE
if (n_elements(orders) eq 0) then begin 		; generate all orders in the NIR
 m0 = floor(min(nirspec_lref/nirspec_lrng))
 m1 = ceil(max(nirspec_lref/nirspec_lrng))
 orders = indgen(m1-m0+1)+m0
endif else orders = [orders]

lamc = nirspec_lref/orders
orders_wave = fltarr(n_elements(orders),2)
for i=0,n_elements(orders)-1 do orders_wave(i,*) = nirspec_lref/orders(i)*(1.+[-1.,1.]*0.5*nirspec_npix/nirspec_res0)


; LARGE LOOP THROUGH ORDERS
for ooo=0,n_elements(orders)-1 do begin
 print, ' Order '+strtrim(string(orders(ooo)),2)

 dl = total(orders_wave(ooo,*))/2./samp_res
 wave_model = findgen(ceil((orders_wave(ooo,1)-orders_wave(ooo,0))/dl)+1)*dl+orders_wave(ooo,0)

; READ IN TRANSMISSION DATA
; read in atmospheric transmission from atlas
 cmrng = 1.e4/reverse(reform(orders_wave(ooo,*)))
 filerng = fix(cmrng/25.)*25
 filerng = filerng(uniq(filerng))
 files = tfold+'wn'+strtrim(string(indgen((filerng(1)-filerng(0))/25+1)*25+filerng(0)),2)
 for i=0,n_elements(files)-1 do begin
  file = file_search(files(i))
  if (file ne '') then begin
   readcol, files(i),cm,sol,atm, format='f,f,f', /silent
   if (n_elements(wave) eq 0) then begin
    wave = 1.e4/cm		; in microns
    trans = atm
   endif else begin
    wave = [wave,1.e4/cm]
    trans = [trans,atm]
   endelse
  endif
 endfor
 if (n_elements(wave) eq 0) then begin
  print, 'No transmission data - assuming unity'
  trans = wave_model*0.+1.
 endif else begin
  s = sort(wave)
  u = uniq(wave(s))
  trans = interpol(trans(s(u)),wave(s(u)),wave_model,/quad)
 endelse


; GENERATE BT-SETTLED (ALLARD ET AL. 2008)
BTSET08: if (flag_settled08 eq 1) then begin
 print, 'Generating models for BT Settled (2008)'
 modelname = model_names(0)
 modelfold = mfolder+'/allard/agss2009/'
 savefile = calfolder+'models_'+modelname+'_order'+strtrim(string(orders(ooo)),2)+'.dat'

; check for presence of models
 chk = file_search(modelfold)
 if (chk eq '') then begin
  print, 'Cannot find model folder '+modelfold+' for '+modelname+' models; skipping'
  goto, SAUMON13 
 endif

 teff = indgen(11)*100+1500
 logg = [4.0,4.5,5.0,5.5]
 z = [0.0]
 cloud = [2.]
 models = dblarr($
  n_elements(wave_model),$
  n_elements(teff),$
  n_elements(logg),$
  n_elements(z),$
  n_elements(cloud),$
  n_elements(model_vsini))*0.

 for iii=0,n_elements(teff)-1 do begin
  for jjj=0,n_elements(logg)-1 do begin
   for kkk=0,n_elements(z)-1 do begin
    for lll=0,n_elements(cloud)-1 do begin

     print, '  '+strtrim(string(teff(iii)),2)+'  '+strtrim(string(logg(jjj)),2)+'  '+strtrim(string(z(kkk)),2)+'  '+strtrim(string(cloud(lll)),2)
     mfile = modelfold+'lte'+strmid(strtrim(string(teff(iii)/100000.),2),2,3)+'-'+strmid(strtrim(string(logg(jjj)),2),0,3)+'-'+strmid(strtrim(string(z(kkk)),2),0,3)+'.BT-Settl.7'
     chk = file_search(mfile+'.gz')
     if (chk ne '') then spawn, 'gunzip '+mfile+'.gz'

     basemodel = wave_model*0.
     chk = file_search(mfile)
     if (chk ne '') then begin
      readcol, mfile, wave, flux, /silent
      spawn, 'gzip '+mfile
      wave = wave/1.e4		; convert Angstroms to microns
      s = sort(wave)
      wave = wave(s)
      flam = 1.e4*10.^(flux(s)-8.d0)		; this is erg/cm2/s/micron at star surface

; fix to keep allow subsampled model

      w = where(wave_model ge min(wave) and wave_model le max(wave),cnt)
      if (cnt gt 0) then begin
       basemodel(w) = interpol(flam,wave,wave_model(w),/quad)

; error checking
       plot, wave_model, basemodel, /xsty, /ysty
       oplot, wave_model, basemodel*trans, color=2
       
; spin model up
       for vvv=0,n_elements(model_vsini)-1 do $
        models(w,iii,jjj,kkk,lll,vvv) = nirspecfit_broaden(wave_model(w),basemodel(w),model_vsini(vvv),/rotate)

      endif 
     endif
    endfor
   endfor
  endfor
 endfor

 str = create_struct($
  'modelname',modelname,$
  'order',orders(ooo),$
  'units',units,$
  'wave_model', wave_model,$
  'trans', trans, $
  'models', models, $
  'teff', teff, $
  'logg', logg, $
  'z',z,$
  'cloud',cloud,$
  'vsini',model_vsini)
  
  stop
  
 save, str, file=savefile
 if (file_search(savefile+'.gz') ne '') then spawn, 'rm '+savefile+'.gz'
 spawn, 'gzip '+savefile
endif

; GENERATE SAUMON 2013
SAUMON13: if (flag_saumon13 eq 1) then begin
 print, 'Generating models for Saumon et al (2013)'
 modelname = model_names(1)
 modelfold = mfolder+'/saumon/20130531/'

 savefile = calfolder+'models_'+modelname+'_order'+strtrim(string(orders(ooo)),2)+'.dat'

; check for presence of models
 chk = file_search(modelfold)
 if (chk eq '') then begin
  print, 'Cannot find model folder '+modelfold+' for '+modelname+' models; skipping'
  goto, ENDLOOP 
 endif

 teff = indgen(10)*100+1500
 teffstr = strtrim(string(teff),2)
 logg = [4.5,5.0,5.5]
 loggstr = ['300','1000','3000'] 
 z = [0.0]
 zstr = ['']
 cloud = [2.]
 cldstr = ['f2']

 models = dblarr($
  n_elements(wave_model),$
  n_elements(teff),$
  n_elements(logg),$
  n_elements(z),$
  n_elements(cloud),$
  n_elements(model_vsini))*0.

 for iii=0,n_elements(teff)-1 do begin
  for jjj=0,n_elements(logg)-1 do begin
   for kkk=0,n_elements(z)-1 do begin
    for lll=0,n_elements(cloud)-1 do begin

     print, '  '+teffstr(iii)+' '+loggstr(jjj)+' '+zstr(kkk)+' '+cldstr(lll)
     mfile = modelfold+'sp_t'+teffstr(iii)+'g'+loggstr(jjj)+cldstr(lll)
     chk = file_search(mfile+'.gz')
     if (chk ne '') then spawn, 'gunzip '+mfile+'.gz'

     basemodel = wave_model*0.
     chk = file_search(mfile)
     if (chk ne '') then begin
      readcol, mfile, wave, flux, skip=3, /silent
      spawn, 'gzip '+mfile
      s = sort(wave)
      wave = wave(s)
      flam = 1.e4*cvel*flux(s)/(wave(s)^2)		; this is erg/cm2/s/micron at star surface

; fix to allow subsampled model

      w = where(wave_model ge min(wave) and wave_model le max(wave),cnt)
      if (cnt gt 0) then begin
       basemodel(w) = interpol(flam,wave,wave_model(w),/quad)

; error checking
       plot, wave_model, basemodel, /xsty, /ysty
       oplot, wave_model, basemodel*trans, color=2
          
; spin model up
       for vvv=0,n_elements(model_vsini)-1 do $
        models(w,iii,jjj,kkk,lll,vvv) = nirspecfit_broaden(wave_model(w),basemodel(w),model_vsini(vvv),/rotate)
      endif
     endif
    endfor
   endfor
  endfor
 endfor

 str = create_struct($
  'modelname',modelname,$
  'order',orders(ooo),$
  'units',units,$
  'wave_model', wave_model,$
  'trans', trans, $
  'models', models, $
  'teff', teff, $
  'logg', logg, $
  'z',z,$
  'cloud',cloud,$
  'vsini',model_vsini)
  
 save, str, file=savefile
 if (file_search(savefile+'.gz') ne '') then spawn, 'rm '+savefile+'.gz'
 spawn, 'gzip '+savefile
endif


ENDLOOP: tmp=0
endfor


return
end
