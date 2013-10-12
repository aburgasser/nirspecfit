; NIRSPECFIT
; ROUTINE TO COMPUTE BEST (MINIMUM CHI2) AND OPTIMAL (WEIGHTED)
; PARAMETERS GIVEN SET OF CHISQ VALUES
; PARAM MUST BE A 1D ARRAY 

PRO nirspecfit_bestparam, param, chisq, dof, best, optimal, optimal_e, prob, flat=flat, reduced=reduced, ftest=ftest

; probabilities
case 1 of
 keyword_set(exp): begin
  prob = exp(-0.5*(chisq-min(chisq,loc)))		; exponential
  if (keyword_set(reduced)) then $
   prob = exp(-0.5*(chisq-min(chisq,loc))*dof)	; exponential for reduced chisq
  end
 keyword_set(flat): prob = chisq*0.+1.		; flat
 keyword_set(ftest): prob = 2.*(1.-f_pdf(chisq/min(chisq,loc),dof,dof))	; f test 
 else: prob = 2.*(1.-f_pdf(chisq/min(chisq,loc),dof,dof))	; f test (default)
endcase

best = param(loc)
optimal =  total(param*prob)/total(prob)
optimal_e =  sqrt(total(param^2*prob)/total(prob)-(optimal)^2)

return
end 
