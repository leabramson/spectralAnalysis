pro plot3res, prefix, pa, OIII = oiii

  if keyword_set(OIII) then oiii = 1 else oiii = 0

  ;; Read results
  if NOT oiii then begin
     innerRes = grabPyspecResults(prefix+'inner'+string(pa, f='(I1)')+'.analyze')
     interRes = grabPyspecResults(prefix+'inter'+string(pa, f='(I1)')+'.analyze')
     outerRes = grabPyspecResults(prefix+'outer'+string(pa, f='(I1)')+'.analyze')
  endif else begin
     innerRes = grabPyspecResults(prefix+'inner'+string(pa, f='(I1)')+'.analyze', /oiii)
     interRes = grabPyspecResults(prefix+'inter'+string(pa, f='(I1)')+'.analyze', /oiii)
     outerRes = grabPyspecResults(prefix+'outer'+string(pa, f='(I1)')+'.analyze', /oiii)
  endelse
   
  ;; Read G102
  readcol, prefix+'inner'+string(pa, f='(I1)')+'.bestspec0', $
           lambdaB, innTraceB, innEerrB, innPolyB, innModelB, innUsedB
  readcol, prefix+'inter'+string(pa, f='(I1)')+'.bestspec0', $
           lambdaB, intTraceB, intEerrB, intPolyB, intModelB, intUsedB
  readcol, prefix+'outer'+string(pa, f='(I1)')+'.bestspec0', $
           lambdaB, outTraceB, outEerrB, outPolyB, outModelB, outUsedB

  ;; Read G141
  readcol, prefix+'inner'+string(pa, f='(I1)')+'.bestspec1', $
           lambdaR, innTraceR, innEerrR, innPolyR, innModelR, innUsedR
  readcol, prefix+'inter'+string(pa, f='(I1)')+'.bestspec1', $
           lambdaR, intTraceR, intEerrR, intPolyR, intModelR, intUsedR
  readcol, prefix+'outer'+string(pa, f='(I1)')+'.bestspec1', $
           lambdaR, outTraceR, outEerrR, outPolyR, outModelR, outUsedR

  ;; Do G102
  innSpecB = innTraceB * innPolyB
  intSpecB = intTraceB * intPolyB
  outSpecB = outTraceB * outPolyB
  innSpecB[where(~innUsedB)] = !values.F_NAN
  intSpecB[where(~intUsedB)] = !values.F_NAN
  outSpecB[where(~outUsedB)] = !values.F_NAN

  innModelB *= innPolyB
  intModelB *= intPolyB
  outModelB *= outPolyB
  innModelB[where(~innUsedB)] = !values.F_NAN
  intModelB[where(~intUsedB)] = !values.F_NAN
  outModelB[where(~outUsedB)] = !values.F_NAN

  ;; Do G141
  innSpecR = innTraceR * innPolyR
  intSpecR = intTraceR * intPolyR
  outSpecR = outTraceR * outPolyR
  innSpecR[where(~innUsedR)] = !values.F_NAN
  intSpecR[where(~intUsedR)] = !values.F_NAN
  outSpecR[where(~outUsedR)] = !values.F_NAN

  innModelR *= innPolyR
  intModelR *= intPolyR
  outModelR *= outPolyR
  innModelR[where(~innUsedR)] = !values.F_NAN
  intModelR[where(~intUsedR)] = !values.F_NAN
  outModelR[where(~outUsedR)] = !values.F_NAN

  ;; Combine
  lambda = [lambdaB, lambdaR] / (1 + innerRes.ZFIT)
  inner  = [innSpecB, innSpecR]
  inter  = [intSpecB, intSpecR]
  outer  = [outSpecB, outSpecR]

  innMod = [innModelB, innModelR]
  intMod = [intModelB, intModelR]
  outMod = [outModelB, outModelR]

  plot, lambda, inner, yran = [0,400], /nodat
  oplot, lambda, inner, col = 255
  oplot, lambda, inter, col = '00a5ff'x
  oplot, lambda, outer, col = 'ffa500'x
  oplot, lambda, innMod, col = 110
  oplot, lambda, intMod, col = '0044ff'x  
  oplot, lambda, outMod, col = 'ff0000'x

  r       = [0, mean([0.25,0.75]), mean([0.75,1.50])]
  ageGrad = [innerRes.AGEFIT, interRes.AGEFIT, outerRes.AGEFIT]
  ageHi   = [innerRes.AGEHI, interRes.AGEHI, outerRes.AGEHI]
  ageLo   = [innerRes.AGELO, interRes.AGELO, outerRes.AGELO]
  ageEHi  = 0.434 * (ageHi - ageGrad) / ageGrad
  ageELo  = 0.434 * (ageLo - ageGrad) / ageGrad * (-1)
  ageGrad = alog10(ageGrad)

  HaGrad  = [innerRes.EWFIT, interRes.EWFIT, outerRes.EWFIT]
  HaHi    = [innerRes.EWHI, interRes.EWHI, outerRes.EWHI]
  HaLo    = [innerRes.EWLO, interRes.EWLO, outerRes.EWLO]
  HaEHi   = HaHi - HaGrad
  HaELo   = HaGrad - HaLo

  OIIIGrad  = [innerRes.EWFIT_OIII, interRes.EWFIT_OIII, outerRes.EWFIT_OIII]
  OIIIHi    = [innerRes.EWHI_OIII, interRes.EWHI_OIII, outerRes.EWHI_OIII]
  OIIILo    = [innerRes.EWLO_OIII, interRes.EWLO_OIII, outerRes.EWLO_OIII]
  OIIIEHi   = OIIIHi - OIIIGrad
  OIIIELo   = OIIIGrad - OIIILo
  
  
  stop
  
end
