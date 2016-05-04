pro plot3res, prefix, OIII = oiii

  if keyword_set(OIII) then oiii = 1 else oiii = 0

  ;; Read results
  if NOT oiii then begin
     innerRes = grabPyspecResults(prefix+'inner.analyze')
     interRes = grabPyspecResults(prefix+'inter.analyze')
     outerRes = grabPyspecResults(prefix+'outer.analyze')
  endif else begin
     innerRes = grabPyspecResults(prefix+'inner.analyze', /oiii)
     interRes = grabPyspecResults(prefix+'inter.analyze', /oiii)
     outerRes = grabPyspecResults(prefix+'outer.analyze', /oiii)
  endelse
   
  ;; Read G102
  readcol, prefix+'inner.bestspec0', $
           lambdaB, innTraceB, innEerrB, innPolyB, innModelB, innUsedB
  readcol, prefix+'inter.bestspec0', $
           lambdaB, intTraceB, intEerrB, intPolyB, intModelB, intUsedB
  readcol, prefix+'outer.bestspec0', $
           lambdaB, outTraceB, outEerrB, outPolyB, outModelB, outUsedB

  ;; Read G141
  readcol, prefix+'inner.bestspec1', $
           lambdaR, innTraceR, innEerrR, innPolyR, innModelR, innUsedR
  readcol, prefix+'inter.bestspec1', $
           lambdaR, intTraceR, intEerrR, intPolyR, intModelR, intUsedR
  readcol, prefix+'outer.bestspec1', $
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
  
  stop
  
end
