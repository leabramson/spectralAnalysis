pro plot1res, blueRes, redRes, redshift

  readcol, blueRes, $
           lambdaB, innTraceB, innEerrB, innPolyB, innModelB, innUsedB
  readcol, redRes, $
           lambdaR, innTraceR, innEerrR, innPolyR, innModelR, innUsedR

  innSpecB = innTraceB * innPolyB
  innSpecB[where(~innUsedB)] = !values.F_NAN

  innModelB *= innPolyB
  innModelB[where(~innUsedB)] = !values.F_NAN

  innSpecR = innTraceR * innPolyR
  innSpecR[where(~innUsedR)] = !values.F_NAN

  innModelR *= innPolyR
  innModelR[where(~innUsedR)] = !values.F_NAN

  lambda = [lambdaB, lambdaR] / (1 + redshift)
  inner  = [innSpecB, innSpecR]
  innMod = [innModelB, innModelR]

  plot, lambda, smooth(inner, 2), yran = [0,400]
  oplot, lambda, innMod, col = 255, thick = 1

  stop
  
end
