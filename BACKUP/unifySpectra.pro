pro unifyForPyspecfit, blueSpecFile, redSpecFile

  bluData = mrdfits(blueSpecFile, 1)
  redData = mrdfits(redSpecFile, 1)

  bl       = bluData.LAMBDA
  binner   = bluData.F_INNER / bluData.SENSITIVITY
  binter   = bluData.F_INTER / bluData.SENSITIVITY
  bouter   = bluData.F_OUTER / bluData.SENSITIVITY
  binn_err = bluData.VAR_INNER / bluData.SENSITIVITY^2
  bint_err = bluData.VAR_INTER / bluData.SENSITIVITY^2
  bout_err = bluData.VAR_OUTER / bluData.SENSITIVITY^2
  
  rl       = redData.LAMBDA
  rinner   = redData.F_INNER / redData.SENSITIVITY
  rinter   = redData.F_INTER / redData.SENSITIVITY
  router   = redData.F_OUTER / redData.SENSITIVITY
  rinn_err = redData.VAR_INNER / redData.SENSITIVITY^2
  rint_err = redData.VAR_INTER / redData.SENSITIVITY^2
  rout_err = redData.VAR_OUTER / redData.SENSITIVITY^2

  dl = 2.
  lrun   = ceil((max(rl) - min(bl)) / dl)
  lambda = findgen(lrun) * dl + min(bl)
  bs     = where(lambda le max(bl))
  rs     = where(lambda ge min(rl))

  binner   = interpol(binner  , bl, lambda)
  binter   = interpol(binter  , bl, lambda)
  bouter   = interpol(bouter  , bl, lambda)
  binn_err = interpol(binn_err, bl, lambda)
  bint_err = interpol(bint_err, bl, lambda)
  bout_err = interpol(bout_err, bl, lambda)
  
  rinner   = interpol(rinner  , rl, lambda)
  rinter   = interpol(rinter  , rl, lambda)
  router   = interpol(router  , rl, lambda)
  rinn_err = interpol(rinn_err, rl, lambda)
  rint_err = interpol(rint_err, rl, lambda)
  rout_err = interpol(rout_err, rl, lambda)

  cut = where(lambda gt max(bl))
  binner[cut] = 0
  binter[cut] = 0
  bouter[cut] = 0
  binn_err[cut] = 1d6
  bint_err[cut] = 1d6
  bout_err[cut] = 1d6

  cut = where(lambda lt min(rl))
  rinner[cut] = 0
  rinter[cut] = 0
  router[cut] = 0
  rinn_err[cut] = 1d6
  rint_err[cut] = 1d6
  rout_err[cut] = 1d6

  inner     = dblarr(lrun)
  inner_var = dblarr(lrun)
  inter     = dblarr(lrun)
  inter_var = dblarr(lrun)
  outer     = dblarr(lrun)
  outer_var = dblarr(lrun)
  for ii = 0, lrun - 1 do begin
     inner[ii] = (binner[ii] / binn_err[ii] + rinner[ii] / rinn_err[ii]) $
                 / (1. / binn_err[ii] + 1. / rinn_err[ii])
     inner_var[ii] = 1. / (1. / binn_err[ii] + 1. / rinn_err[ii])
     inter[ii] = (binter[ii] / bint_err[ii] + rinter[ii] / rint_err[ii]) $
                 / (1. / bint_err[ii] + 1. / rint_err[ii])
     inter_var[ii] = 1. / (1. / bint_err[ii] + 1. / rint_err[ii])
     outer[ii] = (bouter[ii] / bout_err[ii] + router[ii] / rout_err[ii]) $
                 / (1. / bout_err[ii] + 1. / rout_err[ii])
     outer_var[ii] = 1. / (1. / bout_err[ii] + 1. / rout_err[ii])
  endfor
  
  stop
  
end
; unifyForPyspecfit, '753_1_B.fits', '753_1_R.fits'
