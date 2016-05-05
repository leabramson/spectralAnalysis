function grabPyspecResults, filename, $
                            OIII = oiii

  if keyword_set(OIII) then oiii = 1 else oiii = 0

  spawn, 'cat '+filename+' | grep "z " > tmpz.dat'
  spawn, 'cat '+filename+' | grep "age" > tmpAge.dat'
  spawn, 'cat '+filename+' | grep "HpsEW" > tmpEW.dat'
  if oiii then $
     spawn, 'cat '+filename+' | grep "OIIIpsEW" > tmpEW_OIII.dat'
  
  readcol, 'tmpz.dat', zFit, zMin, zMax, f = 'X,F,F,F'
  readcol, 'tmpAge.dat', ageFit, ageMin, ageMax, f = 'X,F,F,F'
  readcol, 'tmpEW.dat', EWFit, EWMin, EWMax, f = 'X,F,F,F'
  if oiii then $
     readcol, 'tmpEW_OIII.dat', EWFit_OIII, EWMin_OIII, EWMax_OIII, f = 'X,F,F,F'

  if NOT oiii then $
     savedata = {ZFIT: zFit[0], ZLO: zMin[0], ZHI: zMax[0], $
                 AGEFIT: ageFit[0], AGELO: ageMin[0], AGEHI: ageMax[0], $
                 EWFIT:EWFit[0], EWLO: EWmin[0], EWHI: EWmax[0]} $
  else $
     savedata = {ZFIT: zFit[0], ZLO: zMin[0], ZHI: zMax[0], $
                 AGEFIT: ageFit[0], AGELO: ageMin[0], AGEHI: ageMax[0], $
                 EWFIT:EWFit[0], EWLO: EWmin[0], EWHI: EWmax[0], $
                 EWFIT_OIII:EWFit_OIII[0], EWLO_OIII: EWmin_OIII[0], EWHI_OIII: EWmax_OIII[0]}
  
  
  RETURN, savedata
end

;;
;;
;;

pro plotPyspecResults, field, objno, PA, $
                       MASS = mass, MAG = mag

  if NOT keyword_set(MASS) then mass = -99
  if NOT keyword_set(MAG) then mag = 1.

  mass /= mag

  pa = string(pa, f = '(I1)')
  objno = strcompress(string(fix(objno)), /rem)
  
  dir = field+'/'+objno+'_results'
  fitsDir = field+'/'
   
  ;; Read results
  innerRes = grabPyspecResults(dir+'/inner'+pa+'.analyze')
  interRes = grabPyspecResults(dir+'/inter'+pa+'.analyze')
  outerRes = grabPyspecResults(dir+'/outer'+pa+'.analyze')

  ;; Read G102
  readcol, dir+'/inner'+pa+'.bestspec0', $
           lambdaB, innTraceB, innErrB, innPolyB, innModelB, innUsedB
  readcol, dir+'/inter'+pa+'.bestspec0', $
           lambdaB, intTraceB, intErrB, intPolyB, intModelB, intUsedB
  readcol, dir+'/outer'+pa+'.bestspec0', $
           lambdaB, outTraceB, outErrB, outPolyB, outModelB, outUsedB

  ;; Read G141
  readcol, dir+'/inner'+pa+'.bestspec1', $
           lambdaR, innTraceR, innErrR, innPolyR, innModelR, innUsedR
  readcol, dir+'/inter'+pa+'.bestspec1', $
           lambdaR, intTraceR, intErrR, intPolyR, intModelR, intUsedR
  readcol, dir+'/outer'+pa+'.bestspec1', $
           lambdaR, outTraceR, outErrR, outPolyR, outModelR, outUsedR

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

  innErrB *= innPolyB
  intErrB *= intPolyB
  outErrB *= outPolyB
  innErrB[where(~innUsedB)] = !values.F_NAN
  intErrB[where(~intUsedB)] = !values.F_NAN
  outErrB[where(~outUsedB)] = !values.F_NAN

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

  innErrR *= innPolyR
  intErrR *= intPolyR
  outErrR *= outPolyR
  innErrR[where(~innUsedR)] = !values.F_NAN
  intErrR[where(~intUsedR)] = !values.F_NAN
  outErrR[where(~outUsedR)] = !values.F_NAN

  ;; Combine
  lambda = [lambdaB, lambdaR] / (1 + innerRes.ZFIT)
  inner  = [innSpecB, innSpecR]
  inter  = [intSpecB, intSpecR]
  outer  = [outSpecB, outSpecR]

  innMod = [innModelB, innModelR]
  intMod = [intModelB, intModelR]
  outMod = [outModelB, outModelR]

  innErr = [innErrB, innErrR]
  intErr = [intErrB, intErrR]
  outErr = [outErrB, outErrR]

  ;; Read raw data
  bluDat = mrdfits(fitsDir+'/'+objNo+'_'+pa+'_B.fits', 1)
  redDat  = mrdfits(fitsDir+'/'+objNo+'_'+pa+'_R.fits', 1)

  bDI = bluDat.DIRECT_IM
  rDI = redDat.DIRECT_IM
  s = size(bDI, /dim)
  pscale = 0.13 ;; WFC3IR "/pix
  run = (findgen(s[0]) - (s[0]+1)/2) * pscale
  run -= run/2
  
  b2D = bluDat.SPEC2D
  r2D = redDat.SPEC2D
  s2D = size(b2D, /dim)

  blam = bluDat.lambda
  rlam = redDat.lambda
  
  innOff = median(bluDat.INNERUP - bluDat.TRACE)
  intOff = median(bluDat.INTERUP - bluDat.TRACE)
  outOff = median(bluDat.OUTERUP - bluDat.TRACE)

  kpc    = 1. / zang(1.0, innerRes.ZFIT) ;; kpc / asec
  kpcPix = kpc * pscale ;; kpc / pix
  hlrKpc = bluDat.RE * kpcPix  ;; re in Kpc
  hlrObs = bluDat.RE * pscale  ;; re in asec
  normRadii = [0, mean([innOff,intOff]),mean([intOff,outOff])] / bluDat.RE
  
  set_plot, 'PS'
  outputName = field+'_'+objno+'_'+PA+'_withBestfit.eps'
  device, filename = outputName, $
          /col, /encap, /decomp, bits_per_pixel = 8, $
          xsize = 7, ysize = 5, /in

;  mass = 10.8 ;; From takahiro + lensing magnification estimated from HFF site @ 2.5x
;  re   = 4.0 * 0.13 / zang(1, 1.086)  ;; In Kpc; from 0753_1_B.fits
;  read_jpeg, 'MACS1149CLS_1096_rgb_rot.jpeg', bdi, true = 1
;  s = size(bdi)
;  w = s[2]
;  tp1 = 20. / 60. * w - 1
;  tp2 = 40. / 60. * w - 1
  
  plot, run[20:40], run[20:40], /nodat, $
        xtitle = greek('delta')+'!18x!X ["]', $
        ytitle = greek('delta')+'!18y!X ["]', /iso, $
        pos = [0.1, 0.75, 0.25, 0.90], $
        xthick = 4, ythick = 4, charthick = 3, charsize = 1, $
        xtickint = 1, ytickint = 1, title = 'direct image'
  cgimage, bdi[20:40,20:40], stretch = 7, /over, /neg
  oplot, run, replicate(innOff, n_elements(run))  * pscale, col = 255      , linesty = 0, thick = 1
  oplot, run, replicate(-innOff, n_elements(run)) * pscale, col = 255      , linesty = 0, thick = 1
  oplot, run, replicate(intOff, n_elements(run))  * pscale, col = '00a500'x, linesty = 0, thick = 1
  oplot, run, replicate(-intOff, n_elements(run)) * pscale, col = '00a500'x, linesty = 0, thick = 1
  oplot, run, replicate(outOff, n_elements(run))  * pscale, col = 'ff5500'x, linesty = 0, thick = 1
  oplot, run, replicate(-outOff, n_elements(run)) * pscale, col = 'ff5500'x, linesty = 0, thick = 1
  plot, run[20:40], run[20:40], /nodat, $
        pos = [0.1, 0.75, 0.25, 0.90], $
        xtickname = replicate(' ', 60), ytickname = replicate(' ',60), $
        xthick = 4, ythick = 4, charthick = 3, charsize = 1, $
        xtickint = 1, ytickint = 1, /noer, /iso, yminor = 2, xminor = 2, yticklen = 0.075, xticklen = 0.075
  
  st = !X.WINDOW[1]+0.05
  plot, blam / (1 + innerRes.ZFIT), findgen(n_elements(b2D[0,*])), /nodat, $
        xtitle = greek('lambda')+'!Drest!N ['+texToIDL('\AA')+']', $
        xthick = 5, ythick = 5, charthick = 3, charsize = 1, $
        pos = [st,!Y.WINDOW[0],0.55,!Y.WINDOW[1]], /noer, $
        ytickname = replicate(' ',60), yran = [20,40], title = 'G102', /xsty
  wid = !X.WINDOW[1] - st
  cgimage, b2D[*,20:40], stretch = 2, /over
;  cgimage, bluDat.EXTRACT_ZONES[*,20:40], stretch = 5, /over
  oplot, blam / (1 + innerRes.ZFIT), bludat.innerUP, col = 255      , linesty = 0, thick = 1
  oplot, blam / (1 + innerRes.ZFIT), bludat.innerDN, col = 255      , linesty = 0, thick = 1
  oplot, blam / (1 + innerRes.ZFIT), bludat.interUp, col = '00a500'x, linesty = 0, thick = 1
  oplot, blam / (1 + innerRes.ZFIT), bludat.interDN, col = '00a500'x, linesty = 0, thick = 1
  oplot, blam / (1 + innerRes.ZFIT), bludat.outerUP, col = 'ff5500'x, linesty = 0, thick = 1
  oplot, blam / (1 + innerRes.ZFIT), bludat.outerDN, col = 'ff5500'x, linesty = 0, thick = 1
  
  plot, rlam / (1 + innerRes.ZFIT), findgen(n_elements(r2D[0,*])), /nodat, $
        xtitle = greek('lambda')+'!Drest!N ['+texToIDL('\AA')+']', $
        xthick = 5, ythick = 5, charthick = 3, charsize = 1, $
        pos = [!X.WINDOW[1]+0.05,!Y.WINDOW[0],!X.WINDOW[1] + 0.05 + wid,!Y.WINDOW[1]], /noer, $
        ytickname = replicate(' ',60), yran = [20,40], title = 'G141', /xsty
  cgimage, r2D[*,20:40], stretch = 2, /over
  oplot, rlam / (1 + innerRes.ZFIT), reddat.innerUP, col = 255      , linesty = 0, thick = 1
  oplot, rlam / (1 + innerRes.ZFIT), reddat.innerDN, col = 255      , linesty = 0, thick = 1
  oplot, rlam / (1 + innerRes.ZFIT), reddat.interUp, col = '00a500'x, linesty = 0, thick = 1
  oplot, rlam / (1 + innerRes.ZFIT), reddat.interDN, col = '00a500'x, linesty = 0, thick = 1
  oplot, rlam / (1 + innerRes.ZFIT), reddat.outerUP, col = 'ff5500'x, linesty = 0, thick = 1
  oplot, rlam / (1 + innerRes.ZFIT), reddat.outerDN, col = 'ff5500'x, linesty = 0, thick = 1
  
  plot, lambda, inner, $
        xran = [3000,8000], yran = [0,400], /xsty, $
        xtitle = greek('lambda')+'!Drest!N ['+texToIDL('\AA')+']', $
        ytitle = '!18f!X', /nodat, pos = [0.1,0.1,0.5,0.6], /noer, $
        xthick = 4, ythick = 4, charthick = 3, charsize = 1, yminor = 2, /ysty
  oplot, lambda, smooth(inner, 2, /nan), col = 255      , thick = 3
  oplot, lambda, smooth(inter, 2, /nan), col = '00a500'x, thick = 3
  oplot, lambda, smooth(outer, 2, /nan), col = 'ffa500'x, thick = 3
  oplot, lambda, innMod, thick = 4, col = 100
  oplot, lambda, intMod, thick = 4, col = '004400'x
  oplot, lambda, outMod, thick = 4, col = 'ff0000'x
  oplot, lambda, innErr, thick = 3, col = '555555'x
  legend, /top, /right, box = 0, $
          ['!18R!X<0.25 !18R!X!De!N', 'fit', $
           '0.25!18R!X!De!N<!18R!X<0.75!18R!X!De!N', 'fit', $
           '0.75!18R!X!De!N<!18R!X<1.50!18R!X!De!N', 'fit'], $
          col = [255,100,'00a500'x,'004400'x,'ffa500'x,'ff0000'x], $
          pspacing = 1, thick = [2,4,2,4,2,4], $
          charsize = 0.7, charthick = 2, linesty = 0

  rgrid = normRadii
  ageTrend = [innerRes.AGEFIT, interRes.AGEFIT, outerRes.AGEFIT]
  ageHi    = [innerRes.AGEHI, interRes.AGEHI, outerRes.AGEHI]
  ageLo    = [innerRes.AGELO, interRes.AGELO, outerRes.AGELO]
  hiBar    = 0.434 * (agehi - ageTrend) / ageTrend
  loBar    = 0.434 * (ageTrend - ageLo) / ageTrend
  
  plotsym, 0, /fill
  plot, rgrid, alog10(ageTrend), /nodat, $
        yran = [min(alog10(ageTrend)) - 0.5,max(alog10(ageTrend)) + 0.5], $
;        ytitle = 'Age [Gyr]', $
        ytickname = replicate(' ',60), $
        xthick = 4, ythick = 4, charthick = 3, $
        pos = [0.55,0.35, 0.875,0.6], $
        /noer, xtickname = replicate(' ', 60), charsize = 1, ysty = 8+1
  axis, yaxis = 1, yran = !Y.CRANGE, ythick = 4, charsize = 1, ytitle = 'log Age [Gyr]', $
        charthick = 3, /ysty
  oplot, !X.CRANGE, replicate(alog10(ageTrend[0]), 2), linesty = 1, col = '555555'x
  oploterror, [rgrid[0]], [alog10(ageTrend[0])], [Hibar[0]], /hibar, psym = 8, symsize = 1, $
              col = cgcolor('red'), errcol = cgcolor('red'), errthick = 3
  oploterror, [rgrid[0]], [alog10(ageTrend[0])], [Lobar[0]], /lobar, psym = 8, symsize = 1, $
              col = cgcolor('red'), errcol = cgcolor('red'), errthick = 3
  oploterror, [rgrid[1]], [alog10(ageTrend[1])], [Hibar[1]], /hibar, psym = 8, symsize = 1, $
              col = '00a500'x, errcol = '00a500'x, errthick = 3
  oploterror, [rgrid[1]], [alog10(ageTrend[1])], [Lobar[1]], /lobar, psym = 8, symsize = 1, $
              col = '00a500'x, errcol = '00a500'x, errthick = 3
  oploterror, [rgrid[2]], [alog10(ageTrend[2])], [Hibar[2]], /hibar, psym = 8, symsize = 1, $
              col = 'ffa500'x, errcol = 'ffa500'x, errthick = 3
  oploterror, [rgrid[2]], [alog10(ageTrend[2])], [Lobar[2]], /lobar, psym = 8, symsize = 1, $
              col = 'ffa500'x, errcol = 'ffa500'x, errthick = 3
  legend, /bottom, /left, box = 0, $ ; pos = [0.17,8.4], /data
          ['!18z!X='+string(bluDat.Z, f = '(F5.3)'), $
           '!18R!X!De!N='+string(hlrObs, f = '(F4.2)')+'"='+string(hlrKpc, f = '(F3.1)')+' kpc'], $
          charsize = 1, charthick = 3
  
  EWTrend = [innerRes.EWFIT, interRes.EWFIT, outerRes.EWFIT]
  EWHi    = [innerRes.EWHI, interRes.EWHI, outerRes.EWHI]
  EWLo    = [innerRes.EWLO, interRes.EWLO, outerRes.EWLO]
  hiBar    = EWhi - EWTrend
  loBar    = EWTrend - EWLo
  
  plot, rgrid, EWTrend, /nodat, $
;        ytitle = 'Age [Gyr]', $
        yran = [min(EWtrend) - 5, max(EWtrend) + 10], $
        ytickname = replicate(' ',60), $
        xthick = 4, ythick = 4, charthick = 3, $
        pos = [0.55,0.1,0.875,0.35], $
        /noer, charsize = 1, ysty = 8+1, $
        xtitle = '!18R/R!X!De!N', yminor = 2
  axis, yaxis = 1, ythick = 4, charsize = 1, ytitle = 'EW(H'+greek('alpha')+') ['+texToIDL('\AA')+']', $
        charthick = 3, yran = !y.crange, yminor = 2
  oplot, !X.CRANGE, replicate(EWTrend[0], 2), linesty = 1, col = '555555'x
  oploterror, [rgrid[0]], [EWTrend[0]], [Hibar[0]], /hibar, psym = 8, symsize = 1, $
              col = cgcolor('red'), errcol = cgcolor('red'), errthick = 3
  oploterror, [rgrid[0]], [EWTrend[0]], [Lobar[0]], /lobar, psym = 8, symsize = 1, $
              col = cgcolor('red'), errcol = cgcolor('red'), errthick = 3
  oploterror, [rgrid[1]], [EWTrend[1]], [Hibar[1]], /hibar, psym = 8, symsize = 1, $
              col = '00a500'x, errcol = '00a500'x, errthick = 3
  oploterror, [rgrid[1]], [EWTrend[1]], [Lobar[1]], /lobar, psym = 8, symsize = 1, $
              col = '00a500'x, errcol = '00a500'x, errthick = 3
  oploterror, [rgrid[2]], [EWTrend[2]], [Hibar[2]], /hibar, psym = 8, symsize = 1, $
              col = 'ffa500'x, errcol = 'ffa500'x, errthick = 3
  oploterror, [rgrid[2]], [EWTrend[2]], [Lobar[2]], /lobar, psym = 8, symsize = 1, $
              col = 'ffa500'x, errcol = 'ffa500'x, errthick = 3
  legend, /bottom, /left, box = 0, $ ; pos = [0.17,max(EWtrend) + 5], /data
          'log!18M!X!D*!N='+string(mass, f = '(F4.1)'), $
          charsize = 1, charthick = 3
  
  device, /close
;  spawn, 'gv '+outputName+' &'

;  stop
  
  
;  plot, lambda, inner, $
;        xran = [3600,8100], yran = [0,400], $
;        xtitle = greek('lambda')+' ['+texToIDL('\AA')+']', $
;        ytitle = '!18f!X', /nodat
;  oplot, lambda, inner, col = 255
;  oplot, lambda, inter, col = '00a500'x
;  oplot, lambda, outer, col = 'ff5500'x
;  oplot, lambda, innMod
;  oplot, lambda, intMod
;  oplot, lambda, outMod
  
end

pro doAll

  plotPyspecResults, 'M0717',   '399', '1', mass = 11.34, mag = 4.29 ;       0.99   FIT 20160504 
  plotPyspecResults, 'M0717',  '1236', '1', mass = 10.73, mag = 7.28 ;      2.08   FIT 20160504
  plotPyspecResults, 'M1149',   '753', '1', mass = 11.23, mag = 2.48 ;      0.94   FIT 20160504
  plotPyspecResults, 'M2129',   '839', '1', mass = 11.23, mag = 3.81 ;     -1.03   FIT 20160504
  plotPyspecResults, 'M2129',   '841', '1', mass = 10.98, mag = 3.81 ;     -2.64   FIT 20160504
  plotPyspecResults, 'M2129',   '843', '1', mass = 11.09, mag = 3.81 ;      0.19   FIT 20160504
  plotPyspecResults, 'M2129',   '845', '1', mass = 11.11, mag = 3.81 ;     -6.32   FIT 20160504
  plotPyspecResults, 'M2129',  '1681', '1', mass = 11.14, mag = 3.81 ;      1.39   FIT 20160504

  plotPyspecResults, 'M0717',   '399', '2', mass = 11.34, mag = 4.29 ;       0.99   FIT 20160504 
  plotPyspecResults, 'M0717',  '1236', '2', mass = 10.73, mag = 7.28 ;      2.08   FIT 20160504
  plotPyspecResults, 'M1149',   '753', '2', mass = 11.23, mag = 2.48 ;      0.94   FIT 20160504
  plotPyspecResults, 'M2129',   '839', '2', mass = 11.23, mag = 3.81 ;     -1.03   FIT 20160504
  plotPyspecResults, 'M2129',   '841', '2', mass = 10.98, mag = 3.81 ;     -2.64   FIT 20160504
  plotPyspecResults, 'M2129',   '843', '2', mass = 11.09, mag = 3.81 ;      0.19   FIT 20160504
  plotPyspecResults, 'M2129',   '845', '2', mass = 11.11, mag = 3.81 ;     -6.32   FIT 20160504
  plotPyspecResults, 'M2129',  '1681', '2', mass = 11.14, mag = 3.81 ;      1.39   FIT 20160504

end
