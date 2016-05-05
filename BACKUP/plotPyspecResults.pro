pro plotPyspecResults, dir

  ;; Read results
  innerRes = grabPyspecResults('to_louis/output_v1/inner.analyze')
  interRes = grabPyspecResults('to_louis/output_v1/inter.analyze')
  outerRes = grabPyspecResults('to_louis/output_v1/outer.analyze')

  ;; Read G102
  readcol, dir+'/inner'+pa+'bestspec0', $
           lambdaB, innTraceB, innErrB, innPolyB, innModelB, innUsedB
  readcol, dir+'/inter'+pa+'bestspec0', $
           lambdaB, intTraceB, intErrB, intPolyB, intModelB, intUsedB
  readcol, dir+'/outer'+pa+'bestspec0', $
           lambdaB, outTraceB, outErrB, outPolyB, outModelB, outUsedB

  ;; Read G141
  readcol, dir+'/inner'+pa+'bestspec1', $
           lambdaR, innTraceR, innErrR, innPolyR, innModelR, innUsedR
  readcol, dir+'/inter'+pa+'bestspec1', $
           lambdaR, intTraceR, intErrR, intPolyR, intModelR, intUsedR
  readcol, dir+'/outer'+pa+'bestspec1', $
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
  bluDat = mrdfits('753_1_B.fits', 1)
  redDat  = mrdfits('753_1_R.fits', 1)

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
  
  set_plot, 'PS'
  device, filename = 'M1149_0753_info.eps', $
          /col, /encap, /decomp, bits_per_pixel = 8, $
          xsize = 7, ysize = 5, /in

  mass = 10.8 ;; From takahiro + lensing magnification estimated from HFF site @ 2.5x
;  re   = 4.0 * 0.13 / zang(1, 1.086)  ;; In Kpc; from 0753_1_B.fits
  read_jpeg, 'MACS1149CLS_1096_rgb_rot.jpeg', bdi, true = 1
  s = size(bdi)
  w = s[2]
  tp1 = 20. / 60. * w - 1
  tp2 = 40. / 60. * w - 1
  
  plot, run[20:40], run[20:40], /nodat, $
        xtitle = greek('delta')+'!18x!X ["]', $
        ytitle = greek('delta')+'!18y!X ["]', /iso, $
        pos = [0.1, 0.75, 0.25, 0.90], $
        xthick = 4, ythick = 4, charthick = 3, charsize = 1, $
        xtickint = 1, ytickint = 1, title = 'direct image'
  cgimage, bdi[*,tp1:tp2,tp1:tp2], /over;, /neg
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
        ytickname = replicate(' ',60), yran = [20,40], title = 'G102'
  wid = !X.WINDOW[1] - st
  cgimage, b2D[*,20:40], stretch = 5, /over
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
        ytickname = replicate(' ',60), yran = [20,40], title = 'G141'
  cgimage, r2D[*,20:40], stretch = 5, /over
  oplot, rlam / (1 + innerRes.ZFIT), reddat.innerUP, col = 255      , linesty = 0, thick = 1
  oplot, rlam / (1 + innerRes.ZFIT), reddat.innerDN, col = 255      , linesty = 0, thick = 1
  oplot, rlam / (1 + innerRes.ZFIT), reddat.interUp, col = '00a500'x, linesty = 0, thick = 1
  oplot, rlam / (1 + innerRes.ZFIT), reddat.interDN, col = '00a500'x, linesty = 0, thick = 1
  oplot, rlam / (1 + innerRes.ZFIT), reddat.outerUP, col = 'ff5500'x, linesty = 0, thick = 1
  oplot, rlam / (1 + innerRes.ZFIT), reddat.outerDN, col = 'ff5500'x, linesty = 0, thick = 1
  
  plot, lambda, inner, $
        xran = [3600,8100], yran = [50,400], /xsty, $
        xtitle = greek('lambda')+'!Drest!N ['+texToIDL('\AA')+']', $
        ytitle = '!18f!X', /nodat, pos = [0.1,0.1,0.5,0.6], /noer, $
        xthick = 4, ythick = 4, charthick = 3, charsize = 1, yminor = 2, /ysty
  oplot, lambda, smooth(inner, 2, /nan), col = 255      , thick = 3
  oplot, lambda, smooth(inter, 2, /nan), col = '00a500'x, thick = 3
  oplot, lambda, smooth(outer, 2, /nan), col = 'ffa500'x, thick = 3
  oplot, lambda, innMod, thick = 4, col = 100
  oplot, lambda, intMod, thick = 4, col = '004400'x
  oplot, lambda, outMod, thick = 4, col = 'ff0000'x
  legend, /top, /right, box = 0, $
          ['!18R!X<0.25 !18R!X!De!N', 'fit', $
           '0.25!18R!X!De!N<!18R!X<0.75!18R!X!De!N', 'fit', $
           '0.75!18R!X!De!N<!18R!X<1.50!18R!X!De!N', 'fit'], $
          col = [255,100,'00a500'x,'004400'x,'ffa500'x,'ff0000'x], $
          pspacing = 1, thick = [2,4,2,4,2,4], $
          charsize = 0.7, charthick = 2, linesty = 0

  rgrid = [mean(bluDat.INNERUP - bluDat.TRACE) / bluDat.RE, $
           mean(bluDat.INTERUP - bluDat.TRACE) / bluDat.RE, $
           mean(bluDat.OUTERUP - bluDat.TRACE) / bluDat.RE]
  ageTrend = [innerRes.AGEFIT, interRes.AGEFIT, outerRes.AGEFIT]
  ageHi    = [innerRes.AGEHI, interRes.AGEHI, outerRes.AGEHI]
  ageLo    = [innerRes.AGELO, interRes.AGELO, outerRes.AGELO]
  hiBar    = 0.434 * (agehi - ageTrend) / ageTrend
  loBar    = 0.434 * (ageTrend - ageLo) / ageTrend
  
  plotsym, 0, /fill
  plot, rgrid, alog10(ageTrend), /nodat, $
        yran = [7.6,10], $
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
  legend, pos = [0.17,8.4], /data, box = 0, $
          ['!18z!X=1.086', $
           '!18R!X!De!N='+string(hlrObs, f = '(F4.2)')+'"='+string(hlrKpc, f = '(F3.1)')+' kpc'], $
          charsize = 1, charthick = 3
  
  EWTrend = [innerRes.EWFIT, interRes.EWFIT, outerRes.EWFIT]
  EWHi    = [innerRes.EWHI, interRes.EWHI, outerRes.EWHI]
  EWLo    = [innerRes.EWLO, interRes.EWLO, outerRes.EWLO]
  hiBar    = EWhi - EWTrend
  loBar    = EWTrend - EWLo
  
  plot, rgrid, EWTrend, /nodat, $
;        ytitle = 'Age [Gyr]', $
        yran = [0,40], $
        ytickname = replicate(' ',60), $
        xthick = 4, ythick = 4, charthick = 3, $
        pos = [0.55,0.1,0.875,0.35], $
        /noer, charsize = 1, ysty = 8+1, $
        xtitle = '!18R/R!X!De!N', yminor = 2
  axis, yaxis = 1, ythick = 4, charsize = 1, ytitle = 'EW(H'+greek('alpha')+') ['+texToIDL('\AA')+']', $
        charthick = 3, yran = [0,40], yminor = 2
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
  legend, pos = [0.17,38], /data, box = 0, $
          'log!18M!X!D*!N='+string(mass, f = '(F4.1)'), $
          charsize = 1, charthick = 3
  
  device, /close
  spawn, 'open M1149_0753_info.eps'
        
  
  stop
  
  
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
  
  stop
end
