pro getFiles, indir, outprefix

  spawn, 'ls '+indir+'/*_1_?.fits > '+outprefix+'_PA1.list'
  spawn, 'ls '+indir+'/*_2_?.fits > '+outprefix+'_PA2.list'

end

pro getUnified, indir, outprefix

  spawn, 'ls '+indir+'/*_1_unified_1D.fits > '+outprefix+'_PA1.list'
  spawn, 'ls '+indir+'/*_2_unified_1D.fits > '+outprefix+'_PA2.list'
  
end

; plotMultiSpec, 'M0717', 'M0717', /unified
; plotMultiSpec, 'M1149', 'M1149', /unified
; plotMultiSpec, 'M2129', 'M2129', /unified

; plotMultiSpec, 'M0717_faint', 'M0717_faint', /unified
; plotMultiSpec, 'M1149_faint', 'M1149_faint', /unified
; plotMultiSpec, 'M1423_faint', 'M1423_faint', /unified
; plotMultiSpec, 'M2129_faint', 'M2129_faint', /unified
; plotMultiSpec, 'R1347_faint', 'R1347_faint', /unified
pro plotMultiSpec, field, outdir, UNIFIED = unified

  if NOT keyword_Set(UNIFIED) then unified = 0 else unified = 1

  getFiles, field, field+'_2Ddata'
  readcol, field+'_2Ddata_PA1.list', Pa1, f ='A'
  readcol, field+'_2Ddata_PA2.list', Pa2, f ='A'

  if unified then begin
     getUnified, field, field+'_U1Ddata'
     readcol, field+'_U1Ddata_PA1.list', uPa1, f ='A'
     readcol, field+'_U1Ddata_PA2.list', uPa2, f ='A'
  endif
  
  pa1files = Pa1
  pa2files = Pa2 
  fieldnames = replicate(field, n_elements(Pa1))
  ids = strarr(n_elements(pa1files))
  for ii = 0, n_elements(pa1files) - 1 do $
     ids[ii] = strmid(pa1files[ii], strpos(pa1files[ii], '/')+1, 4)
  
  npa1 = n_elements(pa1files)
  npa2 = n_elements(pa2files)

;  stop
  
  counter = 0
  while counter le npa1 - 1 do begin

     pa1b = mrdfits(pa1files[counter], 1)
     pa1r = mrdfits(pa1files[counter+1], 1)
     pa2b = mrdfits(pa2files[counter], 1)
     pa2r = mrdfits(pa2files[counter+1], 1)

     di1 = pa1b.DIRECT_IM
     di2 = pa2b.DIRECT_IM

     s = size(di1, /dim)
     pscale = 0.13 ;; WFC3IR "/pix
     run = (findgen(s[0]) - (s[0]+1)/2) * pscale
     run -= run/2
  
     b2D1 = pa1b.SPEC2D
     r2D1 = pa1r.SPEC2D
     b2D2 = pa2b.SPEC2D
     r2D2 = pa2r.SPEC2D
     
     s2D = size(b2D, /dim)
    
     blam1 = pa1b.LAMBDA
     rlam1 = pa1r.LAMBDA
     blam2 = pa2b.LAMBDA
     rlam2 = pa2r.LAMBDA

     if NOT unified then begin
        lambda1 = [blam1, rlam1] / (1 + pa1b.Z)
        inner1 = [pa1b.F_INNER / pa1b.SENSITIVITY, pa1r.F_INNER / pa1r.SENSITIVITY]
        inter1 = [pa1b.F_INTER / pa1b.SENSITIVITY, pa1r.F_INTER / pa1r.SENSITIVITY]
        outer1 = [pa1b.F_OUTER / pa1b.SENSITIVITY, pa1r.F_OUTER / pa1r.SENSITIVITY]
        err1   = sqrt([pa1b.VAR_INNER, pa1r.VAR_INNER]) / [pa1b.SENSITIVITY, pa1r.SENSITIVITY]
        inner1[where(inner1 / err1 lt 2)] = !VALUES.F_NaN
        inter1[where(inter1 / err1 lt 2)] = !VALUES.F_NaN
        outer1[where(outer1 / err1 lt 2)] = !VALUES.F_NaN   
        
        lambda2 = [blam2, rlam2] / (1 + pa1b.Z)
        inner2 = [pa2b.F_INNER / pa2b.SENSITIVITY, pa2r.F_INNER / pa2r.SENSITIVITY]
        inter2 = [pa2b.F_INTER / pa2b.SENSITIVITY, pa2r.F_INTER / pa2r.SENSITIVITY]
        outer2 = [pa2b.F_OUTER / pa2b.SENSITIVITY, pa2r.F_OUTER / pa2r.SENSITIVITY]
        err2   = sqrt([pa2b.VAR_INNER, pa2r.VAR_INNER]) / [pa2b.SENSITIVITY, pa2r.SENSITIVITY]
        inner2[where(inner2 / err2 lt 2)] = !VALUES.F_NaN
        inter2[where(inter2 / err2 lt 2)] = !VALUES.F_NaN
        outer2[where(outer2 / err2 lt 2)] = !VALUES.F_NaN   
     endif else begin
        D1 = mrdfits(uPa1[counter / 2], 1)
        lambda1 = D1.LAMBDA / (1 + D1.REDSHIFT)
        inner1 = D1.INNER
        inter1 = D1.INTER
        outer1 = D1.OUTER
        err1   = sqrt(D1.INNER_VAR)

        D2 = mrdfits(uPa2[counter / 2], 1)
        lambda2 = D2.LAMBDA / (1 + D2.REDSHIFT)
        inner2 = D2.INNER
        inter2 = D2.INTER
        outer2 = D2.OUTER
        err2   = sqrt(D2.INNER_VAR)
     endelse
     
     innOff1 = median(pa1b.INNERUP - pa1b.TRACE)
     intOff1 = median(pa1b.INTERUP - pa1b.TRACE)
     outOff1 = median(pa1b.OUTERUP - pa1b.TRACE)

     innOff2 = median(pa2b.INNERUP - pa2b.TRACE)
     intOff2 = median(pa2b.INTERUP - pa2b.TRACE)
     outOff2 = median(pa2b.OUTERUP - pa2b.TRACE)

     kpc    = 1. / zang(1.0, pa1b.Z)            ;; kpc / asec
     kpcPix = kpc * pscale                      ;; kpc / pix
     hlrObs = sqrt(pa1b.RE * pa2b.RE) * pscale ;; re in asec
     hlrKpc = sqrt(pa1b.RE * pa2b.RE) * kpcPix ;; re in Kpc
     
    
     set_plot, 'PS'
     device, filename = outdir+'/'+fieldnames[counter]+'_'+ids[counter]+'_info.eps', $
             /col, /encap, /decomp, bits_per_pixel = 8, $
             xsize = 7, ysize = 5, /in    
     plot, run[20:40], run[20:40], /nodat, $
           xtitle = greek('delta')+'!18x!X ["]', $
           ytitle = greek('delta')+'!18y!X ["]', /iso, $
           pos = [0.1, 0.8, 0.25, 0.90], $
           xthick = 4, ythick = 4, charthick = 3, charsize = 1, $
           xtickint = 1, ytickint = 1, title = 'direct image'
     cgimage, di1[20:40,20:40], stretch = 7, /over, /neg
     oplot, run, replicate(innOff1, n_elements(run))  * pscale, col = 255      , linesty = 0, thick = 1
     oplot, run, replicate(-innOff1, n_elements(run)) * pscale, col = 255      , linesty = 0, thick = 1
     oplot, run, replicate(intOff1, n_elements(run))  * pscale, col = '00a500'x, linesty = 0, thick = 1
     oplot, run, replicate(-intOff1, n_elements(run)) * pscale, col = '00a500'x, linesty = 0, thick = 1
     oplot, run, replicate(outOff1, n_elements(run))  * pscale, col = 'ff5500'x, linesty = 0, thick = 1
     oplot, run, replicate(-outOff1, n_elements(run)) * pscale, col = 'ff5500'x, linesty = 0, thick = 1
     plot, run[20:40], run[20:40], /nodat, $
           pos = [!X.WINDOW[0], !Y.WINDOW[0], !X.WINDOW[1], !Y.WINDOW[1]], $
           xtickname = replicate(' ', 60), ytickname = replicate(' ',60), $
           xthick = 4, ythick = 4, charthick = 3, charsize = 1, $
           xtickint = 1, ytickint = 1, /noer, /iso, yminor = 2, $
           xminor = 2, yticklen = 0.075, xticklen = 0.075
     
     st = !X.WINDOW[1]+0.05
     plot, blam1 / (1 + pa1b.Z), findgen(n_elements(b2D1[0,*])), /nodat, $
           xtitle = greek('lambda')+'!Drest!N ['+texToIDL('\AA')+']', $
           xthick = 5, ythick = 5, charthick = 3, charsize = 1, $
           pos = [st,!Y.WINDOW[0],0.55,!Y.WINDOW[1]], /noer, $
           ytickname = replicate(' ',60), yran = [20,40], title = 'G102', /xsty
     wid = !X.WINDOW[1] - st
     cgimage, b2D1[*,20:40], stretch = 2, /over
;  cgimage, bluDat.EXTRACT_ZONES[*,20:40], stretch = 2, /over
     oplot, blam1 / (1 + pa1b.Z), pa1b.innerUP, col = 255      , linesty = 0, thick = 1
     oplot, blam1 / (1 + pa1b.Z), pa1b.innerDN, col = 255      , linesty = 0, thick = 1
     oplot, blam1 / (1 + pa1b.Z), pa1b.interUp, col = '00a500'x, linesty = 0, thick = 1
     oplot, blam1 / (1 + pa1b.Z), pa1b.interDN, col = '00a500'x, linesty = 0, thick = 1
     oplot, blam1 / (1 + pa1b.Z), pa1b.outerUP, col = 'ff5500'x, linesty = 0, thick = 1
     oplot, blam1 / (1 + pa1b.Z), pa1b.outerDN, col = 'ff5500'x, linesty = 0, thick = 1
  
     plot, rlam1 / (1 + pa1b.Z), findgen(n_elements(r2D1[0,*])), /nodat, $
        xtitle = greek('lambda')+'!Drest!N ['+texToIDL('\AA')+']', $
        xthick = 5, ythick = 5, charthick = 3, charsize = 1, $
        pos = [!X.WINDOW[1]+0.05,!Y.WINDOW[0],!X.WINDOW[1] + 0.05 + wid,!Y.WINDOW[1]], /noer, $
        ytickname = replicate(' ',60), yran = [20,40], title = 'G141', /xsty
     cgimage, r2D1[*,20:40], stretch = 2, /over
     oplot, rlam1 / (1 + pa1b.Z), pa1r.innerUP, col = 255      , linesty = 0, thick = 1
     oplot, rlam1 / (1 + pa1b.Z), pa1r.innerDN, col = 255      , linesty = 0, thick = 1
     oplot, rlam1 / (1 + pa1b.Z), pa1r.interUp, col = '00a500'x, linesty = 0, thick = 1
     oplot, rlam1 / (1 + pa1b.Z), pa1r.interDN, col = '00a500'x, linesty = 0, thick = 1
     oplot, rlam1 / (1 + pa1b.Z), pa1r.outerUP, col = 'ff5500'x, linesty = 0, thick = 1
     oplot, rlam1 / (1 + pa1b.Z), pa1r.outerDN, col = 'ff5500'x, linesty = 0, thick = 1
     cgtext, 0.95, mean(!Y.WINDOW), /norm, 'PA1', charsize = 1.5, charthick = 4, orien = 270, align = 0.5
     
     plot, run[20:40], run[20:40], /nodat, $
           xtitle = greek('delta')+'!18x!X ["]', $
           ytitle = greek('delta')+'!18y!X ["]', /iso, $
           pos = [0.1, 0.65, 0.25, 0.75], $
           xthick = 4, ythick = 4, charthick = 3, charsize = 1, $
           xtickint = 1, ytickint = 1, /noer
     cgimage, di2[20:40,20:40], stretch = 7, /over, /neg
     oplot, run, replicate(innOff2, n_elements(run))  * pscale, col = 255      , linesty = 0, thick = 1
     oplot, run, replicate(-innOff2, n_elements(run)) * pscale, col = 255      , linesty = 0, thick = 1
     oplot, run, replicate(intOff2, n_elements(run))  * pscale, col = '00a500'x, linesty = 0, thick = 1
     oplot, run, replicate(-intOff2, n_elements(run)) * pscale, col = '00a500'x, linesty = 0, thick = 1
     oplot, run, replicate(outOff2, n_elements(run))  * pscale, col = 'ff5500'x, linesty = 0, thick = 1
     oplot, run, replicate(-outOff2, n_elements(run)) * pscale, col = 'ff5500'x, linesty = 0, thick = 1
     plot, run[20:40], run[20:40], /nodat, $
           pos = [!X.WINDOW[0], !Y.WINDOW[0], !X.WINDOW[1], !Y.WINDOW[1]], $
           xtickname = replicate(' ', 60), ytickname = replicate(' ',60), $
           xthick = 4, ythick = 4, charthick = 3, charsize = 1, $
           xtickint = 1, ytickint = 1, /noer, /iso, $
           yminor = 2, xminor = 2, yticklen = 0.075, xticklen = 0.075
     
     st = !X.WINDOW[1]+0.05
     plot, blam2 / (1 + pa2b.Z), findgen(n_elements(b2D2[0,*])), /nodat, $
           xtitle = greek('lambda')+'!Drest!N ['+texToIDL('\AA')+']', $
           xthick = 5, ythick = 5, charthick = 3, charsize = 1, $
           pos = [st,!Y.WINDOW[0],0.55,!Y.WINDOW[1]], /noer, $
           ytickname = replicate(' ',60), yran = [20,40], /xsty
     wid = !X.WINDOW[1] - st
     cgimage, b2D2[*,20:40], stretch = 2, /over
;  cgimage, bluDat.EXTRACT_ZONES[*,20:40], stretch = 2, /over
     oplot, blam2 / (1 + pa2b.Z), pa2b.innerUP, col = 255      , linesty = 0, thick = 1
     oplot, blam2 / (1 + pa2b.Z), pa2b.innerDN, col = 255      , linesty = 0, thick = 1
     oplot, blam2 / (1 + pa2b.Z), pa2b.interUp, col = '00a500'x, linesty = 0, thick = 1
     oplot, blam2 / (1 + pa2b.Z), pa2b.interDN, col = '00a500'x, linesty = 0, thick = 1
     oplot, blam2 / (1 + pa2b.Z), pa2b.outerUP, col = 'ff5500'x, linesty = 0, thick = 1
     oplot, blam2 / (1 + pa2b.Z), pa2b.outerDN, col = 'ff5500'x, linesty = 0, thick = 1
  
     plot, rlam2 / (1 + pa2b.Z), findgen(n_elements(r2D2[0,*])), /nodat, $
        xtitle = greek('lambda')+'!Drest!N ['+texToIDL('\AA')+']', $
        xthick = 5, ythick = 5, charthick = 3, charsize = 1, $
        pos = [!X.WINDOW[1]+0.05,!Y.WINDOW[0],!X.WINDOW[1] + 0.05 + wid,!Y.WINDOW[1]], /noer, $
        ytickname = replicate(' ',60), yran = [20,40], /xsty
     cgimage, r2D2[*,20:40], stretch = 2, /over
     oplot, rlam2 / (1 + pa2b.Z), pa2r.innerUP, col = 255      , linesty = 0, thick = 1
     oplot, rlam2 / (1 + pa2b.Z), pa2r.innerDN, col = 255      , linesty = 0, thick = 1
     oplot, rlam2 / (1 + pa2b.Z), pa2r.interUp, col = '00a500'x, linesty = 0, thick = 1
     oplot, rlam2 / (1 + pa2b.Z), pa2r.interDN, col = '00a500'x, linesty = 0, thick = 1
     oplot, rlam2 / (1 + pa2b.Z), pa2r.outerUP, col = 'ff5500'x, linesty = 0, thick = 1
     oplot, rlam2 / (1 + pa2b.Z), pa2r.outerDN, col = 'ff5500'x, linesty = 0, thick = 1
     cgtext, 0.95, mean(!Y.WINDOW), /norm, 'PA2', charsize = 1.5, charthick = 4, orien = 270, align = 0.5
     
     plot, lambda1, inner1, $
           yran = [0,500], $
           xtitle = greek('lambda')+'!Drest!N ['+texToIDL('\AA')+']', $
           ytitle = '!18f!X', /nodat, pos = [0.1,0.1,0.5,0.55], /noer, $
           xthick = 4, ythick = 4, charthick = 3, charsize = 1, yminor = 2, /ysty, /xsty
     oplot, lambda1, smooth(outer1, 2, /nan), col = 'ffa500'x, thick = 2
     oplot, lambda1, smooth(inter1, 2, /nan), col = '00a500'x, thick = 3
     oplot, lambda1, smooth(inner1, 2, /nan), col = 255      , thick = 3
     oplot, lambda1, err1, thick = 2, col = 100
     legend, /top, /right, box = 0, $
             ['!18R!X<0.25 !18R!X!De!N', 'fit', $
              '0.25!18R!X!De!N<!18R!X<0.75!18R!X!De!N', 'fit', $
              '0.75!18R!X!De!N<!18R!X<1.50!18R!X!De!N', 'fit'], $
             col = [255,100,'00a500'x,'004400'x,'ffa500'x,'ff0000'x], $
             pspacing = 1, thick = [2,4,2,4,2,4], $
             charsize = 0.7, charthick = 2, linesty = 0
     legend, /top, /left, 'PA1', charsize = 1.1, charthick = 4, box = 0, /clear
     
     plot, lambda2, inner2, $
           yran = [0,500], $
           xtitle = greek('lambda')+'!Drest!N ['+texToIDL('\AA')+']', $
           /nodat, pos = [0.55,0.1,0.95,0.55], /noer, $
           xthick = 4, ythick = 4, charthick = 3, charsize = 1, yminor = 2, /ysty, /xsty
     oplot, lambda2, smooth(outer2, 2, /nan), col = 'ffa500'x, thick = 2
     oplot, lambda2, smooth(inter2, 2, /nan), col = '00a500'x, thick = 3
     oplot, lambda2, smooth(inner2, 2, /nan), col = 255      , thick = 3
     oplot, lambda2, err2, thick = 2, col = 100
     legend, /top, /left, 'PA2', charsize = 1.1, charthick = 4, box = 0, /clear
     legend, /top, /right, box = 0, $
          ['!18z!X='+string(pa1b.Z, f = '(F4.2)'), $
           '!18R!X!De!N='+string(hlrObs, f = '(F4.2)')+'"='+string(hlrKpc, f = '(F3.1)')+' kpc'], $
          charsize = 1, charthick = 3
     device, /close

     device, filename =  outdir+'/'+fieldnames[counter]+'_'+ids[counter]+'_B_2d.eps', $
             /col, /encap, /decomp, bits_per_pixel = 8, $
             xsize = 10, ysize = 8, /in
     !p.multi = [0,2,3]
     plot, blam1 / (1 + pa1b.Z), findgen(60), /xsty, /nodat, $
           ytitle = 'spatial pix [2D spectrum]', charthick = 4, title = 'PA1'
     cgimage, pa1b.SPEC2D, stretch = 2, /over
     plot, blam2 / (1 + pa2b.Z), findgen(60), /xsty, /nodat, $
           charthick = 4, title = 'PA2'
     cgimage, pa2b.SPEC2D, stretch = 2, /over
     plot, blam1 / (1 + pa1b.Z), findgen(60), /xsty, /nodat, $
           ytitle = 'spatial pix [extraction zones]', charthick = 4
     cgimage, pa1b.EXTRACT_ZONES, stretch = 2, /over
     plot, blam2 / (1 + pa2b.Z), findgen(60), /xsty, /nodat, $
           charthick = 4
     cgimage, pa2b.EXTRACT_ZONES, stretch = 2, /over
     plot, blam1 / (1 + pa1b.Z), findgen(60), /xsty, /nodat, $
           ytitle = 'spatial pix [model resid.]', charthick = 4, xtitle = 'rest wavelength'
     cgimage, pa1b.SPEC2D - fan(pa1b.F_TOT * 3 * pa1b.RE, 60) * pa1b.MODEL_TRACE, stretch = 2, /over, $
              minval = min(pa1b.SPEC2D), maxval = max(pa1b.SPEC2D)
     plot, blam2 / (1 + pa2b.Z), findgen(60), /xsty, /nodat, $
           charthick = 4, xtitle = 'rest wavelength'
     cgimage, pa2b.SPEC2D - fan(pa2b.F_TOT * 3 * pa2b.RE, 60) * pa2b.MODEL_TRACE, stretch = 2, /over, $
              minval = min(pa2b.SPEC2D), maxval = max(pa2b.SPEC2D)
     device, /close

     device, filename =  outdir+'/'+fieldnames[counter]+'_'+ids[counter]+'_R_2d.eps', $
             /col, /encap, /decomp, bits_per_pixel = 8, $
             xsize = 10, ysize = 8, /in
     !p.multi = [0,2,3]
     plot, rlam1 / (1 + pa1r.Z), findgen(60), /xsty, /nodat, $
           ytitle = 'spatial pix [2D spectrum]', charthick = 4, title = 'PA1'
     cgimage, pa1r.SPEC2D, stretch = 2, /over
     plot, rlam2 / (1 + pa2r.Z), findgen(60), /xsty, /nodat, $
           charthick = 4, title = 'PA2'
     cgimage, pa2r.SPEC2D, stretch = 2, /over
     plot, rlam1 / (1 + pa1r.Z), findgen(60), /xsty, /nodat, $
           ytitle = 'spatial pix [extraction zones]', charthick = 4
     cgimage, pa1r.EXTRACT_ZONES, stretch = 2, /over
     plot, rlam2 / (1 + pa2r.Z), findgen(60), /xsty, /nodat, $
           charthick = 4
     cgimage, pa2r.EXTRACT_ZONES, stretch = 2, /over
     plot, rlam1 / (1 + pa1r.Z), findgen(60), /xsty, /nodat, $
           ytitle = 'spatial pix [model resid.]', charthick = 4, xtitle = 'rest wavelength'
     cgimage, pa1r.SPEC2D - fan(pa1r.F_TOT * 3 * pa1r.RE, 60) * pa1r.MODEL_TRACE, stretch = 2, /over, $
              minval = min(pa1r.SPEC2D), maxval = max(pa1r.SPEC2D)
     plot, rlam2 / (1 + pa2r.Z), findgen(60), /xsty, /nodat, $
           charthick = 4, xtitle = 'rest wavelength'
     cgimage, pa2r.SPEC2D - fan(pa2r.F_TOT * 3 * pa2r.RE, 60) * pa2r.MODEL_TRACE, stretch = 2, /over, $
              minval = min(pa2r.SPEC2D), maxval = max(pa2r.SPEC2D)
     device, /close

     !p.multi = 0
     
;     spawn, 'open '+strcompress(outdir+'/'+string(fix(counter))+'_info.eps', /rem)
;     stop
     counter += 2
  endwhile
  
  

end
