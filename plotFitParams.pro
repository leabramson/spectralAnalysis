;; Plot output from GETFITPARAMS.pro (FITS FORMAT!)

pro plotFitParams, resultsFile, $
                   PARAM = param, OUTPUT = output

  data   = mrdfits(resultsFile, 1)
  files  = data.FILE
  nfiles = n_elements(files)
  
  ;; Parse the files to get objects
  dir  = []
  targ = []
  pa1targlist = strarr(2, nfiles)  ;;  [Red, Blue] fits files
  pa2targlist = strarr(2, nfiles)
  for ii = 0, nfiles - 1 do begin
     slash = strpos(files[ii], '/')
     dir  = [dir, strmid(files[ii], 0, slash)]
     targ = [targ, strmid(files[ii], slash+1, $
                         strpos(files[ii], '_', /reverse_search) - slash - 1)]
     pa1targlist[*, ii] = [dir[ii]+'/'+targ[ii]+'_1_B.fits', $
                           dir[ii]+'/'+targ[ii]+'_1_R.fits']
     pa2targlist[*, ii] = [dir[ii]+'/'+targ[ii]+'_2_B.fits', $
                           dir[ii]+'/'+targ[ii]+'_2_R.fits']
  endfor

  ;; Get physical parameters
  redshifts = fltarr(nfiles)
  obsRadius  = fltarr(nfiles)
  physRadius = fltarr(nfiles)
  normRadii  = fltarr(3,nfiles,2)
  pscale = 0.13                 ;; Asec/ pix in WFC3 IR Channel
  for ii = 0, nfiles - 1 do begin
     bluDat1 = mrdfits(pa1targlist[0,ii], 1)
     redDat1 = mrdfits(pa1targlist[1,ii], 1)
     bluDat2 = mrdfits(pa2targlist[0,ii], 1)
     redDat2 = mrdfits(pa2targlist[1,ii], 1)
        
     innOff1 = median(bluDat1.INNERUP - bluDat1.TRACE)
     intOff1 = median(bluDat1.INTERUP - bluDat1.TRACE)
     outOff1 = median(bluDat1.OUTERUP - bluDat1.TRACE)

     innOff2 = median(bluDat2.INNERUP - bluDat2.TRACE)
     intOff2 = median(bluDat2.INTERUP - bluDat2.TRACE)
     outOff2 = median(bluDat2.OUTERUP - bluDat2.TRACE)

     redshifts[ii] = bluDat1.Z

     kpc    = 1. / zang(1.0, redshifts[ii]) ;; kpc / asec
     kpcPix = kpc * pscale                  ;; kpc / pix

     reCirc         = sqrt(bluDat1.RE * bluDat2.RE)  ;; circularized Re
     obsRadius[ii]  = reCirc * pscale ;; re in asec 
     physRadius[ii] = reCirc * kpcPix ;; re in kpc
     normRadii[*,ii,0] = [0, mean([innOff1,intOff1]),mean([intOff1,outOff1])] / reCirc
     normRadii[*,ii,1] = [0, mean([innOff2,intOff2]),mean([intOff2,outOff2])] / reCirc
  endfor

  case param of
     'age'  : begin
        yr = [7,10]
        ttl = 'Age Gradients'
     end
     'HpsEW': begin
        yr = [0,100]
        ttl = 'EW(H'+greek('alpha')+') Gradients'
     end
     'OIII' : begin
        yr = [0,100]
        ttl = 'EW([OIII]) Gradients'
     end
  endcase

  set_plot, 'PS'
  device, filename = output, $
          /col, /encap, /decomp, bits_per_pix = 8, $
          xsize = 4, ysize = 10, /in
  multiplot, [1,nfiles], ygap = 0.005
  for ii = 0, nfiles - 1 do begin
     nr1 = reform(normRadii[*,ii,0])
     pr1 = nr1 * physRadius[ii]
     nr2 = reform(normRadii[*,ii,1])
     pr2 = nr2 * physRadius[ii]
     
     mids1 = reform(data[ii].VALUES[*,0])
     mids2 = reform(data[ii].VALUES[*,1])
     his1  = reform(data[ii].HI_ERR[*,0])
     his2  = reform(data[ii].HI_ERR[*,1])
     los1  = reform(data[ii].LO_ERR[*,0])
     los2  = reform(data[ii].LO_ERR[*,1])

     if ii eq 0 then begin
        plot, nr1, mids1, /nodat, $
              xran = [0,1.5], yran = yr, $
              charsize = 1.25, charthick = 4, $
              xthick = 5, ythick = 5, yminor = 2, $
              xtickname = replicate( ' ',60), /ynoz, $
              title = ttl, /xsty
     endif else if ii gt 0 AND ii lt nfiles - 1 then begin
        plot, nr1, mids1, /nodat, $
              xran = [0,1.5], yran = yr, $
              charsize = 1.25, charthick = 4, $
              xthick = 5, ythick = 5, yminor = 2, $
              xtickname = replicate( ' ',60), /ynoz , /xsty
     endif else begin
        plot, nr1, mids1, /nodat, $
              xran = [0,1.5], yran = yr, $
              charsize = 1.25, charthick = 4, $
              xthick = 5, ythick = 5, yminor = 2, $
              xtitle = '!18r / r!X!De!N', /xsty
     endelse
     xxx1 = [nr1, reverse(nr1)] < !X.CRANGE[1]
     xxx2 = [nr2, reverse(nr2)] < !X.CRANGE[1]
     yyy1 = ([mids1 + his1, reverse(mids1 - los1)] < !Y.CRANGE[1]) > !Y.CRANGE[0]
     yyy2 = ([mids2 + his2, reverse(mids2 - los2)] < !Y.CRANGE[1]) > !Y.CRANGE[0]
     polyfill, xxx2, yyy2, /line_fill, $
               col = 255, thick = 1, spacing = 0.025, $
               orien = 45
     polyfill, xxx1, yyy1, /line_fill, $
               col = 'ff0000'x, thick = 1, spacing = 0.025, $
               orien = -45

     legend, /top, /left, box = 0, $
             dir[ii]+' '+targ[ii], $
             charsize = 1.0, charthick = 3
     legend, /top, /right, box = 0, $
             '!18r!X!De, circ!N = '+string(physRadius[ii], $
                                           f = '(F4.1)')+' kpc', $
             charsize = 1.0, charthick = 3
     
     multiplot

  endfor

  device, /close
  spawn, 'gv '+output
  
  stop
  
end
