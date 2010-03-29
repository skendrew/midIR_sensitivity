function metis_ifures, band, x, mode

; this function computes the delta_l for the IFU mode of METIS, according to the optical design
; parameters of the instrument
; INPUTS:
; - band:     band of observation, can be 'l', 'm', 'n' or 'q' (q not supported yet)
; - x   :     wavelength scale in micron
; - xrange:   the wavelength range of the observation band in micron. vector of 2 elements.
; written on 26/08/09 by S. Kendrew, kendrew@strw.leidenuniv.nl

@metis_sens_param

if (band eq 'l') then begin
    infile=profiles_dir+'metis_ifures_lm.dat'
    readcol, infile, f='d,d', wave, res, /silent
    xrange=[3., 4.2]
endif else if (band eq 'm') then begin
    infile=profiles_dir+'metis_ifures_lm.dat'
    readcol, infile, f='d,d', wave, res, /silent
    xrange=[4.5, 5.5]
endif else begin
    infile=profiles_dir+'metis_ifures_n.dat'
    readcol, infile, f='d,d', wave, res, /silent
    xrange=[7., 14.]
endelse

; res is given in units of 1000 so do x 1000.
res=res*1000.

; interpolate onto wavelength grid
resint=interpol(res, wave, x)
dx=x/resint

;plot, x, res
;oplot, x, resint, col=255

if (mode eq 'img') then begin
loc=where((x lt min(wave)) or (x gt max(wave)))
if (loc ne -1) then begin
resint[loc]=mean(res)
dx[loc]=mean(dx[where(x le max(wave)) or (x ge min(wave))])
endif
endif

return, dx


end