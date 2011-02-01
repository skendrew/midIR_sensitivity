function metis_winbgr, x, temp

; returns the entrance window background vs. wavelength in W/cm2/um/arcsec2
@metis_sens_param

n=2.4028  ;real part of refractive index @ 10.6 um
r=(n-1)/(n+1)
win_thick=1.2  ;window thickness [cm]
;read in window absortion coefficient and calculate emissivity
absfile=profiles_dir+'znse_abs.dat' 
print, absfile
readcol, absfile, format='d,d', abs_l, alpha, skipline=4, /silent
win_em=((1-r^2)*(1-exp(-alpha*win_thick)))/(1-(r^2*exp(-alpha*win_thick)))        ;r and d defined in metis_thermal_pre_param
win_em=interpol(win_em, abs_l, x)
;win_em[where(l gt 20.)]=1.                ;set the emissivity to level off after 20 micron at 100%

;calcuate window background
bgr=win_em*planck(x*1d4, temp)*1d-3/!dpi
bgr=bgr/4.25d10      ; convert from sr to arcsec2

return, bgr

end