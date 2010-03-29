function metis_thru, mode, band, x, POL=pol

; this function calculates the throughput of METIS for a given mode of observation using some
; custom functions for particular components, developed for metis_thermal

; the output can feed into the sensitivity model

; INPUTS:
; mode = observation mode, 'img', 'lspec' or 'hspec'
; band = band of observation, 'l', 'm', 'n' or 'q'
; x = the wavelength grid on which to calculate the output

; KEYWORDS:
; pol = sets polarimetry

;************************************
;*CONSTANTS                         *
;************************************
;pre_surf=9      ;no. of reflective optics in the pre-optics
;pre_ref=0.99    ;reflectivity per pre-optics mirror
;pre_thru=pre_ref^pre_surf
;ar_factor=1.3
;thru=1.

;wollaston_thru=0.9      ; throughput of the wollaston prism in METIS

; CONSTANTS ALL READ IN VIA PARAMETER FILES
@metis_sens_param
@metis_pre_param
@metis_thermal_ifu_param
@metis_thermal_img_param

;**********************************************************************
;* WINDOW TRANSMISSION CURVE  also used for prism transmission for IFU*
;**********************************************************************
winfile=profiles_dir+'znse_trans.dat'
readcol, winfile, format='d,d,d', window_l, window_trans4, window_trans12, skipline=3, /silent  ;read in transmission spectrum
window_trans=interpol(window_trans12, window_l, x)          ;interpolate onto l grid
for i=0l, n_elements(x)-1 do begin
  if (x[i] lt 14.) then window_trans[i]=0.72
  if (window_trans[i] lt 0.) then window_trans=0. 
endfor

window_trans=window_trans*ar_factor                       ;apply anti-reflection improvement 

prism_trans=interpol(window_trans4, window_l, x)
for i=0l, n_elements(x)-1 do begin
  if (x[i] lt 14.) then prism_trans[i]=0.72
  if (prism_trans[i] lt 0.) then prism_trans[i]=0.
endfor
prism_trans=prism_trans*ar_factor

;*******************************
;* WFS DICHROIC                *
;*******************************

wfsdic_trans=metis_wfsdic(x)

;*******************************
;* BAND-SPLITTING DICHROIC     *
;*******************************

bandic_profs=metis_bandic(x, band)
bandic_trans=bandic_profs[*,0]

;*******************************
;* IMAGER THROUGHPUT           *
;*******************************

if mode eq 'img' then begin
  
  thru=pre_thru*img_thru*window_trans*wfsdic_trans*bandic_trans
  if keyword_set(pol) then thru=thru*wollaston_thru
;oplot, x, thru, col=255
  
;*******************************
;* LSPEC THROUGHPUT            *
;*******************************
 
  
endif else if mode eq 'lspec' then begin

  grism_trans=window_trans
  thru=pre_thru*img_thru*window_trans*wfsdic_trans*bandic_trans*grism_trans
  ;oplot, x, thru, col=255

;*******************************
;* HSPEC THROUGHPUT            *
;*******************************

endif else if mode eq 'hspec' then begin
  sl_loss=dblarr(n_elements(x))
    for i=0l, n_elements(x)-1 do begin
      if (x[i] le 7.) then begin
        sl_loss[i]=0.85d*(x[i]/3.7d)^0.25
      endif else begin
      if (x[i] gt 15.) then begin
        sl_loss[i]=0.88d*(x[i]/18.d)^0.25
      endif else begin
        sl_loss[i]=0.88d*(x[i]/9.d)^0.25
      endelse
     endelse
endfor

  spo_thru=spo_ref_surf^spo_no_surf*wfsdic_trans*window_trans*bandic_trans
  smo_thru=smo_ref_surf^smo_no_surf*grating_eff*prism_trans*slicer_ref*sl_loss
  thru=spo_thru*smo_thru
  ;oplot, x, thru, col=255 

endif

return, thru



end
 