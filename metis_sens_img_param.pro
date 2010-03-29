
;pro metis_sens_img_param
;parameters file for sensitivity calculations for the IMAGER



no_surf=15				;number of cold surfaces. 4 pre-optics+3 de-rotator+3 collimator+2 folding mirrors+3 camera
ref_surf=0.99			;reflectance per cold surface
window_trans=0.95
no_dichroic=1
dichroic_trans=0.8^no_dichroic
;pix_spat=3.d				;no. pixels per spatial res element. assume nyquist sampled.
thru=ref_surf^no_surf*window_trans*dichroic_trans

pix_scale_lm=8.6		; mas/pixel
pix_scale_nq=17.2		; mas/pixel
;pix_spec=2.				; for testing - REMOVE

print, 'imager throughput = ', thru

;end

