
;pro metis_img_param
;parameters file for sensitivity calculations for the IMAGER

;
;
;no_surf=15				;number of cold surfaces. 4 pre-optics+3 de-rotator+3 collimator+2 folding mirrors+3 camera
;ref_surf=0.99			;reflectance per cold surface
;window_trans=0.95
;no_dichroic=1
;dichroic_trans=0.8^no_dichroic
;;pix_spat=3.d				;no. pixels per spatial res element. assume nyquist sampled.
;thru=ref_surf^no_surf*window_trans*dichroic_trans
;
;pix_scale_lm=8.6		; mas/pixel
;pix_scale_nq=17.2		; mas/pixel
;;pix_spec=2.				; for testing - REMOVE
;
;print, 'imager throughput = ', thru

img_no_surf=7       ;number of cold surfaces after pre-optics
img_ref_surf=0.99     ;reflectance per cold surface
img_thru=img_ref_surf^img_no_surf
img_em=1.-img_ref_surf
img_temp_lm=85.
img_temp_n=25.
;FOR TESTING
;img_pix_scale_lm=10.
;###
img_pix_scale_lm=8.6    ; mas/pixel
img_pix_scale_nq=17.2   ; mas/pixel



;end

