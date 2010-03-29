; Parameter file for the thermal properties of the pre-optics


;temperatures
;win_temp=tel_t  ;entrance window temperature
wfs_em=0.15    ;emissivity of the WFS dichroic
pre_temp=85.    ;pre-optics temperature
wfs_temp=pre_temp  ;dichroic temperature

;parameters for window emissivity calcs
n=2.4028  ;real part of refractive index @ 10.6 um
r=(n-1)/(n+1)
win_thick=1.2  ;window thickness [cm]

;improvement factor in window transmission from anti-reflection coating - TBD
ar_factor=1.3

;re-imaging optics spec
pre_surf=9      ;no. of reflective optics in the pre-optics
pre_ref=0.99    ;reflectivity per pre-optics mirror
;pre_em=1.-pre_ref
;pre_emtot=9*pre_em*(1-(4*pre_em))    ;the emissivity for 9 cold mirrors
pre_thru=pre_ref^pre_surf

