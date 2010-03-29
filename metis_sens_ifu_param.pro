
;pro metis_sens_ifu_param
;parameters file for sensitivity calculations

;provisional temperatures
hspec_temp_lm=65.    ; Kelvin
hspec_temp_n=25.

;SPO path
spo_no_surf=3        ;7      ;# of "normal" reflective optics, e.g. re-imaging & fold mirrors
spo_ref_surf=0.99
;spodic_trans=0.85
slicer_ref=0.99    ;to be multiplied by slicing loss, calculated in code
slices_lm=24.
slices_n=15.

;SMO path
smo_no_surf=7        ;6      ;# "normal" surfaces
smo_ref_surf=0.99
no_grating=1          ;2
grating_eff=0.8^no_grating
pix_spat=2.d        ;no. pixels per spatial res element
pix_spec=2.5        ;dito spectral
