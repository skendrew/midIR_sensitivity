
;pro metis_ifu_param
;parameters file for sensitivity calculations
;
;provisional temperatures
hspec_temp_lm=80.    ; Kelvin
hspec_temp_n=25.

;SPO path
spo_no_surf=7      ;# of "normal" reflective optics, e.g. re-imaging & fold mirrors
spo_ref_surf=0.99
slicer_ref=0.99    ;to be multiplied by slicing loss, calculated in code
slices_lm=24.
slices_n=15.

;SMO path
smo_no_surf=6      ;# "normal" surfaces
smo_ref_surf=0.99
no_grating=2
grating_eff=0.8^no_grating

;sampling parameters
samp_spat_lm=9      ;mas/px
samp_spec_lm=7.2    ;mas/px
samp_spat_n=22      ;mas/px
samp_spec_n=17.6    ;mas/px

