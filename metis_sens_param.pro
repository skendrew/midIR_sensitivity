
;pro metis_sens_param
; COMMON parameters to all modes and bands
;instrument parameters moved to metis_sen_img_param and metis_sens_ifu_param.pro

; file locations
def_dir='c:/metis/code_public/'
atm_dir=def_dir+'atm/'      ; location of atmospehric profiles
ee_dir=def_dir+'ee_files/'  ; location of encircled energy tables
profiles_dir=def_dir+'profiles/'  ; location of miscellaneous input data, mainly optical profiles e.g. throughputs, resolutions etc

;general constants
h=6.63d-34 ;Js
k_b=1.38d-23 ;J/K
c=3d8   ;m/s

;telescope parameters
prim_area=!dpi*(21)^2		;telescope primary area
obs=9.2*prim_area/100		;obscuration
eff_area=(prim_area-obs)*1e4		;effective area in cm^2
mir_ref=0.98d			;mirror reflectivity
no_mir=5				;no. telescope mirrors

;detector description - useful for output reference
lm_det='Teledyne Hawaii-II detector for LM-band, 18 um pitch, 2k x 2k px'
nq_det='Raytheon Aquarius detector for N band, 30 um pitch, 1k x 1k px'

; polarimetry: half-wave plate throughput
;pol_thru=0.64

;detector parameters (per pixel!)
lm_rd_noise=20.			;read noise [e-/frame]
lm_dark=0.5				;dark current [e-/s]
lm_well=6.5e4        ;well depth [e-]
nq_rd_noise=1d3				;read noise [e-/frame]
nq_dark=2d3					;dark current [e-/s]
nq_well=7d6          ;50% well depth (low gain) [e-]

lm_pix_size=0.0018d  ; pixel size [cm]
nq_pix_size=0.003d    


pcg=1.					;photo-conductive gain
gdb=1.					;gain dispersion beta
qe_scale=0.7			;QE scaling factor
emp_fac=sqrt(2)			;empirical factor


;end