function metis_mres_sn, sig, time, band, site, aomode, OVHD=ovhd, TEL_EM=tel_em, $
        EXTENDED=extended, RH=rh

;*************************************************************************************
;* SENSITIVITY CODE FOR MEDIUM-RESOLUTION SLIT SPECTROSCOPY MODE OF METIS             * 
;* METIS IS THE MID-IR INSTRUMENT FOR THE EUROPEAN ELT                                *
;* Written by Sarah Kendrew, Leiden Observatory, 2009                                 * 
;* Questions? Email me: kendrew@strw.leidenuniv.nl                                    *
;*************************************************************************************
; 
; USE OF THIS CODE
; ================
; Anyone is welcome to use this code to look at the likely performance of a mid-IR instrument on a 42-m telescope, and I'm happy to 
; advise on appropriate usage. A paper will be published at the 2010 SPIE conference, so if your work with this code results in a 
; publication, please cite Kendrew et al, 2010 in press. If you have a great idea to make novel use of this code and would like to 
; collaborate, let me know.
;
; DEPENDENCIES
; ============
; To run this code you need the following additional files:
; - Atmospheric profiles paranal_lm_3000.dat, paranal_n_3000.dat, high_and_dry_lm_3000.dat, high_and_dry_n_3000.dat
; - Parameter files metis_pre_param.pro (pre-optics params), metis_sens_param.pro (general telescope params), 
;   metis_img_param.pro
; - Encircled energy lookup tables, EnsquaredEnergy_AO_LMN_V10.dat for SCAO and atlas_eeav.dat for LTAO
; - DQE calculation function: metis_dqe.pro
; - Window background calculation: metis_winbgr.pro
; - Instrument throughput calculation: metis_thru.pro
; - Grism profiles: grisms.dat in profiles/ directory

; INPUTS
; ======
; sig          input signal, in mJy. the code will calculate the S/N of this signal as well as sensitivity of the IFU
; time         total clock time, in seconds. for e.g. sensitivity for S/N in 1 hour, set time=3600.
; band         observation band, options: 'l', 'm' or 'n'. this should correspond to lref.
; site         site of observation. options: 'low' or 'high'
; aomode       mode of adaptive optics used. options: 'scao' for single-conjugate natural guide star AO, 'ltao' for laser tomography AO
; 
; OPTIONAL KEYWORDS
; =================
; OVHD         overhead. if this keyword is set, a 20% overhead will be applied, i.e. time = time*0.8.
; TEL_EM       telescope emissivity. default is 0.1 (10%). this is an important parameter in the calculation, vary it to see the effect.
; EXTENDED     extended source. see the accompanying tech note for a detailed explanation of how this works.
; RH           helf-light radius (only valid if keyword EXTENDED also set!)
; 
; OUTPUT
; ======
; The output is a 2-element vector. First element is the continuum sensitivity in mJy, the second is the line sensitivity
; at a reference wavelength lref specified in the initial code checks, in W/m2. Default reference lines are 
; L-band:   Br-alpha @ 4.05 um
; M-band:   CO       @ 4.66 um
; N-band:   [NeII]   @ 12.81 um


;*******************************************************
; INITIAL CHECKS AND PARAMETERS
;*******************************************************

; check operating system
os = !version.os_family
if os eq 'Windows' then begin
  graphics = 'win'
endif else begin
  graphics = 'x'
endelse

print, 'Your operating system is: ', os

if not keyword_set(tel_em) then tel_em=0.1

if keyword_set(extended) then begin
    if not keyword_set(rh) then begin
       print, 'Source type: EXTENDED, CONSTANT SURFACE BRIGHTNESS (mJy/arcsec^2)'
    endif else begin
      print, 'Source type: EXTENDED (mJy)
    endelse
endif else begin
  print, 'Source type: POINT (mJy)'
endelse

if keyword_set(rh) and not keyword_set(extended) then begin
  print, 'Cant have half-light radius for a point source!'
  return, 0
endif


; read in common parameters
@metis_sens_param
; read in pre-optics parameters
@metis_pre_param
;read in imager-specific parameters
@metis_img_param


if site eq 'low' then tel_t=286.3 else tel_t=270.
win_temp=tel_t

; set clock
t=systime(1)

; read in filters file. for now these are single-band filters, not a full set of "scientifically interesting"
; ones. just read in effective filter width for now. given in um.

grismfile=profiles_dir+'grisms.dat'
readcol, grismfile, f='d,d,d', lmin, lmax, grismrp, /silent

; BAND CHECKS
case band of
  'l': begin
    f1=10.286d
    f2=10.286d
    rd_noise=lm_rd_noise    ; as specified in parameter file [e-/pixel]
    dark=lm_dark            ; pick dark current value for LM band
    pix_scale=img_pix_scale_lm    ; dito pixel scale
    lmin=lmin[0]
    lmax=lmax[0]
    res=grismrp[0]
    lref=4.05
    lslit=3.            ; wavelength for slit width reference - from optical design
    pix_size=lm_pix_size
    well=lm_well
    fileband='lm'       ; for the atmosphere filename
    end

  'm': begin
    f1=10.286d
    f2=10.286d
    rd_noise=lm_rd_noise    ; as specified in parameter file [e-/pixel]
    dark=lm_dark
    pix_scale=img_pix_scale_lm
    lmin=lmin[0]
    lmax=lmax[0]
    lref=4.65
    lslit=3.            ; wavelength for slit width reference
    res=grismrp[0]
    pix_size=lm_pix_size
    well=lm_well
    fileband='lm'       ; for the atmosphere filename....
    end

  'n': begin
    f1=8.571d 
    f2=8.571d
    rd_noise=nq_rd_noise    ; as specified in parameter file
    dark=nq_dark
    pix_scale=img_pix_scale_nq
    lmin=lmin[3]
    lmax=lmax[3]
    res=grismrp[3]
    lref=12.81
    lslit=7.              ; wavelength for slit width reference
    pix_size=nq_pix_size
    well=nq_well
    ;dit=ndit
    fileband=band       ; for the atmsophere filename....
    end


  else: print, "Band not recognised"

endcase


; READ IN ATMOSPHERIC PROFILES
;------------------------------
lowfile=atm_dir+'paranal_'+fileband+'_3000.dat'
highfile=atm_dir+'high_and_dry_'+fileband+'_3000.dat'

if site EQ 'low' then begin
  readcol, lowfile, F='D,D,x,D',l, trans, em, /silent
  tel_t=286.3d
  name='paranal'
endif else if site EQ 'high' then begin
  readcol, highfile,F='D,D,x,D', l, trans, em, /silent
  tel_t=270.0d
  name='high_and_dry'
endif
print, 'atmosphere data read'


;perform some unit conversions on the emissivity values (from photons/s/um/m^2/arcsec^2 to W/cm^2/sr/um)
;em=em*(h*c/l)*4.25d10. UNITS WOMAN, UNITS! PAY ATTENTION!
em_l=em*(h*c/(l*1d-6))*1d-4*4.25d10

; clip these arrays to the region of interest
em_l=em_l[where((l gt lmin) and (l le lmax))]
trans=trans[where((l gt lmin) and (l le lmax))]
l=l[where((l gt lmin) and (l le lmax))]
;plot, l, trans

; calculate delta_l
print, 'resolving power = ', res
delta_l=l/res


;# pixels in S/N area: spectral sampling is 2px/dl, spatial calculated from pixel scale
;npix=2.*lref*1d-6*206265./42./(pix_scale*1d-3)

if keyword_set(extended) then begin
  if keyword_set(rh) then begin
    npix=(2*lslit*1d-6*206265/42./(pix_scale*1d-3))*(2*rh/(pix_scale*1d-3))
  endif else begin
    npix=(2*lslit*1d-6*206265/42./(pix_scale*1d-3))*(1./(pix_scale*1d-3))
  endelse
endif else begin
  npix=(lref*1d-6*206265./42./(pix_scale*1d-3))*(2*lslit*1d-6*206265/42./(pix_scale*1d-3))
endelse    

npix=ceil(npix)

print, 'npix = ', npix

;calculate S/N reference area, in arcsec^2. Is equal to width of slit (2l/D) x psf FWHM (1.22l/D)
sn_size=npix*(pix_scale*1d-3)^2

print, 'sn_size = ', sn_size, ' arcsec^2'

diff=l-lref
tmp=min(diff, line, /absolute)
print, 'wavelengh at position lines = ', l[line]

;look up EE using appropriate AO mode
if keyword_set(extended) then begin
  if keyword_set(rh) then ee=0.5 else ee=1. 
endif else begin
  if aomode eq 'ltao' then begin
    ee_file='ee_files/atlas_eeav.dat'
    readcol, ee_file, f='f,d,d,d', ee_snsize, ee_l, ee_m, ee_n, /silent
  endif else begin    
    ee_file='ee_files/EnsquaredEnergy_AO_LMN_V10.dat'
    readcol, ee_file, format='f,d,d,d', ee_snsize, ee_l, ee_m, ee_n, /silent
  endelse
  
  ee_list=dblarr(n_elements(ee_snsize))

  case band of
    'l': ee_list=ee_l
    'm': ee_list=ee_m
    'n': ee_list=ee_n
    ;  'q': ee_list=ee_q
    end

  ;but now sn_size and ee_snsize are vecors of different length -> need to lookup element by element
  diff=ee_snsize-sqrt(sn_size*1d6)
  tmp=min(diff, loc, /absolute)
  if aomode eq 'ltao' then ee=ee_list[loc]/100. else ee=ee_list[loc]
endelse

print, 'ref = ', sqrt(sn_size*1d6), ' mas'

print, 'ee = ', ee

; calculate telescope transmission
if (tel_em gt 0.1) then tel_thru=1.-tel_em else tel_thru=mir_ref^no_mir
print, 'telescope transmission = ', tel_thru

; calculate throughput
obsmode='img'
thru=metis_thru(obsmode, band, l)
print, 'mean throughput of instrument = ', mean(thru)
plot, l, thru, /ynozero, xtitle='micron', title=band


; calculate dqe
dqe=metis_dqe(l, band)
print, 'mean DQE in filter band = ', mean(dqe)

; calculate efficiency
eff=thru*tel_thru*mean(dqe)*pcg
print, 'mean total efficiency including DQE, pcg, telescope, and instrument = ', mean(eff)

; total conversion factor
conv=eff*delta_l*eff_area*lref/1.985d-19

; flux from object
fobj=sig*trans*3.*1d-19/lref^2    ; convert flux from mJy to [W/cm2/um]

; BACKGROUNDS
; telescope
tel_bgrsr=planck(l*1d4, tel_t)*1d-3/!dpi    ; in W/cm2/sr/um
tel_bgr=tel_bgrsr*tel_em/4.25d10             ; in W/cm2/arcsec2/um 

;plot, l, tel_bgr, /ylog, yrange=[1d-18, 1d-13]

; atmosphere
em_l=em_l/4.25d10*tel_thru

;window background
win_bgr=metis_winbgr(l, tel_t)

; add telescope and sky
sky=tel_bgr+em_l+win_bgr                         ; in W/cm2/um/arcsec2

; add in overhead to the clock time if specified
if keyword_set(ovhd) then time=time*0.8

; # detected electrons from sky & telescope
nskypers=sky*sn_size*conv      ; electrons per s per resolution element

dit=well/max(nskypers/npix)
frames=round(time/dit)
nsky=nskypers*dit              ; electrons per DIT in sn area
print, 'optimal dit = ', dit, ' s'
print, '# frames = ', frames
print, 'max # electrons from sky in sn area = ', max(nsky), ' e-/DIT'

; # detected electrons from object in s/n reference area
nobj=fobj*ee*conv*dit
sn=sqrt(frames)*nobj/sqrt(nobj+nsky+(npix*rd_noise^2)+(npix*dark*dit))


;calculate minimum detectable flux in 1 hr to s/n of 10.
minflx=10.*emp_fac*sqrt(nsky+(npix*rd_noise^2)+(npix*dit*dark))/(ee*trans*conv*dit*sqrt(frames))    ; in W/cm2/um
minsig=minflx*l^2/3./1d-19
minsig_line=minflx*delta_l*1d4

if keyword_set(extended) then begin
  if keyword_set(rh) then begin
    print, 'sensitivity at s/n 10 in 1 hr = ', minsig[line]*1000., ' microJy in half-light radius of ', rh, 'arsec'
    print, 'line sensitivity at s/n in 1 hr = ', minsig_line[line]*1d19, '*1d-19 W/m2'
  endif else begin
    print, 'sensitivity at s/n 10 in 1 hr = ', minsig*1000., ' microJy/arcsec^2'
    print, 'line sensitivity at s/n in 1 hr = ', minsig_line[line]*1d19, '*1d-19 W/m2/arcsec^2'
  endelse
endif else begin
  print, 'sensitivity at s/n 10 in 1 hr = ', minsig[line]*1000., ' microJy'
   print, 'line sensitivity at s/n 10 in 1 hr = ', minsig_line[line]*1d19, '*1d-19 W/m2'
endelse


out=[minsig[line], minsig_line[line]]

return, out
end


