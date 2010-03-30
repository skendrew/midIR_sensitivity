function metis_img_sn, sig, time, band, site, aomode, OVHD=ovhd, TEL_EM=tel_em, $
          EXTENDED=extended, RH=rh

;*****************************************************************************************
;* SENSITIVITY CODE FOR IMAGER MODE OF METIS, THE MID-IR INSTRUMENT FOR THE EUROPEAN ELT *
;* Written by Sarah Kendrew, Leiden Observatory, 2009                                    * 
;* Questions? Email me: kendrew@strw.leidenuniv.nl                                       *
;*****************************************************************************************
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
; - Atmospheric profiles paranal_lm_100.dat, paranal_n_100.dat, high_and_dry_lm_100.dat, high_and_dry_n_100.dat
; - Parameter files metis_pre_param.pro (pre-optics params), metis_sens_param.pro (general telescope params), 
;   metis_imgparam.pro (imager-specific parameters)
; - Encircled energy lookup tables, EnsquaredEnergy_AO_LMN_V10.dat for SCAO and atlas_eeav.dat for LTAO
; - DQE calculation function: metis_dqe.pro
; - Window background calculation: metis_winbgr.pro
; - Instrument throughput calculation: metis_thru.pro
; - Filter specifications: filters.dat in profiles/ directory. These are based on ISAAC and VISIR filters on VLT.

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
; RH           half-light radius of an extended object (in arcsec)
; 
; OUTPUT
; ======
; - the minimum detectable signal with a S/N of 10. within the time you've specified in the input in mJy


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

; read in common parameters
@metis_sens_param
; read in pre-optics parameters
@metis_pre_param
;read in imager-specific parameters
@metis_img_param

if site eq 'low' then tel_t=286.3 else tel_t=270.
win_temp=tel_t

if not keyword_set(tel_em) then tel_em=0.10

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

; set clock
t=systime(1)

; read in filters file. for now these are single-band filters, not a full set of "scientifically interesting"
; ones. just read in effective filter width for now. given in um.

filtersfile=profiles_dir+'filters.dat'
readcol, filtersfile, format='x,d,d', lc, dl, skipline=3, /silent


; BAND CHECKS
case band of
  'l': begin
    f1=10.286d          
    rd_noise=lm_rd_noise    ; as specified in parameter file [e-/pixel]
    dark=lm_dark        ; pick dark current value for LM band
    pix_scale=img_pix_scale_lm    ; dito pixel scale
    delta_l=dl[0]       ; dito band width
    lc=lc[0]
    pix_size=lm_pix_size
    well=lm_well
    fileband='lm'       ; for the atmosphere filename...
    end

  'm': begin
    f1=10.286d          ; formerly 9.73
    rd_noise=lm_rd_noise    ; as specified in parameter file [e-/pixel]
    dark=lm_dark
    pix_scale=img_pix_scale_lm
    delta_l=dl[1]       ; dito band width
    lc=lc[1]
    pix_size=lm_pix_size
    well=lm_well
    fileband='lm'       ; for the atmosphere filename....
    end

  'n': begin
    f1=8.571d         ; formerly 6.67
    rd_noise=nq_rd_noise    ; as specified in parameter file
    dark=nq_dark
    pix_scale=img_pix_scale_nq
    delta_l=dl[2]     
    lc=lc[2]
    pix_size=nq_pix_size
    well=nq_well
    fileband=band       ; for the atmsophere filename....
    end

  else: print, "Band not recognised"

endcase


; READ IN ATMOSPHERIC PROFILES
;------------------------------
resfile=100
lowfile=atm_dir+'paranal_'+fileband+'_'+strcompress(string(resfile),/remove_all)+'.dat'
highfile=atm_dir+'high_and_dry_'+fileband+'_'+strcompress(string(resfile),/remove_all)+'.dat'

if site EQ 'low' then begin
  readcol, lowfile, F='D,D,x,D',l, trans, em, /silent
  tel_t=286.d
  name='paranal'
endif else if site EQ 'high' then begin
  readcol, highfile,F='D,D,x,D', l, trans, em, /silent
  tel_t=270.0d
  name='high_and_dry'
endif
if not keyword_set(silent) then begin
  print, 'atmosphere data read'
endif

;perform some unit conversions on the emissivity values (from photons/s/um/m^2/arcsec^2 to W/cm^2/sr/um)
em_l=em*(h*c/(l*1d-6))*1d-4*4.25d10

; clip these arrays to where trans is non-zero
  em_l=em_l[where((l gt (lc-(delta_l/2))) and (l le lc+(delta_l/2)))]
  trans=trans[where((l gt (lc-(delta_l/2))) and (l le lc+(delta_l/2)))]
  l=l[where((l gt (lc-(delta_l/2))) and (l le lc+(delta_l/2)))]


; effective transmission in this filter band is:
transeff=mean(trans)
if not keyword_set(silent) then begin
print, 'effective atmospheric transmission in filter = ', transeff
endif

; extended sources: calculate the surface flux per arcsec^2
if keyword_set(extended) then begin
  sig=sig/!dpi/rh^2    
  if not keyword_set(silent) then begin
    print, 'Extended source: flux/area = ', sig, ' mJy/arcsec^2'
  endif
endif

;calculate S/N reference area, in arcsec^2.
;option 1 extended and rh set: define as square of side 2*half-light radius
;option 2 extended, const surf brightness: 1 arcsec^2
;option 3 point source, define as square of side lambda/D

;calculate # pixels in S/N area using pixel scale
;for the point source case we want to enforce a strictly square geometry, as 
;here the EE is most important, hence slightly complicated-looking formula  
if keyword_set (extended) then begin
  if keyword_set(rh) then npix=ceil(!dpi*rh^2/(pix_scale*1d-3)^2) else npix = ceil(1./(pix_scale*1d-3)^2) 
endif else npix=(ceil((lc*1d-6)/42.*206265*1d3/pix_scale))^2

print, 'number of px in s/n area = ', npix

if keyword_set(extended) then begin
  if keyword_set(rh) then sn_size=!dpi*rh^2 else sn_size=1. 
endif else sn_size=npix*(pix_scale*1d-3)^2

print, 'size of s/n area = ', sn_size, ' arcsec^2'

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
if keyword_set(pol) then tel_thru=tel_thru*pol_thru
print, 'telescope transmission = ', tel_thru

; calculate throughput
obsmode='img'
if keyword_set(pol) then thru=metis_thru(obsmode, band, l, /pol) else thru=metis_thru(obsmode, band, l)
print, 'mean throughput of instrument = ', mean(thru)

; calculate dqe
dqe=metis_dqe(l, band)
print, 'mean DQE in filter band = ', mean(dqe)

; calculate efficiency
eff=mean(thru)*tel_thru*mean(dqe)*pcg
print, 'total efficiency including DQE, pcg, telescope, and instrument = ', eff

; total conversion factor
conv=eff*delta_l*eff_area*lc/1.985d-19

print, ' conv = ', conv

; flux from object
fobj=sig*transeff*3.*1d-19/lc^2    ; convert flux from mJy to [W/cm2/um] or .../arcsec^2 for extended

; for extended sources we need to calculate the total flux over the sn area as is resolved
if keyword_set(extended) then fobj2=fobj*sn_size


; BACKGROUNDS
; telescope
tel_bgrsr=planck(l*1d4, tel_t)*1d-3/!dpi    ; in W/cm2/sr/um
tel_bgr=mean(tel_bgrsr/4.25d10)             ; in W/cm2/arcsec2/um 
tel_bgr=tel_bgr*tel_em

; atmosphere
emeff=mean(em_l/4.25d10)*tel_thru

;window
win_bgr=metis_winbgr(l, tel_t)
win_bgr=mean(win_bgr)

; add telescope and sky
sky=tel_bgr+emeff+win_bgr                         ; in W/cm2/um/arcsec2

; add in overhead to the clock time if specified
if keyword_set(ovhd) then begin
  if keyword_set(pol) then time=time*0.65 else time=time*0.8
endif

; # detected electrons from sky & telescope
nskypers=sky*sn_size*conv      ; electrons per s in sn area (1 arcsec^2 if extended source)
dit=well/(nskypers/npix)
frames=round(time/dit)
nsky=nskypers*dit              ; electrons per DIT in sn area
print, 'optimal dit = ', dit, ' s'
print, '# frames = ', frames

print, '# electrons from sky in sn area = ', nsky, ' e-/DIT'

if keyword_set(ditcheck) then begin
  print, 'Optimal DIT to avoid detector saturation is:     ', dit, ' s'
  return, dit
endif

; # detected electrons from object in s/n reference area
nobj=fobj*ee*conv*dit
if keyword_set(extended) then begin
    print, '# electrons from object per DIT = ', nobj, ' e-';/arcsec^2'
    print, 'flux from object in sn area = ', fobj, ' W/cm2/um';/arcsec^2'
  endif else begin
    print, '# electrons from object per DIT = ', nobj, ' e-'
    print, 'flux from object in sn area = ', fobj, ' W/cm2/um'
endelse  

sn=sqrt(frames)*nobj/sqrt(nsky+(npix*rd_noise^2)+(npix*dark*dit))

print, '# electrons from sky in sn area = ', nsky, ' e-/DIT'
print, 'read noise from sn area         = ', npix*rd_noise^2, ' e-/DIT'
print, 'dark current from sn area       = ', npix*dark*dit, ' e-/DIT'
if keyword_set(extended) then begin
    print, 's/n = ', sn
  endif else begin
    print, 's/n =', sn
endelse

;calculate minimum detectable flux in 1 hr to s/n of 10.
minflx=10.*sqrt(nsky+(npix*rd_noise^2)+(npix*dit*dark))/(ee*conv*dit*sqrt(frames)*transeff)    ; in W/cm2/um
minsig=minflx*emp_fac*lc^2/3./1d-19
if keyword_set(extended) then begin
  if keyword_set(rh) then begin
    print, 'sensitivity at s/n 10 in 1 hr = ', minsig*1000., ' microJy in half-light radius of ', rh, 'arsec'
  endif else begin
    print, 'sensitivity at s/n 10 in 1 hr = ', minsig*1000., ' microJy/arcsec^2'
  endelse
endif else begin
  print, 'sensitivity at s/n 10 in 1 hr = ', minsig*1000., ' microJy'
endelse


out=minsig

return, out
end