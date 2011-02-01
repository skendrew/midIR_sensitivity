function metis_ifu_sn, sig, time, lref, band, site, aomode, OVHD=ovhd, TEL_EM=tel_em,$
        EXTENDED=extended

;*************************************************************************************
;* SENSITIVITY CODE FOR IFU MODE OF METIS, THE MID-IR INSTRUMENT FOR THE EUROPEAN ELT *
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
; - Atmospheric profiles paranal_lm_100000.dat, paranal_n_100000.dat, high_and_dry_lm_100000.dat, high_and_dry_n_100000.dat
; - Parameter files metis_pre_param.pro (pre-optics params), metis_sens_param.pro (general telescope params), 
;   metis_ifu_param.pro
; - Encircled energy lookup tables, EnsquaredEnergy_AO_LMN_V10.dat for SCAO and atlas_eeav.dat for LTAO
; - DQE calculation function: metis_dqe.pro
; - IFU resolving power calculation: metis_ifures.pro
; - Window background calculation: metis_winbgr.pro
; - Instrument throughput calculation: metis_thru.pro
; 
; INPUTS
; ======
; sig          input signal, in mJy. the code will calculate the S/N of this signal as well as sensitivity of the IFU
; time         total clock time, in seconds. for e.g. sensitivity for S/N in 1 hour, set time=3600.
; lref         reference wavelength, in micron. must lie within selected band
; band         observation band, options: 'l', 'm' or 'n'. this should correspond to lref.
; site         site of observation. options: 'low' or 'high'
; aomode       mode of adaptive optics used. options: 'scao' for single-conjugate natural guide star AO, 'ltao' for laser tomography AO
; 
; OPTIONAL KEYWORDS
; =================
; OVHD         overhead. if this keyword is set, a 20% overhead will be applied, i.e. time = time*0.8.
; TEL_EM       telescope emissivity. default is 0.1 (10%). this is an important parameter in the calculation, vary it to see the effect.
; EXTENDED     extended source. see the accompanying tech note for a detailed explanation of how this works.
; 
; OUTPUT
; ======
; There are two possible outputs coded into the program.
; 1. line sensitivity @ reference wavelength (lref) -> this gives a single number at the chosen line position, in W/m2
; 2. a 2D array of wavelength [um] vs. line sensitivity [W/m2] for the typical wavelength coverage range around the reference line.
;    this is useful for plotting.



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


; read in common parameters (including file locations)
@metis_sens_param
; read in pre-optics parameters
@metis_pre_param
;read in IFU-specific parameters
@metis_ifu_param


if site eq 'low' then tel_t=286.3 else tel_t=270.
win_temp=tel_t

if not keyword_set(tel_em) then tel_em=0.1

; set clock
t=systime(1)

; BAND CHECKS
case band of
  'l': begin
    fspat=9.7d
    fspec=12.2d
    rp=100000
    rd_noise=lm_rd_noise    ; as specified in parameter file [e-/pixel]
    dark=lm_dark        ; pick dark current value for LM band
    pix_size=lm_pix_size
    scale_spat=samp_spat_lm
    scale_spec=samp_spec_lm
    well=lm_well
    fileband='lm'       ; for the atmosphere filename...
    lrange=[3.,4.2]
    end

  'm': begin
    fspat=9.7d          ; formerly 9.73
    fspec=12.2d
    rp=100000
    rd_noise=lm_rd_noise    ; as specified in parameter file [e-/pixel]
    dark=lm_dark
    pix_size=lm_pix_size
    scale_spat=samp_spat_lm
    scale_spec=samp_spec_lm
    well=lm_well
    fileband='lm'       ; for the atmosphere filename....
    lrange=[4.5, 5.5]
    end

  'n': begin
    fspat=6.65d
    fspec=8.38d
    rp=50000
    rd_noise=nq_rd_noise    ; as specified in parameter file
    dark=nq_dark
    pix_size=nq_pix_size
    scale_spat=samp_spat_n
    scale_spec=samp_spec_n
    well=nq_well
    fileband=band       ; for the atmsophere filename....
    lrange=[7.5, 14.]
    end

  else: print, "Band not recognised"

endcase

if (band eq 'n') then pltx=[lref-0.06, lref+0.06] else pltx=[lref-0.03, lref+0.03]


; READ IN ATMOSPHERIC PROFILES
;------------------------------

if site EQ 'low' then begin
    atfile=atm_dir+'paranal_'+fileband+'_'+strcompress(string(rp), /remove_all)+'.dat'
    readcol, atfile, F='D,D,x,D',l, trans, em, /silent
    name='paranal'
endif else if site EQ 'high' then begin
  atfile=atm_dir+'high_and_dry_'+fileband+'_'+strcompress(string(rp), /remove_all)+'.dat'
  readcol, atfile ,F='D,D,x,D', l, trans, em, /silent
  name='high_and_dry'
endif
print, 'Atmosphere data read'

;perform some unit conversions on the emissivity values (from photons/s/um/m^2/arcsec^2 to W/cm^2/sr/um)
em_l=em*(h*c/(l*1d-6))*1d-4*4.25d10

; clip these arrays to wavelength range of interest
em_l=em_l[where((l gt lrange[0]) and (l le lrange[1]))]
trans=trans[where((l gt lrange[0]) and (l le lrange[1]))]
l=l[where((l gt lrange[0]) and (l le lrange[1]))]


; calculate resolution element
delta_l=metis_ifures(band, l, 'ifu')

; find reference wavelength location
diff=l-lref
tmp=min(diff,ref, /absolute)
print, 'delta_l @ lref = ', delta_l[ref]*1000.,   ' nm'
print, 'atmospheric transmission @ lref = ', trans[ref]

!p.charsize=2.
;plot, l, delta_l

;spectral sampling: 2.5 pixels per delta_l optimised at 4.7 and 12.8 um
pix_spec=2.5        ; amend this later

;spatial pixels: sampling optimised at 3.7 um for LM-band, 9 um for N-band 
; for extended sources: calculate per spaxel
if (band eq 'n') then pix_spat=2.*(l/9.) else pix_spat=2.*(l/3.7)

;calculate # pixels per diffraction element using pixel scale from design and diffraction limited core
npix=ceil(pix_spat*pix_spec)

;calculate S/N reference area, in mas^2. if seeing limited then size of seeing disk, else 1.22l/D
;if aomode eq 'noao' then sn_area=1000*0.8*(lc/0.5)^(-0.2) else sn_area=pix_scale*pix_spat*pix_spec
sn_size=npix*scale_spat*1d-3*scale_spec*1d-3

print, '# spatial pixels in sn area @ lref = ', pix_spat[ref]
print, 'sn size @ lref = ', sn_size[ref], ' arcsec^2'


;look up EE using appropriate AO mode
;but now sn_size and ee_snsize are vectors of different length -> need to lookup element by element
if keyword_set(extended) then begin
  if keyword_set(rh) then ee=0.5 else ee=1. 
endif else begin
  if aomode eq 'ltao' then begin
    ee_file=ee_dir+'atlas_eeav.dat' 
    readcol, ee_file, f='f,d,d,d', ee_snsize, ee_l, ee_m, ee_n, /silent
  endif else begin    
    ee_file=ee_dir+'EnsquaredEnergy_AO_LMN_V10.dat'
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
  diff=ee_snsize-sqrt(sn_size[ref]*1d6)
  tmp=min(diff, loc, /absolute)
  if aomode eq 'ltao' then ee=ee_list[loc]/100. else ee=ee_list[loc]
endelse

print, 'ee = ', ee

; calculate telescope transmission
if (tel_em gt 0.1) then tel_thru=1.-tel_em else tel_thru=mir_ref^no_mir
print, 'telescope transmission = ', tel_thru

; calculate throughput
obsmode='hspec'
thru=metis_thru(obsmode, band, l)
print, 'throughput of instrument @ lref = ', thru[ref]


; calculate dqe
dqe=metis_dqe(l, band)
print, 'DQE in filter band @ lref = ', dqe[ref]

; calculate efficiency
eff=thru*tel_thru*dqe*pcg

print, 'total efficiency including DQE, pcg, telescope, and instrument @ lref = ', eff[ref]

; total conversion factor [e-/s per W/cm2/um] (eff_area is in cm2)
conv=eff*delta_l*eff_area*l/1.985d-19

; flux from object
fobj=sig*3.*1d-19*trans/l^2    ; convert flux from mJy to [W/cm2/um]

; BACKGROUNDS
; telescope background = blackbody ftn @ ambient T * telescope emissivity
tel_bgrsr=planck(l*1d4, tel_t)*1d-3/!dpi    ; in W/cm2/sr/um
tel_bgr=tel_bgrsr/4.25d10             ; in W/cm2/arcsec2/um 
tel_bgr=tel_bgr*tel_em

; atmosphere background = emissivity from the atm file * telescope throughput
emeff=em_l*tel_thru/4.25d10          ; in W/cm2/arcsec2/um

; entrance window background
win_bgr=metis_winbgr(l, tel_t)      ; in W/cm2/arcsec2/um

; add telescope and sky
sky=tel_bgr+emeff+win_bgr                     ; in W/cm2/um/arcsec2

; add in overhead to the clock time if specified
if keyword_set(ovhd) then time=time*0.8

nskypers=sky*sn_size*conv           ; electrons per s in sn area

; # detected electrons from object in sn area
if keyword_set(extended) then begin
  nobjpers=fobj*sn_size*conv
endif else begin
  nobjpers=fobj*ee*conv              ; electrons per s in sn area
endelse

; calculate integration time to ensure background limited (i.e. not read noise limited)
dit=2*rd_noise^2/((nskypers[ref]/npix[ref])-dark)

frames=round(time/dit)
nsky=nskypers*dit              ; electrons per DIT in sn area
print, 'optimal dit = ', dit, ' s'
print, '# frames = ', frames


; # detected electrons from object in s/n reference area
if keyword_set(extended) then nobj=fobj*sn_size*conv*dit else nobj=fobj*ee*conv*dit
print, 'flux from object in sn area at ', lref, ' um = ', fobj[ref], ' W/cm2/um'
print, '# electrons from object in sn area at ', lref, ' um = ', nobj[ref], ' e-/DIT'

; calculate signal to noise
sn=sqrt(frames)*nobj/sqrt(nsky+(npix*rd_noise^2)+(npix*dark*dit))
print, 'S/N at ', lref, ' um = ', sn[ref]

trange=[dit-dit*0.9, dit+dit*0.9]
exptime=findgen(1000)+10.
 

;calculate minimum detectable flux in 1 hr to s/n of 10.
if keyword_set(extended) then minflx=10.*emp_fac*sqrt(nsky+(npix*rd_noise^2)+(npix*dit*dark))/(sn_size*trans*conv*dit*sqrt(frames)) $
    else minflx=10.*emp_fac*sqrt(nsky+(npix*rd_noise^2)+(npix*dit*dark))/(ee*trans*conv*dit*sqrt(frames))    ; in W/cm2/um

minsig=minflx*lref^2/3./1d-19        ; CONTINUUM SENSITIVITY in mJy
minsig_line=minflx*delta_l*1d4        ; LINE SENSITIVITY in W/m2


print, 'continuum sensitivity at s/n 10 in 1 hr at ', lref, ' um', minsig[ref], ' mJy'
print, 'line sensitivity at s/n 10 in 1 hr at ', lref, ' um = ', minsig_line[ref], ' W/m2'

!p.thick=2.
!p.charsize=1.
!p.charthick=1.5

plottitle='Line sensitivity, '+ band+ '-band'
plot, l, minsig_line, /ylog, yrange=[1d-22, 1d-18],title=plottitle, xtitle='micron', ytitle='w/m2'


;out=[[l], [minsig_line]] ;IF WANT SENSITIVITY OVER FULL WAVELENGTH RANGE
out=minsig_line[ref]      ;IF WANT SENSITIVITY AT REFERENCE WAVELENGTH ONLY 
return, out
end