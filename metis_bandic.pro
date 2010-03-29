function metis_bandic, x, band

;this is a transmission curve for the dichroic splitting the L/M from the N band. Based on the 
; MIRI MRS dichroic for channel 1a, obtained from Alistair
; output is an array of [transmission, reflectivity, emissivity] for the relevant band


@metis_sens_param ;-> so it knows the file locations

bandicfile=profiles_dir+'bandic_prof.dat'
readcol, bandicfile, format='d,d,d', skipline=1, wave, ref, trans, /silent

;convert from percentage to decimal
ref=ref/100.
trans=trans/100.

;ref+trans+em=1
em=1.-ref-trans

;interpolate these onto the wavelength grid
ref=interpol(ref, wave, x)
trans=interpol(trans, wave, x)
em=interpol(em, wave, x)


for i=0l, n_elements(x)-1 do begin
  if (x[i] lt 1.5) then begin
    ref[i]=0.
    trans[i]=0.
    em[i]=0.
    endif
endfor 

; now need to switch trans and ref for LM bands
if ((band eq 'l') or (band eq 'm')) then profs=[[ref], [trans], [em]] else profs=[[trans], [ref], [em]]
 
 

return, profs

end