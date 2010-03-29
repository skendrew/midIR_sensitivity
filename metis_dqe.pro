function metis_dqe, l, band

; this function calculates the quantum efficiency of METIS detectors. 
; based on experimental curves for the detector types.
; arguments: 
; l - wavelength grid in micron
; band - band of observation - 'l', 'm' or 'n'

dqe=dblarr(n_elements(l))

if (band eq 'l') or (band eq 'm') then begin
for i=0l, n_elements(l)-1 do begin
    if ((l[i] gt 1.) and (l[i] le 5.3)) then dqe[i] = 0.7 else begin
      if ((l[i] gt 5.3) and (l[i] le 5.5)) then dqe[i] = 0.5 else dqe[i] = 0.0
    endelse
endfor
endif else begin
for i=0l, n_elements(l)-1 do begin
    if (l[i] lt 8.) then dqe[i] = (0.06*(l[i]-3.))+0.3 else begin
       if ((l[i] gt 8.) and (l[i] le 14.)) then dqe[i]=0.6 else begin
          if ((l[i] gt 14.) and (l[i] le 19.)) then dqe[i]=0.7 else begin
            if ((l[i] gt 19.) and (l[i] le 24.)) then dqe[i]=0.6 else dqe[i]=(-0.092*(l[i]-24.))+0.6
            endelse
          endelse
        endelse    
endfor
endelse          

dqe=dqe*1.2

return, dqe
 
end