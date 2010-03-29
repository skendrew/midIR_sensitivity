function metis_wfsdic, x

; this calculates the transmission profile of the METIS dichroic that splits the light off to the WFS
; it is located just behind the entrance window in the cryostat

trans=dblarr(n_elements(x))

for i=0l, n_elements(x)-1 do begin
  if (x[i] le 2.5) then trans[i]=0. else begin
    if ((x[i] gt 2.5) and (x[i] le 3.5)) then trans[i]=0.85*(x[i]-2.5) else trans[i]=0.85
    endelse
endfor
     

return, trans
end