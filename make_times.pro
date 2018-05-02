PRO make_times,zmax=zmax

  ;Get cosmic time at each z.
  
  IF (n_elements(zmax EQ 0)) THEN zmax=3.
  
  z=(findgen(zmax*1000.+1.)/1000.)
  tv=fltarr(n_elements(z))
  
  FOR i=0L,n_elements(z)-1 DO BEGIN
     tmp=cosmocalc(z[i])
     tv[i]=tmp.t_l
  ENDFOR
  save,tv,file='times_fine.sav'

  stop
END
