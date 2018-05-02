FUNCTION generate_weights,logm,rest=rest

  ;Generate an array of weights for halo masses of logm
  ;IF /rest keyword set, restores .sav file with weights already
  ;generated in halo_weights.sav.
  ;Assumes you already have P(k) in file camb_matterpower_0.00000.dat'
  ;generated from CAMB.

  IF ~keyword_set(rest) THEN BEGIN
     ;MAD Generate weights based on HMF
     dndm=halo_mass_function(logm,0.,power_spec='camb_matterpower_0.00000.dat')
     weight=dndm/total(dndm)
     zv = reverse(findgen(3001)/1000.)

     ;MAD Copy weights for every z (must be a faster way...)
     weights=weight
     FOR i=0L,n_elements(zv)-2 DO BEGIN
        weights=[[weights],[weight]]
     ENDFOR
     weights=transpose(weights)
     save,weights,filename='halo_weights.sav'
  ENDIF ELSE BEGIN
     restore,'halo_weights.sav'
  ENDELSE
  return,weights

END
