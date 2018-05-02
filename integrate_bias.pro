FUNCTION integrate_bias,halomasses,z,h=h,omega_m=omega_m,omega_l=omega_l,omega_b=omega_b,$
                   power_spec=power_spec,weights=weights

;Take an array of halo masses, fit a spline to get dNdM, then get
;the bias of that distribution.
  
IF (n_elements(weights) EQ 0) THEN weights=fltarr(n_elements(halomasses))+1.

IF ~keyword_set(h) THEN h=0.702
IF ~keyword_set(omega_m) THEN omega_m=0.275
IF ~keyword_set(omega_b) THEN omega_b=0.046
IF ~keyword_set(omega_l) THEN omega_l=0.725
IF ~keyword_set(spec_ind) THEN spec_ind=0.968

plothist_weight,halomasses,weights,xvals,yvals,bin=0.1,/noplot

IF (n_elements(where(yvals GT 0)) GE 3) THEN BEGIN
	;MAD Generate list of log(masses), b(M)
	logm=findgen(1000)/100.+8
	m=10.^logm
	bofm=mhalo2bias(logm,z,power_spec=power_spec,Delta=200.,$
    	            h=h,omega_m=omega_m,omega_l=omega_l,omega_b=omega_b,$
    	            spec_ind=spec_ind,/silent)

	fit_masses,halomasses,logm,dndm,binsize=0.1,weights=weights

	num=int_tabulated(logm,bofm*dndm,/double)
	den=int_tabulated(logm,dndm,/double)
	bias=num/den
ENDIF ELSE BEGIN
	usebins=where(yvals GT 0)
	centers=xvals[usebins]
	natcenter=yvals[usebins]
	FOR i=0L,n_elements(centers)-1 DO BEGIN
		b=mhalo2bias(centers[i],z,power_spec=power_spec,Delta=200.,$
    	             h=h,omega_m=omega_m,omega_l=omega_l,omega_b=omega_b,$
    	             spec_ind=spec_ind,/silent)
		IF (n_elements(bvals) EQ 0) THEN bvals=b[0] ELSE $
			bvals=[bvals,b[0]]
	ENDFOR
	bias=wmean(bvals,weights=natcenter)
ENDELSE

return,bias
END
