FUNCTION generate_samples,obsc_life_frac=obsc_life_frac,$
                          lcut=lcut,ratio_lim=ratio_lim,$
                          zlims=zlims
;+
; After evoloving your system, apply some cuts to it to get an 'observed'
; sample. Returns a structure with a bunch of properties of the sample.
;
;  Optional Inputs:
;    obsc_life_frac - the fraction of the AGN phase spent in obscured
;                     phase. Default: 0.5
;    lcut - Luminosity cut for detection (log erg/s). Default: 45.
;    ratio_lim - parameter controlling how far the BH mass can fall
;                behind the "expected" value. Default: 0.6
;    zlims - 2 element array with lower and upper z limits. Default: [0.5,1.5]
;-  
  ;MAD restore masses
  restore,'masses_grown.sav'

  ;MAD Get halo weights
  test=file_search('halo_weights.sav')
  IF (test EQ '') THEN weights=generate_weights(mhalo_z[0,*]) ELSE $
     weights=generate_weights(mhalo_z[0,*],/rest)

  IF ~keyword_set(ratio_lim) THEN ratio_lim = 0.3
  ;MAD When a quasar is on, what fraction of time is in obscured phase
  IF ~keyword_set(obsc_life_frac) THEN obsc_life_frac=0.5
  obsc_ratio_lim = ratio_lim*(1.-obsc_life_frac)
  ;MAD Set Luminosity limit
  IF ~keyword_set(lcut) THEN lcut=45.
  ;MAD set redshift range
  IF ~keyword_set(zlims) THEN zlims=[0.5,1.5]
  
  ;MAD Limit to z slice
  usez = where((zv GE zlims[0]) AND (zv LE zlims[1]))
  tv_usez=tv[usez]
  mhalo_usez = mhalo_z[usez,*]
  lbol_usez = lbol[usez,*]
  grow_usez=grow[usez,*]
  mbh_expect_usez=mbh_expect[usez,*]
  mbh_usez=mbh[usez,*]
  mstar_usez=mstar[usez,*]
  fedd_usez=fedd[usez,*]
  weights_usez=weights[usez,*]

  ;MAD Select actively growing things (that we could see as AGN)
  ;MAD above luminosity cut and separate into the initial
  ;MAD (obscured) and later (unobscured) phase
  obsc = WHERE_XYZ((grow_usez EQ 1) AND (mbh_expect_usez - mbh_usez GE obsc_ratio_lim) $
  					AND (lbol_usez GE lcut),xind=obscz,yind=obscm)
  unob = WHERE_XYZ((grow_usez EQ 1) AND (mbh_expect_usez - mbh_usez LT obsc_ratio_lim) $
  					AND (lbol_usez GE lcut),xind=unobz,yind=unobm)

  IF (obscz[0] EQ -1 OR obscm[0] EQ -1 OR unobz[0] EQ -1 OR unobm[0] EQ -1) THEN BEGIN
     IF (obscz[0] EQ -1 OR obscm[0] EQ -1) THEN n_obsc=0 ELSE $
     	n_obsc=total(weights_usez[obscz,obscm])
     IF (unobz[0] EQ -1 OR unobm[0] EQ -1) THEN n_unob=0 ELSE $
     	n_unob=total(weights_usez[unobz,unobm])
     samp={mhalo_obsc_mean:-99., $
          mhalo_obsc_med:-99., $
          mhalo_unob_mean:-99., $
          mhalo_unob_med:-99., $
          mstar_obsc_mean:-99., $
          mstar_obsc_med:-99., $
          mstar_unob_mean:-99., $
          mstar_unob_med:-99., $
          mbh_obsc_mean:-99., $
          mbh_obsc_med:-99., $
          mbh_unob_mean:-99., $
          mbh_unob_med:-99., $
          lbol_obsc_mean:-99., $
          lbol_obsc_med:-99., $
          lbol_unob_mean:-99., $
          lbol_unob_med:-99., $
          fedd_obsc_mean:-99., $
          fedd_obsc_med:-99., $
          fedd_unob_mean:-99., $
          fedd_unob_med:-99., $
          z_obsc_mean:-99., $
          z_obsc_med:-99., $
          z_unob_mean:-99., $
          z_unob_med:-99., $
          b_obsc:-99., $
          b_unob:-99., $
          n_obsc:n_obsc, $
          n_unob:n_unob, $
          max_nspurts:-99}
  ENDIF ELSE BEGIN
     ;MAD Get effective bias of obscured/unobscured 
	 obschalos=mhalo_usez[obscz,obscm]
	 unobhalos=mhalo_usez[unobz,unobm]
	 obschalos_uniq=obschalos[rem_dup(obschalos)]
	 unobhalos_uniq=unobhalos[rem_dup(unobhalos)]
	 IF (n_elements(obschalos_uniq) GT 1) THEN $
	     b_obsc=integrate_bias(mhalo_usez[obscz,obscm],1.,weights=weights_usez[obscz,obscm],$
	     	power_spec='camb_matterpower_1.00000.dat') ELSE $
	     	b_obsc=mhalo2bias(obschalos_uniq,1.,power_spec='camb_matterpower_1.00000.dat',/silent)
	 IF (n_elements(unobhalos_uniq) GT 1) THEN $
	     b_unob=integrate_bias(mhalo_usez[unobz,unobm],1.,weights=weights_usez[unobz,unobm],$
	     	power_spec='camb_matterpower_1.00000.dat') ELSE $
	     	b_unob=mhalo2bias(unobhalos_uniq,1.,power_spec='camb_matterpower_1.00000.dat',/silent)
     
    samp={mhalo_obsc_mean:wmean(mhalo_usez[obscz,obscm],weights=weights_usez[obscz,obscm]), $
          mhalo_obsc_med:wmedian(mhalo_usez[obscz,obscm],weights=weights_usez[obscz,obscm]), $
          mhalo_unob_mean:wmean(mhalo_usez[unobz,unobm],weights=weights_usez[unobm,unobm]), $
          mhalo_unob_med:wmedian(mhalo_usez[unobz,unobm],weights=weights_usez[unobz,unobm]), $
          mstar_obsc_mean:wmean(mstar_usez[obscz,obscm],weights=weights_usez[obscz,obscm]), $
          mstar_obsc_med:wmedian(mstar_usez[obscz,obscm],weights=weights_usez[obscz,obscm]), $
          mstar_unob_mean:wmean(mstar_usez[unobz,unobm],weights=weights_usez[unobz,unobm]), $
          mstar_unob_med:wmedian(mstar_usez[unobz,unobm],weights=weights_usez[unobz,unobm]), $
          mbh_obsc_mean:wmean(mbh_usez[obscz,obscm],weights=weights_usez[obscz,obscm]), $
          mbh_obsc_med:wmedian(mbh_usez[obscz,obscm],weights=weights_usez[obscz,obscm]), $
          mbh_unob_mean:wmean(mbh_usez[unobz,unobm],weights=weights_usez[unobz,unobm]), $
          mbh_unob_med:wmedian(mbh_usez[unobz,unobm],weights=weights_usez[unobz,unobm]), $
          lbol_obsc_mean:wmean(lbol_usez[obscz,obscm],weights=weights_usez[obscz,obscm]), $
          lbol_obsc_med:wmedian(lbol_usez[obscz,obscm],weights=weights_usez[obscz,obscm]), $
          lbol_unob_mean:wmean(lbol_usez[unobz,unobm],weights=weights_usez[unobz,unobm]), $
          lbol_unob_med:wmedian(lbol_usez[unobz,unobm],weights=weights_usez[unobz,unobm]), $
          fedd_obsc_mean:wmean(fedd_usez[obscz,obscm],weights=weights_usez[obscz,obscm]), $
          fedd_obsc_med:wmedian(fedd_usez[obscz,obscm],weights=weights_usez[obscz,obscm]), $
          fedd_unob_mean:wmean(fedd_usez[unobz,unobm],weights=weights_usez[unobz,unobm]), $
          fedd_unob_med:wmedian(fedd_usez[unobz,unobm],weights=weights_usez[unobz,unobm]), $
          z_obsc_mean:wmean(zv[usez[obscz]],weights=weights_usez[obscz,obscm]), $
          z_obsc_med:wmedian(zv[usez[obscz]],weights=weights_usez[obscz,obscm]), $
          z_unob_mean:wmean(zv[usez[unobz]],weights=weights_usez[unobz,unobm]), $
          z_unob_med:wmedian(zv[usez[unobz]],weights=weights_usez[unobz,unobm]), $
          b_obsc:b_obsc, $
          b_unob:b_unob, $
          n_obsc:total(weights_usez[obscz,obscm]), $
          n_unob:total(weights_usez[unobz,unobm]), $
          max_nspurts:max(nspurts)}
   ENDELSE
   return,samp
END
