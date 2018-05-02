FUNCTION halo_bh_growth,rest=rest,ratio_lim=ratio_lim,$
                        mean_fedd=mean_fedd,vary_fedd=vary_fedd,nhalos=nhalos,$
                        maxz=maxz
;+
;  Get an initial (flat) distribution of halos, get their 
;  corresponding stellar and BH masses, and grow them from
;  zmax to z=0
;
;  Optional keywords:
;    ratio_lim - parameter controlling how far a BH can lag behind
;                expected value before growing.
;    mean_fedd - the mean eddington ratio of the full AGN population
;    vary_fedd - set to have a Schechter function distribution of
;                f_edd, otherwise all AGN have mean value
;    nhalos - number of sources to start with at max z.
;    maxz - the redshift that things start at.
;  
;-                         
  
;  st=timer()
  
  ;MAD Set some parameters
  h=0.702

  ;MAD Set the minimum ratio of BH mass to expected mass for when
  ;MAD a quasar turns on
  IF ~keyword_set(ratio_lim) THEN ratio_lim = 0.3
  ;MAD Set the Eddington ratio
  IF ~keyword_set(mean_fedd) THEN mean_fedd = 0.1
  
  ;MAD Set fine array of z 
  zv = reverse(findgen(maxz*1000+1.)/1000.)

  ;Read in the times these zs correspond to (from make_times.pro)
  restore,'times_fine.sav'
  tv=reverse(tv)

  IF ~keyword_set(rest) THEN BEGIN
     ;MAD Generate flat distribution of halo masses
	 mhalo_0=randomu(713,nhalos)*4.5+10.5
     
     ;MAD Generate array of halo masses at each z, for each initial halo mass
     mhalo_z = dblarr(n_elements(zv),n_elements(mhalo_0))
     ;MAD Project each halo mass from z=3 to z=0
     FOR i=0, n_elements(mhalo_0)-1 DO BEGIN
        tmp=halo_growth(mhalo_0[i],3.)
        mhalo_z[*,i]=reverse(tmp.mhalo)
     ENDFOR

     ;MAD Calculate expected stellar mass for halo (Guo et al. 2010)
     term1=(((10.^mhalo_z)/h)/(10.^11.4))^(-0.926)
     term2=(((10.^mhalo_z)/h)/(10.^11.4))^(0.261)
     mstar=ALOG10(((10.^mhalo_z)/h)*0.129*((term1+term2)^(-2.44)))

     ;MAD Calculate expected BH mass for stellar mass (Haring & Rix 2004)
     mbh_expect = 8.2 + 1.12*alog10((10.^mstar)/(10.^11.)) 

     save,mhalo_0,mhalo_z,mstar,mbh_expect,filename='masses.sav'
  ENDIF ELSE BEGIN
     restore,'masses.sav'
  ENDELSE

  ;MAD Initialize BH masses
  mbh = mbh_expect*0.
  ;MAD Generate initial BH masses as a fraction of their expected values
  mbh[0,*]=mbh_expect[0,*]+ALOG10((randomu(637,n_elements(mbh_expect[0,*]))*0.9+0.1))

  IF keyword_set(vary_fedd) THEN BEGIN
	 feddtmp=generate_fedd(10000,mean_fedd,seed1=341,seed2=267)
	 tmp=floor(randomu(122,n_elements(zv)*n_elements(mhalo_0))*(n_elements(feddtmp)-1.))
	 fedd=reform(feddtmp[tmp],n_elements(zv),n_elements(mhalo_0))
  ENDIF ELSE BEGIN
     fedd=fltarr(n_elements(zv),n_elements(mhalo_0))+(mean_fedd)
  ENDELSE

  ;MAD Initialize Lbol array
  lbol = mbh*0.

  ;MAD Calculate delta t values for each z step
  dt = diff(tv)
  dt = [0,dt]

  ;MAD The grow flag marks if the BH is growing at each z
  grow=mhalo_z*0

  ;MAD Count how many times each BH goes through an AGN phase
  nspurts=intarr(n_elements(mhalo_0))
  
  ;MAD Loop over time, turning on quasars when
  ;MAD they fall outside of expected range by ratio_lim
  FOR i=1L,n_elements(zv)-1 DO BEGIN
     ;MAD First just keep previous BH mass
     mbh[i,*] = mbh[i-1,*]

     ;MAD Determine if need to start growing - either previous 
     ;MAD mass is outside ratio, or grow flag is still on
     igrow=WHERE((mbh_expect[i-1,*]-mbh[i-1,*] GT ratio_lim) OR (grow[i-1,*] eq 1),ngrow)

     ;MAD Grow the ones that need to grow
     IF (ngrow GT 0) THEN BEGIN        
        ;MAD Grow by the amount allowed in delta t by Salpeter time
        tsalp = 4.5d7*(fedd[i,igrow]^(-1))
        mbh[i,igrow] = mbh[i-1,igrow]+alog10(exp(abs(dt[i])*1d9/tsalp))

        ;MAD Turn Mbh into Lbol based on Eddington fraction
        lbol[i,igrow] = alog10(1.51d38*fedd[i,igrow])+mbh[i,igrow]

        ;MAD Mark it as growing at this step
        grow[i,igrow] = 1

        ;MAD Mark as a unique start to growth if needed
        xx=where(grow[i-1,igrow] EQ 0 AND grow[i,igrow] EQ 1,cnt)
        IF (cnt NE 0) THEN nspurts[igrow[xx]]=nspurts[igrow[xx]]+1
        
        ;MAD Determine which ones are done growing (over expected mass)
        idone=WHERE(mbh[i,*] GE mbh_expect[i,*],ndone)
        IF (ndone GT 0) THEN grow[i,idone] = 0
     ENDIF
  ENDFOR
  
  save,zv,tv,mhalo_z,lbol,grow,mbh_expect,mbh,mstar,fedd,nspurts,filename='masses_grown.sav'
     
;  et=timer(/fin,st=st,unit='s')
end
