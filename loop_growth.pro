PRO loop_growth

  ;Basically just runs generate_samples.pro many times
  ;for all combinations of parameters defined below.
  
  Lcuts=dindgen(30)/10.+44.5
  ratio_lims=dindgen(10)/10.+0.1
  obsc_life_fracs=dindgen(7)/10.+0.2
  mean_fedds=dindgen(10.)/10.+0.1

  v=0.
  
  nruns=n_elements(lcuts)*n_elements(ratio_lims)*n_elements(obsc_life_fracs)*n_elements(mean_fedds)
  out={lcut:0., ratio_lim:0., obsc_life_frac:0., mean_fedd:0., $
       mhalo_obsc_mean:-99., $
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
       n_obsc:-99., $
       n_unob:-99., $
       max_nspurts:-99.}
  out=replicate(out,nruns)

  FOR l=0L,n_elements(mean_fedds)-1 DO BEGIN  
     FOR j=0L,n_elements(ratio_lims)-1 DO BEGIN
     	runit=halo_bh_growth(nhalos=5000,ratio_lim=ratio_lims[j],$
                             mean_fedd=mean_fedds[l],/rest,/vary_fedd)
        FOR k=0L,n_elements(obsc_life_fracs)-1 DO BEGIN
           FOR i=0L,n_elements(lcuts)-1 DO BEGIN           
              counter,v,nruns

              x=generate_samples(lcut=lcuts[i],ratio_lim=ratio_lims[j],$
                                 obsc_life_frac=obsc_life_fracs[k])

              out[v].lcut=lcuts[i]
              out[v].ratio_lim=ratio_lims[j]
              out[v].obsc_life_frac=obsc_life_fracs[k]
              out[v].mean_fedd=mean_fedds[l]
              out[v].mhalo_obsc_mean  = x.mhalo_obsc_mean    
              out[v].mhalo_obsc_med   = x.mhalo_obsc_med     
              out[v].mhalo_unob_mean  = x.mhalo_unob_mean    
              out[v].mhalo_unob_med   = x.mhalo_unob_med     
              out[v].mstar_obsc_mean  = x.mstar_obsc_mean    
              out[v].mstar_obsc_med   = x.mstar_obsc_med     
              out[v].mstar_unob_mean  = x.mstar_unob_mean    
              out[v].mstar_unob_med   = x.mstar_unob_med     
              out[v].mbh_obsc_mean    = x.mbh_obsc_mean      
              out[v].mbh_obsc_med     = x.mbh_obsc_med       
              out[v].mbh_unob_mean    = x.mbh_unob_mean      
              out[v].mbh_unob_med     = x.mbh_unob_med       
              out[v].lbol_obsc_mean   = x.lbol_obsc_mean     
              out[v].lbol_obsc_med    = x.lbol_obsc_med      
              out[v].lbol_unob_mean   = x.lbol_unob_mean     
              out[v].lbol_unob_med    = x.lbol_unob_med      
              out[v].fedd_obsc_mean   = x.fedd_obsc_mean
              out[v].fedd_obsc_med    = x.fedd_obsc_med
              out[v].fedd_unob_mean   = x.fedd_unob_mean
              out[v].fedd_unob_med    = x.fedd_unob_med
              out[v].z_obsc_mean      = x.z_obsc_mean        
              out[v].z_obsc_med       = x.z_obsc_med         
              out[v].z_unob_mean      = x.z_unob_mean        
              out[v].z_unob_med       = x.z_unob_med         
              out[v].b_obsc           = x.b_obsc             
              out[v].b_unob           = x.b_unob             
              out[v].n_obsc           = x.n_obsc             
              out[v].n_unob           = x.n_unob             
              out[v].max_nspurts      = x.max_nspurts        
              v=v+1
           ENDFOR
        ENDFOR
     ENDFOR
  ENDFOR

  mwrfits,out,'looped_results_varyfedd.fits',/create
  
  stop
END
