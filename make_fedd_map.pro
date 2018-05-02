PRO go

  fedds=dindgen(999)/1000.+0.001
;  fedds=findgen(9)/10.+0.1
  modfedds=(1.0842100*fedds)-0.07324078d
  
  meanlogfedd=fltarr(n_elements(fedds))
  meanfedd=fltarr(n_elements(fedds))
  FOR i=0L,n_elements(fedds)-1 DO BEGIN
     fedddist=generate_fedd_tmp(10000,((fedds[i])),seed1=98,seed2=156)
     meanlogfedd[i]=mean(fedddist)
     meanfedd[i]=mean(10.^fedddist)
  ENDFOR

  meanfedd_out=meanfedd
  meanfedd_in=fedds
  save,meanfedd_in,meanfedd_out,filename='fedd_map.sav'
  

  
  stop
END
