FUNCTION generate_fedd,n,meanfedd,alpha=alpha,maxfedd=maxfedd,seed1=seed1,seed2=seed2
 
  ;Generates a Schechter function eddington ratio distribution.
  
  IF (n_elements(alpha) EQ 0) THEN alpha=0.4
  IF (n_elements(maxfedd) EQ 0) THEN maxfedd=ALOG10(1.)

  restore,'fedd_map.sav'
  xx=closest(meanfedd_out,meanfedd)
  meanfedd_use=meanfedd_in[xx]

  logfedd=(dindgen(6000)/1000)-5
  logfedd=logfedd[WHERE(logfedd LE maxfedd)]
  fedd=10.^logfedd
  
  dndf=((fedd)^(alpha*(-1.)))*exp((-1.)*(fedd))
  
;  logfeddcut=dindgen(600)/100-5.
;  logfeddcut=logfeddcut[where(logfeddcut LE ALOG10(maxfedd))]
  logfeddcut=logfedd
  FOR i=0L,n_elements(logfeddcut)-1 DO BEGIN
     xx=where(logfedd GE logfeddcut[i])
     IF (n_elements(meanfedds) EQ 0) THEN $
        meanfedds=wmean(logfedd[xx],weights=dndf[xx]) ELSE $
           meanfedds=[meanfedds,wmean(logfedd[xx],weights=dndf[xx])]
  ENDFOR

  usecut=closest(meanfedds,alog10(meanfedd_use))
  logfeddcut=logfeddcut[usecut]
  dndf=dndf[where(logfedd GE logfeddcut)]
  logfedd=logfedd[where(logfedd GE logfeddcut)]
  dndf=dndf/max(dndf)
  range=maxfedd-logfeddcut

  WHILE (n_elements(usefedd) LT n) DO BEGIN
	 IF ~keyword_set(seed1) THEN testfedd=randomu(systime_seed,1,/double)*range+logfeddcut ELSE $
	 	testfedd=randomu(seed1,1,/double)*range+logfeddcut
	 IF ~keyword_set(seed2) THEN testnum=randomu(systime_seed,/double) ELSE $
	     testnum=randomu(seed2,/double)
         indx=closest(logfedd,testfedd[0])
     IF (testnum LE dndf[indx]) THEN BEGIN
        IF (n_elements(usefedd)) EQ 0 THEN usefedd=testfedd ELSE $
           usefedd=[usefedd,testfedd]
     ENDIF
  ENDWHILE
  usefedd=10.^usefedd
return,usefedd
END
