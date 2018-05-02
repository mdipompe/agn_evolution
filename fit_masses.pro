;+
;  NAME:
;    fit_masses
;  PURPOSE:
;    Take a mass distribution and fit a model dndm
;
;  USE:
;    fit_masses,masses,sampled_m,dndm,binsize=binsize
;
;  INPUT:
;    masses - the mass values of the real data set to fit
;    sampled_m - the z values you want to know the fit at
;
;  Optional Inputs:
;    binsize - the bin size of the histogram to fit (default 0.2)
;
;  OUTPUT:
;    dndm - the values of the fit at sampled_m
;    dndm_fit.txt - a text file with m and the fit, so you can check
;                   that things worked if you want.
;
;  HISTORY:
;    11-7-14 - Written - MAD (UWyo)
;-
PRO fit_masses,masses,sampled_m,dndm,binsize=binsize,weights=weights

;MAD Set defaults
IF ~keyword_set(binsize) THEN binsize=0.1

;MAD Bin the real data
h=histogram_weight(masses,weight=weights,binsize=binsize,min=0,max=max(masses))

;MAD Make array of bin centers
x=fltarr(ceil((max(masses)/binsize))+1)
x[0]=binsize/2.
i=1
WHILE (max(x) LT max(masses)) DO BEGIN
 x[i]=x[i-1]+binsize
 i=i+1
ENDWHILE

;MAD Smooth over gaps in the histogram
xx=where(h NE 0)
h=h[xx]
x=x[xx]

;MAD Force fit to N=0 at z=0
x=[0,x]
h=[0,h]

;MAD Fit cubic spline, normalize
fit=spline(x,h,sampled_m)
fit[where((sampled_m GT max(masses)) OR (sampled_m LT min(masses)))]=0
area=0.
FOR i=0L,n_elements(sampled_m)-1 DO BEGIN
 area=area+(fit[i]*(sampled_m[1]-sampled_m[0]))
ENDFOR
dndm=fit/area

openw,1,'dndm_fit.txt'
printf,1,';z     fit    dndm (normalized fit)'
FOR i=0,n_elements(sampled_m)-1 DO BEGIN
  printf,1,sampled_m[i],fit[i],dndm[i],format='(F,1x,F,1x,F)'
ENDFOR
close,1

return
END
