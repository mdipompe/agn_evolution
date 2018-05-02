FUNCTION wmedian,data,weights=weights

  IF (n_elements(weights) EQ 0) THEN weights=fltarr(n_elements(data))+1.

  sweights=weights[bsort(data)]
  sdata=data[bsort(data)]
  totw=total(sweights)
  cumtot=total(sweights,/cum)
  xx=where(cumtot GT total(sweights)/2.)
  meddata=sdata[xx[0]]

  return,meddata
END
