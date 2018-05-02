FUNCTION wmean,data,weights=weights

  IF (n_elements(weights) EQ 0) THEN weights=fltarr(n_elements(data))+1.

  wmean=total(data*weights,/double)/total(weights,/double)

  return,wmean
END
