PRO P2_Calc_20150731, X,Y,Z, p2
;KS: take the ACF of the rank-2 Legendre polynomial of the cosine of angel between two vectors (that of time t and time 0)
;x,y,z are 1D inputs...
;we want to calculate angle change relative to the angle at time 0
;
;vectors generated from tkh_rot... are all unit vectors, so....

dthet = fltarr(n_elements(x))
FOR i=1L, n_elements(x)-1 DO dthet[i] = x[0]*x[i] + y[0]*y[i] + z[0]*z[i]
p2 = (3*dthet^2 - 1)/2

PRINT, 'END OF PROGRAM P2_CALC'
END
