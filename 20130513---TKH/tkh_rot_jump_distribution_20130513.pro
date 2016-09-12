

;Generates a random jump distance sampled from the probability distribution:
;
;P(r) = (r/(D)*exp(-r^2/(2D))
;
;Uses the rejection method to generate this distribution froom 
;uniformly distibuted numbers
;
;This distribution has a maximum of P(r) = exp(-0.5)/((D)^0.5) 
;which occurs at r = (D)^0.5


FUNCTION TKH_Rot_Jump_Distribution_20130513, D_Rot, N_Points=N_Points

IF N_Elements(D_Rot) EQ 0 THEN Print, 'You must include a diffusion constant'


;If no number of points specified, assume only 1
IF N_Elements(N_Points) EQ 0 THEN N_Points = 1

Output_array = dblarr(N_Points + 10000)

counter = 0

;Generate 10,000 candidate 'r's at a time, keep the acceptable
;ones, then repeat until enough have been generated.
WHILE counter LT N_Points DO BEGIN

  ;Generate a pair of uniform random numbers for each point, (r,x)
  ;If x < P(r)  then accept that value of x.  If not, then throw it away

  R  = RandomU(seed, /double, 10000)
  X  = RandomU(seed, /double, 10000)

  ;In order to minimize the number of rejected values, while still allowing
  ;rare large jump events, distribute the R values between 0 and (10 x D^0.5)
  ;where P(r) < 1x10^(-8) if D > 0.0001 (i.e. tau<1000).

  R    = R * 10 *(D_Rot)^0.5
  X    = X * exp(-0.5)/((2*D_Rot)^0.5)

  P_R  = R * exp(-R^2/(4*D_Rot)) / (2*D_Rot)

  Good = where(X LT P_R, count)
  IF count GT 0 THEN BEGIN
     R_Good = R[Good]
     Output_Array[counter] = R_Good
     counter = counter + count
  ENDIF

ENDWHILE
results = (N_Points EQ 1)? Output_array[0] : Output_Array[0:N_Points-1]
RETURN, results

R = 0
X = 0
END
