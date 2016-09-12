;Calculates the degree of depolarization caused by a large aperature
;microscope objective

FUNCTION obj_depol_Terms, NA, n=Refractive_Index
IF N_Elements(refractive_index) NE 1 THEN refractive_index = 1D

Alpha = ASin(NA/Refractive_Index) ; Angle of light cone

;Calculating formulae from Fourkas paper...
    A = 1.0D/12 * (2 - 3*cos(alpha) + (cos(alpha))^3)

    B = 1.0D/8 * (cos(alpha)-(cos(alpha))^3)

    C = 1.0D/48 * (7 -  3 * cos(alpha) - 3 * (cos(alpha))^2 - (cos(alpha))^3)

Results = {a:a, b:b, c:c}

RETURN, Results

END
