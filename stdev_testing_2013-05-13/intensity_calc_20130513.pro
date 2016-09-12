PRO Intensity_Calc_20130513, X,Y,Z, NA_detect, NA_excit, I_left_anal, I_right_anal, refractive_index=refractive_index

;calculate the intensity for the left and right side based on the x,y,z coordinates from simulation
;the intensitiy is returned with values from 0 -2 with mean 1.

fourkas = obj_depol_terms(double(NA_detect), n=refractive_index)
A = fourkas.a
B = fourkas.b
C = fourkas.c

if NA_excit EQ 0 then NA_excit = 0.001D ; value cannot be excactly 0.000, because it is not defined
fourkas = obj_depol_terms(double(NA_excit),n=refractive_index)
A_excit = fourkas.a
B_excit = fourkas.b
C_excit = fourkas.c

Excitation  =  2.0*A_excit + 2.0*B_excit *(1-Z^2)
Det_Left  =  (A + B * (1-Z^2) + C * (X^2 - Y^2))
Det_Right =  (A + B * (1-Z^2) - C * (X^2 - Y^2))
I_Left_Anal  = abs(Excitation * Det_Left)
I_Right_Anal = abs(Excitation * Det_Right)
;I_ave = mean([I_left_anal,I_Right_anal], /nan)
;mean intensity should be one and range should be 0 - 2, so we need to scale with average
I_left_anal /= mean(I_left_anal,/nan)
I_right_anal /= mean(I_right_anal,/nan)

PRINT, 'END OF PROGRAM INTENSITY_CALC'
END