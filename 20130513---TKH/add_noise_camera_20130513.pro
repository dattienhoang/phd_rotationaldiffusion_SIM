PRO ADD_NOISE_CAMERA_20130513, I_left_anal, I_right_anal, LD, noise_factor = noise_factor, thresholding = thresholding

;this program adds noise to intensity traces for our experimental setup.
;in a first step we add camera noise and then background noise
;the number of photons and the gain, as well as the level of background noise is
;based on experimental traces. 

;I_left_anal, I_right_anal are the analytical signal with mean 1 and range 0-2. This need to be supplied from the
;main program. LD should be an array with !values.f_nan and have the same number of elements as I_left_anal(I_right_anal)
; It needs to be supplied from main program as well. When returning to main program, LD will contain the LD values
;based on noisy intensities and NaN where the total intensity was below 2*standard deviation of noise.

;now we have the analytical signal and we start adding the various noises.


if n_elements(noise_factor) EQ 0 then noise_factor = 0.30

n_points = double(n_elements(I_left_anal)) ; point in trajectory

;print, 'PRO ADD_NOISE CAMERA number of points =', n_points
;1.) set number of photons per frame, should be 200 in our case
I_left_anal  *= 200.00
I_right_anal *= 200.00

;2.) obtain number of electron by multiplying photons with Gain, here about 300, in paper actually M^2

I_left_anal *= 300.00
I_right_anal *=300.00


;3.) calculate the average electron number 
I_L_AVE = mean(I_left_anal, /nan)
I_R_AVE = mean(I_right_anal,/nan)

;I_left_anal += I_L_AVE*0.30
;I_right_anal += I_R_AVE*0.30
;4.) Add statistics of detection: gaussian random number with mean = # of electrons
;at a particular time and standard deviation = 2*mean (see Andor manual page 184).

;Original reference: Transaction on electronic devices, Vol50, 2003, 1227-1232:
;The noise perfomance of electron Multiplying charge coupled devices by Mark Stanford Robbins and BJ Hadwen

;Need to be done for each step seperately. 
;produce gaussian distribution with
;array has mean =0, standard deviation (sigma) = 1 from randomN
;change distibution to mean of current intensity (in electrons now) with standard deviation (sigma) = 2
;in formula below: 
;sigma = 2.0*sqrt(2*I_left_anal) ; standard deviation twice the number of detected electrons (Andor Manual p.184)

random_array = randomN(seed,n_points)
I_left_anal += random_array * sqrt(I_left_anal*2.0)
random_array = randomN(seed,n_points)
I_right_anal += random_array * sqrt(I_right_anal * 2.0) 


;step4: After we added noise from the detector we need to add background noise
;We add background noise to each intensity. From experimental data we now
; background noise is about 30% 
I_Left_anal += randomn(seed,n_points)*I_L_ave *noise_factor
I_right_anal += randomn(seed, n_points)*I_R_ave *noise_factor
;determine the standard deviation of background in terms of the some of both channels
;background does not just add up but instead with sqrt(2)
noise_stdev = sqrt(2.0)*(noise_factor*I_L_AVE + noise_factor*I_R_AVE)/2.0 ;standard deviation of noise in I_total

;threshold is 2 time standard deviation of noise
if thresholding EQ 1 then begin
print, 'thresholding data'


;both = where(I_left_anal LE 0.00000 AND I_right_anal LE 0.0000, count_both)
;if count_both GT 0 then begin
;; too avoid that LD is NAN in case both intensities are zero we set the higher one to 1 even when negative.
;bigger=where(I_left_anal[both] GE I_right_anal[both], complement=Nbigger)
;I_left_anal[both[bigger]] = 1.0
;I_right_anal[both[Nbigger]] = 1.0
;endif

bad_L=where(I_Left_anal LE 0.00000, count_L)
bad_R=where(I_right_anal LE 0.00000, count_R)
if count_L Gt 0 then I_left_anal[bad_L] = 0.0
if count_R Gt 0 then I_right_anal[bad_R] = 0.0
LD = (I_left_anal - I_right_anal) / (I_left_anal + I_right_anal)
endif else begin
;print, 'no thresholding for data'

LD = (I_left_anal - I_right_anal) / (I_left_anal + I_right_anal)
endelse

;print, 'END OF PROGRAM ADD_NOISE_CAMERA'
END
