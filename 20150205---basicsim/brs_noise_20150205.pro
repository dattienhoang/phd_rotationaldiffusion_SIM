;brs_noise_20150205.pro -- 
;Routine that adds noise to simulated rotational diffusion trajectories
;
;Dat Tien Hoang, Last major edit: 2015-02-05
;
;This program adds noise to intensity traces after they have been calculated
;from Fourkas' relations. The analytical signal that is inputted is expressed
;in arbitrary units whose mean is 1.
;
;First, photobleaching is added as a percentage of the trajectory completed.
;Then a constant background is added (also expressed in arbitrary units). Next,
;the trajectory is scaled by photons and EM gain and then camera noise is added.
;
;Most of the inputted parameters need to be calculated from experiment. For
;example, the PSF size needs to be empirically determined and to determine how
;many photons are collected on that pixel. We assume, as per TKH's calculations,
;that the entire molecule emits 200 photons per frame at 1Hz. The value for 
;background needs to be back-calculated to attain the value in arbitrary units.
;
;CAVEATS:
;1.) Remember to keep units consistent! Camera noise is expressed for one pixel!

PRO brs_noise_20150205, I_left_anal, I_right_anal, I_tot_sm, LD, threshold, thresh, $
                     noise_factor = noise_factor, $
                     blch=blch

;print, mean(i_left_anal)
;------------------------------------------------------------------------------
;SET PARAMETERS
  ;Signal processing parameters:
  cnoise = 1 ;[0 = do not add camera noise
              ;1 = add camera noise
  corrINT = 1;[0 = do not 
             ; 1 = remove intensities lower than 0
  ;
  ;various parameters about generated background
  t_blchs = 0.25 ;the % of the trajectory before the molecule photobleaches
  B = 0.03       ;the background level before photons and EM gain are factored in
                 ;back-calculate from experiment using parameters in the next 
                 ;chunkto get this value
  SNB = 1.25     ;the SBR desired
  ;
  ;various scaling factors... you need to know a lot about what you are 
  ;measuring (feature size in image, photons emitted, EM gain characteristics,
  ;how to extract intensity the best way, etc), but just simply making these 
  ;following numbers to your experiment works fine too
  pix = 5.0    ;1D particle size
    ;How many photons are we collecting in a pixel? Follow experimental 
    ;psf = psf_gaussian(npixel=pix, fwhm=[pix/2.0,pix/2.0], /normalize) * 200 ;assumes 200 photons from molecule
  phtn_on = 20;200/10;median(psf);photons per frame per pixel
  EM = 300.0   ;EM gain
  ;
  ;various instructions on how to treat background
  sub_bk = 2 ;[0 = do not subtract background
             ; 1 = add "noise due to background subtraction" (30%; Mackowiak+Kaufman, J. Chem. Phys)
             ; 2 = subtract a background (just background noise)
             ; 3 = subtract a smoothed intensity background] 
  ;various instructions on how to threshold
  thresh = 0;[0 = no threshold
            ; 1 = flat threshold
            ; 2 = smooth intensity threshold]
;------------------------------------------------------------------------------
;SIGNAL PREPARATION
  n_points = double(n_elements(I_left_anal)) ; points in trajectory
  
  ;choose a bleaching time. this is distince from the trajectory length set in the main level routine!
  IF blch EQ 1 THEN BEGIN
    ;I_left_anal *= snr_fac & I_right_anal *= snr_fac
    t_blch = t_blchs * n_elements(I_left_anal)
    ;t_blch = 0.25 * n_elements(I_left_anal)
    print, '...bleached @', t_blch, 'of', n_points
    IF t_blchs NE 1 THEN I_left_anal[t_blch:*] = 0
    IF t_blchs NE 1 THEN I_right_anal[t_blch:*] = 0
    I_left_anal *= B*(SNB-1) & I_right_anal *= B*(SNB-1)
  ENDIF
  
  I_left_anal += B & I_right_anal += B

  ;1.) Set number of photons per frame. (Determined empirically above)
  I_left_anal *= phtn_on & I_right_anal *= phtn_on
  ;print, mean(i_left_anal)

  ;4.) obtain number of electron by multiplying photons with EMGAIN, 
  ;    here about 300, in paper actually M^2
  I_left_anal *= EM & I_right_anal *= EM
  ;print, mean(i_left_anal)
  
  ;3.) Add statistics of detection: gaussian random number with mean = # of 
  ;    electrons at a particular time and standard deviation = 2*mean (see 
  ;    Andor manual page 184).
  ;Original Ref: Transaction on Electronic Devices, Vol50, 2003, 1227-1232:
  ;The Noise Perfomance of Electron Multiplying Charge Coupled Devices 
  ;by Mark Stanford Robbins and BJ Hadwen
  ;
  ;Need to be done for each point seperately. Produce Gaussian distribution with
  ;array has mean = 0, standard deviation (sigma) = 1 using randomN.
  ;Change distibution to mean of current intensity (in electrons now) with 
  ;standard deviation (sigma) = 2 using the formula below: 
  ;  sigma = 2.0*sqrt(2*I_left_anal) 
  ;^the standard deviation of twice the number of detected e- (Andor Manual p.184)
  IF cnoise EQ 1 THEN BEGIN
    print,'...adding camera noise.'
    ;All Gaussian
    random_array = randomN(seed,n_points)
    I_left_anal += random_array * sqrt(I_left_anal * 2.0)
    random_array = randomN(seed,n_points)
    I_right_anal += random_array * sqrt(I_right_anal * 2.0)
  ENDIF

;------------------------------------------------------------------------------
;NOW START SIGNAL PROCESSING: BACKGROUND SUBTRACTION
  IF sub_bk EQ 0.0 THEN BEGIN
    print, '...no further background-related treatment.'
  ENDIF
  
  ;as in original stdev_testing.pro
  IF sub_bk EQ  1.0 THEN BEGIN
    IF n_elements(noise_factor) EQ 0 THEN noise_factor = 0.30
    print,'...adding NABS (CR Rubrene Paper) (v1.0) @ NF', noise_factor, '.'
    ;step4: After we added noise from the detector we need to add background noise
    ;5.) calculate the average electron number 
    IF t_blch NE 1 THEN BEGIN
      I_L_AVE = mean(I_left_anal[0:t_blch-1], /nan)
      I_R_AVE = mean(I_right_anal[0:t_blch-1], /nan)
      ENDIF ELSE BEGIN
      I_L_AVE = mean(I_left_anal, /nan)
      I_R_AVE = mean(I_right_anal, /nan)
    ENDELSE
    ;We add background noise to each intensity. From experimental data we now
    ; background noise is about 30% 
    I_Left_anal += randomn(seed, n_points)*I_L_ave *noise_factor
    I_right_anal += randomn(seed, n_points)*I_R_ave *noise_factor
    ;determine the standard deviation of background in terms of the some of both channels
    ;background does not just add up but instead with sqrt(2)
    noise_stdev = sqrt(2.0)*(noise_factor*I_L_AVE + noise_factor*I_R_AVE)/2.0 ;standard deviation of noise in I_total
    ;print, noise_stdev
  ENDIF

  IF sub_bk EQ 2.0 THEN BEGIN
    print, '...subtracting a background (v2.0)'; @ NF', noise_factor, '.'
    random_array = randomN(seed, n_points)
    bk_left_sub = make_array(n_points, value=B)
      bk_left_sub *= phtn_on & bk_left_sub *= EM
      bk_left_sub += random_array * sqrt(bk_left_sub * 2.0)
    random_array = randomN(seed, n_points)
    bk_right_sub = make_array(n_points, value=B)
      bk_right_sub *= phtn_on & bk_right_sub *= EM
      bk_right_sub += random_array * sqrt(bk_right_sub * 2.0)
    I_left_anal -= bk_left_sub & I_right_anal -= bk_right_sub
  ENDIF
  
  IF sub_bk EQ 3.0 THEN BEGIN
    smth = 500
    print, '...subtracting a smooth background (v3.0)'; @ NF', noise_factor, '.'
    random_array = randomN(seed, n_points)
    bk_left_sub = make_array(n_points, value=B)
      bk_left_sub *= phtn_on & bk_left_sub *= EM
      bk_left_sub += random_array * sqrt(bk_left_sub * 2.0)
    random_array = randomN(seed, n_points)
    bk_right_sub = make_array(n_points, value=B)
      bk_right_sub *= phtn_on & bk_right_sub *= EM
      bk_right_sub += temporary(random_array) * sqrt(bk_right_sub * 2.0)
    bk_left_sub = ts_smooth(bk_left_sub, smth) & bk_right_sub = ts_smooth(bk_right_sub, smth)
    I_left_anal -= bk_left_sub & I_right_anal -= bk_right_sub
  ENDIF
  
  ;convert to counts? electron-to-counts relation had been previously calculated by TKH
  I_left_anal *= 10 & I_right_anal *= 10
  
  plot, I_left_anal[0:n_elements(I_left_anal)-1]
  oplot, I_right_anal[0:n_elements(I_right_anal)-1], COLOR='FFCC66'x
  ;oplot, bk_left_sub[0:n_elements(I_left_anal)-1] - mean(bk_left_sub), COLOR='FF0000'x
  ;oplot, bk_left_sub[0:n_elements(I_left_anal)-1], COLOR='FF0000'x]

;------------------------------------------------------------------------------
;NOW START SIGNAL PROCESSING: LINEAR DICHROISM
  
  IF corrINT EQ 1 THEN BEGIN
    print, '...correcting intensity data for LD'
    bad_L=where(I_Left_anal LE 0.00000, count_L)
    bad_R=where(I_right_anal LE 0.00000, count_R)
    IF count_L GT 0 THEN I_left_anal[bad_L] = 0.0
    IF count_R GT 0 THEN I_right_anal[bad_R] = 0.0
    LD = (I_left_anal - I_right_anal) / (I_left_anal + I_right_anal) 
  ENDIF ELSE print, '...not correcting intensity data for LD'

  IF thresh EQ 0 THEN BEGIN
    LD = (I_left_anal - I_right_anal) / (I_left_anal + I_right_anal)
    print, '...no threshold...'
    ENDIF ELSE BEGIN

    LD = (I_left_anal - I_right_anal) / (I_left_anal + I_right_anal)
    I_tot_anal = (I_left_anal + I_right_anal) / 2.0
  
    IF thresh EQ 1 THEN BEGIN
      threshold = max(I_tot_anal[t_blch+1:*])-10
      print, '...flat threshold @', threshold
    ENDIF
  
    IF thresh EQ 2 THEN BEGIN
      ;equiv to setting thresh to smoothed intensity, mapping to actual intensity
      i_tot_sm = ts_smooth(i_tot_anal, 350)
      threshold = max(I_tot_sm[t_blch+350:*]) + 10
      print, '...smoothed intensity threshold @', threshold 
    ENDIF
  
    ;20140430...just to test something HZ vs SI THR...how gap-preserving is each technique?
    ;get rid of this when i am done...
    ;i_tot_anal[t_blch:*]=!values.f_nan
    ;LD[t_blch:*]=!values.f_nan
    ;end test...
  
    noise=where(I_tot_anal lt threshold, n_count)
    IF n_count NE 0 THEN LD[noise] = !values.f_nan
  ENDELSE

  ;LD = LD[0:T_BLCH-5]
  ;LD[T_BLCH:*]=!values.f_nan
  
  ;print, mean(I_left_anal[0:t_blch-1]), mean(I_right_anal[0:t_blch-1])
  ;print, stddev(I_left_anal[t_blch+1:*]), stddev(I_right_anal[t_blch+1:*])
  
  ;free various things...speeds up program!
  B=0 & t_blchs=0 & SNBs =0 & adds=0 & subs=0 & threshs=0

END
