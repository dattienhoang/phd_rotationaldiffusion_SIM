PRO SM_Simulations_tau_beta_20130513

new_traj = 0
add_noise = 0
n_molecules = 1000
N_steps = 1000L;400000L 
windowstep = 10
window_width = 20 ;width of window in terms of tau
stand_dev = fltarr(n_molecules) ; storage of standard deviation for each molecule
log_stand_dev = fltarr(n_molecules) ; storage of standard deviation for each molecule log(taus)
SKEW = FLTARR(n_molecules); storage of skewness for each molecule
kurt = fltarr(n_molecules); storage of kurtosis for each molecule
frame_time = 1.0
tau = 100L
tau_set = tau
d_rot = 1.5/(tau*6.0)
if new_traj EQ 1 then begin
tkh_rot_20130513, d_rot = d_rot, N_steps = N_steps, n_traj = n_molecules

tau_log_normal_20130513, 2.0, 0.0 , n_molecules, tau_distrib 

print, 'Calculation of trajectories finished'
endif

compression = 1;number of frames skipped from original trajectory

;result storage file


file =dialog_pickfile()
;file = '/home/stephan/stdev_test_tau_1000/trajectory1.dat'
file_name = file_basename(file,'.dat')
file_dir = file_dirname(file)
v_file = file_info(file)
print, file_dir

;this are adjustable parameters
;N_points = floor(n_steps/compression)
N_Points = floor(n_steps/compression);1000 ; total number of points
switch_point = 10
tau_ratio = 1.0
finish =n_points

x=fltarr(n_points)
y=fltarr(n_points)
z=fltarr(n_points)

time=findgen(floor(n_points))

acf_tau = fltarr(n_molecules)
median_tau = fltarr(n_molecules)
tau_array = fltarr(n_molecules,ceil(FINISH/windowstep)+1)
tau_cutoff = fltarr(n_molecules)
scale = 5 ;tile increment factor
tile = fltarr(n_molecules, scale)



;@@@@@@@@@@@@@@@@@@@THIS IS JUST FOR PRINTING LDs, ACFS, AND FIT RESULTS@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
tau = fltarr(n_molecules)
beta = fltarr(n_molecules)
betastretched = fltarr(n_molecules)
trajectorylength = fltarr(n_molecules)
acf_array = fltarr(n_molecules,n_elements(time))
taufit = fltarr(n_molecules)
normalizedacf = fltarr(n_molecules,n_elements(time))
prefactor = fltarr(n_molecules)
perror = fltarr(n_molecules,3)
;@@@@@@@@@@@@@@@@@@@THIS IS JUST FOR PRINTING LDs, ACFS, AND FIT RESULTS@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


;filename from tkh_rot increasing automatically by number trajectory1.dat, trajectory2.dat ....
for mol = 0L, n_molecules-1 do begin; n_molecules - 1 do begin; testnum -1 do begin;n_molecules -1 do begin; n_molecules - 1 do begin  ;loop over all molecules
print, 'molecule', mol

starttime = systime(/seconds)

file_full_path = strcompress(file_dir + '/trajectory' + string(mol+1) + '.dat',/remove_all)

;full filepath + filename is required
v_file = file_info(file_full_path)
steps = double(v_file.size/12) ; 4 bytes per point, 3 points per step --> 12 bytes per step
OpenR, lun, file_full_path, /get_lun 
array = fltarr(3,steps)
ReadU, lun, array
close, lun
free_lun, lun
print, 'end of read in file'

array = temporary(array[*,0:*:compression]) ;now you have an trajectory loaded into the array
                                            ;into the procedure and you can do whattever you want with it.

x[0:switch_point] = array[0,0:switch_point]
y[0:switch_point] = array[1,0:switch_point]
z[0:switch_point] = array[2,0:switch_point] 
for k=0L, n_points - switch_point - 1 DO BEGIN
x[switch_point + k] = array[0, switch_point + k*tau_ratio]
y[switch_point + k] = array[1, switch_point + k*tau_ratio]
z[switch_point + k] = array[2, switch_point + k*tau_ratio]
ENDFOR


intensity_calc_20130513, x,y,z, 0.75, 0.001, I_left_anal, I_right_anal
;x, y, and z coordinates of vector, NA = 0.75, 0.001 (no focusing bc of widefield, 0.75 for confocal), 
;program outputs intensity in right and left channels

if add_noise EQ 1 then begin  
print, 'Adding noise'   
add_noise_camera_20130513, I_left_anal, I_right_anal, LD, thresholding = 1 ;thresholding forces -1<LD<1
endif else LD = (I_left_anal - I_right_anal)/(I_left_anal + I_right_anal)


;Calculate average tau to pick window size
result=smf_plot_autocorrelation(time, LD,LD, maxframe=300, stretched='stretched')
;result=smf_plot_autocorrelation(time, LD,LD, maxframe=300, stretched='stretched')
;perror(mol,*) = result.gen_error
prefactor[mol] = result.fit_parameters[0]
trajec_length= n_points/(kww_tau(result.fit_parameters[1],result.fit_parameters[2])/frame_time)   
trajectorylength[mol] = trajec_length


taufit[mol] = result.fit_parameters[1]
print, "Tauc for entire traj for molecule", mol, kww_tau(result.fit_parameters[1],result.fit_parameters[2])
betastretched[mol] = result.fit_parameters[2]
print, 'forced stretched beta', betastretched[mol]
acf_array(mol,0:n_elements(result.acf)-1) = result.acf


filename = strcompress('/home/user/Lindsay/SimTraj_LDACF_trajlengtheq'+string(n_steps/tau_set)+'_mol'+string(mol+1)+'_022713.txt', /remove_all)
get_lun, lun
openW, lun, filename
for k=0L, 299 DO BEGIN ; starting this loop allows the arrary of numbers to be printed line by line
printf, lun, acf_array(mol,k), ','
endfor ; closes file writing for loop
close, lun ; closes file (need to do this or no copy will exist)
free_lun, lun ; frees existing lun assignment

;************************************
; smf_plot_autocorrelation.pro is a function which calculates the 
; autocorrelation function of each  molecules within a file(which contains many frames of data)
; pairs is a structure within the analysis_info structure.
; stretched is B not equal to 1
;dummy=where(finite(LD, points))
trajec_length= n_points/(kww_tau(result.fit_parameters[1],result.fit_parameters[2])/frame_time)   
trajectorylength[mol] = trajec_length
print,'trajectorylength', mol, trajectorylength[mol]


;if trajec_length LT 50 AND $
;  (kww_tau(result.fit_parameters[1],result.fit_parameters[2])/frame_time) GT 20 then begin
;  print, 'linear fit'
;  
;  result=SMF_Plot_Autocorrelation($
;            time,$
;            LD, $
;            LD,$
;            autocorrelation = result.acf,$
;            acf_err = result.acf_err,$
;           stretched = 'stretched',$
;           Method='DIETER_LD')            ; Calculation method     
;  result.fit_parameters[2] = 1.0000
;  print, 'Linear fit imposed'
;  trajec_length= N_points/(kww_tau(result.fit_parameters[1],result.fit_parameters[2])/frame_time)   
;  trajectorylength[mol] = trajec_length
;  print,'trajectorylength', mol, trajectorylength[mol]
;  
;endif

tau[mol] = kww_tau(result.fit_parameters[1], result.fit_parameters[2])
print, 'tau[mol]', tau[mol]
beta[mol]= result.fit_parameters[2]
print, 'beta[mol]', beta[mol]

endfor;endloop over molecules


good = where(finite(tau))
print, 'good molecules', n_elements(good)
median_tau = median(tau[good])
median_lifetime = median(trajectorylength(good))
median_beta = median(betastretched[good])

if add_noise eq 0 then begin
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
filename = strcompress('/home/user/Lindsay/SimHomoRotatingTraj_FORCEDSTRETCHED_trajlengtheq'+string(n_steps/tau_set)+$
          '_taueq'+string(tau_set)+'_compeq'+string(compression)+'NONOISE_030613.txt', /remove_all)
get_lun, lun
openW, lun, filename
printf, lun, 'Input tau = ', tau_set
printf, lun, 'Input fwhm = ', 0.0
printf, lun, 'Number of molecules = ', N_molecules
printf, lun, 'Output:'
printf, lun, 'median tau', median_tau
printf, lun, 'median beta (forced stretched)', median_beta
printf, lun, 'med_lifetime', median_lifetime
printf, lun, 'frames per tau', tau_set/compression/frame_time
printf, lun, 'tau, beta, beta str, lifetime/tau, logtau';
for x = 0, n_molecules-1 do begin
printf, lun, tau[x], ',', beta[x], ',', betastretched[x], ',', trajectorylength[x], ',', alog10(tau[x]), ','
endfor
close, lun ; closes file (need to do this or no copy will exist)
free_lun, lun ; frees existing lun assignment
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
ENDIF

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
filename = '/home/user/Lindsay/egPDI_simtraj_sametrajlength_120512.txt'
get_lun, lun
openW, lun, filename
printf, lun, 'No. good molecules', n_elements(good)
printf, lun, 'median tau (s) =', median_tau
printf, lun, 'median_lifetime (sec) =', median_lifetime
printf, lun, '  taufit    tau     str beta  lifetime/tau'
for k=0L, n_molecules - 1 DO BEGIN ; starting this loop allows the arrary of numbers to be printed line by line
printf, lun, taufit[k], ',',tau[k], ',',$
betastretched[k], ',', trajectorylength[k], ','; indicates array 1x (k+1) molecules
endfor ; closes file writing for loop
close, lun ; closes file (need to do this or no copy will exist)
free_lun, lun ; frees existing lun assignment
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
filename = '/home/user/Lindsay/egPDI_simtraj_varyingtrajlength_error_120512.txt'
good = where(beta ne 1, count)
print, 'good', good
help, tau
help, betastretched
get_lun, lun
openW, lun, filename
printf, lun, 'No. good molecules', n_elements(good)
printf, lun, 'median tau (s) =', median_tau
printf, lun, 'median_lifetime (sec) =', median_lifetime
printf, lun, '  taufit    tau     str beta  lifetime/tau'
for k=0L, count -1 DO BEGIN ; starting this loop allows the arrary of numbers to be printed line by line
print, 'good[k]', good[k]
printf, lun, good[k], ',',tau[good[k]], ',',$
betastretched[good[k]], ',', trajectorylength[good[k]], ','; indicates array 1x (k+1) molecules
endfor ; closes file writing for loop
close, lun ; closes file (need to do this or no copy will exist)
free_lun, lun ; frees existing lun assignment

filename = '/home/user/Lindsay/egPDI_simtraj_varyingtrajlength_PREFACTORerror_120512.txt'
good = where(beta ne 1, count)
get_lun, lun
openW, lun, filename
printf, lun, 'No. good molecules', n_elements(good)
printf,lun, 'molecule, prefactor, error'
for k=0L, count - 1 DO BEGIN ; starting this loop allows the arrary of numbers to be printed line by line
printf, lun, good[k], ',',prefactor[good[k]], ',', perror[good[k],0], ','
endfor ; closes file writing for loop
close, lun ; closes file (need to do this or no copy will exist)
free_lun, lun ; frees existing lun assignment

filename = '/home/user/Lindsay/egPDI_simtraj_varyingtrajlength_TAUerror_120512.txt'
good = where(beta ne 1, count)
get_lun, lun
openW, lun, filename
printf, lun, 'No. good molecules', n_elements(good)
printf, lun, 'molecule, tau, error'
for k=0L, count - 1 DO BEGIN ; starting this loop allows the arrary of numbers to be printed line by line
printf, lun, good[k], ',',tau[good[k]], ',', perror[good[k],1], ','
endfor ; closes file writing for loop
close, lun ; closes file (need to do this or no copy will exist)
free_lun, lun ; frees existing lun assignment

filename = '/home/user/Lindsay/egPDI_simtraj_varyingtrajlength_BETAerror_120512.txt'
good = where(beta ne 1, count)
get_lun, lun
openW, lun, filename
printf, lun, 'No. good molecules', n_elements(good)
printf, lun, 'molecule, beta, error'
for k=0L, count - 1 DO BEGIN ; starting this loop allows the arrary of numbers to be printed line by line
printf, lun, good[k], ',',beta[good[k]], ',', perror[good[k],2], ','
endfor ; closes file writing for loop
close, lun ; closes file (need to do this or no copy will exist)
free_lun, lun ; frees existing lun assignment
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;@@@@@@@@@@@@@@ BETA QE PARY @@@@@@@@@@@@@@@@@@@

help, acf_array
;print, 'acf_array(x,*)', acf_array[1,0:10]

for x = 0, n_molecules-1 do begin ;n_elements(acf_array[*,0])-1 do begin

  normacf = fltarr(n_elements(acf_array(x,*)))
  ;print, 'acf_array(x,*)', acf_array[x,0:10]
  ;print, 'acf_array[x,1]', acf_array[x,1]
  normacf = acf_array[x,*]/acf_array[x,1]
  normalizedacf[x, 0:n_elements(normacf)-2] = normacf[1:n_elements(normacf)-1]
 ;print, ' normalizedacf[0:4,j]',  normalizedacf[x,0:4]
  endfor
  


  str_ensemble_acf = fltarr(n_elements(time)-1)
  for frame = 0L, n_elements(time) - 2 DO BEGIN
  str_ensemble_acf[frame]=mean(acf_array(*,frame));mean((*info.recalc)[realstrtaubeta].acf1[frame], /nan)
  ENDFOR


print,  'mean(normalizedacf(0,frame))', mean(normalizedacf(*,0))
norm_ensemble_acf = fltarr(n_elements(time)-1) 
for frame = 0L, n_elements(normalizedacf(0,*)) - 2 DO BEGIN
norm_ensemble_acf[frame] = mean(normalizedacf[*,frame])
ENDFOR 

 ; Get rid of NaNs and zeros in ensemble acf
  goodacf = where(finite(str_ensemble_acf))
  str_ensemble_acf = str_ensemble_acf[goodacf]
  goodacf2 = where(str_ensemble_acf ne 0)
  str_ensemble_acf = str_ensemble_acf[goodacf2] 
  help, str_ensemble_acf
  ;------------------------------------------------------------------------------------------------
   
   acf_err1 = fltarr(n_elements(str_ensemble_acf))
   
    
  result=SMF_Plot_Autocorrelation_error($
            time,$
            LD, $
            LD,$
           autocorrelation = str_ensemble_acf,$
           
           acf_err = acf_err1, $
           stretched = 'stretched',$
           Method='BASIC_ACF') 

    ensemble_acf1 = result.acf
    ensemble_acf_fit = result.acf_fit
    str_ensemble_prefactor = result.fit_parameters[0]
    str_ensemble_tau = kww_tau(result.fit_parameters[1], result.fit_parameters[2])
    str_ensemble_beta = result.fit_parameters[2]  
    print, 'NO NORMALIZED beta QE', str_ensemble_beta
    perror_qe = fltarr(3)
    perror_qe = result.gen_error
   acf_err2 = fltarr(n_elements(norm_ensemble_acf))
   
    
  result=SMF_Plot_Autocorrelation_error($
            time,$
            LD, $
            LD,$
           autocorrelation = norm_ensemble_acf,$
           
           acf_err = acf_err2, $
           stretched = 'stretched',$
           Method='BASIC_ACF') 

    norm_ensemble_acf = result.acf
    norm_ensemble_acf_fit = result.acf_fit
    norm_str_ensemble_prefactor = result.fit_parameters[0]
    norm_str_ensemble_tau = kww_tau(result.fit_parameters[1], result.fit_parameters[2])
    norm_str_ensemble_beta = result.fit_parameters[2]  
    perror_norm = fltarr(3)
    perror_norm = result.gen_error

print, 'str acf beta', str_ensemble_beta
print, 'error in str', perror_qe
print, 'norm str acf bet', norm_str_ensemble_beta
print, 'error in norm', perror_norm

;filename = '/home/user/Lindsay/dpPDI_distro_normalized_ACFQE_112712.txt'
filename = strcompress('/home/user/Lindsay/SimulatedHomoRotatingTraj_BETAQEinfo_trajlengtheq'+string(n_steps/tau_set)+$
          '_taueq'+string(tau_set)+'_compeq'+string(compression)+'_021813.txt', /remove_all)
get_lun, lun
openW, lun, filename
printf, lun, 'NON NORMALIZED'
printf, lun, 'ensemble prefactor', str_ensemble_prefactor
printf, lun, 'ensemble tau', str_ensemble_tau
printf, lun, 'ensemble beta',str_ensemble_beta
printf, lun, ''
printf, lun, 'NORMALIZED'
printf, lun, 'ensemble prefactor', norm_str_ensemble_prefactor
printf, lun, 'ensemble tau', norm_str_ensemble_tau
printf, lun, 'ensemble beta', norm_str_ensemble_beta
printf, lun, ''
printf, lun, 'non normalized'
printf, lun, 'acf,    fit'
FOR q = 0, N_elements(ensemble_acf1) -1 DO BEGIN
printf, lun, ensemble_acf1[q], ',', ensemble_acf_fit[q], ','
ENDFOR
printf, lun, 'normalized'
printf, lun, 'acf,    fit'
FOR m = 0, n_elements(norm_ensemble_acf) -1 DO BEGIN
printf, lun, norm_ensemble_acf[m], ',', norm_ensemble_acf_fit[m], ','
ENDFOR
close, lun 
free_lun, lun

;filename = '/home/user/Lindsay/dpPDI_distro_str_ensemble_ACFQE_112712.txt'
;get_lun, lun
;openW, lun, filename
printf, lun, 'ensemble prefactor', str_ensemble_prefactor
printf, lun, 'ensemble tau', str_ensemble_tau
printf, lun, 'ensemble beta', str_ensemble_beta
printf, lun, 'acf,    fit'
FOR q = 0, N_elements(ensemble_acf1) -1 DO BEGIN
printf, lun, ensemble_acf1[q], ',', ensemble_acf_fit[q], ','
ENDFOR
close, lun 
free_lun, lun

;@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
;@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
;@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
if mol eq 0 then begin
filename = '/home/user/Lindsay/dpPDI_Simulatedhomorotatingtraj_acf_and_fit_molecule11.txt'
;filename = info.file + '_tau_fit1.txt
get_lun, lun
openW, lun, filename
printf, lun, 'Molecule number', mol
printf, lun, 'tau = ', kww_tau(result.fit_parameters[1], result.fit_parameters[2])
printf, lun, 'prefactor', result.fit_parameters[0]
printf, lun, 'taufit = ', result.fit_parameters[1]
printf, lun, 'beta = ', result.fit_parameters[2]
printf, lun, 'lifetime/tau', trajectorylength[mol]
printf, lun, 'time, acf, acf fit'
for point = 1, 299 do begin
printf, lun, time[point], ',', result.acf[point], ',', result.acf_fit[point], ','
endfor ; closes file writing for loop
close, lun ; closes file (need to do this or no copy will exist)
free_lun, lun ; frees existing lun assignment
endif
;@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
;@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
;@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@









;@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


for mol = 0L, 0 do begin;n_molecules - 1 DO BEGIN ;loop over all molecules
;filename from tkh_rot increasing automatically by number trajectory1.dat, trajectory2.dat ....

starttime = systime(/seconds)

file_full_path = strcompress(file_dir + '/trajectory' + string(mol) + '.dat',/remove_all)
;full filepath + filename is required
v_file = file_info(file_full_path)
steps = double(v_file.size/12) ; 4 bytes per point, 3 points per step --> 12 bytes per step

OpenR, lun, file_full_path, /get_lun
array = fltarr(3,steps)
ReadU, lun, array
close, lun
free_lun, lun
print, 'end of read in file'

array = temporary(array[*,0:*:compression])

x[0:switch_point] = array[0,0:switch_point]
y[0:switch_point] = array[1,0:switch_point]
z[0:switch_point] = array[2,0:switch_point]
 
for k=0L, n_points - switch_point - 1 DO BEGIN
x[switch_point + k] = array[0, switch_point + k*tau_ratio]
y[switch_point + k] = array[1, switch_point + k*tau_ratio]
z[switch_point + k] = array[2, switch_point + k*tau_ratio]
ENDFOR


intensity_calc_20130513, x,y,z, 0.75, 0.001, I_left_anal, I_right_anal

;LD = (I_left_anal - I_right_anal)/(I_left_anal + I_right_anal)
  
if add_noise EQ 1 then begin  
print, 'Adding noise'   
add_noise_camera_20130513, I_left_anal, I_right_anal, LD, thresholding = 1
endif

;Calculate average tau to pick window size
result=smf_plot_autocorrelation(time, LD,LD, maxframe=300, stretched='stretched')
print, 'molecule =', mol, ', tau of entire trajectory:', result.fit_parameters
start = 0L
acf_tau[mol] = result.fit_parameters[1]
print, acf_tau[mol]
acf_max = floor(window_width * result.fit_parameters[1]) ;this is the windowsize 
 ;sliding window step

if finite(result.fit_parameters[1]) EQ 1 then begin 
maxframe = floor(2.0*result.fit_parameters[1])      

N_segment=ceil(FINISH/windowstep)+1 ;windowstep and finish are indices for arrays
temptau=make_array(N_segment,/float, value=!values.F_NAN) ;here the result is save for every window step
tempbeta=make_array(N_segment,/float, value=!values.F_NAN)
tempacf=make_array(N_segment,maxframe,/float, value=!values.F_NAN)
temp_acf_fit=make_array(N_segment,maxframe,/float, value=!values.F_NAN)
j=0L


FOR i=START, FINISH - acf_max-1, windowstep DO BEGIN; START AND FINISH DETERMINE THE START AND FINISH POINTS OF USED
                                                     ; trajectory for calculation

temp = smf_plot_autocorrelation(time[0:acf_max],$
                                LD[i:i+acf_max],$
                                LD[i:i+acf_max],$
                                maxframe = maxframe,$
                                method = 'BASIC_ACF',$
                                stretched='stretched')

Result = SMF_Plot_Autocorrelation(time[0:acf_max],$
          LD[i:(i+acf_max)], $
          LD[i:(i+acf_max)],$
          stretched = 'stretched',$
          Method='Dieter_LD', $            ; Calculation method
          autocorrelation=temp.acf,$
          acf_err = temp.acf_err)
dummy=where(finite(result.acf), count1); Check that ACF was calculated
dummy=where(finite(result.acf_fit), count2); ensure that acf was fitted

;plot, temp.acf, psym = 1
;oplot, temp.acf_fit
;oplot, result.acf_fit
;print, 'stretched', temp.fit_parameters, result.fit_parameters[1]
;wait, 3.0


if count1 GT 1 THEN BEGIN
tempacf[j,0:maxframe-1]=result.acf[0:maxframe-1]
  if count2 GT 1 THEN BEGIN; ensure that acf was fitted
  temp_acf_fit[j,0:maxframe-1] = result.acf_fit[0:maxframe-1]
  temptau[j]=result.fit_parameters[1] ;save result for every window
  
 
;exclude values where the number of available points is less then 30% of the theoretically possible number. 
;on average available points is 50% of theoretically possible.
dummy = where(LD[i:i+acf_max-1], count)
if count LT 0.33 * acf_max then temptau[j] = !values.f_NaN
  ;the linear fit cannot be 3 times bigger than the stretched exponential fit
  if result.fit_parameters[1] GT (3.0 * temp.fit_parameters[1]) then temptau[j] = !values.f_NAN
  ENDIF
ENDIF
;ENDIF; end of if loop to check whether acf was fitted by exponential function
j++

ENDFOR ;for loop over start finsih trajectory

print, 'j'
goodtemptau = where(finite(temptau))
if n_elements(goodtemptau) gt 1 then begin
help, temptau
help, goodtemptau
temptaugood = fltarr(n_elements(goodtemptau))
temptaugood = temptau[goodtemptau]
endif 

;get_lun, lun
;openw, lun, strcompress('/home/user/temptautest_'$
;                +string(mol) +'.txt')
;for tau = 0L, n_elements(temptaugood) -1 DO PRINTF, lun, tau, ',', temptaugood[tau], ','
;close, lun
;free_lun, lun


good=where(finite(temptau), count)
if count GT 1 then begin
median_tau[mol] = median(temptau[good])

mom_result = moment(temptau, /nan)
stand_dev[mol] = sqrt(mom_result[1])

log_mom_result = moment(alog10(temptau), /nan)
log_stand_dev[mol] = sqrt(log_mom_result[1])

SKEW[MOL] = mom_result[2]
KURT[MOL] = mom_result[3]
tau_array[mol,*] = temptau

sort_taus=sort(temptau)
tau_cut = temptau[sort_taus[ceil(0.95*count-1)]] ;95% OF ALL TAUS ARE BELOW TAU_CUTOFF
tau_cutoff[mol] = tau_cut/median_tau[mol]
;plot, temptau
;print, 'moment', mom_result
;print, 'taucut', tau_cut
;print, 'median_tau', median_tau[mol]
;SKEWNESS AND kURTOSIS

endif ;check count

print, '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
print, 'molecule', mol
print, 'median tau', median_tau[mol]
print, 'tau cutoff (95%)', tau_cutoff[mol]
print, 'standard devation of taus:', stand_dev[mol]
print, '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'

ENDIF; check whether fit_parameter finite
endtime = systime(/seconds)
print, 'time =', endtime - starttime

N_segment=0
temptau=0
tempbeta=0
tempacf=0
temp_acf_fit=0



ENDFOR; end of loop over molecules

;find 95% confidence for tau cutoff and standard deviation
sortcutoff = sort(tau_cutoff)
taucut_confidence = tau_cutoff(sortcutoff[ceil(0.95*n_elements(tau_cutoff)-1)])
print, 'taucut_confidence', taucut_confidence
sort_logstdev = sort(log_stand_dev)
logstdev_confidence = log_stand_dev(sort_logstdev[ceil(0.95*n_elements(log_stand_dev)-1)])
sort_stdev = sort(stand_dev)
stdev_confidence = stand_dev(sort_stdev[ceil(0.95*n_elements(stand_dev)-1)])

; normalize all stand deviations to the median tau
norm_stdev = stand_dev/median_tau 
sort_norm_stdev = sort(norm_stdev)
norm_stdev_confidence = norm_stdev(sort_norm_stdev[ceil(0.95*n_elements(norm_stdev)-1)])



get_lun, lun
openw, lun, 'stdev_tau_corr_061112.txt'
printf, lun, 'tau cut off confidence', taucut_confidence
printf, lun, '95% confidence for st. dev normalized:', norm_stdev_confidence
printf, lun, 'stdev_confidence', stdev_confidence
printf, lun, 'logstdev_confidence', logstdev_confidence
printf, lun, 'med tau, total traj tau, stdev, stdev for tau traj, tau cut off value'
for mol = 0L, n_molecules - 1 DO printf, lun, median_tau[mol], ',', acf_tau[mol],',', stand_dev[mol], ',', tau_cutoff[mol]
printf, lun, 'log stand dev of tau traj for individual molecules'
printf, lun, 'normalized st devs to the median of each tau trajectory:'
for mol2 = 0L, n_molecules -1 DO PRINTF, lun, mol2, ',', stand_dev[mol2]/median_tau[mol2], ','
close, lun
free_lun, lun




save, tau_array, filename = 'sim_save.sav'
print, 'PRO stdev_testing finished'
END