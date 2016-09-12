;+
; :Author: Toby
;-
;This program is designed to simulate the rotational diffusion of a vector
;in 3-D.  It uses quaternions to simulate each rotation, and records the 
;position of the vector at each point in its trajectory.

;The quaternion rotations are accomplished using the library written by 
;Craig Markwardt.

;Since each trajectory is quite long, the data is stored directly into 
;files as it is generated, in 100,000 step chunks (to limit RAM usage).  
;The user selects a destination folder and the filenames automtically 
;iterate within that folder.  Note that there is currently no provision to
;check for the existence of previous files/trajectories, so be careful
;where you tell the program to save its data

;Vector information is stored in the form [Vx, Vy, Vz] where these three 
;components define the direction of a unit vector.  This is redundant (it
;could be more briefly be stored as angles [theta, phi]) but is the form
;of the raw data produced by the QTROT.pro procedure.  

;The results are stored as single precision floating point data, although
;all computations within the program are performed using double precision.
;During computation I am very concerned about rounding errors, but there 
;is no need to store the final positions with the same precision since no
;additional computations will require the same degree of precision as the 
;serial rotation calculations used in this procedure.  In addition, to save
;space, one could just save a subset of the position vectors but I haven't
;added that yet.  I suspect that for 2 degree rotations, every tenth point 
;wouldbe more than sufficient

;The random numbers used to select the axis of rotation are incredibly 
;important, and I am not wholly happy relying on the built-in IDL random
;number generator.  It should be okay for the simulations here as long as
;I use the '/DOUBLE' flag,  but to be certain I should probably replace 
;it with a more up-to-date RNG.


;N_TRAJECORIES
  ;The number of simulated trajectories.  Each one is stored in a seperate
  ;file whose name includes the number of that trajectory.  The beginning
  ;orientation of the n+1 trajectory is the last orientation of the n 
  ;trajectory.  If undefined, it is assumed N_Trajectories=1.
  
;SEED
  ;The seed used for the RNG.  It can be left undefined or read in as a 
  ;SAV file from the previous run.  By saving the seed afer each run I
  ;hope to minimize problems with repeating a 'random' sequence, and 
  ;screwing up the uncorrelated nature of my random jumps.
  
;ROTATION_ANGLE  
  ;The angle of an individual jump, in radians.  The default value corresponds
  ;to about 2 degree jumps.  Either the ROTATION_ANGLE or D_ROTATION parameter  
  ;should be used, but not both.
  
;D_ROTATION
  ;The rotational diffusion constant used to generate the distribution of
  ;rotational step angles.  Either the ROTATION_ANGLE or D_ROTATION parameter
  ;should be used, but not both.
  
;N_STEPS
  ;The number of steps in a single trajectory.  The data file for each 
  ;trajectory will have 3*N_Steps values sored in a float array.  This should
  ;be a multiple of 100,000 -- if it isn't it will be rounded up to the next
  ;multiple of 100,000 steps.
  
;WAITTIME
  ;This is an optional wait command that interrupts the code every 100,000
  ;points to wait for WAITTIME seconds.  The purpose of this is to allow my 
  ;computer time to process other things going on, since its resource allocation
  ;sucks and I don't want to bother making it better.
  
;START_NUMBER
  ;If a number of trajectories already exist, and you wish to continue the 
  ;calculations to add more, this parameter allws you to enter the number 
  ;corresponding to the first trajectory file the program will save.  Thus, if
  ;you have already calculated fifty trajectories, setting START_NUMBER=51 will
  ;ensure that the first trajectory calculated will be saved as TRAJECTORY51.dat
  
  
PRO TKH_Rot_20130513, N_trajectories = N_trajectories, $ ;Number of output files/trajectories
              Seed=Seed, $                      ;seed for RNG -- read in as SAV file
              rotation_angle= rotation_angle, $  ;rotation step in radians
              D_rotation= D_rotation, $          ;Rotational diffusion constant
              N_Steps = N_Steps, $              ;steps in a single trajectory
              Waittime = waittime, $            ;Delay per 100,000 steps, in seconds
              Start_Number = Start_Number;,$      ;trajecory start number
              ;file = file, $ ;filename where to save trajectories
              ;V_FILE = V_FILE,$ ;filename for last trajectory
              ;seedfile = seedfile ;seed for RNG -- read in as SAV file

;print, 'd_rot (tkh_rot.pro) d_rot = ', d_rot
print, 'd_rotational (tkh_rot.pro) d_rotational = ', d_rotation
;print, 'n_traj (tkh_rot.pro) n_trajectories =', n_traj
print, 'n_trajectories (tkh_rot.pro) n_trajectories =', n_trajectories

;If number of trajectories is not indicated, then assume a value of 1
IF N_Elements(N_trajectories) EQ 0 THEN begin
N_trajectories=1
print, 'n_trajectories was not specified!'
endif

;If number of steps is not indicated, then assume a value of 100,000
IF N_Elements(N_Steps) EQ 0 THEN N_Steps=500000

;If no rotation angle or diffusion constant given, assume a 
;jump of about 2 degrees
IF (N_Elements(rotation_angle)+ N_Elements(D_Rotation)) EQ 0 $
     THEN begin
     rotation_angle = 1.2*!DPI/180 ;0.034906585D
     print, 'd_rotation was not specified!'
     endif

;If no "Start_Number" is given, assume the user wishes to start at 
;trajectory 1
IF N_Elements(Start_Number) EQ 0 THEN Start_Number = 1

;Choose a folder to hold the data files
file = dialog_pickfile(/directory, /write, $
      Title='Please choose a folder for the data file(s)')
      
      
;If the user wishes to seed the RNG with a given seed, then locate the save file
;Otherwise leave SEED undefined and allow IDL to use the generic state to seed
;the RNG.
seedfile = dialog_pickfile(/read, $
 title='If you have previously saved the seed variable, please select the SAV file.',$
 filter='*.sav')
IF seedfile NE '' THEN RESTORE, seedfile 



;If user wishes to continue from the last know position of the orientation vector,
;then allow user to read in the last known position from a previous trajecory file
V_File = dialog_pickfile(/read, filter = '*.dat', $
 title = 'If you wish to continue from a previous trajecory, select the V_Out file.')

;Reading in the vector file -- check that it exiss, then find the number of points in 
;it by dividing its total size by 12 (4Bytes/#, 3#/timepoint)
IF V_File NE '' THEN BEGIN 
  V_File_Info = File_Info(V_File)
  IF V_File_Info.exists EQ 1 THEN BEGIN
     V_Size = V_File_Info.size/12 ; Find # of points in file
     ;To avoid calculation errors, use the fact that it must be in the form
     ;100 000 x n +1
     V_Size = round(V_Size/100000.0, /L64)* 100000L + 1
     OpenR, lun, V_File, /get_lun
     V_Old = Assoc(lun, fltarr(3,1))
     V_Init = V_Old[V_Size-1]     ;Take the last point as the new starting point
     Close, lun
     Free_lun, lun
     V_Old = 0 ;clear variable
  ENDIF  ;of file exists
ENDIF  ; of file selected
  
;If no file is read in, then set initial vector position to x=y=z=0.57735027
IF N_Elements(V_Init) EQ 0 THEN V_init = [0.57735027,0.57735027,0.57735027]  
print, V_Init
     
;Initialize array to hold the orientation vector information     
Vout = dblarr(3,100001)
Vout[*,0] = V_Init

;Introduce a counter to keep track of the trajectory currently being simulated
n=0


;************************************************
;BEGINNING OF LOOP STRUCTURES  -- THE HEART OF THE CODE HAPPENS BELOW

WHILE n LT N_trajectories DO BEGIN

;Record systime for probing calculation speed
starttime = systime(1)


;Set up file for the trajectory
 filename = strcompress(file + 'trajectory' + string(n+round(Start_Number)) + '.dat',/remove_all)
;filename = strcompress(file + 'trajectory' + string(start_number) + '.dat',/remove_all)

;Open data file for this trajectory
OPENW, lun, filename, /get_lun


;*********
;Now calculate trajectory in 100,000 point chunks
FOR step = 0, fix((N_Steps-1.)/100000) DO BEGIN
;Generate one hundred thousand random number pairs with double precision 
;for each trajectory.  These will be used to define rotation axis.
;randnumbers = randomu(seed, 2, 100000, /double, /uniform)
randnumbers = randomu(seed, 2, 100000, /double)
;Define a matrix to hold the angles for the rotation axes
angles = dblarr(2,100000)

;Define masks for the theta and phi angles
phi_mask = findgen(2,100000)/2
phi_mask = (phi_mask - long(phi_mask)) LT 0.2
theta_mask = 1 - phi_mask

;print, phi_mask
;Compute angles for rotation axes using random numbers
;Phi angles are uniformly distributed over 0 < phi < 2*pi
;Theta angles have a probablility ditribution of sin(theta)

;angles= 2*!dpi*randnumbers*phi_mask + (acos(2*randnumbers - 1))*theta_mask
angles[0,*] = 2*!dpi*randnumbers[0,*]
angles[1,*] = (acos(2*randnumbers[1,*] - 1))
;angles[0,*] = randnumbers[0,*] * 360
;angles[1,*] = !RADEG*asin(2*randnumbers[1,*] -1)
;map_set, 0.5*!DPI,0,-0.5*!DPI,/lambert,/isotropic
;plots, angles[0,*], angles[1,*], psym = 3
;wait, 300
;Calculate the size of the rotational jump.  
;If the ROTATION_ANGLE keyword was used, this will be a constant value
;if the D_ROTATION keyword was used, this will be a distribution of
;jump sizes.
IF N_Elements(D_ROTATION) NE 0 THEN $
     Rotation_Jumps = TKH_ROT_JUMP_DISTRIBUTION_20130513(D_ROTATION, N_POINTS=100000) $
   ELSE Rotation_Jumps = make_array(100000,value=rotation_angle)
;help, rotation_jumps
;print,'rotation_jumps[0,10]*180/3.14', rotation_jumps[0:10]*180/3.14

;filename = '/home/user/Lindsay/RotationJumps_taueq100_compressioneq1_38xtaulength.txt'
;filename = info.file + '_tau_fit1.txt
;get_lun, lun
;openW, lun, filename
;for x = 0, 3800 do begin ; 3800/5=760
;printf, lun, rotation_jumps[x]*180/3.14
;endfor
;close, lun ; closes file (need to do this or no copy will exist)
;free_lun, lun ; frees existing lun assignment



;Calculate quaternion axis vectors from rotation axis angles then calculate
;quaternion
Quaternionaxis = dblarr(3,100000)
Q = dblarr(4,100000)

Quaternionaxis[0,*] = cos(angles[0,*])*sin(angles[1,*])
Quaternionaxis[1,*] = sin(angles[0,*])*sin(angles[1,*])
Quaternionaxis[2,*] = cos(angles[1,*])
Q =  QTCOMPOSE(Quaternionaxis, Rotation_Jumps)


;Apply quaternion roation and record position one step at a time

FOR i=0L,99999 DO BEGIN ;100,000 iterations
  ;Calculate new orientation
  Vout[*,i+1] = QTVROT(Vout[*,i], Q[*,i])
  ;Renormalize to remove rounding errors
  IF (i+1) mod 100 EQ 0 THEN Vout[*,i+1] = Vout[*,i+1]/(total(Vout[*,i+1]^2))
ENDFOR

IF N_Elements(waittime) NE 0 THEN wait, waittime

;print, 'Loop #', step, ' took ', systime(1)-starttime, ' seconds.'
WriteU, lun, float(Vout)

;move final orientation up to initial position for next iteration
Vout[*,0] = Vout[*,100000]

ENDFOR ;(N_Steps)
;**********************
;
;Close data file for this trajectory
Close, lun & free_lun, lun
print, filename, systime(1)-starttime , ' seconds'

txtfile = strcompress(file + 'trajectory_coords' + string(n+round(Start_Number)) + '.txt',/remove_all)
get_lun, lun
openW, lun, txtfile
for f = 0L, N_Steps/4 -1 do begin
printf, lun, VOut[*,f]
endfor
for f = N_Steps/4L, N_Steps/2-1 do begin
printf, lun, VOut[*,f]
endfor
for f = N_Steps/2L, (N_Steps*3/4)-1 do begin
printf, lun, VOut[*,f]
endfor
for f = N_Steps*3/4L, N_Steps-1 do begin
printf, lun, VOut[*,f]
endfor
close, lun ; closes file (need to do this or no copy will exist)
free_lun, lun 

;Save final state of the seed variable in a .SAV file.  Do this after 
;each trajectory in case the program crashes while running
seedfile = strcompress(file + 'seedfile.sav', /remove_all)
Save, seed, filename = seedfile

n = n+1
ENDWHILE  ;For trajecories


END
