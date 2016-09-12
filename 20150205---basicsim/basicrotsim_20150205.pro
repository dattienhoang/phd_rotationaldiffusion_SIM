;basicrotsim_20150205.pro -- 
;main routine for generating and then manipulating rotational diffusion
;data
;
;Dat Tien Hoang, Last major edit: 2015-02-05


PRO basicrotsim_20150205

print, 'BEGIN basicrotsim_20150205'

  ;basic parameters for the simulation
  new_traj = 0 ;[0 = load previously generated trajectories from *.txt files
               ; 1 = generate new trajectories using parameters defined below]
    n_molecules = 500  ;the amount of trajectores to load or generate
    N_steps = 100000   ;number of points in trajectory
    ;we typically do not change the following parameters for several reasons...
    ;  d_rot is a physical constant, and so is defined as such
    ;  tau is kept at 100-steps because it is a convenient number
    ;  frame_time is 1-sec just for unit purposes
    tau = 100L
    d_rot = 1.5/(tau*6.0)
    frame_time = 1.0
  ;post-processing parameters for the simulation
  add_noise = 1 ;[0 = do not add noise, do not assign B or S
                ; 1 = add noise, assign B and S (defined in subprocedure)]
  compression = 1 ;number of frames to skip from original trajectory
                  ;useful for testing frame-rate effects (ie sampling rate)
                  ;if not skipping any, then compression=1, else it is the 
                  ;number of points you want to skip
  blch = 1 ;[0 = do not add photobleaching
           ; 1 = photobleach...parameters set in brs_noise_20150205.pro]

  IF new_traj EQ 1 THEN BEGIN
    ;you will be prompted about several things:
    ;1.) The folder to place the new trajectories
    ;2.) A "previously saved seed file". IDL uses seed files to generate 
    ;    random numbers, such as those required for a random walk simulation. 
    ;    If you don't have one, then select cancel and a seed will self-generate.
    ;3.) A previous location...all new trajectories begin from the same coordinates. If you want to begin somewhere random, then a good way is to generate f
    print, '...generating new trajectories'
    ;tkh_rot, d_rot = d_rot, N_steps = N_steps, n_traj = n_molecules
    print, '......d_rot, n_steps, n_molecules', d_rot, n_steps, n_molecules
    tkh_rot_20130513,  d_rot = d_rot, N_steps = N_steps, n_traj = n_molecules
    ;tau_log_normal_20130513, 2.0, 0.0, n_molecules, tau_distrib 
    print, '...finished generating new trajectories'
  ENDIF
  
  ;this creates an IDL widget to prompt you to pick a *.dat file containing 
  ;a trajectory. if you kept the filenames as they were generated, and also 
  ;kept the files in the same folder, then the next FOR loop later will know 
  ;how to handle the files and then read all the files on order...so it is
  ;IMPORTANT TO NOT CHANGE THE FILENAMES OF THE TRAJECTORY *.DAT FILES and to
  ;KEEP ALL THE TRAJECTORIES IN THE SAME FOLDER.
  file =dialog_pickfile()
  file_name = file_basename(file,'.dat')
  file_dir = file_dirname(file)
  v_file = file_info(file)
  print, '...loading trajectories from folder:', file_dir

  ;do not change these parameters...they are here to help with reading
  ;large trajectories fast!
  N_Points = floor(n_steps/compression)
  switch_point = 10
  tau_ratio = 1.0
  ;finish =n_points
  
  ;arrays to temporarily hold the data as it is being read in
  x=fltarr(n_points)
  y=fltarr(n_points)
  z=fltarr(n_points)
  ;time=findgen(floor(n_points))

  FOR mol = 1L, n_molecules DO BEGIN
    starttime = systime(/seconds)
    
    IF mol MOD 10 EQ 0 THEN print, '......reading file ', $
      strcompress(string(mol) + '/' + string(n_molecules), /remove_all)
    file_full_path = strcompress(file_dir + '/trajectory' + string(mol) + '.dat', /remove_all)
    v_file = file_info(file_full_path)
    steps = double(v_file.size/12) ; 4 bytes per point, 3 points per step --> 12 bytes per step
    openr, lun, file_full_path, /get_lun 
    array = fltarr(3,steps)
    readu, lun, array
    close, lun & free_lun, lun
    
    array = temporary(array[*,0:*:compression]) ;now you have an trajectory loaded into the array
                                                ;into the procedure and you can do whattever you want with it.
    
    ;reformat the coordinates in the large array into individual arrays for easier handling
    x[0:switch_point] = array[0,0:switch_point]
    y[0:switch_point] = array[1,0:switch_point]
    z[0:switch_point] = array[2,0:switch_point] 
    FOR k=0L, n_points - switch_point - 1 DO BEGIN
      x[switch_point + k] = array[0, switch_point + k*tau_ratio]
      y[switch_point + k] = array[1, switch_point + k*tau_ratio]
      z[switch_point + k] = array[2, switch_point + k*tau_ratio]
    ENDFOR
    ;this subroutine performs the calculation (according to Fourkas) to 
    ;determine relative intensities in orthogonal polarization channels given 
    ;moleular dipole position. Here the NA is set to 0.75. The focusing is set
    ;to 0.001 since we are using wide-field, but should be set to 0.75
    ;if you are using confocal
    print, '......convert positions to intensities via Fourkas'
    brs_intensity_20150205, x,y,z, 0.75, 0.001, I_left_anal, I_right_anal
    
    ;this is where we add noise and stuff if desired. you set parameters 
    ;in this subroutine rather than here
    IF add_noise EQ 1 THEN BEGIN
      print, '......adding noise and assigning SBR'   
      ;add_noise_camera_20130513, I_left_anal, I_right_anal, LD, thresholding = 1
      brs_noise_20150205, I_left_anal, I_right_anal, I_tot_sm, LD, threshold, $
                          thresh, blch=blch
    ENDIF ELSE LD = (I_left_anal - I_right_anal)/(I_left_anal + I_right_anal)
    
    ;this part here is where I would usually do the part that is relevant to 
    ;my work, such as printing the trajectory coordinates to a *.txt file so I 
    ;can plot it elsewhere, or calculate an autocorrelation function. I removed
    ;this part because I am not sure about your application
    ;however, I templated a file printing code right here:
    fileparts = strcompress(file_dir + 'summaryfile_' + string(mol) + '.txt')
    openw, lun, fileparts, /get_lun
    printf, lun, 'rotation simulation of:'
    printf, lun, '...', file_dir
    printf, lun, 'x(t), y(t), z(t)'
    FOR i=0L, n_points-1 DO printf, lun, x[i], ',', y[i], ',', z[i]
    close, lun & free_lun, lun
  ENDFOR
  
  ;this part here is where I would print out some kind of summary file 
  fileparts = strcompress(file_dir + 'summaryfile' + '.txt')
  openw, lun, fileparts, /get_lun
  printf, lun, 'summarized data here...'
  close, lun & free_lun, lun

print, 'END of basicrotsim_20150205'
END
