!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!! NEMO/OPA  :  1 - run manager      (namrun)
!! namelists    2 - Domain           (namcfg, namzgr, namzgr_sco, namdom, namtsd)
!!              3 - Surface boundary (namsbc, namsbc_ana, namsbc_flx, namsbc_clio, namsbc_core, namsbc_sas
!!                                    namsbc_cpl, namtra_qsr, namsbc_rnf,
!!                                    namsbc_apr, namsbc_ssr, namsbc_alb)
!!              4 - lateral boundary (namlbc, namcla, namagrif, nambdy, nambdy_tide)
!!              5 - bottom  boundary (nambfr, nambbc, nambbl)
!!              6 - Tracer           (nameos, namtra_adv, namtra_ldf, namtra_dmp)
!!              7 - dynamics         (namdyn_adv, namdyn_vor, namdyn_hpg, namdyn_spg, namdyn_ldf)
!!              8 - Verical physics  (namzdf, namzdf_ric, namzdf_tke, namzdf_kpp, namzdf_ddm, namzdf_tmx, namzdf_tmx_new)
!!              9 - diagnostics      (namnc4, namtrd, namspr, namflo, namhsb, namsto)
!!             10 - miscellaneous    (namsol, nammpp, namctl)
!!             11 - Obs & Assim      (namobs, nam_asminc)
!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!!======================================================================
!!                   ***  Run management namelists  ***
!!======================================================================
!!   namrun       parameters of the run
!!======================================================================
!
!-----------------------------------------------------------------------
&namrun        !   parameters of the run
!-----------------------------------------------------------------------
   nn_no       =       0   !  job number (no more used...)
   cn_exp      =  'eNATL60-BLBT02Y'
   nn_it000    =    1393201 !  first time step
   nn_itend    =    1447200 !  last  time step (std 5475)
   nn_date0    =    20101008 !  date at nit_0000 (format yyyymmdd) used if ln_rstart=F or (ln_rstart=T and nn_rstctl=0 or 1)
   nn_leapy    =       1    !  Leap year calendar (1) or not (0)
   ln_rstart   = .true.  !  start from rest (F) or from a restart file (T)
   nn_euler    =       1    !  = 0 : start with forward time step if ln_rstart=T  #lulu
   nn_rstctl   =  2   !  restart control => activated only if ln_rstart = T
                            !    = 0 nn_date0 read in namelist ; nn_it000 : read in namelist
                            !    = 1 nn_date0 read in namelist ; nn_it000 : check consistancy between namelist and restart
                            !    = 2 nn_date0 read in restart  ; nn_it000 : check consistancy between namelist and restart
   cn_ocerst_in    = 'eNATL60-BLBT02Y_01393200_restart_oce'
   cn_ocerst_indir = '/gpfs/scratch/pr1egh00/pr1egh01/eNATL60/eNATL60-BLBT02Y-R/01393200'   !  directory from which to read input ocean restarts
   cn_ocerst_out = 'restart_oce'
   cn_ocerst_outdir = '/gpfs/scratch/pr1egh00/pr1egh01/eNATL60/eNATL60-BLBT02Y-R/VRAC_OCE'      !  directory in which to write output ocean restarts
   nn_istate   =    0    !  output the initial state (1) or not (0)
!!
!! 1D is 2160 time step at dt=40s !
   nn_stock    =   10800   !  5d at dt=40s  frequency of creation of a restart file (modulo referenced to 1)
   nn_write    =   54000   ! 25d at dt=40s  frequency of write in the output file   (modulo referenced to nn_it000)
!!   nn_stock    =   12960   !  6d at dt=40s  frequency of creation of a restart file (modulo referenced to 1)
!!   nn_write    =   12960   !  6d at dt=40s  frequency of write in the output file   (modulo referenced to nn_it000)
!!   
!!   nn_write    = 43200   ! 20d at dt=40s  frequency of write in the output file   (modulo referenced to nn_it000)
!!   nn_write    =   32400   ! 15d at dt=40s  frequency of write in the output file   (modulo referenced to nn_it000)
!!
!! nn_write    =   64800   ! 30d at dt=40s  frequency of write in the output file   (modulo referenced to nn_it000)
!!   nn_write    =   10800   ! 5d at dt=40s  frequency of write in the output file   (modulo referenced to nn_it000)
   ln_dimgnnn  = .false.   !  DIMG file format: 1 file for all processors (F) or by processor (T)
   ln_mskland  = .false.   !  mask land points in NetCDF outputs (costly: + ~15%)
   ln_clobber  = .true.    ! #LOLO!  clobber (overwrite) an existing file
   nn_chunksz  =       0   !  chunksize (bytes) for NetCDF file (works only with iom_nf90 routines)
/
!
!!======================================================================
!!                      ***  Domain namelists  ***
!!======================================================================
!!   namcfg       parameters of the configuration
!!   namzgr       vertical coordinate
!!   namzgr_sco   s-coordinate or hybrid z-s-coordinate
!!   namdom       space and time domain (bathymetry, mesh, timestep)
!!   namtsd       data: temperature & salinity
!!======================================================================
!
!-----------------------------------------------------------------------
&namcfg     !   parameters of the configuration
!-----------------------------------------------------------------------
   cp_cfg      =  'enatl'              !  name of the configuration
   jp_cfg      =       60              !  resolution of the configuration
   jpidta      =     8354              !  1st lateral dimension ( >= jpi )
   jpjdta      =     4729              !  2nd    "         "    ( >= jpj )
   jpkdta      =     300               !  number of levels      ( >= jpk )
   jpiglo      =     8354              !  1st dimension of global domain --> i =jpidta
   jpjglo      =     4729              !  2nd    -                  -    --> j =jpjdta
   jpizoom     =       1               !  left bottom (i,j) indices of the zoom
   jpjzoom     =       1               !  in data domain indices
   jperio      =       0               !  lateral cond. type (between 0 and 6)
/
!-----------------------------------------------------------------------
&namzgr        !   vertical coordinate
!-----------------------------------------------------------------------
   ln_zps      = .true.    !  z-coordinate - partial steps   (T/F)
/
!-----------------------------------------------------------------------
&namzgr_sco    !   s-coordinate or hybrid z-s-coordinate
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&namdom        !   space and time domain (bathymetry, mesh, timestep)
!-----------------------------------------------------------------------
   nn_closea   =    0      !  remove (=0) or keep (=1) closed seas and lakes (ORCA)
   nn_msh      =    0      !  create (=1) a mesh file or not (=0) #lulu
   rn_hmin     =   -16.    ! #LOLO!  min depth of the ocean (>0) or min number of ocean level (<0) #lulu
   rn_e3zps_min=   25.     !  partial step thickness is set larger than the minimum of
   rn_e3zps_rat=    0.2    !  rn_e3zps_min and rn_e3zps_rat*e3t, with 0<rn_e3zps_rat<1
                           !
   rn_rdt      =   40.     !  time step for the dynamics (and tracer if nn_acc=0) #lulu
   rn_atfp     =    0.1    ! #LOLO! asselin time filter parameter
   ppglam0     =  999999.d0            !  longitude of first raw and column T-point (jphgr_msh = 1)
   ppgphi0     =  999999.d0            ! latitude  of first raw and column T-point (jphgr_msh = 1)
   ppe1_deg    =  999999.d0            !  zonal      grid-spacing (degrees)
   ppe2_deg    =  999999.d0            !  meridional grid-spacing (degrees)
   ppe1_m      =  999999.d0            !  zonal      grid-spacing (degrees)
   ppe2_m      =  999999.d0            !  meridional grid-spacing (degrees)
   ppsur       =  -12539.3187051d0     !  ORCA r4, r2 and r05 coefficients
   ppa0        =    48.7693273601d0    ! (default coefficients)
   ppa1        =    42.8722426445d0    !
   ppkth       =    317.238187329d0    !
   ppacr       =    38.4998412469d0    !
   ppdzmin     =        999999.d0      !  Minimum vertical spacing
   pphmax      =        999999.d0      !  Maximum depth
   ldbletanh   =           .true.      !  Use/do not use double tanf function for vertical coordinates
   ppa2        =  20.1768193931d0      !  Double tanh function parameters
   ppkth2      =  31.5541059316d0      !
   ppacr2      =  121.356963487d0      !
/
!-----------------------------------------------------------------------
&namsplit      !   time splitting parameters                            ("key_dynspg_ts")
!-----------------------------------------------------------------------
   ln_bt_fw      =    .false.          ! #LOLO? <needs u2b et v2b in restarts...>  Forward integration of barotropic equations
   ln_bt_av      =    .true.           !  Time filtering of barotropic variables
   ln_bt_nn_auto =    .true.           !  Set nn_baro automatically to be just below
                                       !  a user defined maximum courant number (rn_bt_cmax)
   nn_baro       =    30               ! #LOLO!  Number of iterations of barotropic mode
                                       !  during rn_rdt seconds. Only used if ln_bt_nn_auto=F
   rn_bt_cmax    =    0.65             ! #LOLO! Maximum courant number allowed if ln_bt_nn_auto=T  #lulu
   nn_bt_flt     =    2                ! #LOLO! Time filter choice ( => Jouano p.c. ?)
                                       !  = 0 None
                                       !  = 1 Boxcar over   nn_baro barotropic steps
                                       !  = 2 Boxcar over 2*nn_baro     "        "
/
!-----------------------------------------------------------------------
&namcrs        !   Grid coarsening for dynamics output and/or
               !   passive tracer coarsened online simulations
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&namtsd    !   data : Temperature  & Salinity
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!          !  file name                            ! frequency (hours) ! variable  ! time interp. !  clim  ! 'yearly'/ ! weights  ! rotation ! land/sea mask !
!          !                                       !  (if <0  months)  !   name    !   (logical)  !  (T/F) ! 'monthly' ! filename ! pairing  ! filename      !
   sn_tem  = 'votemper_GLORYS12V1-eNATL60_2008.nc', -12     ,'votemper',   .true.    , .true. , 'yearly'  ,    ''    ,    ''    ,    ''
   sn_sal  = 'vosaline_GLORYS12V1-eNATL60_2008.nc', -12     ,'vosaline',   .true.    , .true. ,'yearly'  ,    ''    ,    ''    ,    ''
   !
   ln_tsd_init   = .false.  ! Initialisation of ocean T & S with T &S input data (T) or not (F)
   ln_tsd_tradmp = .false.  !  damping of ocean T & S toward T & S input data (T) or not (F)
/
!!======================================================================
!!            ***  Surface Boundary Condition namelists  ***
!!======================================================================
!!   namsbc          surface boundary condition
!!   namsbc_ana      analytical         formulation
!!   namsbc_flx      flux               formulation
!!   namsbc_clio     CLIO bulk formulae formulation
!!   namsbc_core     CORE bulk formulae formulation
!!   namsbc_mfs      MFS  bulk formulae formulation
!!   namsbc_cpl      CouPLed            formulation                     ("key_oasis3")
!!   namsbc_sas      StAndalone Surface module
!!   namtra_qsr      penetrative solar radiation
!!   namsbc_rnf      river runoffs
!!   namsbc_isf      ice shelf melting/freezing
!!   namsbc_apr      Atmospheric Pressure
!!   namsbc_ssr      sea surface restoring term (for T and/or S)
!!   namsbc_alb      albedo parameters
!!======================================================================
!
!-----------------------------------------------------------------------
&namsbc        !   Surface Boundary Condition (surface module)
!-----------------------------------------------------------------------
   nn_fsbc     = 6         ! #LOLO! (4 min @ dt=40s)  frequency of surface boundary condition computation #lulu
   ln_blk_core = .true.    !  CORE bulk formulation                     (T => fill namsbc_core)
   ln_apr_dyn  = .true.    ! #LOLO?  Patm gradient added in ocean & ice Eqs.   (T => fill namsbc_apr ) #lulu
   nn_ice_embd = 1         !  =0 levitating ice (no mass exchange, concentration/dilution effect)
                           !  =1 levitating ice with mass and salt exchange but no presure effect
                           !  =2 embedded sea-ice (full salt and mass exchanges and pressure)
   ln_dm2dc    = .true.    !  daily mean to diurnal cycle on short wave
   ln_rnf      = .true.    !  runoffs                                   (T   => fill namsbc_rnf)
   nn_isf      = 0         !  ice shelf melting/freezing                (/=0 => fill namsbc_isf)
                           !  0 =no isf                  1 = presence of ISF
                           !  2 = bg03 parametrisation   3 = rnf file for isf
                           !  4 = ISF fwf specified
                           !  option 1 and 4 need ln_isfcav = .true. (domzgr)
   ln_ssr      = .true.    !  Sea Surface Restoring on T and/or S       (T => fill namsbc_ssr)
   nn_fwb      = 0         ! #LOLO! FreshWater Budget: =0 unchecked
                           !     =1 global mean of e-p-r set to zero at each time step
                           !     =2 annual global mean of e-p-r set to zero
/
!-----------------------------------------------------------------------
&namsbc_core   !   namsbc_core  CORE bulk formulae
!-----------------------------------------------------------------------
!              !  file name                    ! frequency (hours) ! variable  ! time interp. !  clim  ! 'yearly'/ ! weights                               ! rotation ! land/sea mask !
!              !                               !  (if <0  months)  !   name    !   (logical)  !  (T/F) ! 'monthly' ! filename                              ! pairing  ! filename      !
   sn_wndi     = 'drowned_u10_DFS5.2'          ,         3         , 'u10'     ,   .false.  , .false., 'yearly'  , 'weight_bicubic_512x256-eNATL60.nc' , 'Uwnd'   , ''
   sn_wndj     = 'drowned_v10_DFS5.2'          ,         3         , 'v10'     ,   .false.  , .false., 'yearly'  , 'weight_bicubic_512x256-eNATL60.nc' , 'Vwnd'   , ''
   sn_qsr      = 'drowned_radsw_DFS5.2'        ,        24         , 'radsw'   ,   .false.  , .false., 'yearly'  , 'weight_bilinear_512x256-eNATL60.nc' , ''       , ''
   sn_qlw      = 'drowned_radlw_DFS5.2'        ,        24         , 'radlw'   ,   .false.  , .false., 'yearly'  , 'weight_bilinear_512x256-eNATL60.nc' , ''       , ''
   sn_tair     = 'drowned_t2_DFS5.2'           ,         3         , 't2'      ,   .false.  , .false., 'yearly'  , 'weight_bilinear_512x256-eNATL60.nc' , ''       , ''
   sn_humi     = 'drowned_q2_DFS5.2'           ,         3         , 'q2'      ,   .false.  , .false., 'yearly'  , 'weight_bilinear_512x256-eNATL60.nc' , ''       , ''
   sn_prec     = 'drowned_precip_DFS5.2'       ,        24         , 'precip'  ,   .false.  , .false., 'yearly'  , 'weight_bilinear_512x256-eNATL60.nc' , ''       , ''
   sn_snow     = 'drowned_snow_DFS5.2'         ,        24         , 'snow'    ,   .false.  , .false., 'yearly'  , 'weight_bilinear_512x256-eNATL60.nc' , ''       , ''
   cn_dir      = './FATM/'  !  root directory for the location of the bulk files
   rn_zqt      =  2.        !  Air temperature and humidity reference height (m)
   rn_zu       = 10.        !  Wind vector reference height (m)
   !!
   rn_vfac     = 0.5       !#LOLO Julien Jouano! multiplicative factor for ocean/ice velocity
                           !  in the calculation of the wind stress (0.=absolute winds or 1.=relative winds)
/
!-----------------------------------------------------------------------
&namtra_qsr    !   penetrative solar radiation
!-----------------------------------------------------------------------
!              !  file name  ! frequency (hours) ! variable  ! time interp. !  clim  ! 'yearly'/ ! weights  ! rotation ! land/sea mask !
!              !             !  (if <0  months)  !   name    !   (logical)  !  (T/F) ! 'monthly' ! filename ! pairing  ! filename      !
   sn_chl      ='CHLA_ESACCI-OC-L3S-4km-eNATL60_y2010', -1    , 'CHLA'    ,   .true.     , .true. , 'yearly'  , ''       , ''       , ''
   ln_traqsr   = .true.    !  Light penetration (T) or not (F)
   ln_qsr_rgb  = .true.    !  RGB (Red-Green-Blue) light penetration
   ln_qsr_2bd  = .false.   !  2 bands              light penetration
   nn_chldta   =      0    ! #LOLO!  RGB : 2D Chl data (=1), 3D Chl data (=2) or cst value (=0)
/
!-----------------------------------------------------------------------
&namsbc_rnf    !   runoffs namelist surface boundary condition
!-----------------------------------------------------------------------
!              !  file name           ! frequency (hours) ! variable  ! time interp. !  clim  ! 'yearly'/ ! weights  ! rotation ! land/sea mask !
!              !                      !  (if <0  months)  !   name    !   (logical)  !  (T/F) ! 'monthly' ! filename ! pairing  ! filename      !
   sn_rnf      = 'runoff_3.1.3.nc'    ,-1                 ,'sorunoff' ,.true.        , .true. , 'yearly'  , ''       , ''       , ''
   sn_cnf      = 'runoff_3.1.3.nc'    , 0                 ,'socoefr'  , .false.      , .true. , 'yearly'  , ''       , ''       , ''
   ln_rnf_mouth = .true.     ! #LOLO!  specific treatment at rivers mouths
   rn_hrnf      =  15.e0     !  depth over which enhanced vertical mixing is used
   rn_avt_rnf   =   2.e-3    ! #LOLO! value of the additional vertical mixing coef. [m2/s]
   ln_rnf_depth = .false.    !  read in depth information for runoff
   ln_rnf_tem   = .false.    !  read in temperature information for runoff
   ln_rnf_sal   = .false.    !  read in salinity information for runoff
   ln_rnf_depth_ini = .true. !  compute depth at initialisation from runoff file
   rn_rnf_max   = 0.05       !  max value of the runoff climatology over global domain ( if ln_rnf_depth_ini = .true )
   rn_dep_max = 150.         !  depth over which runoffs is spread ( ln_rnf_depth_ini = .true )
   nn_rnf_depth_file = 0     ! create (=1) a runoff depth file or not (=0)
/
!-----------------------------------------------------------------------
&namsbc_apr    !   Atmospheric pressure used as ocean forcing or in bulk
!-----------------------------------------------------------------------
!              !  file name  ! frequency (hours) ! variable  ! time interp. !  clim  ! 'yearly'/ ! weights  ! rotation ! land/sea mask !
!              !             !  (if <0  months)  !   name    !   (logical)  !  (T/F) ! 'monthly' ! filename ! pairing  ! filename      !
   sn_apr      = 'drowned_msl_ERAinterim', 3     ,   'MSL'   ,    .true.    ,.false. , 'yearly'  ,'weight_bilinear_512x256-eNATL60.nc','',''
   cn_dir      = './FATM/'  !  root directory for the location of the bulk files
   rn_pref     = 101000.    !  reference atmospheric pressure   [N/m2]/
   ln_ref_apr  = .false.    !  ref. pressure: global mean Patm (T) or a constant (F)
   ln_apr_obc  = .false.    !  inverse barometer added to OBC ssh data
/
!-----------------------------------------------------------------------
&namsbc_ssr    !   surface boundary condition : sea surface restoring
!-----------------------------------------------------------------------
!              !  file name  ! frequency (hours) ! variable  ! time interp. !  clim  ! 'yearly'/ ! weights  ! rotation ! land/sea mask !
!              !             !  (if <0  months)  !   name    !   (logical)  !  (T/F) ! 'monthly' ! filename ! pairing  ! filename      !
   sn_sss      = 'sss_WOA2013-1440x720-eNATL60_mean_2005-2012.nc',-1,  'sss'   ,    .true.   , .true. , 'yearly'  ,    ''    ,    ''    ,     ''
   nn_sssr     =     2     !  add a damping     term in the surface freshwater flux (=2)
   rn_deds     =  -166.67  ! #LOLO?  magnitude of the damping on salinity   [mm/day]
   rn_sssr_bnd =   4.e0    ! #LOLO?  ABS(Max/Min) value of the damping erp term [mm/day]
/
!-----------------------------------------------------------------------
&namsbc_alb    !   albedo parameters
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&namberg       !   iceberg parameters
!-----------------------------------------------------------------------
      ln_icebergs              = .false.
/
!-----------------------------------------------------------------------
&namlbc        !   lateral momentum boundary condition
!-----------------------------------------------------------------------
   rn_shlat    =    2.     ! #LOLO? shlat = 0  !  0 < shlat < 2  !  shlat = 2  !  2 < shlat #lulu
                           !       free slip  !   partial slip  !   no slip   ! strong slip
   ln_vorlat   = .false.   !  consistency of vorticity boundary condition with analytical eqs.
/
!-----------------------------------------------------------------------
&namcla        !   cross land advection
!-----------------------------------------------------------------------
   nn_cla      =    0      !  advection between 2 ocean pts separates by land
/
!-----------------------------------------------------------------------
&nam_tide      !   tide parameters (#ifdef key_tide)
!-----------------------------------------------------------------------
   clname(1)  =  'M2'
   clname(2)  =  'S2'
   clname(3)  =  'N2'
   clname(4)  =  'K1'
   clname(5)  =  'O1'
/
!-----------------------------------------------------------------------
&nambdy        !  unstructured open boundaries                          ("key_bdy")
!-----------------------------------------------------------------------
   nb_bdy         = 2        !  number of open boundary sets
   ln_coords_file = .false.,.false.      !  =T : read bdy coordinates from file
   cn_dyn2d       = 'flather','flather'  ! #LOLO!
   nn_dyn2d_dta   = 3,3                  !  = 0, bdy data are equal to the initial state
      !                       !  = 1, bdy data are read in 'bdydata   .nc' files
      !                       !  = 2, use tidal harmonic forcing data from files
      !                       !  = 3, use external data AND tidal harmonic forcing
   cn_dyn3d      =  'frs','frs'        !
   nn_dyn3d_dta  =  1,1                !  = 0, bdy data are equal to the initial state
   !                                   !  = 1, bdy data are read in 'bdydata   .nc' files
   cn_tra        =  'frs','frs'        !
   nn_tra_dta    =  1,1                !  = 0, bdy data are equal to the initial state
   !                                   !  = 1, bdy data are read in 'bdydata   .nc' files
   cn_ice_lim    =  'none','frs'       ! #LOLO!
   nn_ice_lim_dta=  0,1                !  = 0, bdy data are equal to the initial state
   !                                   !  = 1, bdy data are read in 'bdydata   .nc' files
   rn_ice_tem    = 270.,270.           !  lim3 only: arbitrary temperature of incoming sea ice
   rn_ice_sal    = 10. ,10.            !  lim3 only:      --   salinity           --
   rn_ice_age    = 30. ,30.            !  lim3 only:      --   age                --
   !
   ln_tra_dmp    =.false.,.false.       !  open boudaries conditions for tracers
   ln_dyn3d_dmp  =.false.,.false.        !  open boundary condition for baroclinic velocities
   rn_time_dmp   =  1.,1.     ! Damping time scale in days
   rn_time_dmp_out =  1.,1.   ! Outflow damping time scale
   nn_rimwidth   = 10,10        !  width of the relaxation zone
   ln_vol        = .true.     !  total volume correction (see nn_volctl parameter)
   nn_volctl     =  1         !  = 0, the total water flux across open boundaries is zero
/
!-----------------------------------------------------------------------
&nambdy_index  !  structured open boundaries definition     
!-----------------------------------------------------------------------
    ctypebdy ='S'                   ! Open boundary type (W,E,S or N)
    nbdyind  = 2                    ! indice of velocity row or column
                                    ! if ==-1, set obc at the domain boundary
                                    !        , discard start and end indices
    nbdybeg  = 2412                 ! indice of segment start
    nbdyend  = 5226                 ! indice of segment end
/
!-----------------------------------------------------------------------
&nambdy_index  !  structured open boundaries definition
!-----------------------------------------------------------------------
    ctypebdy ='N'                   ! Open boundary type (W,E,S or N)
    nbdyind  = 4727                 ! indice of velocity row or column
                                    ! if ==-1, set obc at the domain boundary
                                    !        , discard start and end indices
    nbdybeg  = 1860                 ! indice of segment start
    nbdyend  = 6026                 ! indice of segment end
/
!-----------------------------------------------------------------------
&nambdy_dta      !  open boundaries SOUTH - external data
!-----------------------------------------------------------------------
!              !  file name      ! frequency (hours) ! variable   ! time interp.   !  clim   ! 'yearly'/ ! weights  ! rotation ! land/sea mask !
!              !                 !  (if <0  months)  !   name     !   (logical)    !  (T/F ) ! 'monthly' ! filename ! pairing  ! filename      !
   bn_ssh      = 'sossheig_GLORYS12V1-BDY_t_S_eNATL60_3.6',24, 'sossheig',    .true.   , .false.,  'yearly'  ,    ''            ,   ''     ,     ''
   bn_u2d      = 'boo.nc'        ,         24        , 'vobtcrtx',    .true.   , .false.,  'daily'  ,    ''            ,   ''     ,     ''
   bn_v2d      = 'boo.nc'        ,         24        , 'vobtcrty',    .true.   , .false.,  'daily'  ,    ''            ,   ''     ,     ''
   bn_u3d      = 'vozocrtx_GLORYS12V1-BDY_u_S_eNATL60_3.6',24, 'vozocrtx',    .true.   , .false.,  'yearly'  ,    ''            ,   ''     ,     ''
   bn_v3d      = 'vomecrty_GLORYS12V1-BDY_v_S_eNATL60_3.6',24, 'vomecrty',    .true.   , .false.,  'yearly'  ,    ''            ,   ''     ,     ''
   bn_tem      = 'votemper_GLORYS12V1-BDY_t_S_eNATL60_3.6',24, 'votemper',    .true.   , .false.,  'yearly'  ,    ''            ,   ''     ,     ''
   bn_sal      = 'vosaline_GLORYS12V1-BDY_t_S_eNATL60_3.6',24, 'vosaline',    .true.   , .false.,  'yearly'  ,    ''            ,   ''     ,     ''
   cn_dir  =    './BDY/'
   ln_full_vel = .true.
/
!
!!-----------------------------------------------------------------------
&nambdy_dta    !  open boundaries NORTH - external data                       
!-----------------------------------------------------------------------
!              !  file name      ! frequency (hours) ! variable   ! time interp.   !  clim   ! 'yearly'/ ! weights  ! rotation ! land/sea mask !
!              !                 !  (if <0  months)  !   name     !   (logical)    !  (T/F ) ! 'monthly' ! filename ! pairing  ! filename      !
   bn_ssh      = 'sossheig_GLORYS12V1-BDY_t_N_eNATL60_3.6',24, 'sossheig',    .true.   , .false.,  'yearly'  ,    ''            ,   ''     ,     ''
   bn_u2d      = 'boo.nc'        ,         24        , 'vobtcrtx',    .true.   , .false.,  'daily'  ,    ''            ,   ''     ,     ''
   bn_v2d      = 'boo.nc'        ,         24        , 'vobtcrty',    .true.   , .false.,  'daily'  ,    ''            ,   ''     ,     ''
   bn_u3d      = 'vozocrtx_GLORYS12V1-BDY_u_N_eNATL60_3.6',24, 'vozocrtx',    .true.   , .false.,  'yearly'  ,    ''            ,   ''     ,     ''
   bn_v3d      = 'vomecrty_GLORYS12V1-BDY_v_N_eNATL60_3.6',24, 'vomecrty',    .true.   , .false.,  'yearly'  ,    ''            ,   ''     ,     ''
   bn_tem      = 'votemper_GLORYS12V1-BDY_t_N_eNATL60_3.6',24, 'votemper',    .true.   , .false.,  'yearly'  ,    ''            ,   ''     ,     ''
   bn_sal      = 'vosaline_GLORYS12V1-BDY_t_N_eNATL60_3.6',24, 'vosaline',    .true.   , .false.,  'yearly'  ,    ''            ,   ''     ,     ''
!
   bn_a_i      = 'ileadfra_GLORYS2V4-BDY_t_N_eNATL60_3.6',-1, 'ileadfra',    .true.   , .false. , 'yearly'  ,    ''    ,   ''     , '' ! #LOLO: MONTHLY of GLORYS2V4 !!! (not daily GLORYS12 as above)
   bn_ht_i     = 'iicethic_GLORYS2V4-BDY_t_N_eNATL60_3.6',-1, 'iicethic',    .true.   , .false. , 'yearly'  ,    ''    ,   ''     , '' ! #LOLO: MONTHLY of GLORYS2V4 !!! (not daily GLORYS12 as above)
   bn_ht_s     = 'isnowthi_GLORYS2V4-BDY_t_N_eNATL60_3.6',-1, 'isnowthi',    .true.   , .false. , 'yearly'  ,    ''    ,   ''     , '' ! #LOLO: MONTHLY of GLORYS2V4 !!! (not daily GLORYS12 as above)
! 
   cn_dir  =    './BDY/'
   ln_full_vel = .true.
/
!
!-----------------------------------------------------------------------
&nambdy_tide     ! tidal forcing at open boundaries - SOUTH
!-----------------------------------------------------------------------
   filtide      = './BDY/bdytide_FES2014b-eNATL60_S_'         !  file name root of tidal forcing files
   ln_bdytide_2ddta = .false.
   ln_bdytide_conj  = .false.
/
!-----------------------------------------------------------------------
&nambdy_tide     ! tidal forcing at open boundaries - NORTH
!-----------------------------------------------------------------------
   filtide      = './BDY/bdytide_FES2014b-eNATL60_N_'         !  file name root of tidal forcing files
   ln_bdytide_2ddta = .false.
   ln_bdytide_conj  = .false.
/
!
!
!
!!======================================================================
!!                 ***  Bottom boundary condition  ***
!!======================================================================
!!   nambfr        bottom friction
!!   nambbc        bottom temperature boundary condition
!!   nambbl        bottom boundary layer scheme                         ("key_trabbl")
!!======================================================================
!
!-----------------------------------------------------------------------
&nambfr        !   bottom friction
!-----------------------------------------------------------------------
   nn_bfr      =    2      !  type of bottom friction :   = 0 : free slip,  = 1 : linear friction
                           !                              = 2 : nonlinear friction
   rn_bfri2    =    2.5e-3 ! #LOLO! bottom drag coefficient (non linear case). Minimum coeft if ln_loglayer=T !! P.Mattiot => 3.e-3
   rn_bfri2_max =   1.e-1  ! max. bottom drag coefficient (non linear case and ln_loglayer=T)
   rn_bfeb2    =    0.0    ! #LOLO!  bottom turbulent kinetic energy background  (m2/s2)
   rn_bfrz0    =    3.e-3  !  bottom roughness [m] if ln_loglayer=T
   ln_bfr2d    = .true.    ! #LOLO! horizontal variation of the bottom friction coef (read a 2D mask file )
   rn_bfrien   =    50.    ! #LOLO! local multiplying factor of bfr (ln_bfr2d=T)
   ln_bfrimp   = .true.    !  implicit bottom friction (requires ln_zdfexp = .false. if true)
   ln_loglayer = .true.    ! #LOLO!  logarithmic formulation (non linear case)
/
!-----------------------------------------------------------------------
&nambbc        !   bottom temperature boundary condition
!-----------------------------------------------------------------------
!              !                              !  (if <0  months)  !  
!              !  file name      ! frequency (hours) ! variable   ! time interp.   !  clim   ! 'yearly'/ ! weights  ! rotation ! land/sea mask !
!              !                 !  (if <0  months)  !   name     !   (logical)    !  (T/F ) ! 'monthly' ! filename ! pairing  ! filename      !
   sn_qgh      ='geothermal_flux_Goutorbe',  -12.  , 'gh_flux'    ,   .false.     , .true. , 'yearly'  ,   ''     ,   ''     ,   ''
   !
   cn_dir      = './'      !  root directory for the location of the runoff files
   ln_trabbc   = .false.   ! #LOLO? Apply a geothermal heating at the ocean bottom
   nn_geoflx   =    0      !    geothermal heat flux: = 0 no flux
                           !     = 1 constant flux
                           !     = 2 variable flux (read in geothermal_heating.nc in mW/m2)
   rn_geoflx_cst = 86.4e-3 !  Constant value of geothermal heat flux [W/m2]
/
!-----------------------------------------------------------------------
&nambbl        !   bottom boundary layer scheme
!-----------------------------------------------------------------------
   nn_bbl_ldf  =  1      ! #LOLO?  diffusive bbl (=1)   or not (=0)
   nn_bbl_adv  =  1      !  advective bbl (=1/2) or not (=0)
   rn_ahtbbl   =  1000.  ! #LOLO?  lateral mixing coefficient in the bbl  [m2/s]
   rn_gambbl   =  10.    !  advective bbl coefficient                 [s]
/

!!======================================================================
!!                        Tracer (T & S ) namelists
!!======================================================================
!!   nameos        equation of state
!!   namtra_adv    advection scheme
!!   namtra_adv_mle   mixed layer eddy param. (Fox-Kemper param.)
!!   namtra_ldf    lateral diffusion scheme
!!   namtra_dmp    T & S newtonian damping
!!======================================================================
!
!-----------------------------------------------------------------------
&nameos        !   ocean physical parameters
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&namtra_adv    !   advection scheme for tracer
!-----------------------------------------------------------------------
!!   ln_traadv_tvd    =  .false.  !  TVD scheme
!!   ln_traadv_ubs    =  .true.   ! #LOLO!  UBS scheme
   ln_traadv_tvd    =  .true.  !  TVD scheme
   ln_traadv_ubs    =  .false.   ! #LOLO!  UBS scheme
!! LOLO: try if instabilities!!!   ln_traadv_tvd_zts=  .false.  !  TVD scheme with sub-timestepping of vertical tracer advection
/
!-----------------------------------------------------------------------
&namtra_adv_mle !   mixed layer eddy parametrisation (Fox-Kemper param)
!-----------------------------------------------------------------------
   ln_mle    = .true.      ! #LOLO! (T) use the Mixed Layer Eddy (MLE) parameterisation
/
!----------------------------------------------------------------------------------
&namtra_ldf    !   lateral diffusion scheme for tracers
!----------------------------------------------------------------------------------
   !                       !  Operator type:
   ln_traldf_lap    =  .true.  ! #LOLO!  laplacian operator
   ln_traldf_bilap  =  .false. ! #LOLO!  bilaplacian operator    !LULU
   ln_traldf_iso    =  .true.  ! #LOLO! iso-neutral                 (needs "key_ldfslp") !LULU
   ln_traldf_hor    =  .false. !  horizontal (geopotential)   (needs "key_ldfslp" when ln_sco=T) !LULU
   rn_aeiv_0        =    0.    !  eddy induced velocity coefficient [m2/s]
   rn_aht_0         =   20.    ! [m^2/s] Tweaked traldf_lap.F90 !!!! LOLO!!!
!!
!! LOLO add:
!   rthr_grad_sst    = 0.0007  ! [K/m] | threshold value from which |grad(SST)| is considered too extreme!!
!   ramp_aht         =  10.0   ! multiplicative factor to ahtu,ahtv to apply where |grad(SST)| == rthr_grad_sst !
!   !!                         !    => can become larger where  |grad(SST)| > rthr_grad_sst !
!   rmax_msk         = 15.0    !    => but never larger than rmax_msk
!   rmin_msk         = 0.01    ! minimum value of the mask (background for the rest of the domain non problematic)
!   !!
!   !!nt_often         =  45     !  !!! 45 @40s => 30 minutes!
!   nt_often         =  15     !  !!! 15 @40s => 10 minutes!
!   !!nb_smooth_sst    = 8
!   nb_smooth_sst    = 2
!   nb_smooth_mask   = 10
!
/
!-----------------------------------------------------------------------
&namtra_dmp    !   tracer: T & S newtonian damping
!-----------------------------------------------------------------------
   ln_tradmp   =  .false.  !  add a damping termn (T) or not (F)
/

!!======================================================================
!!                      ***  Dynamics namelists  ***
!!======================================================================
!!   namdyn_adv    formulation of the momentum advection
!!   namdyn_vor    advection scheme
!!   namdyn_hpg    hydrostatic pressure gradient
!!   namdyn_spg    surface pressure gradient                            (CPP key only)
!!   namdyn_ldf    lateral diffusion scheme
!!======================================================================
!
!-----------------------------------------------------------------------
&namdyn_adv    !   formulation of the momentum advection
!-----------------------------------------------------------------------
   ln_dynadv_vec = .false.  ! #LOLO?  vector form (T) or flux form (F)
   nn_dynkeg     = 1        ! #LOLO? scheme for grad(KE): =0   C2  ;  =1   Hollingsworth correction
   ln_dynadv_cen2= .false.  !  flux form - 2nd order centered scheme
   ln_dynadv_ubs = .true.   ! #LOLO!  flux form - 3rd order UBS      scheme
   ln_dynzad_zts = .false.  ! <#LOLO? [requires vector form!]> Use (T) sub timestepping for vertical momentum advection
/
!-----------------------------------------------------------------------
&nam_vvl    !   vertical coordinate options
!-----------------------------------------------------------------------
   ln_vvl_zstar  = .true.           !  zstar vertical coordinate
/
!-----------------------------------------------------------------------
&namdyn_vor    !   option of physics/algorithm (not control by CPP keys)
!-----------------------------------------------------------------------
   ln_dynvor_ene = .false. !  energy conserving scheme GM! #LOLO!
   ln_dynvor_ens = .false. !  enstrophy conserving scheme GM!
   ln_dynvor_een = .true.     ! #LOLO!!!!  energy & enstrophy scheme  
   ln_dynvor_een_old = .false.  ! #LOLO?  energy & enstrophy scheme - original formulation
/
!-----------------------------------------------------------------------
&namdyn_hpg    !   Hydrostatic pressure gradient option
!-----------------------------------------------------------------------
   ln_hpg_zco  = .false.   !  z-coordinate - full steps  
   ln_hpg_zps  = .false.   !  z-coordinate - partial steps (interpolation)
   ln_hpg_sco  = .true.    ! #LOLO! <because VVL on> s-coordinate (standard jacobian formulation)
/
!-----------------------------------------------------------------------
!namdyn_spg    !   surface pressure gradient   (CPP key only)
!-----------------------------------------------------------------------
!                          !  explicit free surface                     ("key_dynspg_exp")
!                          !  filtered free surface                     ("key_dynspg_flt")
!                          !  split-explicit free surface               ("key_dynspg_ts")

!-----------------------------------------------------------------------
&namdyn_ldf    !   lateral diffusion on momentum
!-----------------------------------------------------------------------
   !                       !  Type of the operator :
   ln_dynldf_lap    = .false.   ! #LOLO! laplacian operator     #lulu
   ln_dynldf_bilap  = .false.   ! #LOLO!  bilaplacian operator  #lulu
   !                       !  Direction of action  :
   ln_dynldf_level  =  .false.  !  iso-level
   ln_dynldf_hor    =  .true.   !  horizontal (geopotential)            (require "key_ldfslp" in s-coord.)
   ln_dynldf_iso    =  .false.  !  iso-neutral                          (require "key_ldfslp")
   !                       !  Coefficient
   rn_ahm_0_lap     =     0.    !  horizontal laplacian eddy viscosity   [m2/s]
   rn_ahmb_0        =     0.    !  background eddy viscosity for ldf_iso [m2/s]
   rn_ahm_0_blp     = -0.25e8   !  horizontal bilaplacian eddy viscosity [m4/s]
   rn_ahm_m_blp     = -2.00e8   !  upper limit for bilap  abs(ahm) < min( dx^4/128rdt, rn_ahm_m_blp)
/

!!======================================================================
!!             Tracers & Dynamics vertical physics namelists
!!======================================================================
!!    namzdf            vertical physics
!!    namzdf_ric        richardson number dependent vertical mixing     ("key_zdfric")
!!    namzdf_tke        TKE dependent vertical mixing                   ("key_zdftke")
!!    namzdf_kpp        KPP dependent vertical mixing                   ("key_zdfkpp")
!!    namzdf_ddm        double diffusive mixing parameterization        ("key_zdfddm")
!!    namzdf_tmx        tidal mixing parameterization                   ("key_zdftmx")
!!    namzdf_tmx_new    new tidal mixing parameterization               ("key_zdftmx_new")
!!======================================================================
!
!-----------------------------------------------------------------------
&namzdf        !   vertical physics
!-----------------------------------------------------------------------
   rn_avm0     =   1.e-4   !  vertical eddy viscosity   [m2/s]          (background Kz if not "key_zdfcst")
   rn_avt0     =   1.e-5   !  vertical eddy diffusivity [m2/s]          (background Kz if not "key_zdfcst")
   nn_avb      =    0      !  profile for background avt & avm (=1) or not (=0)
   nn_havtb    =    1      !  horizontal shape for avtb (=1) or not (=0)
   ln_zdfevd   = .true.    !LOLO!  enhanced vertical diffusion (evd) (T) or not (F)
   nn_evdm     =    0      ! #LOLO?  evd apply on tracer (=0) or on tracer and momentum (=1)
   rn_avevd    =  10.      ! #LOLO?   evd mixing coefficient [m2/s]
   ln_zdfnpc   = .false.   !LOLO!!!  Non-Penetrative Convective algorithm (T) or not (F)
   nn_npc      =    1      !  frequency of application of npc
   nn_npcp     =  365      !  npc control print frequency
   ln_zdfexp   = .false.   !  time-stepping: split-explicit (T) or implicit (F) time stepping
   nn_zdfexp   =    3      !  number of sub-timestep for ln_zdfexp=T
/
!-----------------------------------------------------------------------
&namzdf_tke    !   turbulent eddy kinetic dependent vertical diffusion  ("key_zdftke")
!-----------------------------------------------------------------------
   rn_ediff    =   0.1     !  coef. for vertical eddy coef. (avt=rn_ediff*mxl*sqrt(e) )
   rn_ediss    =   0.7     !  coef. of the Kolmogoroff dissipation
   rn_ebb      =   67.83   !  coef. of the surface input of tke (=67.83 suggested when ln_mxl0=T)
   nn_mxl      =   3       !#LOLO: same as CJM165  mixing length: = 0 bounded by the distance to surface and bottom
                           !                 = 1 bounded by the local vertical scale factor
                           !                 = 2 first vertical derivative of mixing length bounded by 1
                           !                 = 3 as =2 with distinct disspipative an mixing length scale
   rn_mxl0     =   0.01    !#LOLO: same as CJM165  surface  buoyancy lenght scale minimum value (0.04)
   nn_etau     =   1       ! #LOLO? (will help get shalower LabSea MLDs) penetration of tke below the mixed layer (ML) due to internal & intertial waves
                           !        = 0 no penetration
                           !        = 1 add a tke source below the ML
                           !        = 2 add a tke source just at the base of the ML
                           !        = 3 as = 1 applied on HF part of the stress    ("key_oasis3")
/
!-----------------------------------------------------------------------
&namzdf_ddm    !   double diffusive mixing parameterization             ("key_zdfddm")
!-----------------------------------------------------------------------
!!! #LOLO: not used  (key_zdfddm not used!)
!   rn_avts     = 1.e-4     !  maximum avs (vertical mixing on salinity)
!   rn_hsbfr    = 1.6       !  heat/salt buoyancy flux ratio
/
!-----------------------------------------------------------------------
&namzdf_tmx    !   tidal mixing parameterization                        ("key_zdftmx")
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&namzdf_tmx_new    !   new tidal mixing parameterization                ("key_zdftmx_new")
!-----------------------------------------------------------------------
/
!!======================================================================
!!                  ***  Miscellaneous namelists  ***
!!======================================================================
!!   namsol            elliptic solver / island / free surface
!!   nammpp            Massively Parallel Processing                    ("key_mpp_mpi)
!!   namctl            Control prints & Benchmark
!!   namc1d            1D configuration options                         ("key_c1d")
!!   namc1d_uvd        data: U & V currents                             ("key_c1d")
!!   namc1d_dyndmp     U & V newtonian damping                          ("key_c1d")
!!   namsto            Stochastic parametrization of EOS
!!======================================================================
!
!-----------------------------------------------------------------------
&namsol        !   elliptic solver / island / free surface
!-----------------------------------------------------------------------
   nn_nmax     =   2800    ! #LOLO?  maximum of iterations for the SOR solver
   rn_sor      =  1.973    ! #LOLO?  optimal coefficient for SOR solver (to be adjusted with the domain)
/
!-----------------------------------------------------------------------
&nammpp        !   Massively Parallel Processing                        ("key_mpp_mpi)
!-----------------------------------------------------------------------
jpni        =   251  !  jpni   number of processors following i (set automatically if < 1)
jpnj        =   121  !  jpnj   number of processors following j (set automatically if < 1)
jpnij       =   16215 !  jpnij  number of local domains (set automatically if < 1)
/
!-----------------------------------------------------------------------
&namctl        !   Control prints & Benchmark
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&namptr       !   Poleward Transport Diagnostic
!-----------------------------------------------------------------------
   ln_diaptr  = .false.    !  Poleward heat and salt transport (T) or not (F)
   ln_subbas  = .false.     !  Atlantic/Pacific/Indian basins computation (T) or not
/
!-----------------------------------------------------------------------
&namhsb       !  Heat and salt budgets
!-----------------------------------------------------------------------
   ln_diahsb  = .false.    !  check the heat and salt budgets (T) or not (F)
/
!-----------------------------------------------------------------------
&namdyn_nept  !   Neptune effect (simplified: lateral and vertical diffusions removed)
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&namzdf_gls                !   GLS vertical diffusion                   ("key_zdfgls")
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&namtrd        !   diagnostics on dynamics and/or tracer trends
!              !       and/or mixed-layer trends and/or barotropic vorticity
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&namsto       ! Stochastic parametrization of EOS
!-----------------------------------------------------------------------
/
