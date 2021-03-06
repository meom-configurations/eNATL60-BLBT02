!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!! LIM3 namelist :  
!!              1 - Generic parameters                 (namicerun)
!!              2 - Ice initialization                 (namiceini)
!!              3 - Ice discretization                 (namiceitd)
!!              4 - Ice dynamics and transport         (namicedyn)
!!              5 - Ice thermodynamics                 (namicethd)
!!              6 - Ice salinity                       (namicesal)
!!              7 - Ice mechanical redistribution      (namiceitdme)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!------------------------------------------------------------------------------
&namicerun     !   Generic parameters
!------------------------------------------------------------------------------
   jpl              =    1                 ! #LOLO!  number of ice  categories
   cn_icerst_in     = 'eNATL60-BLBT02Y_01393200_restart_ice'     !  suffix of ice restart name (input)
   cn_icerst_indir  = '/gpfs/scratch/pr1egh00/pr1egh01/eNATL60/eNATL60-BLBT02Y-R/01393200'  !  directory from which to read input ice restarts
   cn_icerst_outdir = '/gpfs/scratch/pr1egh00/pr1egh01/eNATL60/eNATL60-BLBT02Y-R/VRAC_ICE' !  directory in which to write output ice restarts
/
!------------------------------------------------------------------------------
&namiceini     !   Ice initialization
!------------------------------------------------------------------------------
/
!------------------------------------------------------------------------------
&namiceitd     !   Ice discretization
!------------------------------------------------------------------------------
/
!------------------------------------------------------------------------------
&namicedyn     !   Ice dynamics and transport
!------------------------------------------------------------------------------
/
!------------------------------------------------------------------------------
&namicehdf     !   Ice horizontal diffusion
!------------------------------------------------------------------------------
   nn_ahi0        =  -1            ! #LOLO! Clem Rousset! => coute cher pour rien! horizontal diffusivity computation
/
!------------------------------------------------------------------------------
&namicethd     !   Ice thermodynamics
!------------------------------------------------------------------------------
   nn_monocat  = 1                 ! #LOLO! Clem Rousset! virtual ITD mono-category parameterizations (1, jpl = 1 only) or not (0)
                                   !     2: simple piling instead of ridging --- temporary option
                                   !     3: activate G(he) only              --- temporary option
                                   !     4: activate lateral melting only    --- temporary option
/
!------------------------------------------------------------------------------
&namicesal     !   Ice salinity
!------------------------------------------------------------------------------
/
!------------------------------------------------------------------------------
&namiceitdme   !   Ice mechanical redistribution (ridging and rafting)
!------------------------------------------------------------------------------
   rn_astar    = 0.03             ! #LOLO! Clem Rousset! exponential measure of ridging ice fraction (nn_partfun = 1)
   rn_hstar    = 25.0             ! #LOLO! Clem Rousset!  determines the maximum thickness of ridged ice (m) (Hibler, 1980)
/

