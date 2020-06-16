MODULE ldfdyn
!!======================================================================
!!                       ***  MODULE  ldfdyn  ***
!! Ocean physics:  lateral viscosity coefficient
!!=====================================================================
!! History :  OPA  ! 1997-07  (G. Madec)  multi dimensional coefficients
!!   NEMO     1.0  ! 2002-09  (G. Madec)  F90: Free form and module
!!----------------------------------------------------------------------

!!----------------------------------------------------------------------
!!   ldf_dyn_init : initialization, namelist read, and parameters control
!!   ldf_dyn_c3d   : 3D eddy viscosity coefficient initialization
!!   ldf_dyn_c2d   : 2D eddy viscosity coefficient initialization
!!   ldf_dyn_c1d   : 1D eddy viscosity coefficient initialization
!!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE ldfdyn_oce      ! ocean dynamics lateral physics
   USE phycst          ! physical constants
   USE ldfslp          ! ???
   USE ioipsl
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! distribued memory computing library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE wrk_nemo        ! Memory Allocation

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ldf_dyn_init   ! called by opa.F90

  INTERFACE ldf_zpf
     MODULE PROCEDURE ldf_zpf_1d, ldf_zpf_1d_3d, ldf_zpf_3d
  END INTERFACE

!! * Substitutions
!!----------------------------------------------------------------------
!!                    ***  domzgr_substitute.h90   ***
!!----------------------------------------------------------------------
!! ** purpose :   substitute fsdep. and fse.., the vert. depth and scale
!!      factors depending on the vertical coord. used, using CPP macro.
!!----------------------------------------------------------------------
!! History :  1.0  !  2005-10  (A. Beckmann, G. Madec) generalisation to all coord.
!!            3.1  !  2009-02  (G. Madec, M. Leclair)  pure z* coordinate
!!----------------------------------------------------------------------


! s* or z*-coordinate (3D + time dependency) + use of additional now arrays (..._n)





































! This part should be removed one day ...
! ... In that case all occurence of the above statement functions
!     have to be replaced in the code by xxx_n












!!----------------------------------------------------------------------
!! NEMO/OPA 3.3 , NEMO Consortium (2010)
!! $Id: domzgr_substitute.h90 4488 2014-02-06 10:43:09Z rfurner $
!! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
!!----------------------------------------------------------------------
!!----------------------------------------------------------------------
!! NEMO/OPA 3.3 , NEMO Consortium (2010)
!! $Id: ldfdyn.F90 4624 2014-04-28 12:09:03Z acc $
!! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
!!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ldf_dyn_init
!!----------------------------------------------------------------------
!!                  ***  ROUTINE ldf_dyn_init  ***
!!
!! ** Purpose :   set the horizontal ocean dynamics physics
!!
!! ** Method  :
!!      -  default option : ahm = constant coef. = rn_ahm_0 (namelist)
!!      - 'key_dynldf_c1d': ahm = F(depth)                     see ldf_dyn_c1d.h90
!!      - 'key_dynldf_c2d': ahm = F(latitude,longitude)        see ldf_dyn_c2d.h90
!!      - 'key_dynldf_c3d': ahm = F(latitude,longitude,depth)  see ldf_dyn_c3d.h90
!!
!!      N.B. User defined include files.  By default, 3d and 2d coef.
!!      are set to a constant value given in the namelist and the 1d
!!      coefficients are initialized to a hyperbolic tangent vertical
!!      profile.
!!
!! Reference :   Madec, G. and M. Imbard, 1996: Climate Dynamics, 12, 381-388.
!!----------------------------------------------------------------------
      INTEGER ::   ioptio         ! ???
      INTEGER ::   ios            ! Local : output status for namelist read
      LOGICAL ::   ll_print = .FALSE.    ! Logical flag for printing viscosity coef.
!!
      NAMELIST/namdyn_ldf/ ln_dynldf_lap  , ln_dynldf_bilap,                  &
         &                 ln_dynldf_level, ln_dynldf_hor  , ln_dynldf_iso,   &
         &                 rn_ahm_0_lap   , rn_ahmb_0      , rn_ahm_0_blp ,   &
         &                 rn_cmsmag_1    , rn_cmsmag_2    , rn_cmsh,         &
         &                 rn_ahm_m_lap   , rn_ahm_m_blp

!!----------------------------------------------------------------------

      REWIND( numnam_ref )              ! Namelist namdyn_ldf in reference namelist : Lateral physics
      READ  ( numnam_ref, namdyn_ldf, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namdyn_ldf in reference namelist', lwp )

      REWIND( numnam_cfg )              ! Namelist namdyn_ldf in configuration namelist : Lateral physics
      READ  ( numnam_cfg, namdyn_ldf, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namdyn_ldf in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, namdyn_ldf )

      IF(lwp) THEN                      ! Parameter print
         WRITE(numout,*)
         WRITE(numout,*) 'ldf_dyn : lateral momentum physics'
         WRITE(numout,*) '~~~~~~~'
         WRITE(numout,*) '   Namelist namdyn_ldf : set lateral mixing parameters'
         WRITE(numout,*) '      laplacian operator                      ln_dynldf_lap   = ', ln_dynldf_lap
         WRITE(numout,*) '      bilaplacian operator                    ln_dynldf_bilap = ', ln_dynldf_bilap
         WRITE(numout,*) '      iso-level                               ln_dynldf_level = ', ln_dynldf_level
         WRITE(numout,*) '      horizontal (geopotential)               ln_dynldf_hor   = ', ln_dynldf_hor
         WRITE(numout,*) '      iso-neutral                             ln_dynldf_iso   = ', ln_dynldf_iso
         WRITE(numout,*) '      horizontal laplacian eddy viscosity     rn_ahm_0_lap    = ', rn_ahm_0_lap
         WRITE(numout,*) '      background viscosity                    rn_ahmb_0       = ', rn_ahmb_0
         WRITE(numout,*) '      horizontal bilaplacian eddy viscosity   rn_ahm_0_blp    = ', rn_ahm_0_blp
         WRITE(numout,*) '      upper limit for laplacian eddy visc     rn_ahm_m_lap    = ', rn_ahm_m_lap
         WRITE(numout,*) '      upper limit for bilap eddy viscosity    rn_ahm_m_blp    = ', rn_ahm_m_blp

      ENDIF

      ahm0     = rn_ahm_0_lap              ! OLD namelist variables defined from DOCTOR namelist variables
      ahmb0    = rn_ahmb_0
      ahm0_blp = rn_ahm_0_blp

! ... check of lateral diffusive operator on tracers
!           ==> will be done in trazdf module

! ... Space variation of eddy coefficients
      ioptio = 0


      IF(lwp) WRITE(numout,*) '   momentum mixing coef. = F( latitude, longitude)'
      ioptio = ioptio+1


      IF( ioptio == 0 ) THEN
          IF(lwp) WRITE(numout,*) '   momentum mixing coef. = constant  (default option)'
        ELSEIF( ioptio > 1 ) THEN
           CALL ctl_stop( 'use only one of the following keys: key_dynldf_c3d, key_dynldf_c2d, key_dynldf_c1d' )
      ENDIF


      IF( ln_dynldf_bilap ) THEN
         IF(lwp) WRITE(numout,*) '   biharmonic momentum diffusion'
         IF( .NOT. ln_dynldf_lap ) ahm0 = ahm0_blp   ! Allow spatially varying coefs, which use ahm0 as input
         IF( ahm0_blp > 0 .AND. .NOT. lk_esopa )   CALL ctl_stop( 'The horizontal viscosity coef. ahm0 must be negative' )
      ELSE
         IF(lwp) WRITE(numout,*) '   harmonic momentum diff. (default)'
         IF( ahm0 < 0 .AND. .NOT. lk_esopa )   CALL ctl_stop( 'The horizontal viscosity coef. ahm0 must be positive' )
      ENDIF


! Lateral eddy viscosity
! ======================

      CALL ldf_dyn_c2d( ll_print )   ! ahm = 1D coef. = F( longitude, latitude )

     nkahm_smag = 0


!
   END SUBROUTINE ldf_dyn_init


!!----------------------------------------------------------------------
!!                      ***  ldfdyn_c2d.h90  ***
!!----------------------------------------------------------------------
!!   ldf_dyn_c2d  : set the lateral viscosity coefficients
!!   ldf_dyn_c2d_orca : specific case for orca r2 and r4
!!----------------------------------------------------------------------

!!----------------------------------------------------------------------
!! NEMO/OPA 3.3 , NEMO Consortium (2010)
!! $Id: ldfdyn_c2d.h90 5400 2015-06-10 15:29:08Z cbricaud $
!! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
!!----------------------------------------------------------------------

   SUBROUTINE ldf_dyn_c2d( ld_print )
!!----------------------------------------------------------------------
!!                 ***  ROUTINE ldf_dyn_c2d  ***
!!
!! ** Purpose :   initializations of the horizontal ocean physics
!!
!! ** Method :
!!      2D eddy viscosity coefficients ( longitude, latitude )
!!
!!       harmonic operator   : ahm1 is defined at t-point
!!                             ahm2 is defined at f-point
!!           + isopycnal     : ahm3 is defined at u-point
!!           or geopotential   ahm4 is defined at v-point
!!           iso-model level : ahm3, ahm4 not used
!!
!!       biharmonic operator : ahm3 is defined at u-point
!!                             ahm4 is defined at v-point
!!                           : ahm1, ahm2 not used
!!
!!----------------------------------------------------------------------
      LOGICAL, INTENT (in) :: ld_print   ! If true, output arrays on numout
!
      INTEGER  ::   ji, jj
      REAL(wp) ::   za00, zd_max, zetmax, zeumax, zefmax, zevmax
!!----------------------------------------------------------------------

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'ldf_dyn_c2d : 2d lateral eddy viscosity coefficient'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'

! harmonic operator (ahm1, ahm2) : ( T- and F- points) (used for laplacian operators
! ===============================                       whatever its orientation is)
      IF( ln_dynldf_lap ) THEN
! define ahm1 and ahm2 at the right grid point position
! (USER: modify ahm1 and ahm2 following your desiderata)

         zd_max = MAX( MAXVAL( e1t(:,:) ), MAXVAL( e2t(:,:) ) )
         IF( lk_mpp )   CALL mpp_max( zd_max )   ! max over the global domain

         IF(lwp) WRITE(numout,*) '              laplacian operator: ahm proportional to e1'
         IF(lwp) WRITE(numout,*) '              maximum grid-spacing = ', zd_max, ' maximum value for ahm = ', ahm0

         za00 = ahm0 / zd_max
         DO jj = 1, jpj
            DO ji = 1, jpi
               zetmax = MAX( e1t(ji,jj), e2t(ji,jj) )
               zefmax = MAX( e1f(ji,jj), e2f(ji,jj) )
               ahm1(ji,jj) = za00 * zetmax
               ahm2(ji,jj) = za00 * zefmax
            END DO
         END DO

         IF( ln_dynldf_iso ) THEN
            IF(lwp) WRITE(numout,*) '              Caution, as implemented now, the isopycnal part of momentum'
            IF(lwp) WRITE(numout,*) '                 mixing use aht0 as eddy viscosity coefficient. Thus, it is'
            IF(lwp) WRITE(numout,*) '                 uniform and you must be sure that your ahm is greater than'
            IF(lwp) WRITE(numout,*) '                 aht0 everywhere in the model domain.'
         ENDIF

! Special case for ORCA R1, R2 and R4 configurations (overwrite the value of ahm1 ahm2)
! ==============================================
         IF( cp_cfg == "orca" .AND. ( jp_cfg == 2 .OR. jp_cfg == 4 ) )   CALL ldf_dyn_c2d_orca( ld_print )
         IF( cp_cfg == "orca" .AND.   jp_cfg == 1)                       CALL ldf_dyn_c2d_orca_R1( ld_print )

! Control print
         IF( lwp .AND. ld_print ) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'inildf: 2D ahm1 array'
            CALL prihre(ahm1,jpi,jpj,1,jpi,1,1,jpj,1,1.e-3,numout)
            WRITE(numout,*)
            WRITE(numout,*) 'inildf: 2D ahm2 array'
            CALL prihre(ahm2,jpi,jpj,1,jpi,1,1,jpj,1,1.e-3,numout)
         ENDIF
      ENDIF

! biharmonic operator (ahm3, ahm4) : at U- and V-points (used for bilaplacian operator
! =================================                      whatever its orientation is)
      IF( ln_dynldf_bilap ) THEN
! (USER: modify ahm3 and ahm4 following your desiderata)
! Here: ahm is proportional to the cube of the maximum of the gridspacing
!       in the to horizontal direction

         zd_max = MAX( MAXVAL( e1u(:,:) ), MAXVAL( e2u(:,:) ) )
         IF( lk_mpp )   CALL mpp_max( zd_max )   ! max over the global domain

         IF(lwp) WRITE(numout,*) '              bi-laplacian operator: ahm proportional to e1**3 '
         IF(lwp) WRITE(numout,*) '              maximum grid-spacing = ', zd_max, ' maximum value for ahm = ', ahm0

         za00 = ahm0_blp / ( zd_max * zd_max * zd_max )
         DO jj = 1, jpj
            DO ji = 1, jpi
               zeumax = MAX( e1u(ji,jj), e2u(ji,jj) )
               zevmax = MAX( e1v(ji,jj), e2v(ji,jj) )
               ahm3(ji,jj) = za00 * zeumax * zeumax * zeumax
               ahm4(ji,jj) = za00 * zevmax * zevmax * zevmax
            END DO
         END DO

! Control print
         IF( lwp .AND. ld_print ) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'inildf: ahm3 array'
            CALL prihre(ahm3,jpi,jpj,1,jpi,1,1,jpj,1,1.e-3,numout)
            WRITE(numout,*)
            WRITE(numout,*) 'inildf: ahm4 array'
            CALL prihre(ahm4,jpi,jpj,1,jpi,1,1,jpj,1,1.e-3,numout)
         ENDIF
      ENDIF
!
   END SUBROUTINE ldf_dyn_c2d


   SUBROUTINE ldf_dyn_c2d_orca( ld_print )
!!----------------------------------------------------------------------
!!                 ***  ROUTINE ldf_dyn_c2d  ***
!!
!!                   **** W A R N I N G ****
!!
!!                ORCA R2 and R4 configurations
!!
!!                   **** W A R N I N G ****
!!
!! ** Purpose :   initializations of the lateral viscosity for orca R2
!!
!! ** Method  :   blah blah blah...
!!
!!----------------------------------------------------------------------
      USE ldftra_oce, ONLY:   aht0
      USE iom
!
      LOGICAL, INTENT (in) ::   ld_print   ! If true, output arrays on numout
!
      INTEGER  ::   ji, jj, jn   ! dummy loop indices
      INTEGER  ::   inum, iim, ijm            ! local integers
      INTEGER  ::   ifreq, il1, il2, ij, ii
      INTEGER  ::   ijpt0,ijpt1, ierror
      REAL(wp) ::   zahmeq, zcoft, zcoff, zmsk
      CHARACTER (len=15) ::   clexp
      INTEGER,     POINTER, DIMENSION(:,:)  :: icof
      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: ztemp2d  ! temporary array to read ahmcoef file
!!----------------------------------------------------------------------
!
      CALL wrk_alloc( jpi   , jpj   , icof  )
!
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'inildf: 2d eddy viscosity coefficient'
      IF(lwp) WRITE(numout,*) '~~~~~~  --'
      IF(lwp) WRITE(numout,*) '        orca ocean configuration'

      IF( cp_cfg == "orca" .AND. cp_cfz == "antarctic" ) THEN
!
! 1.2 Modify ahm
! --------------
         IF(lwp)WRITE(numout,*) ' inildf: Antarctic ocean'
         IF(lwp)WRITE(numout,*) '         no tropics, no reduction of ahm'
         IF(lwp)WRITE(numout,*) '         north boundary increase'

         ahm1(:,:) = ahm0
         ahm2(:,:) = ahm0

         ijpt0=max(1,min(49 -njmpp+1,jpj))
         ijpt1=max(0,min(49-njmpp+1,jpj-1))
         DO jj=ijpt0,ijpt1
            ahm2(:,jj)=ahm0*2.
            ahm1(:,jj)=ahm0*2.
         END DO
         ijpt0=max(1,min(48 -njmpp+1,jpj))
         ijpt1=max(0,min(48-njmpp+1,jpj-1))
         DO jj=ijpt0,ijpt1
            ahm2(:,jj)=ahm0*1.9
            ahm1(:,jj)=ahm0*1.75
         END DO
         ijpt0=max(1,min(47 -njmpp+1,jpj))
         ijpt1=max(0,min(47-njmpp+1,jpj-1))
         DO jj=ijpt0,ijpt1
            ahm2(:,jj)=ahm0*1.5
            ahm1(:,jj)=ahm0*1.25
         END DO
         ijpt0=max(1,min(46 -njmpp+1,jpj))
         ijpt1=max(0,min(46-njmpp+1,jpj-1))
         DO jj=ijpt0,ijpt1
            ahm2(:,jj)=ahm0*1.1
         END DO

      ELSE IF( cp_cfg == "orca" .AND. cp_cfz == "arctic" ) THEN
! 1.2 Modify ahm
! --------------
         IF(lwp)WRITE(numout,*) ' inildf: Arctic ocean'
         IF(lwp)WRITE(numout,*) '         no tropics, no reduction of ahm'
         IF(lwp)WRITE(numout,*) '         south and west boundary increase'


         ahm1(:,:) = ahm0
         ahm2(:,:) = ahm0

         ijpt0=max(1,min(98-jpjzoom+1-njmpp+1,jpj))
         ijpt1=max(0,min(98-jpjzoom+1-njmpp+1,jpj-1))
         DO jj=ijpt0,ijpt1
            ahm2(:,jj)=ahm0*2.
            ahm1(:,jj)=ahm0*2.
         END DO
         ijpt0=max(1,min(99-jpjzoom+1-njmpp+1,jpj))
         ijpt1=max(0,min(99-jpjzoom+1-njmpp+1,jpj-1))
         DO jj=ijpt0,ijpt1
            ahm2(:,jj)=ahm0*1.9
            ahm1(:,jj)=ahm0*1.75
         END DO
         ijpt0=max(1,min(100-jpjzoom+1-njmpp+1,jpj))
         ijpt1=max(0,min(100-jpjzoom+1-njmpp+1,jpj-1))
         DO jj=ijpt0,ijpt1
            ahm2(:,jj)=ahm0*1.5
            ahm1(:,jj)=ahm0*1.25
         END DO
         ijpt0=max(1,min(101-jpjzoom+1-njmpp+1,jpj))
         ijpt1=max(0,min(101-jpjzoom+1-njmpp+1,jpj-1))
         DO jj=ijpt0,ijpt1
            ahm2(:,jj)=ahm0*1.1
         END DO
      ELSE
! Read 2d integer array to specify western boundary increase in the
! ===================== equatorial strip (20N-20S) defined at t-points
!
         ALLOCATE( ztemp2d(jpi,jpj) )
         ztemp2d(:,:) = 0.
         CALL iom_open ( 'ahmcoef.nc', inum )
         CALL iom_get  ( inum, jpdom_data, 'icof', ztemp2d)
         icof(:,:)  = NINT(ztemp2d(:,:))
         CALL iom_close( inum )
         DEALLOCATE(ztemp2d)

! Set ahm1 and ahm2  ( T- and F- points) (used for laplacian operator)
! =================
! define ahm1 and ahm2 at the right grid point position
! (USER: modify ahm1 and ahm2 following your desiderata)


! Decrease ahm to zahmeq m2/s in the tropics
! (from 90 to 20 degre: ahm = constant
! from 20 to  2.5 degre: ahm = decrease in (1-cos)/2
! from  2.5 to  0 degre: ahm = constant
! symmetric in the south hemisphere)

         zahmeq = aht0

         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( ABS( gphif(ji,jj) ) >= 20. ) THEN
                  ahm2(ji,jj) =  ahm0
               ELSEIF( ABS( gphif(ji,jj) ) <= 2.5 ) THEN
                  ahm2(ji,jj) =  zahmeq
               ELSE
                  ahm2(ji,jj) = zahmeq + (ahm0-zahmeq)/2.   &
                     * ( 1. - COS( rad * ( ABS(gphif(ji,jj))-2.5 ) * 180. / 17.5 ) )
               ENDIF
               IF( ABS( gphit(ji,jj) ) >= 20. ) THEN
                  ahm1(ji,jj) =  ahm0
               ELSEIF( ABS( gphit(ji,jj) ) <= 2.5 ) THEN
                  ahm1(ji,jj) =  zahmeq
               ELSE
                  ahm1(ji,jj) = zahmeq + (ahm0-zahmeq)/2.   &
                     * ( 1. - COS( rad * ( ABS(gphit(ji,jj))-2.5 ) * 180. / 17.5 ) )
               ENDIF
            END DO
         END DO

! increase along western boundaries of equatorial strip
! t-point
         DO jj = 1, jpjm1
            DO ji = 1, jpim1
               zcoft = FLOAT( icof(ji,jj) ) / 100.
               ahm1(ji,jj) = zcoft * ahm0 + (1.-zcoft) * ahm1(ji,jj) 
            END DO
         END DO
! f-point
         icof(:,:) = icof(:,:) * tmask(:,:,1)
         DO jj = 1, jpjm1
            DO ji = 1, jpim1   ! NO vector opt.
               zmsk = tmask(ji,jj+1,1) + tmask(ji+1,jj+1,1) + tmask(ji,jj,1) + tmask(ji,jj+1,1)
               IF( zmsk == 0. ) THEN
                  zcoff = 1.
               ELSE
                  zcoff = FLOAT( icof(ji,jj+1) + icof(ji+1,jj+1) + icof(ji,jj) + icof(ji,jj+1) )   &
                     / (zmsk * 100.)
               ENDIF
               ahm2(ji,jj) = zcoff * ahm0 + (1.-zcoff) * ahm2(ji,jj)
            END DO
         END DO
      ENDIF
      
! Lateral boundary conditions on ( ahm1, ahm2 )
!                                ==============
      CALL lbc_lnk( ahm1, 'T', 1. )   ! T-point, unchanged sign
      CALL lbc_lnk( ahm2, 'F', 1. )   ! F-point, unchanged sign

! Control print
      IF( lwp .AND. ld_print ) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'inildf: 2D ahm1 array'
         CALL prihre(ahm1,jpi,jpj,1,jpi,1,1,jpj,1,1.e-3,numout)
         WRITE(numout,*)
         WRITE(numout,*) 'inildf: 2D ahm2 array'
         CALL prihre(ahm2,jpi,jpj,1,jpi,1,1,jpj,1,1.e-3,numout)
      ENDIF
!
      CALL wrk_dealloc( jpi   , jpj   , icof  )
!
   END SUBROUTINE ldf_dyn_c2d_orca


   SUBROUTINE ldf_dyn_c2d_orca_R1( ld_print )
!!----------------------------------------------------------------------
!!                 ***  ROUTINE ldf_dyn_c2d  ***
!!
!!                   **** W A R N I N G ****
!!
!!                ORCA R1 configuration
!!
!!                   **** W A R N I N G ****
!!
!! ** Purpose :   initializations of the lateral viscosity for orca R1
!!
!! ** Method  :   blah blah blah...
!!
!!----------------------------------------------------------------------
      USE ldftra_oce, ONLY:   aht0
      USE iom
!
      LOGICAL, INTENT (in) ::   ld_print   ! If true, output arrays on numout
!
      INTEGER ::   ji, jj, jn      ! dummy loop indices
      INTEGER ::   inum            ! temporary logical unit
      INTEGER ::   iim, ijm
      INTEGER ::   ifreq, il1, il2, ij, ii
      INTEGER ::   ijpt0,ijpt1, ierror
      REAL(wp) ::   zahmeq, zcoft, zcoff, zmsk, zam20s
      CHARACTER (len=15) ::   clexp
      INTEGER,     POINTER, DIMENSION(:,:)  :: icof
      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: ztemp2d  ! temporary array to read ahmcoef file
!!----------------------------------------------------------------------
!
      CALL wrk_alloc( jpi   , jpj   , icof  )
!
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'inildf: 2d eddy viscosity coefficient'
      IF(lwp) WRITE(numout,*) '~~~~~~  --'
      IF(lwp) WRITE(numout,*) '        orca_r1 configuration'

      IF( cp_cfg == "orca" .AND. cp_cfz == "antarctic" ) THEN
!
! 1.2 Modify ahm
! --------------
         IF(lwp)WRITE(numout,*) ' inildf: Antarctic ocean'
         IF(lwp)WRITE(numout,*) '         no tropics, no reduction of ahm'
         IF(lwp)WRITE(numout,*) '         north boundary increase'

         ahm1(:,:) = ahm0
         ahm2(:,:) = ahm0

         ijpt0=max(1,min(49 -njmpp+1,jpj))
         ijpt1=max(0,min(49-njmpp+1,jpj-1))
         DO jj=ijpt0,ijpt1
            ahm2(:,jj)=ahm0*2.
            ahm1(:,jj)=ahm0*2.
         END DO
         ijpt0=max(1,min(48 -njmpp+1,jpj))
         ijpt1=max(0,min(48-njmpp+1,jpj-1))
         DO jj=ijpt0,ijpt1
            ahm2(:,jj)=ahm0*1.9
            ahm1(:,jj)=ahm0*1.75
         END DO
         ijpt0=max(1,min(47 -njmpp+1,jpj))
         ijpt1=max(0,min(47-njmpp+1,jpj-1))
         DO jj=ijpt0,ijpt1
            ahm2(:,jj)=ahm0*1.5
            ahm1(:,jj)=ahm0*1.25
         END DO
         ijpt0=max(1,min(46 -njmpp+1,jpj))
         ijpt1=max(0,min(46-njmpp+1,jpj-1))
         DO jj=ijpt0,ijpt1
            ahm2(:,jj)=ahm0*1.1
         END DO

      ELSE IF( cp_cfg == "orca" .AND. cp_cfz == "arctic" ) THEN
! 1.2 Modify ahm
! --------------
         IF(lwp)WRITE(numout,*) ' inildf: Arctic ocean'
         IF(lwp)WRITE(numout,*) '         no tropics, no reduction of ahm'
         IF(lwp)WRITE(numout,*) '         south and west boundary increase'


         ahm1(:,:) = ahm0
         ahm2(:,:) = ahm0

         ijpt0=max(1,min(98-jpjzoom+1-njmpp+1,jpj))
         ijpt1=max(0,min(98-jpjzoom+1-njmpp+1,jpj-1))
         DO jj=ijpt0,ijpt1
            ahm2(:,jj)=ahm0*2.
            ahm1(:,jj)=ahm0*2.
         END DO
         ijpt0=max(1,min(99-jpjzoom+1-njmpp+1,jpj))
         ijpt1=max(0,min(99-jpjzoom+1-njmpp+1,jpj-1))
         DO jj=ijpt0,ijpt1
            ahm2(:,jj)=ahm0*1.9
            ahm1(:,jj)=ahm0*1.75
         END DO
         ijpt0=max(1,min(100-jpjzoom+1-njmpp+1,jpj))
         ijpt1=max(0,min(100-jpjzoom+1-njmpp+1,jpj-1))
         DO jj=ijpt0,ijpt1
            ahm2(:,jj)=ahm0*1.5
            ahm1(:,jj)=ahm0*1.25
         END DO
         ijpt0=max(1,min(101-jpjzoom+1-njmpp+1,jpj))
         ijpt1=max(0,min(101-jpjzoom+1-njmpp+1,jpj-1))
         DO jj=ijpt0,ijpt1
            ahm2(:,jj)=ahm0*1.1
         END DO
      ELSE
         
! Read 2d integer array to specify western boundary increase in the
! ===================== equatorial strip (20N-20S) defined at t-points
         ALLOCATE( ztemp2d(jpi,jpj) )
         ztemp2d(:,:) = 0.
         CALL iom_open ( 'ahmcoef.nc', inum )
         CALL iom_get  ( inum, jpdom_data, 'icof', ztemp2d)
         icof(:,:)  = NINT(ztemp2d(:,:))
         CALL iom_close( inum )
         DEALLOCATE(ztemp2d)

! Set ahm1 and ahm2  ( T- and F- points) (used for laplacian operator)
! =================
! define ahm1 and ahm2 at the right grid point position
! (USER: modify ahm1 and ahm2 following your desiderata)


! Decrease ahm to zahmeq m2/s in the tropics
! (from 90   to 20   degrees: ahm = scaled by local metrics
!  from 20   to  2.5 degrees: ahm = decrease in (1-cos)/2
!  from  2.5 to  0   degrees: ahm = constant
! symmetric in the south hemisphere)

         zahmeq = aht0
         zam20s = ahm0*COS( rad * 20. )

         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( ABS( gphif(ji,jj) ) >= 20. ) THEN
!              leave as set in ldf_dyn_c2d
               ELSEIF( ABS( gphif(ji,jj) ) <= 2.5 ) THEN
                  ahm2(ji,jj) =  zahmeq
               ELSE
                  ahm2(ji,jj) =  zahmeq + (zam20s-zahmeq)/2.   &
                     * ( 1. - COS( rad * ( ABS(gphif(ji,jj))-2.5 ) * 180. / 17.5 ) )
               ENDIF
               IF( ABS( gphit(ji,jj) ) >= 20. ) THEN
!             leave as set in ldf_dyn_c2d
               ELSEIF( ABS( gphit(ji,jj) ) <= 2.5 ) THEN
                  ahm1(ji,jj) =  zahmeq
               ELSE
                  ahm1(ji,jj) =  zahmeq + (zam20s-zahmeq)/2.   &
                     * ( 1. - COS( rad * ( ABS(gphit(ji,jj))-2.5 ) * 180. / 17.5 ) )
               ENDIF
            END DO
         END DO

! increase along western boundaries of equatorial strip
! t-point
         DO jj = 1, jpjm1
            DO ji = 1, jpim1
               IF( ABS( gphit(ji,jj) ) < 20. ) THEN
                  zcoft = FLOAT( icof(ji,jj) ) / 100.
                  ahm1(ji,jj) = zcoft * ahm0 + (1.-zcoft) * ahm1(ji,jj) 
               ENDIF
            END DO
         END DO
! f-point
         icof(:,:) = icof(:,:) * tmask(:,:,1)
         DO jj = 1, jpjm1
            DO ji = 1, jpim1
               IF( ABS( gphif(ji,jj) ) < 20. ) THEN
                  zmsk = tmask(ji,jj+1,1) + tmask(ji+1,jj+1,1) + tmask(ji,jj,1) + tmask(ji,jj+1,1)
                  IF( zmsk == 0. ) THEN
                     zcoff = 1.
                  ELSE
                     zcoff = FLOAT( icof(ji,jj+1) + icof(ji+1,jj+1) + icof(ji,jj) + icof(ji,jj+1) )   &
                        / (zmsk * 100.)
                  ENDIF
                  ahm2(ji,jj) = zcoff * ahm0 + (1.-zcoff) * ahm2(ji,jj)
               ENDIF
            END DO
         END DO
      ENDIF
      
! Lateral boundary conditions on ( ahm1, ahm2 )
!                                ==============
      CALL lbc_lnk( ahm1, 'T', 1. )   ! T-point, unchanged sign
      CALL lbc_lnk( ahm2, 'F', 1. )   ! F-point, unchanged sign

! Control print
      IF( lwp .AND. ld_print ) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'inildf: 2D ahm1 array'
         CALL prihre(ahm1,jpi,jpj,1,jpi,1,1,jpj,1,1.e-3,numout)
         WRITE(numout,*)
         WRITE(numout,*) 'inildf: 2D ahm2 array'
         CALL prihre(ahm2,jpi,jpj,1,jpi,1,1,jpj,1,1.e-3,numout)
      ENDIF
!
      CALL wrk_dealloc( jpi   , jpj   , icof  )
!
   END SUBROUTINE ldf_dyn_c2d_orca_R1



   SUBROUTINE ldf_zpf_1d( ld_print, pdam, pwam, pbot, pdep, pah )
!!----------------------------------------------------------------------
!!                  ***  ROUTINE ldf_zpf  ***
!!
!! ** Purpose :   vertical adimensional profile for eddy coefficient
!!
!! ** Method  :   1D eddy viscosity coefficients ( depth )
!!----------------------------------------------------------------------
      LOGICAL , INTENT(in   )                 ::   ld_print   ! If true, output arrays on numout
      REAL(wp), INTENT(in   )                 ::   pdam       ! depth of the inflection point
      REAL(wp), INTENT(in   )                 ::   pwam       ! width of inflection
      REAL(wp), INTENT(in   )                 ::   pbot       ! bottom value (0<pbot<= 1)
      REAL(wp), INTENT(in   ), DIMENSION(jpk) ::   pdep       ! depth of the gridpoint (T, U, V, F)
      REAL(wp), INTENT(inout), DIMENSION(jpk) ::   pah        ! adimensional vertical profile
!!
      INTEGER  ::   jk           ! dummy loop indices
      REAL(wp) ::   zm00, zm01, zmhb, zmhs       ! temporary scalars
!!----------------------------------------------------------------------

      zm00 = TANH( ( pdam - gdept_1d(1    ) ) / pwam )
      zm01 = TANH( ( pdam - gdept_1d(jpkm1) ) / pwam )
      zmhs = zm00 / zm01
      zmhb = ( 1.e0 - pbot ) / ( 1.e0 - zmhs ) / zm01

      DO jk = 1, jpk
         pah(jk) = 1.e0 + zmhb * ( zm00 - TANH( ( pdam - pdep(jk) ) / pwam )  )
      END DO

      IF(lwp .AND. ld_print ) THEN      ! Control print
         WRITE(numout,*)
         WRITE(numout,*) '         ahm profile : '
         WRITE(numout,*)
         WRITE(numout,'("  jk      ahm       ","  depth t-level " )')
         DO jk = 1, jpk
            WRITE(numout,'(i6,2f12.4,3x,2f12.4)') jk, pah(jk), pdep(jk)
         END DO
      ENDIF
!
   END SUBROUTINE ldf_zpf_1d


   SUBROUTINE ldf_zpf_1d_3d( ld_print, pdam, pwam, pbot, pdep, pah )
!!----------------------------------------------------------------------
!!                  ***  ROUTINE ldf_zpf  ***
!!
!! ** Purpose :   vertical adimensional profile for eddy coefficient
!!
!! ** Method  :   1D eddy viscosity coefficients ( depth )
!!----------------------------------------------------------------------
      LOGICAL , INTENT(in   )                         ::   ld_print   ! If true, output arrays on numout
      REAL(wp), INTENT(in   )                         ::   pdam       ! depth of the inflection point
      REAL(wp), INTENT(in   )                         ::   pwam       ! width of inflection
      REAL(wp), INTENT(in   )                         ::   pbot       ! bottom value (0<pbot<= 1)
      REAL(wp), INTENT(in   ), DIMENSION          (:) ::   pdep       ! depth of the gridpoint (T, U, V, F)
      REAL(wp), INTENT(inout), DIMENSION      (:,:,:) ::   pah        ! adimensional vertical profile
!!
      INTEGER  ::   jk           ! dummy loop indices
      REAL(wp) ::   zm00, zm01, zmhb, zmhs, zcf  ! temporary scalars
!!----------------------------------------------------------------------

      zm00 = TANH( ( pdam - gdept_1d(1    ) ) / pwam )
      zm01 = TANH( ( pdam - gdept_1d(jpkm1) ) / pwam )
      zmhs = zm00 / zm01
      zmhb = ( 1.e0 - pbot ) / ( 1.e0 - zmhs ) / zm01

      DO jk = 1, jpk
         zcf = 1.e0 + zmhb * ( zm00 - TANH( ( pdam - pdep(jk) ) / pwam )  )
         pah(:,:,jk) = zcf
      END DO

      IF(lwp .AND. ld_print ) THEN      ! Control print
         WRITE(numout,*)
         WRITE(numout,*) '         ahm profile : '
         WRITE(numout,*)
         WRITE(numout,'("  jk      ahm       ","  depth t-level " )')
         DO jk = 1, jpk
            WRITE(numout,'(i6,2f12.4,3x,2f12.4)') jk, pah(1,1,jk), pdep(jk)
         END DO
      ENDIF
!
   END SUBROUTINE ldf_zpf_1d_3d


   SUBROUTINE ldf_zpf_3d( ld_print, pdam, pwam, pbot, pdep, pah )
!!----------------------------------------------------------------------
!!                  ***  ROUTINE ldf_zpf  ***
!!
!! ** Purpose :   vertical adimensional profile for eddy coefficient
!!
!! ** Method  :   3D for partial step or s-coordinate
!!----------------------------------------------------------------------
      LOGICAL , INTENT(in   )                         ::   ld_print   ! If true, output arrays on numout
      REAL(wp), INTENT(in   )                         ::   pdam       ! depth of the inflection point
      REAL(wp), INTENT(in   )                         ::   pwam       ! width of inflection
      REAL(wp), INTENT(in   )                         ::   pbot       ! bottom value (0<pbot<= 1)
      REAL(wp), INTENT(in   ), DIMENSION      (:,:,:) ::   pdep       ! dep of the gridpoint (T, U, V, F)
      REAL(wp), INTENT(inout), DIMENSION      (:,:,:) ::   pah        ! adimensional vertical profile
!!
      INTEGER  ::   jk           ! dummy loop indices
      REAL(wp) ::   zm00, zm01, zmhb, zmhs       ! temporary scalars
!!----------------------------------------------------------------------

      zm00 = TANH( ( pdam - gdept_1d(1    ) ) / pwam )   
      zm01 = TANH( ( pdam - gdept_1d(jpkm1) ) / pwam )
      zmhs = zm00 / zm01
      zmhb = ( 1.e0 - pbot ) / ( 1.e0 - zmhs ) / zm01

      DO jk = 1, jpk
         pah(:,:,jk) = 1.e0 + zmhb * ( zm00 - TANH( ( pdam - pdep(:,:,jk) ) / pwam )  )
      END DO

      IF(lwp .AND. ld_print ) THEN      ! Control print
         WRITE(numout,*)
         WRITE(numout,*) '         ahm profile : '
         WRITE(numout,*)
         WRITE(numout,'("  jk      ahm       ","  depth t-level " )')
         DO jk = 1, jpk
            WRITE(numout,'(i6,2f12.4,3x,2f12.4)') jk, pah(1,1,jk), pdep(1,1,jk)
         END DO
      ENDIF
!
   END SUBROUTINE ldf_zpf_3d

!!======================================================================
END MODULE ldfdyn
