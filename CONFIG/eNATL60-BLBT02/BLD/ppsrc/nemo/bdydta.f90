MODULE bdydta
!!======================================================================
!!                       ***  MODULE bdydta  ***
!! Open boundary data : read the data for the unstructured open boundaries.
!!======================================================================
!! History :  1.0  !  2005-01  (J. Chanut, A. Sellar)  Original code
!!             -   !  2007-01  (D. Storkey) Update to use IOM module
!!             -   !  2007-07  (D. Storkey) add bdy_dta_fla
!!            3.0  !  2008-04  (NEMO team)  add in the reference version
!!            3.3  !  2010-09  (E.O'Dea) modifications for Shelf configurations
!!            3.3  !  2010-09  (D.Storkey) add ice boundary conditions
!!            3.4  !  2011     (D. Storkey) rewrite in preparation for OBC-BDY merge
!!            3.6  !  2012-01  (C. Rousset) add ice boundary conditions for lim3
!!----------------------------------------------------------------------

!!----------------------------------------------------------------------
!!   'key_bdy'                     Open Boundary Conditions
!!----------------------------------------------------------------------
!!    bdy_dta        : read external data along open boundaries from file
!!    bdy_dta_init   : initialise arrays etc for reading of external data
!!----------------------------------------------------------------------
   USE timing          ! Timing
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE bdy_oce         ! ocean open boundary conditions
   USE bdytides        ! tidal forcing at boundaries
   USE fldread         ! read input fields
   USE iom             ! IOM library
   USE in_out_manager  ! I/O logical units
   USE dynspg_oce, ONLY: lk_dynspg_ts ! Split-explicit free surface flag

   USE ice
   USE limvar          ! redistribute ice input into categories

   USE sbc_oce
   USE sbcapr

   IMPLICIT NONE
   PRIVATE

   PUBLIC   bdy_dta          ! routine called by step.F90 and dynspg_ts.F90
   PUBLIC   bdy_dta_init     ! routine called by nemogcm.F90

   INTEGER, ALLOCATABLE, DIMENSION(:)   ::   nb_bdy_fld        ! Number of fields to update for each boundary set.
   INTEGER                              ::   nb_bdy_fld_sum    ! Total number of fields to update for all boundary sets.

   LOGICAL,           DIMENSION(jp_bdy) ::   ln_full_vel_array ! =T => full velocities in 3D boundary conditions
! =F => baroclinic velocities in 3D boundary conditions
!$AGRIF_DO_NOT_TREAT
   TYPE(FLD), PUBLIC, ALLOCATABLE, DIMENSION(:), TARGET ::   bf        ! structure of input fields (file informations, fields read)
!$AGRIF_END_DO_NOT_TREAT
   TYPE(MAP_POINTER), ALLOCATABLE, DIMENSION(:) :: nbmap_ptr   ! array of pointers to nbmap


   LOGICAL :: ll_bdylim3                  ! determine whether ice input is lim2 (F) or lim3 (T) type
   INTEGER :: jfld_hti, jfld_hts, jfld_ai ! indices of ice thickness, snow thickness and concentration in bf structure
   INTEGER, DIMENSION(jp_bdy) :: jfld_htit, jfld_htst, jfld_ait ! indices of ice thickness, snow thickness and concentration in bf structure


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
!! $Id: bdydta.F90 9811 2018-06-19 15:15:12Z clem $
!! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
!!----------------------------------------------------------------------
CONTAINS

      SUBROUTINE bdy_dta( kt, jit, time_offset )
!!----------------------------------------------------------------------
!!                   ***  SUBROUTINE bdy_dta  ***
!!
!! ** Purpose :   Update external data for open boundary conditions
!!
!! ** Method  :   Use fldread.F90
!!
!!----------------------------------------------------------------------
!!
      INTEGER, INTENT( in )           ::   kt    ! ocean time-step index
      INTEGER, INTENT( in ), OPTIONAL ::   jit   ! subcycle time-step index (for timesplitting option)
      INTEGER, INTENT( in ), OPTIONAL ::   time_offset  ! time offset in units of timesteps. NB. if jit
! is present then units = subcycle timesteps.
! time_offset = 0 => get data at "now" time level
! time_offset = -1 => get data at "before" time level
! time_offset = +1 => get data at "after" time level
! etc.
!!
      INTEGER     ::  ib_bdy, jfld, jstart, jend, ib, ii, ij, ik, igrd, jl  ! local indices
      INTEGER,          DIMENSION(jpbgrd) ::   ilen1 
      INTEGER, POINTER, DIMENSION(:)      ::   nblen, nblenrim  ! short cuts
      TYPE(OBC_DATA), POINTER             ::   dta              ! short cut
!!
!!---------------------------------------------------------------------------
!!
      IF( nn_timing == 1 ) CALL timing_start('bdy_dta')

! Initialise data arrays once for all from initial conditions where required
!---------------------------------------------------------------------------
      IF( kt .eq. nit000 .and. .not. PRESENT(jit) ) THEN

! Calculate depth-mean currents
!-----------------------------
         
         DO ib_bdy = 1, nb_bdy

            nblen => idx_bdy(ib_bdy)%nblen
            nblenrim => idx_bdy(ib_bdy)%nblenrim
            dta => dta_bdy(ib_bdy)

            IF( nn_dyn2d_dta(ib_bdy) .eq. 0 ) THEN 
               ilen1(:) = nblen(:)
               IF( dta%ll_ssh ) THEN 
                  igrd = 1
                  DO ib = 1, ilen1(igrd)
                     ii = idx_bdy(ib_bdy)%nbi(ib,igrd)
                     ij = idx_bdy(ib_bdy)%nbj(ib,igrd)
                     dta_bdy(ib_bdy)%ssh(ib) = sshn(ii,ij) * tmask(ii,ij,1)         
                  END DO 
               END IF
               IF( dta%ll_u2d ) THEN 
                  igrd = 2
                  DO ib = 1, ilen1(igrd)
                     ii = idx_bdy(ib_bdy)%nbi(ib,igrd)
                     ij = idx_bdy(ib_bdy)%nbj(ib,igrd)
                     dta_bdy(ib_bdy)%u2d(ib) = un_b(ii,ij) * umask(ii,ij,1)         
                  END DO 
               END IF
               IF( dta%ll_v2d ) THEN 
                  igrd = 3
                  DO ib = 1, ilen1(igrd)
                     ii = idx_bdy(ib_bdy)%nbi(ib,igrd)
                     ij = idx_bdy(ib_bdy)%nbj(ib,igrd)
                     dta_bdy(ib_bdy)%v2d(ib) = vn_b(ii,ij) * vmask(ii,ij,1)         
                  END DO 
               END IF
            ENDIF

            IF( nn_dyn3d_dta(ib_bdy) .eq. 0 ) THEN 
               ilen1(:) = nblen(:)
               IF( dta%ll_u3d ) THEN 
                  igrd = 2 
                  DO ib = 1, ilen1(igrd)
                     DO ik = 1, jpkm1
                        ii = idx_bdy(ib_bdy)%nbi(ib,igrd)
                        ij = idx_bdy(ib_bdy)%nbj(ib,igrd)
                        dta_bdy(ib_bdy)%u3d(ib,ik) =  ( un(ii,ij,ik) - un_b(ii,ij) ) * umask(ii,ij,ik)         
                     END DO
                  END DO 
               END IF
               IF( dta%ll_v3d ) THEN 
                  igrd = 3 
                  DO ib = 1, ilen1(igrd)
                     DO ik = 1, jpkm1
                        ii = idx_bdy(ib_bdy)%nbi(ib,igrd)
                        ij = idx_bdy(ib_bdy)%nbj(ib,igrd)
                        dta_bdy(ib_bdy)%v3d(ib,ik) =  ( vn(ii,ij,ik) - vn_b(ii,ij) ) * vmask(ii,ij,ik)         
                        END DO
                  END DO 
               END IF
            ENDIF

            IF( nn_tra_dta(ib_bdy) .eq. 0 ) THEN 
               ilen1(:) = nblen(:)
               IF( dta%ll_tem ) THEN
                  igrd = 1 
                  DO ib = 1, ilen1(igrd)
                     DO ik = 1, jpkm1
                        ii = idx_bdy(ib_bdy)%nbi(ib,igrd)
                        ij = idx_bdy(ib_bdy)%nbj(ib,igrd)
                        dta_bdy(ib_bdy)%tem(ib,ik) = tsn(ii,ij,ik,jp_tem) * tmask(ii,ij,ik)         
                     END DO
                  END DO 
               END IF
               IF( dta%ll_sal ) THEN
                  igrd = 1 
                  DO ib = 1, ilen1(igrd)
                     DO ik = 1, jpkm1
                        ii = idx_bdy(ib_bdy)%nbi(ib,igrd)
                        ij = idx_bdy(ib_bdy)%nbj(ib,igrd)
                        dta_bdy(ib_bdy)%sal(ib,ik) = tsn(ii,ij,ik,jp_sal) * tmask(ii,ij,ik)         
                     END DO
                  END DO 
               END IF
            ENDIF


            IF( nn_ice_lim_dta(ib_bdy) .eq. 0 ) THEN 
               ilen1(:) = nblen(:)
               IF( dta%ll_a_i ) THEN
                  igrd = 1   
                  DO jl = 1, jpl
                     DO ib = 1, ilen1(igrd)
                        ii = idx_bdy(ib_bdy)%nbi(ib,igrd)
                        ij = idx_bdy(ib_bdy)%nbj(ib,igrd)
                        dta_bdy(ib_bdy)%a_i (ib,jl) =  a_i(ii,ij,jl) * tmask(ii,ij,1) 
                     END DO
                  END DO
               ENDIF
               IF( dta%ll_ht_i ) THEN
                  igrd = 1   
                  DO jl = 1, jpl
                     DO ib = 1, ilen1(igrd)
                        ii = idx_bdy(ib_bdy)%nbi(ib,igrd)
                        ij = idx_bdy(ib_bdy)%nbj(ib,igrd)
                        dta_bdy(ib_bdy)%ht_i (ib,jl) =  ht_i(ii,ij,jl) * tmask(ii,ij,1) 
                     END DO
                  END DO
               ENDIF
               IF( dta%ll_ht_s ) THEN
                  igrd = 1   
                  DO jl = 1, jpl
                     DO ib = 1, ilen1(igrd)
                        ii = idx_bdy(ib_bdy)%nbi(ib,igrd)
                        ij = idx_bdy(ib_bdy)%nbj(ib,igrd)
                        dta_bdy(ib_bdy)%ht_s (ib,jl) =  ht_s(ii,ij,jl) * tmask(ii,ij,1) 
                     END DO
                  END DO
               ENDIF
            ENDIF


         ENDDO ! ib_bdy


      ENDIF ! kt .eq. nit000

! update external data from files
!--------------------------------
     
      jstart = 1
      DO ib_bdy = 1, nb_bdy   
         dta => dta_bdy(ib_bdy)
         IF( nn_dta(ib_bdy) .eq. 1 ) THEN ! skip this bit if no external data required
      
            IF( PRESENT(jit) ) THEN
! Update barotropic boundary conditions only
! jit is optional argument for fld_read and bdytide_update
               IF( cn_dyn2d(ib_bdy) /= 'none' ) THEN
                  IF( nn_dyn2d_dta(ib_bdy) .eq. 2 ) THEN ! tidal harmonic forcing ONLY: initialise arrays
                     IF( dta%ll_ssh ) dta%ssh(:) = 0.0
                     IF( dta%ll_u2d ) dta%u2d(:) = 0.0
                     IF( dta%ll_u3d ) dta%v2d(:) = 0.0
                  ENDIF
                  IF (cn_tra(ib_bdy) /= 'runoff') THEN
                     IF( nn_dyn2d_dta(ib_bdy) .EQ. 1 .OR. nn_dyn2d_dta(ib_bdy) .EQ. 3 ) THEN

                        jend = jstart + dta%nread(2) - 1
                        CALL fld_read( kt=kt, kn_fsbc=1, sd=bf(jstart:jend), map=nbmap_ptr(jstart:jend),  &
                                     & kit=jit, kt_offset=time_offset )

! If full velocities in boundary data then extract barotropic velocities from 3D fields
                        IF( ln_full_vel_array(ib_bdy) .AND.                                             &
                          &    ( nn_dyn2d_dta(ib_bdy) .EQ. 1 .OR. nn_dyn2d_dta(ib_bdy) .EQ. 3 .OR.  &
                          &      nn_dyn3d_dta(ib_bdy) .EQ. 1 ) )THEN

                           igrd = 2                      ! zonal velocity
                           dta%u2d(:) = 0.0
                           DO ib = 1, idx_bdy(ib_bdy)%nblen(igrd)
                              ii   = idx_bdy(ib_bdy)%nbi(ib,igrd)
                              ij   = idx_bdy(ib_bdy)%nbj(ib,igrd)
                              DO ik = 1, jpkm1
                                 dta%u2d(ib) = dta%u2d(ib) &
                       &                          + e3u_n(ii,ij,ik) * umask(ii,ij,ik) * dta%u3d(ib,ik)
                              END DO
                              dta%u2d(ib) =  dta%u2d(ib) * hur(ii,ij)
                           END DO
                           igrd = 3                      ! meridional velocity
                           dta%v2d(:) = 0.0
                           DO ib = 1, idx_bdy(ib_bdy)%nblen(igrd)
                              ii   = idx_bdy(ib_bdy)%nbi(ib,igrd)
                              ij   = idx_bdy(ib_bdy)%nbj(ib,igrd)
                              DO ik = 1, jpkm1
                                 dta%v2d(ib) = dta%v2d(ib) &
                       &                       + e3v_n(ii,ij,ik) * vmask(ii,ij,ik) * dta%v3d(ib,ik)
                              END DO
                              dta%v2d(ib) =  dta%v2d(ib) * hvr(ii,ij)
                           END DO
                        ENDIF                    
                     ENDIF
                     IF( nn_dyn2d_dta(ib_bdy) .ge. 2 ) THEN ! update tidal harmonic forcing
                        CALL bdytide_update( kt=kt, idx=idx_bdy(ib_bdy), dta=dta, td=tides(ib_bdy),   & 
                          &                 jit=jit, time_offset=time_offset )
                     ENDIF
                  ENDIF
               ENDIF
            ELSE
               IF (cn_tra(ib_bdy) == 'runoff') then      ! runoff condition
                  jend = nb_bdy_fld(ib_bdy)
                  CALL fld_read( kt=kt, kn_fsbc=1, sd=bf(jstart:jend),  &
                               & map=nbmap_ptr(jstart:jend), kt_offset=time_offset )
!
                  igrd = 2                      ! zonal velocity
                  DO ib = 1, idx_bdy(ib_bdy)%nblen(igrd)
                     ii   = idx_bdy(ib_bdy)%nbi(ib,igrd)
                     ij   = idx_bdy(ib_bdy)%nbj(ib,igrd)
                     dta%u2d(ib) = dta%u2d(ib) / ( e2u(ii,ij) * hu_0(ii,ij) )
                  END DO
!
                  igrd = 3                      ! meridional velocity
                  DO ib = 1, idx_bdy(ib_bdy)%nblen(igrd)
                     ii   = idx_bdy(ib_bdy)%nbi(ib,igrd)
                     ij   = idx_bdy(ib_bdy)%nbj(ib,igrd)
                     dta%v2d(ib) = dta%v2d(ib) / ( e1v(ii,ij) * hv_0(ii,ij) )
                  END DO
               ELSE
                  IF( nn_dyn2d_dta(ib_bdy) .eq. 2 ) THEN ! tidal harmonic forcing ONLY: initialise arrays
                     IF( dta%ll_ssh ) dta%ssh(:) = 0.0
                     IF( dta%ll_u2d ) dta%u2d(:) = 0.0
                     IF( dta%ll_v2d ) dta%v2d(:) = 0.0
                  ENDIF
                  IF( dta%nread(1) .gt. 0 ) THEN ! update external data
                     jend = jstart + dta%nread(1) - 1
                     CALL fld_read( kt=kt, kn_fsbc=1, sd=bf(jstart:jend), &
                                  & map=nbmap_ptr(jstart:jend), kt_offset=time_offset )
                  ENDIF
! If full velocities in boundary data then split into barotropic and baroclinic data
                  IF( ln_full_vel_array(ib_bdy) .and.                                             &
                    & ( nn_dyn2d_dta(ib_bdy) .EQ. 1 .OR. nn_dyn2d_dta(ib_bdy) .EQ. 3 .OR. &
                    &   nn_dyn3d_dta(ib_bdy) .EQ. 1 ) ) THEN
                     igrd = 2                      ! zonal velocity
                     dta%u2d(:) = 0.0
                     DO ib = 1, idx_bdy(ib_bdy)%nblen(igrd)
                        ii   = idx_bdy(ib_bdy)%nbi(ib,igrd)
                        ij   = idx_bdy(ib_bdy)%nbj(ib,igrd)
                        DO ik = 1, jpkm1
                           dta%u2d(ib) = dta%u2d(ib) &
                 &                       + e3u_n(ii,ij,ik) * umask(ii,ij,ik) * dta%u3d(ib,ik)
                        END DO
                        dta%u2d(ib) =  dta%u2d(ib) * hur(ii,ij)
                        DO ik = 1, jpkm1
                           dta%u3d(ib,ik) = dta%u3d(ib,ik) - dta%u2d(ib)
                        END DO
                     END DO
                     igrd = 3                      ! meridional velocity
                     dta%v2d(:) = 0.0
                     DO ib = 1, idx_bdy(ib_bdy)%nblen(igrd)
                        ii   = idx_bdy(ib_bdy)%nbi(ib,igrd)
                        ij   = idx_bdy(ib_bdy)%nbj(ib,igrd)
                        DO ik = 1, jpkm1
                           dta%v2d(ib) = dta%v2d(ib) &
                 &                       + e3v_n(ii,ij,ik) * vmask(ii,ij,ik) * dta%v3d(ib,ik)
                        END DO
                        dta%v2d(ib) =  dta%v2d(ib) * hvr(ii,ij)
                        DO ik = 1, jpkm1
                           dta%v3d(ib,ik) = dta%v3d(ib,ik) - dta%v2d(ib)
                        END DO
                     END DO
                  ENDIF

               ENDIF

               IF( .NOT. ll_bdylim3 .AND. cn_ice_lim(ib_bdy) /= 'none' .AND. nn_ice_lim_dta(ib_bdy) == 1 ) THEN ! bdy ice input (case input is lim2 type)
                  jfld_hti = jfld_htit(ib_bdy)
                  jfld_hts = jfld_htst(ib_bdy)
                  jfld_ai  = jfld_ait(ib_bdy)
     	          CALL lim_var_itd ( bf(jfld_hti)%fnow(:,1,1), bf(jfld_hts)%fnow(:,1,1), bf(jfld_ai)%fnow(:,1,1), &
                                  & dta_bdy(ib_bdy)%ht_i,     dta_bdy(ib_bdy)%ht_s,     dta_bdy(ib_bdy)%a_i     )
               ENDIF

            ENDIF
            jstart = jstart + dta%nread(1)
         END IF ! nn_dta(ib_bdy) = 1
      END DO  ! ib_bdy

! bg jchanut tschanges

! Add tides if not split-explicit free surface else this is done in ts loop
      IF (.NOT.lk_dynspg_ts) CALL bdy_dta_tides( kt=kt, time_offset=time_offset )

! end jchanut tschanges

      IF( ln_apr_dyn )THEN
         IF( ln_apr_obc ) THEN
            DO ib_bdy = 1, nb_bdy
               IF (cn_tra(ib_bdy) /= 'runoff')THEN
                  igrd = 1                      ! meridional velocity
                  DO ib = 1, idx_bdy(ib_bdy)%nblenrim(igrd)
                     ii   = idx_bdy(ib_bdy)%nbi(ib,igrd)
                     ij   = idx_bdy(ib_bdy)%nbj(ib,igrd)
                     dta_bdy(ib_bdy)%ssh(ib) = dta_bdy(ib_bdy)%ssh(ib) + ssh_ib(ii,ij)
                  ENDDO
               ENDIF
            ENDDO
         ENDIF
      ENDIF

      IF( nn_timing == 1 ) CALL timing_stop('bdy_dta')

      END SUBROUTINE bdy_dta


      SUBROUTINE bdy_dta_init
!!----------------------------------------------------------------------
!!                   ***  SUBROUTINE bdy_dta_init  ***
!!
!! ** Purpose :   Initialise arrays for reading of external data
!!                for open boundary conditions
!!
!! ** Method  :
!!
!!----------------------------------------------------------------------
      USE dynspg_oce, ONLY: lk_dynspg_ts
!!
      INTEGER     ::  ib_bdy, jfld, jstart, jend, ierror  ! local indices
      INTEGER      ::   ios                               ! Local integer output status for namelist read
!!
      CHARACTER(len=100)                     ::   cn_dir        ! Root directory for location of data files
      CHARACTER(len=100), DIMENSION(nb_bdy)  ::   cn_dir_array  ! Root directory for location of data files
      CHARACTER(len = 256)::   clname                           ! temporary file name
      LOGICAL                                ::   ln_full_vel   ! =T => full velocities in 3D boundary data
! =F => baroclinic velocities in 3D boundary data
      INTEGER                                ::   ilen_global   ! Max length required for global bdy dta arrays
      INTEGER, ALLOCATABLE, DIMENSION(:)     ::   ilen1, ilen3  ! size of 1st and 3rd dimensions of local arrays
      INTEGER, ALLOCATABLE, DIMENSION(:)     ::   ibdy           ! bdy set for a particular jfld
      INTEGER, ALLOCATABLE, DIMENSION(:)     ::   igrid         ! index for grid type (1,2,3 = T,U,V)
      INTEGER, POINTER, DIMENSION(:)         ::   nblen, nblenrim  ! short cuts
      TYPE(OBC_DATA), POINTER                ::   dta           ! short cut

      INTEGER               ::   zndims   ! number of dimensions in an array (i.e. 3 = wo ice cat; 4 = w ice cat)
      INTEGER               ::   inum,id1 ! local integer

      TYPE(FLD_N), ALLOCATABLE, DIMENSION(:) ::   blf_i         !  array of namelist information structures
      TYPE(FLD_N) ::   bn_tem, bn_sal, bn_u3d, bn_v3d   !
      TYPE(FLD_N) ::   bn_ssh, bn_u2d, bn_v2d           ! informations about the fields to be read

      TYPE(FLD_N) ::   bn_a_i, bn_ht_i, bn_ht_s      

      NAMELIST/nambdy_dta/ cn_dir, bn_tem, bn_sal, bn_u3d, bn_v3d, bn_ssh, bn_u2d, bn_v2d 

      NAMELIST/nambdy_dta/ bn_a_i, bn_ht_i, bn_ht_s

      NAMELIST/nambdy_dta/ ln_full_vel
!!---------------------------------------------------------------------------

      IF( nn_timing == 1 ) CALL timing_start('bdy_dta_init')

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'bdy_dta_ini : initialization of data at the open boundaries'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~'
      IF(lwp) WRITE(numout,*) ''

! Set nn_dta
      DO ib_bdy = 1, nb_bdy
         nn_dta(ib_bdy) = MAX(  nn_dyn2d_dta(ib_bdy)       &
                               ,nn_dyn3d_dta(ib_bdy)       &
                               ,nn_tra_dta(ib_bdy)         &

                              ,nn_ice_lim_dta(ib_bdy)    &

                              )
         IF(nn_dta(ib_bdy) .gt. 1) nn_dta(ib_bdy) = 1
      END DO

! Work out upper bound of how many fields there are to read in and allocate arrays
! ---------------------------------------------------------------------------
      ALLOCATE( nb_bdy_fld(nb_bdy) )
      nb_bdy_fld(:) = 0
      DO ib_bdy = 1, nb_bdy         
         IF( cn_dyn2d(ib_bdy) /= 'none' .and. ( nn_dyn2d_dta(ib_bdy) .eq. 1 .or. nn_dyn2d_dta(ib_bdy) .eq. 3 ) ) THEN
            nb_bdy_fld(ib_bdy) = nb_bdy_fld(ib_bdy) + 3
         ENDIF
         IF( cn_dyn3d(ib_bdy) /= 'none' .and. nn_dyn3d_dta(ib_bdy) .eq. 1 ) THEN
            nb_bdy_fld(ib_bdy) = nb_bdy_fld(ib_bdy) + 2
         ENDIF
         IF( cn_tra(ib_bdy) /= 'none' .and. nn_tra_dta(ib_bdy) .eq. 1  ) THEN
            nb_bdy_fld(ib_bdy) = nb_bdy_fld(ib_bdy) + 2
         ENDIF

         IF( cn_ice_lim(ib_bdy) /= 'none' .and. nn_ice_lim_dta(ib_bdy) .eq. 1  ) THEN
            nb_bdy_fld(ib_bdy) = nb_bdy_fld(ib_bdy) + 3
         ENDIF

         IF(lwp) WRITE(numout,*) 'Maximum number of files to open =',nb_bdy_fld(ib_bdy)
      ENDDO            

      nb_bdy_fld_sum = SUM( nb_bdy_fld )

      ALLOCATE( bf(nb_bdy_fld_sum), STAT=ierror )
      IF( ierror > 0 ) THEN   
         CALL ctl_stop( 'bdy_dta: unable to allocate bf structure' )   ;   RETURN  
      ENDIF
      ALLOCATE( blf_i(nb_bdy_fld_sum), STAT=ierror )
      IF( ierror > 0 ) THEN   
         CALL ctl_stop( 'bdy_dta: unable to allocate blf_i structure' )   ;   RETURN  
      ENDIF
      ALLOCATE( nbmap_ptr(nb_bdy_fld_sum), STAT=ierror )
      IF( ierror > 0 ) THEN   
         CALL ctl_stop( 'bdy_dta: unable to allocate nbmap_ptr structure' )   ;   RETURN  
      ENDIF
      ALLOCATE( ilen1(nb_bdy_fld_sum), ilen3(nb_bdy_fld_sum) ) 
      ALLOCATE( ibdy(nb_bdy_fld_sum) ) 
      ALLOCATE( igrid(nb_bdy_fld_sum) ) 

! Read namelists
! --------------
      REWIND(numnam_ref)
      REWIND(numnam_cfg)
      jfld = 0 
      DO ib_bdy = 1, nb_bdy         
         IF( nn_dta(ib_bdy) .eq. 1 ) THEN
            READ  ( numnam_ref, nambdy_dta, IOSTAT = ios, ERR = 901)
901         IF( ios /= 0 ) CALL ctl_nam ( ios , 'nambdy_dta in reference namelist', lwp )

            READ  ( numnam_cfg, nambdy_dta, IOSTAT = ios, ERR = 902 )
902         IF( ios /= 0 ) CALL ctl_nam ( ios , 'nambdy_dta in configuration namelist', lwp )
            IF(lwm) WRITE ( numond, nambdy_dta )

            cn_dir_array(ib_bdy) = cn_dir
            ln_full_vel_array(ib_bdy) = ln_full_vel

            nblen => idx_bdy(ib_bdy)%nblen
            nblenrim => idx_bdy(ib_bdy)%nblenrim
            dta => dta_bdy(ib_bdy)
            dta%nread(2) = 0

! Only read in necessary fields for this set.
! Important that barotropic variables come first.
            IF( nn_dyn2d_dta(ib_bdy) .eq. 1 .or. nn_dyn2d_dta(ib_bdy) .eq. 3 ) THEN 

               IF( dta%ll_ssh ) THEN 
                  if(lwp) write(numout,*) '++++++ reading in ssh field'
                  jfld = jfld + 1
                  blf_i(jfld) = bn_ssh
                  ibdy(jfld) = ib_bdy
                  igrid(jfld) = 1
                  ilen1(jfld) = nblen(igrid(jfld))
                  ilen3(jfld) = 1
                  dta%nread(2) = dta%nread(2) + 1
               ENDIF

               IF( dta%ll_u2d .and. .not. ln_full_vel_array(ib_bdy) ) THEN
                  if(lwp) write(numout,*) '++++++ reading in u2d field'
                  jfld = jfld + 1
                  blf_i(jfld) = bn_u2d
                  ibdy(jfld) = ib_bdy
                  igrid(jfld) = 2
                  ilen1(jfld) = nblen(igrid(jfld))
                  ilen3(jfld) = 1
                  dta%nread(2) = dta%nread(2) + 1
               ENDIF

               IF( dta%ll_v2d .and. .not. ln_full_vel_array(ib_bdy) ) THEN
                  if(lwp) write(numout,*) '++++++ reading in v2d field'
                  jfld = jfld + 1
                  blf_i(jfld) = bn_v2d
                  ibdy(jfld) = ib_bdy
                  igrid(jfld) = 3
                  ilen1(jfld) = nblen(igrid(jfld))
                  ilen3(jfld) = 1
                  dta%nread(2) = dta%nread(2) + 1
               ENDIF

            ENDIF

! read 3D velocities if baroclinic velocities require OR if
! barotropic velocities required and ln_full_vel set to .true.
            IF( nn_dyn3d_dta(ib_bdy) .eq. 1 .or. &
           &  ( ln_full_vel_array(ib_bdy) .and. ( nn_dyn2d_dta(ib_bdy) .eq. 1 .or. nn_dyn2d_dta(ib_bdy) .eq. 3 ) ) ) THEN

               IF( dta%ll_u3d .or. ( ln_full_vel_array(ib_bdy) .and. dta%ll_u2d ) ) THEN
                  if(lwp) write(numout,*) '++++++ reading in u3d field'
                  jfld = jfld + 1
                  blf_i(jfld) = bn_u3d
                  ibdy(jfld) = ib_bdy
                  igrid(jfld) = 2
                  ilen1(jfld) = nblen(igrid(jfld))
                  ilen3(jfld) = jpk
                  IF( ln_full_vel_array(ib_bdy) .and. dta%ll_u2d ) dta%nread(2) = dta%nread(2) + 1
               ENDIF

               IF( dta%ll_v3d .or. ( ln_full_vel_array(ib_bdy) .and. dta%ll_v2d ) ) THEN
                  if(lwp) write(numout,*) '++++++ reading in v3d field'
                  jfld = jfld + 1
                  blf_i(jfld) = bn_v3d
                  ibdy(jfld) = ib_bdy
                  igrid(jfld) = 3
                  ilen1(jfld) = nblen(igrid(jfld))
                  ilen3(jfld) = jpk
                  IF( ln_full_vel_array(ib_bdy) .and. dta%ll_v2d ) dta%nread(2) = dta%nread(2) + 1
               ENDIF

            ENDIF

! temperature and salinity
            IF( nn_tra_dta(ib_bdy) .eq. 1 ) THEN

               IF( dta%ll_tem ) THEN
                  if(lwp) write(numout,*) '++++++ reading in tem field'
                  jfld = jfld + 1
                  blf_i(jfld) = bn_tem
                  ibdy(jfld) = ib_bdy
                  igrid(jfld) = 1
                  ilen1(jfld) = nblen(igrid(jfld))
                  ilen3(jfld) = jpk
               ENDIF

               IF( dta%ll_sal ) THEN
                  if(lwp) write(numout,*) '++++++ reading in sal field'
                  jfld = jfld + 1
                  blf_i(jfld) = bn_sal
                  ibdy(jfld) = ib_bdy
                  igrid(jfld) = 1
                  ilen1(jfld) = nblen(igrid(jfld))
                  ilen3(jfld) = jpk
               ENDIF

            ENDIF


! sea ice
            IF( nn_ice_lim_dta(ib_bdy) .eq. 1 ) THEN
! Test for types of ice input (lim2 or lim3)
! Build file name to find dimensions
               clname=TRIM( cn_dir )//TRIM(bn_a_i%clname)
               IF( .NOT. bn_a_i%ln_clim ) THEN   
                                                  WRITE(clname, '(a,"_y",i4.4)' ) TRIM( clname ), nyear    ! add year
                  IF( bn_a_i%cltype /= 'yearly' ) WRITE(clname, '(a,"m" ,i2.2)' ) TRIM( clname ), nmonth   ! add month
               ELSE
                  IF( bn_a_i%cltype /= 'yearly' ) WRITE(clname, '(a,"_m",i2.2)' ) TRIM( clname ), nmonth   ! add month
               ENDIF
               IF( bn_a_i%cltype == 'daily' .OR. bn_a_i%cltype(1:4) == 'week' ) &
               &                                  WRITE(clname, '(a,"d" ,i2.2)' ) TRIM( clname ), nday     ! add day
!
               CALL iom_open  ( clname, inum )
               id1 = iom_varid( inum, bn_a_i%clvar, kndims=zndims, ldstop = .FALSE. )
               CALL iom_close ( inum )

      	       IF ( zndims == 4 ) THEN
                 ll_bdylim3 = .TRUE.   ! lim3 input
               ELSE
                 ll_bdylim3 = .FALSE.  ! lim2 input
               ENDIF
! End test

               IF( dta%ll_a_i ) THEN
                  jfld = jfld + 1
                  blf_i(jfld) = bn_a_i
                  ibdy(jfld) = ib_bdy
                  igrid(jfld) = 1
                  ilen1(jfld) = nblen(igrid(jfld))
                  IF ( ll_bdylim3 ) THEN ; ilen3(jfld)=jpl ; ELSE ; ilen3(jfld)=1 ; ENDIF
               ENDIF

               IF( dta%ll_ht_i ) THEN
                  jfld = jfld + 1
                  blf_i(jfld) = bn_ht_i
                  ibdy(jfld) = ib_bdy
                  igrid(jfld) = 1
                  ilen1(jfld) = nblen(igrid(jfld))
                  IF ( ll_bdylim3 ) THEN ; ilen3(jfld)=jpl ; ELSE ; ilen3(jfld)=1 ; ENDIF
               ENDIF

               IF( dta%ll_ht_s ) THEN
                  jfld = jfld + 1
                   blf_i(jfld) = bn_ht_s
                  ibdy(jfld) = ib_bdy
                  igrid(jfld) = 1
                  ilen1(jfld) = nblen(igrid(jfld))
                  IF ( ll_bdylim3 ) THEN ; ilen3(jfld)=jpl ; ELSE ; ilen3(jfld)=1 ; ENDIF
               ENDIF

            ENDIF

! Recalculate field counts
!-------------------------
            IF( ib_bdy .eq. 1 ) THEN 
               nb_bdy_fld_sum = 0
               nb_bdy_fld(ib_bdy) = jfld
               nb_bdy_fld_sum     = jfld              
            ELSE
               nb_bdy_fld(ib_bdy) = jfld - nb_bdy_fld_sum
               nb_bdy_fld_sum = nb_bdy_fld_sum + nb_bdy_fld(ib_bdy)
            ENDIF

            dta%nread(1) = nb_bdy_fld(ib_bdy)

         ENDIF ! nn_dta .eq. 1
      ENDDO ! ib_bdy

      DO jfld = 1, nb_bdy_fld_sum
         ALLOCATE( bf(jfld)%fnow(ilen1(jfld),1,ilen3(jfld)) )
         IF( blf_i(jfld)%ln_tint ) ALLOCATE( bf(jfld)%fdta(ilen1(jfld),1,ilen3(jfld),2) )
         nbmap_ptr(jfld)%ptr => idx_bdy(ibdy(jfld))%nbmap(:,igrid(jfld))
         nbmap_ptr(jfld)%ll_unstruc = ln_coords_file(ibdy(jfld))
      ENDDO

! fill bf with blf_i and control print
!-------------------------------------
      jstart = 1
      DO ib_bdy = 1, nb_bdy
         jend = jstart - 1 + nb_bdy_fld(ib_bdy) 
         CALL fld_fill( bf(jstart:jend), blf_i(jstart:jend), cn_dir_array(ib_bdy), 'bdy_dta',   &
         &              'open boundary conditions', 'nambdy_dta' )
         jstart = jend + 1
      ENDDO

! Initialise local boundary data arrays
! nn_xxx_dta=0 : allocate space - will be filled from initial conditions later
! nn_xxx_dta=1 : point to "fnow" arrays
!-------------------------------------

      jfld = 0
      DO ib_bdy=1, nb_bdy

         nblen => idx_bdy(ib_bdy)%nblen
         dta => dta_bdy(ib_bdy)

         if(lwp) then
            write(numout,*) '++++++ dta%ll_ssh = ',dta%ll_ssh
            write(numout,*) '++++++ dta%ll_u2d = ',dta%ll_u2d
            write(numout,*) '++++++ dta%ll_v2d = ',dta%ll_v2d
            write(numout,*) '++++++ dta%ll_u3d = ',dta%ll_u3d
            write(numout,*) '++++++ dta%ll_v3d = ',dta%ll_v3d
            write(numout,*) '++++++ dta%ll_tem = ',dta%ll_tem
            write(numout,*) '++++++ dta%ll_sal = ',dta%ll_sal
         endif

         IF ( nn_dyn2d_dta(ib_bdy) .eq. 0 .or. nn_dyn2d_dta(ib_bdy) .eq. 2 ) THEN
            if(lwp) write(numout,*) '++++++ dta%ssh/u2d/u3d allocated space'
            IF( dta%ll_ssh ) ALLOCATE( dta%ssh(nblen(1)) )
            IF( dta%ll_u2d ) ALLOCATE( dta%u2d(nblen(2)) )
            IF( dta%ll_v2d ) ALLOCATE( dta%v2d(nblen(3)) )
         ENDIF
         IF ( nn_dyn2d_dta(ib_bdy) .eq. 1 .or. nn_dyn2d_dta(ib_bdy) .eq. 3 ) THEN
            IF( dta%ll_ssh ) THEN
               if(lwp) write(numout,*) '++++++ dta%ssh pointing to fnow'
               jfld = jfld + 1
               dta%ssh => bf(jfld)%fnow(:,1,1)
            ENDIF
            IF ( dta%ll_u2d ) THEN
               IF ( ln_full_vel_array(ib_bdy) ) THEN
                  if(lwp) write(numout,*) '++++++ dta%u2d allocated space'
                  ALLOCATE( dta%u2d(nblen(2)) )
               ELSE
                  if(lwp) write(numout,*) '++++++ dta%u2d pointing to fnow'
                  jfld = jfld + 1
                  dta%u2d => bf(jfld)%fnow(:,1,1)
               ENDIF
            ENDIF
            IF ( dta%ll_v2d ) THEN
               IF ( ln_full_vel_array(ib_bdy) ) THEN
                  if(lwp) write(numout,*) '++++++ dta%v2d allocated space'
                  ALLOCATE( dta%v2d(nblen(3)) )
               ELSE
                  if(lwp) write(numout,*) '++++++ dta%v2d pointing to fnow'
                  jfld = jfld + 1
                  dta%v2d => bf(jfld)%fnow(:,1,1)
               ENDIF
            ENDIF
         ENDIF

         IF ( nn_dyn3d_dta(ib_bdy) .eq. 0 ) THEN
            if(lwp) write(numout,*) '++++++ dta%u3d/v3d allocated space'
            IF( dta%ll_u3d ) ALLOCATE( dta_bdy(ib_bdy)%u3d(nblen(2),jpk) )
            IF( dta%ll_v3d ) ALLOCATE( dta_bdy(ib_bdy)%v3d(nblen(3),jpk) )
         ENDIF
         IF ( nn_dyn3d_dta(ib_bdy) .eq. 1 .or. &
           &  ( ln_full_vel_array(ib_bdy) .and. ( nn_dyn2d_dta(ib_bdy) .eq. 1 .or. nn_dyn2d_dta(ib_bdy) .eq. 3 ) ) ) THEN
            IF ( dta%ll_u3d .or. ( ln_full_vel_array(ib_bdy) .and. dta%ll_u2d ) ) THEN
               if(lwp) write(numout,*) '++++++ dta%u3d pointing to fnow'
               jfld = jfld + 1
               dta_bdy(ib_bdy)%u3d => bf(jfld)%fnow(:,1,:)
            ENDIF
            IF ( dta%ll_v3d .or. ( ln_full_vel_array(ib_bdy) .and. dta%ll_v2d ) ) THEN
               if(lwp) write(numout,*) '++++++ dta%v3d pointing to fnow'
               jfld = jfld + 1
               dta_bdy(ib_bdy)%v3d => bf(jfld)%fnow(:,1,:)
            ENDIF
         ENDIF

         IF( nn_tra_dta(ib_bdy) .eq. 0 ) THEN
            if(lwp) write(numout,*) '++++++ dta%tem/sal allocated space'
            IF( dta%ll_tem ) ALLOCATE( dta_bdy(ib_bdy)%tem(nblen(1),jpk) )
            IF( dta%ll_sal ) ALLOCATE( dta_bdy(ib_bdy)%sal(nblen(1),jpk) )
         ELSE
            IF( dta%ll_tem ) THEN
               if(lwp) write(numout,*) '++++++ dta%tem pointing to fnow'
               jfld = jfld + 1
               dta_bdy(ib_bdy)%tem => bf(jfld)%fnow(:,1,:)
            ENDIF
            IF( dta%ll_sal ) THEN 
               if(lwp) write(numout,*) '++++++ dta%sal pointing to fnow'
               jfld = jfld + 1
               dta_bdy(ib_bdy)%sal => bf(jfld)%fnow(:,1,:)
            ENDIF
         ENDIF


         IF (cn_ice_lim(ib_bdy) /= 'none') THEN
            IF( nn_ice_lim_dta(ib_bdy) .eq. 0 ) THEN
               ALLOCATE( dta_bdy(ib_bdy)%a_i (nblen(1),jpl) )
               ALLOCATE( dta_bdy(ib_bdy)%ht_i(nblen(1),jpl) )
               ALLOCATE( dta_bdy(ib_bdy)%ht_s(nblen(1),jpl) )
            ELSE
               IF ( ll_bdylim3 ) THEN ! case input is lim3 type
                  jfld = jfld + 1
                  dta_bdy(ib_bdy)%a_i  => bf(jfld)%fnow(:,1,:)
                  jfld = jfld + 1
                  dta_bdy(ib_bdy)%ht_i => bf(jfld)%fnow(:,1,:)
                  jfld = jfld + 1
                  dta_bdy(ib_bdy)%ht_s => bf(jfld)%fnow(:,1,:)
               ELSE ! case input is lim2 type
                  jfld_ait(ib_bdy)  = jfld + 1
                  jfld_htit(ib_bdy) = jfld + 2
                  jfld_htst(ib_bdy) = jfld + 3
                  jfld     = jfld + 3
                  ALLOCATE( dta_bdy(ib_bdy)%a_i (nblen(1),jpl) )
                  ALLOCATE( dta_bdy(ib_bdy)%ht_i(nblen(1),jpl) )
                  ALLOCATE( dta_bdy(ib_bdy)%ht_s(nblen(1),jpl) )
                  dta_bdy(ib_bdy)%a_i (:,:) = 0.0
                  dta_bdy(ib_bdy)%ht_i(:,:) = 0.0
                  dta_bdy(ib_bdy)%ht_s(:,:) = 0.0
               ENDIF

            ENDIF
         ENDIF


      ENDDO ! ib_bdy

      IF( nn_timing == 1 ) CALL timing_stop('bdy_dta_init')

      END SUBROUTINE bdy_dta_init



!!==============================================================================
END MODULE bdydta
