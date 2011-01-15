module run_parameters
  implicit none

  public :: init_run_parameters, finish_run_parameters

  public :: beta, zeff, tite
  public :: fphi, fapar, fbpar
!  public :: delt, delt_max, wunits, woutunits, tunits, funits, tnorm
  public :: code_delt_max, wunits, woutunits, tunits, funits, tnorm
  public :: nstep, wstar_units, eqzip, margin
  public :: secondary, tertiary, harris
  public :: ieqzip
  public :: k0
  public :: vnm_init
  public :: avail_cpu_time

  private

  real :: beta, zeff, tite
  real :: fphi, fapar, fbpar, faperp
  real :: delt, code_delt_max, user_delt_max, funits, tnorm, margin
  real, dimension (:), allocatable :: wunits, woutunits, tunits
  real, dimension (2) :: vnm_init
  real :: avail_cpu_time
  integer :: nstep
  logical :: wstar_units, eqzip
  logical :: secondary, tertiary, harris
  real :: k0
  integer :: delt_option_switch
  integer, parameter :: delt_option_hand = 1, delt_option_auto = 2
  logical :: initialized = .false.

  integer, allocatable :: ieqzip(:,:)
  integer :: eqzip_option_switch
  integer, parameter :: &
       eqzip_option_none = 1, &
       eqzip_option_secondary = 2, &
       eqzip_option_tertiary = 3, &
       eqzip_option_equilibrium = 4

contains

  subroutine init_run_parameters
    use kt_grids, only: init_kt_grids, naky, nakx => ntheta0
    use gs2_time, only: init_delt, user2code
    
    implicit none
!    logical, save :: initialized = .false.

    if (initialized) return
    initialized = .true.

    call read_parameters

    call init_kt_grids (tnorm)
    call init_delt (delt, tnorm)
    call user2code (user_delt_max, code_delt_max)

!    delt = delt * tnorm

    allocate (wunits(naky))
    allocate (woutunits(naky))
    allocate (tunits(naky))

! omega_* normalization of time: 
    call adjust_time_norm

    if(.not.allocated(ieqzip)) allocate(ieqzip(nakx,naky))
    ieqzip(1:nakx,1:naky)=1
    select case (eqzip_option_switch)
    case (eqzip_option_secondary)
       ! suppress evolution of secondary mode
       ieqzip(1,2) = 0
    case (eqzip_option_tertiary)
       ! suppress evolution of tertiary mode
       ieqzip(2,1) = 0
       ieqzip(nakx,1) = 0
    case (eqzip_option_equilibrium)
       ! suppress evolution of 1D equilibrium (x dependent)
       ieqzip(1:nakx,1) = 0
    end select
  end subroutine init_run_parameters

  subroutine read_parameters
    use file_utils, only: input_unit, error_unit, input_unit_exist
    use mp, only: proc0, broadcast
    use gs2_save, only: init_dt, init_vnm
    use text_options, only: text_option, get_option_value
    implicit none
    type (text_option), dimension (4), parameter :: eqzipopts = &
         (/ text_option('none', eqzip_option_none), &
            text_option('secondary', eqzip_option_secondary), &
            text_option('tertiary', eqzip_option_tertiary), &
            text_option('equilibrium', eqzip_option_equilibrium) /)
    character (len=20) :: eqzip_option
    type (text_option), dimension (3), parameter :: deltopts = &
         (/ text_option('default', delt_option_hand), &
            text_option('set_by_hand', delt_option_hand), &
            text_option('check_restart', delt_option_auto) /)
    character(20) :: delt_option
    integer :: ierr, istatus, in_file
    real :: delt_saved
    real, dimension (2) :: vnm_saved

    real :: teti  ! for back-compatibility
    logical :: exist
    namelist /parameters/ beta, zeff, tite, teti, k0
    namelist /knobs/ fphi, fapar, fbpar, delt, nstep, wstar_units, eqzip, &
         delt_option, margin, secondary, tertiary, faperp, harris, &
         avail_cpu_time, eqzip_option

    if (proc0) then
       fbpar = -1.0
       faperp = 0.0
       beta = 0.0
       zeff = 1.0
       tite = 1.0
       teti = -100.0
       wstar_units = .false.
       eqzip_option = 'none'
       eqzip = .false.
       secondary = .true.
       tertiary = .false.
       harris = .false.
       k0 = 1.
       delt_option = 'default'
       margin = 0.05
       avail_cpu_time = 1.e10

       in_file = input_unit_exist("parameters", exist)
!       if (exist) read (unit=input_unit("parameters"), nml=parameters)
       if (exist) read (unit=in_file,nml=parameters)

       in_file = input_unit_exist("knobs", exist)
!       if (exist) read (unit=input_unit("knobs"), nml=knobs)
       if (exist) read (unit=in_file, nml=knobs)

       if (teti /= -100.0) tite = teti

! Allow faperp-style initialization for backwards compatibility.
! Only fbpar is used outside of this subroutine.
       if (fbpar == -1.) then
          fbpar = faperp
       end if

       if (eqzip) then
          if (secondary .and. tertiary) then
             ierr = error_unit()
             write (ierr, *) 'Forcing secondary = FALSE'
             write (ierr, *) 'because you have chosen tertiary = TRUE'
             secondary = .false.
          end if
          if (secondary .and. harris) then
             ierr = error_unit()
             write (ierr, *) 'Forcing secondary = FALSE'
             write (ierr, *) 'because you have chosen harris = TRUE'
             secondary = .false.
          end if
          if (tertiary .and. harris) then
             ierr = error_unit()
             write (ierr, *) 'Forcing tertiary = FALSE'
             write (ierr, *) 'because you have chosen harris = TRUE'
             tertiary = .false.
          end if
       end if

       ierr = error_unit()
       call get_option_value &
            (delt_option, deltopts, delt_option_switch, ierr, &
            "delt_option in knobs")

       call get_option_value ( &
            eqzip_option, eqzipopts, eqzip_option_switch, error_unit(), &
            "eqzip_option in knobs")

!!$       ! eqzip_option replaces eqzip, secondary, tertiary, harris
!!$       if (eqzip .and. eqzip_option_switch == eqzip_option_none) then
!!$          if (harris) then
!!$             eqzip_option_switch = eqzip_option_equilibrium
!!$             write(error_unit(),*) 'eqzip_option is set to equilibrium'
!!$          else
!!$             if (tertiary) then
!!$                eqzip_option_switch = eqzip_option_tertiary
!!$                write(error_unit(),*) 'eqzip_option is set to tertiary'
!!$             else
!!$                eqzip_option_switch = eqzip_option_secondary
!!$                write(error_unit(),*) 'eqzip_option is set to secondary'
!!$             end if
!!$          end if
!!$       else if (eqzip .and. eqzip_option_switch /= eqzip_option_none) then
!!$          write(error_unit(),*) 'eqzip, secondary, tertiary, harris are ignored'
!!$          write(error_unit(),*) 'because eqzip_option exists'
!!$       end if

    end if

    call broadcast (delt_option_switch)
    call broadcast (delt)
    call broadcast (beta)
    call broadcast (zeff)
    call broadcast (tite)
    call broadcast (fphi)
    call broadcast (fapar)
    call broadcast (fbpar)
    call broadcast (nstep)
    call broadcast (wstar_units)
    call broadcast (eqzip)
    call broadcast (secondary)
    call broadcast (tertiary)
    call broadcast (harris)
    call broadcast (margin)
    call broadcast (k0)
    call broadcast (avail_cpu_time)
    call broadcast (eqzip_option_switch)
    
    user_delt_max = delt

    delt_saved = delt
    if (delt_option_switch == delt_option_auto) then
       vnm_init = 1.0
       call init_vnm (vnm_saved, istatus)
       if (istatus == 0) vnm_init = vnm_saved
       call init_dt (delt_saved, istatus)
       if (istatus == 0) delt  = delt_saved
    endif

  end subroutine read_parameters

  subroutine adjust_time_norm
    use file_utils, only: error_unit
    use mp, only: proc0
    use kt_grids, only: aky
    implicit none
!CMR: Sep 2010
! Attempt to understand time normalisation variables, which are arrays(naky)
!    TUNITS: DT(KY)=TUNITS(KY).CODE_DT
!            This is a generally very useful variable to store ky dependent 
!            timestep in the code time normalisation.
!            Used to multiply ky independent source terms on RHS of GKE.
!    WUNITS: WUNITS(KY)=AKY(KY)*TUNITS(KY)/2
!            Auxiliary variable.  Used to save compute operations when 
!            evaluating source terms on RHS of GKE that are proportional to ky.
!            !! The Mysterious factor 1/2 Explained !!
!            The factor of 1/2 arises because those source terms were first
!            specified using the normalisation Tref=mref vtref^2 
![old note by BD and MK on "Microinstabilities in Axisymmetric Configurations"]
!            Nowadays, internally gs2 uses Tref=(1/2) mref vtref^2 everywhere.
!            WUNITS restores the consistency of the normalisations in gs2.
! [R Numata et al, "AstroGK: Astrophysical gyrokinetics code", JCP, 2010].
!    WOUTUNITS: convert output frequencies to appear in USER v_t normalisation
!    FUNITS: convert output fluxes to appear in USER v_t normalisation        
!CMRend
    if (wstar_units) then
       wunits = 1.0
       where (aky /= 0.0)
          tunits = 2.0/aky
       elsewhere
          tunits = 0.0
       end where
       if (any(tunits == 0.0) .and. proc0) then
          write (error_unit(), *) &
               "WARNING: wstar_units=.true. and aky=0.0: garbage results"
          print *, &
               "WARNING: wstar_units=.true. and aky=0.0: garbage results"
       end if
!CMR: Sep 2010
!  Changes to allow wstar_units to be used consistently with either 
!  v_t normalisation option.   (Wasn't quite right before)
!CMRend
    else
       tunits = 1.0
       wunits = aky/2.0
    end if
    funits = tnorm
    woutunits = tnorm/tunits
  end subroutine adjust_time_norm

  subroutine finish_run_parameters

    implicit none

    if (allocated(wunits)) deallocate (wunits, woutunits, tunits)

    initialized = .false.

  end subroutine finish_run_parameters

end module run_parameters
