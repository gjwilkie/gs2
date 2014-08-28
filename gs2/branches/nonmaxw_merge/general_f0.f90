!> A module which allows the simulation of species with a
!! arbitrary background distribution function (e.g. alphas, 
!! beam ions).
!! 
!! Depending on the species type and the input parameters,
!! it will return either a Maxwellian or some other analytic
!! function or a solution for F_0 returned from an external
!! solver. 
!!
!! This is free software released under the GPLv3
!! Written by
!!    Edmund Highcock (edmundhighcock@sourceforge.net)
!!    George Wilkie  (gwilkie@umd.edu)
!!  

module general_f0

  !> Takes in the gridpoints in energy, as calculated 
  !! for the Maxwellian species (which will be in the first
  !! column of the input array, epoints(:,1)), 
  !! calculates any other grids required
  !! and calculates values of f0
  !! and all other arrays on those grids
  public :: calculate_f0_arrays

  !> A single local function that calculates f0. Takes one argument: 
  !! this feature will be useful with generalized quadrature 
  public :: eval_f0

  !> To facilitate the need to only call one function with one argument
  !! in eval_f0, user must first call this function to set the species
  public :: set_current_f0_species

  !> Initialises the module, chiefly reading the parameters.
  !! NB does not allocate arrays, as negrid must be provided
  !! by the egrid module.
  public :: init_general_f0

  !> Deallocate arrays and close output files. 
  public :: finish_general_f0

  !> Grid of f0 as function of energy and species
  public :: f0_values
  real, dimension (:,:), allocatable :: f0_values

  !> Grid of dF0/dE
  !! dF0/dE = (1/F0) * F0'(E) * E_ref
  !! where E_ref is the normalizing energy: temperature for Maxwellians,
  !! injection energy for alphas, or whichever convention is used.
  !! But E_N / T_ref = spec(is)%temp always!
  !! "generalised temperature" = -1/df0dE
  public :: df0dE
  real, dimension(:,:), allocatable:: df0dE

  !> Grid of generalised temperature 
  !! = -1/T^*_sN d(F_s/d energy_s)_N
  !! = - T_r/T^*_s d(F_s/d energy_s) T*_s/F_s  
  !! (where T*_s is just temperature for Maxwellian species)
  !!  as function of energy and species.
  !! For Maxwellian species this grid is just equal to T_r/T_s
  !! For alphas, T^*_s = E_alpha, the injection energy
  !! and this grid is calculated in this module
!  public :: generalised_temperature
!  real, dimension (:,:), allocatable :: generalised_temperature

  ! These arrays below replace the species quantites of the same name

  !> Equal to sqrt(generalised_temperature/mass) as a 
  !! function of energy and species
  public :: stm
  real, dimension (:,:), allocatable :: stm
  
  !> Generalized dF0/drho  for Maxwellian speices
  public :: f0prim
  real,dimension(:,:), allocatable:: f0prim

  ! get rid of this!
  ! GW: why?
  !> Equal to Z/sqrt(generalised_temperature/mass) as a 
  !! function of energy and species
  public :: zstm
  real, dimension (:,:), allocatable :: zstm

  !> Equal to generalised_temperature/Z as a 
  !! function of energy and species
  public :: gtempoz
  real, dimension (:,:), allocatable :: gtempoz

  !> Equal to Z/generalised_temperature as a 
  !! function of energy and species
  public :: zogtemp
  real, dimension (:,:), allocatable :: zogtemp

  ! probably get rid of this!
  ! GW: why?
  !> Equal to abs(sqrt(generalised_temperature/mass)/Z) as a 
  !! function of energy and species
  public :: smz
  real, dimension (:,:), allocatable :: smz

  !> Arrays are initially only calculated from proc0 as 
  !! calculate_f0_arrays is called from setvgrid. Later
  !! this function is called to put them on all procs
  public :: broadcast_arrays

  !> Unit tests
  !public :: general_f0_unit_test_init_general_f0
  !public :: general_f0_unit_test_calculate_f0_arrays

  private

  !> genquad not added yet
  !logical :: genquad_flag = .false.

  !> Sets the species that this module is concerned with.
  integer :: is_local = -1

  integer :: f0_opt_switch
  
  !> Flag that controls whether or not an externally-supplied f0 
  !! is rescaled to fit the given species parameters.
  !! Either F0 or dens must be changed so that they are self-consistent.
  !! f0_rescale = T -- F0 is rescaled to fit spec(is)%dens (recommended)
  !!            = F -- F0 is taken literally, spec(is)%dens is changed                     
  !! This option does not rescale according to spec(is)%temp!
  logical :: rescale_f0 = .true.

  integer, parameter :: f0_maxwellian = 1, &           !< Regular Maxwellian species
                        f0_external = 2, &             !< f0 determined from table in external input file
                        f0_sdanalytic = 3, &           !< Analytic Gaffey slowing-down distribution
                        f0_sdsemianalytic = 4, &       !< Abel's semi-analytic form of the slowing-down distribution, including ash (not currently implemented)
                        f0_kappa = 5                   !< Kappa distribution

  integer :: gen_f0_output_file
  
  integer :: negrid

  !real :: vcut
  !real :: energy_min

  real, dimension(:,:), allocatable :: egrid
  real, dimension(:,:), allocatable :: weights

  !real, dimension(:), allocatable :: egrid_maxwell
  !real, dimension(:), allocatable :: weights_maxwell

  !> Which species to use as the main ions for some options of falpha
  integer :: main_ion_species

  !> Options for calculating v_crit. Especially when electrons are adiabatic,
  !! their properties need to be inferred somehow because it affects the SD dist.
  !! vcrit_opt = 1: Electron temperature gradient a/LTe is taken to be the same as the ions
  !!           = 2: Use the a/LT input for the non-Maxwellian species as a proxy for the 
  !!                electron temperature gradient: a/LTe = spec(is)%tprim
  integer :: vcrit_opt = 1

  !> "Critical speed" that parameterizes the analytic slowing-down distribution 
  !! If negative or unspecified, vc is calculated from other species
  !! If positive, use this value for vc/v_ref,s
  real:: vcrit = -1.0

  !> Input parameter that gives a simple gradient in F0alpha (not a function of velocity - simply 1/Lnalpha)
  !! Not physical, but could diagnose problems.
  !! Obsoltete
  !! logical :: simple_gradient

  !> For the kappa distribution, the value of kappa
  real :: kappa

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Module admin
!!
!! This subsection contains routines for initializing/finalizing 
!! the module and reading parameters. 

  subroutine init_general_f0
    use file_utils, only: open_output_file
    use mp, only: proc0
    use species, only: has_nonmaxw_species, spec
    logical, save :: initialized = .false.
    !write (*,*) "initializing parameter_scan"

    if (initialized) return
    initialized = .true.

    if (has_nonmaxw_species(spec)) then
       call read_parameters
       if (proc0) call open_output_file(gen_f0_output_file, ".general_f0")
    end if
      
  end subroutine init_general_f0

  subroutine finish_general_f0
    use file_utils, only: close_output_file

    call close_output_file(gen_f0_output_file)
  end subroutine finish_general_f0

  subroutine read_parameters
    use file_utils, only: input_unit, error_unit, input_unit_exist
    use text_options, only: text_option, get_option_value
    use mp, only: proc0, broadcast
    implicit none
    namelist /general_f0_knobs/ rescale_f0, main_ion_species, &
            vcrit, kappa, vcrit_opt
    integer :: ierr, in_file
    logical :: exist

    if (proc0) then

       main_ion_species = -1 
       vcrit = -1.0
       kappa = -1.0
       vcrit_opt = 1
       rescale_f0 = .true.

       in_file = input_unit_exist ("general_f0_knobs", exist)
       if (exist) read (unit=in_file, nml=general_f0_knobs)

    end if

    call broadcast (main_ion_species)
    call broadcast (vcrit)
    call broadcast (kappa)
    call broadcast (vcrit_opt)
    call broadcast (rescale_f0)

  end subroutine read_parameters

  subroutine set_current_f0_species(is)
    implicit none
    integer,intent(in):: is
    is_local = is
  end subroutine

  real function eval_f0(v)
    use mp, only: mp_abort
    use species, only: ion_species, electron_species, nonmaxw_species, spec 
    use species, only: f0_maxwellian,f0_sdanalytic,f0_sdsemianalytic,f0_kappa, f0_external
    use constants, only: pi
    implicit none
    real,intent(in):: v
    integer:: is    
    real:: f0, dummy
  
    if ( is_local .LT. 0 ) then 
       write(*,*) "ERROR: Species needs to be set with set_current_f0_species before calling eval_f0" 
       call mp_abort('')
    end if

    is = is_local

    select case (spec(is)%type)
      case (ion_species)
        f0 = exp(-v**2)/(pi**1.5)
      case (electron_species)
        f0 = exp(-v**2)/(pi**1.5)
      case (nonmaxw_species)
        select case (spec(is)%f0_type)
        case (f0_maxwellian)
          f0 = exp(-v**2)/(pi**1.5)
        case (f0_sdanalytic)
          call eval_f0_sdanalytic(is,v,f0,dummy,dummy) 
        case (f0_sdsemianalytic)
          write(*,*) "ERROR: eval_f0 cannot yet handle semianalytic option."
          call mp_abort('')
        case (f0_kappa)
          call eval_f0_kappa(is,v,f0,dummy,dummy) 
        case (f0_external)
          write(*,*) "ERROR: eval_f0 cannot yet handle external option."
          call mp_abort('')
        end select
    end select

    eval_f0 = f0
    return
  end function eval_f0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! General functions
!  Functions which apply to all species

  !> NB this only gets called by proc0
  !! Arrays are broadcast later
  subroutine calculate_f0_arrays(epoints, wgts)
    use species, only: nspec, spec
    use species, only: ion_species, electron_species, nonmaxw_species
    real, dimension(:,:), intent(inout) :: epoints
    real, dimension(:,:), intent(inout) :: wgts
    !real, intent(in) :: vcut_in
    !logical,intent(in):: genquad_flag_in
    integer :: is

!    genquad_flag = genquad_flag_in
    negrid = size(epoints(:,1))

    call allocate_arrays

    ! Store input epoints and wgts into module-defined egrid and weights
    egrid = epoints
    weights = wgts

!    egrid_maxwell(:) = epoints(:,1)
!    weights_maxwell(:) = wgts(:,1)

    do is = 1,nspec
      select case (spec(is)%type)
      case (ion_species)
        call calculate_f0_arrays_maxwellian(is)
      case (electron_species)
        call calculate_f0_arrays_maxwellian(is)
      case (nonmaxw_species)
        select case (spec(is)%f0_type)
        case (f0_maxwellian)
          call calculate_f0_arrays_maxwellian(is)
        case (f0_sdanalytic)
          call calculate_f0_arrays_sdanalytic(is)
        case (f0_kappa)
          call calculate_f0_arrays_kappa(is)
!        case (f0_semianalytic)
!          call calculate_f0_arrays_sdsemianalytic(is)
        case (f0_external)
          call calculate_f0_arrays_external(is)
        end select
      end select
    end do
    epoints = egrid
    wgts = weights
    
  end subroutine calculate_f0_arrays

  subroutine allocate_arrays
    use species, only: nspec

    allocate(egrid(negrid,nspec))
    allocate(weights(negrid,nspec))
    !allocate(egrid_maxwell(negrid))
    !allocate(weights_maxwell(negrid))
    allocate(f0_values(negrid,nspec))
    allocate(df0dE(negrid,nspec))
    allocate(stm(negrid,nspec))
    allocate(zstm(negrid,nspec))
    allocate(zogtemp(negrid,nspec))
    allocate(gtempoz(negrid,nspec))
    allocate(smz(negrid,nspec))
    allocate(f0prim(negrid,nspec))

  end subroutine allocate_arrays

  subroutine broadcast_arrays
    use mp, only: proc0, broadcast
    
    call broadcast(negrid)
    if (.not. proc0) call allocate_arrays
    call broadcast(egrid)
    call broadcast(weights)
    !call broadcast(egrid_maxwell)
    !call broadcast(weights_maxwell)
    call broadcast(f0_values)
    call broadcast(df0dE)
    call broadcast(stm)
    call broadcast(zstm)
    call broadcast(zogtemp)
    call broadcast(gtempoz)
    call broadcast(smz)
    call broadcast(f0prim)
  end subroutine broadcast_arrays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Maxwellian
!! 
!! This subsection contains routines for handling 
!! Maxwellian distributions. 

  subroutine calculate_f0_arrays_maxwellian(is)
    use constants, only: pi
    use species, only: spec
    integer, intent(in) :: is
    integer :: ie
    ! Energy grid now set for every species in module egrid
    !egrid(:,is) = egrid_maxwell(:)
    f0_values(:, is) = exp(-egrid(:,is))/(pi**1.5)
    !if (.NOT. genquad_flag) 
    weights(:,is) = weights(:,is) * f0_values(:,is)
    do ie = 1,negrid
       df0dE(:,is) = -1.0/spec(is)%temp
       write(gen_f0_output_file,*) ie, egrid(ie,is), f0_values(ie,is), & 
                                   df0dE(ie,is)
    end do
    gtempoz(:,is) = -1.0 / (df0dE(:,is) * spec(is)%z)
    zogtemp(:,is) = -1.0 * df0dE(:,is) * spec(is)%z 
    
    f0prim(:,is) = -( spec(is)%fprim + (egrid(:,is) - 1.5)*spec(is)%tprim)
  end subroutine calculate_f0_arrays_maxwellian

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! External
!! 
!! This subsection contains routines for handling 
!! distributions that are given by an external
!! input file.

!> This subroutine calculates f0 on the grid from an external
!! input file. The grid on the input file can differ from that 
!! of gs2. A cubic spline is used to interpolate between the two.
!! The user can either specify df0/dE or it can be estimated internally.
!! The radial derivative of f0 as a function of energy should also be given.
  subroutine calculate_f0_arrays_external(is)
    !
    ! The egrid is already given by get_legendre_grids. 
    !
    ! The input file can take two forms, as controlled by the
    ! control parameter mode, the first integer read from the file
    !  mode = 1: E points in first columns, F0(E) in second. To calculate
    !            generalized temperature, a cubic spline is used to 
    !            interpolate and the slope is taken from that.
    !       = 2: First column is E, second is F0(E), the third is (1/F0)*dF0/dE
    !       = 3: E, F0(E), -(1/F0)*dF0/drho 
    !       = 4: E, F0(E),dF0/dE, -(1/F0)*dF0/drho 
    use mp, only: broadcast
    use constants, only: pi
    use species, only: spec, nspec
    use file_utils, only: run_name
    use splines, only: fitp_curvd, fitp_curv1, fitp_curv2
    implicit none
    integer, intent(in) :: is
    integer:: f0in_unit = 21, mode, ie, ierr, numdat, il, it
    logical:: acceptable_spline
    real:: test, moment0, moment2, sigma, sigma_fac
    real, dimension(:), allocatable:: f0_values_dat, df0dE_dat, egrid_dat, &
                                      yp, temp, df0drho_dat
    
    ! Open file and read column option
    open(unit=f0in_unit,file=trim(run_name)//'.f0in',status='old',action='read')
    read(f0in_unit,*) mode
    read(f0in_unit,*) numdat

    allocate(f0_values_dat(numdat))
    allocate(df0dE_dat(numdat))
    allocate(df0drho_dat(numdat))
    allocate(egrid_dat(numdat))
    allocate(yp(numdat))
    allocate(temp(numdat))

    sigma = 1.0
    sigma_fac = 2.0
    acceptable_spline = .false.

    if (mode .EQ. 1) then
       ! Read f0 values
       do ie = 1,numdat
          read(f0in_unit,*) egrid_dat(ie), f0_values_dat(ie)
       end do
       close(f0in_unit)

       ! Perform cubic spline to get F0 and its slope at grid points
       ! If F0 ends up being negative or zero anywhere on the grid, double the tension
       ! factor and try again

       do while (.NOT. acceptable_spline)
          ! Generate spline parameters
          call fitp_curv1(numdat,egrid_dat,f0_values_dat,0.0,0.0,3,yp,temp,sigma,ierr)
          if (ierr .NE. 0) then
             write(*,*) "fitp_curv1 returned error code ", ierr
             stop 1
          end if

          do ie = 1,negrid
             ! Interpolate to get f0 at grid points
             f0_values(ie,is) = fitp_curv2(egrid(ie,is),numdat,egrid_dat, &
                             f0_values_dat,yp,sigma)

             ! Calculate d/dE lnF0 to get generalised temperature
             df0dE(ie,is) =  fitp_curvd(egrid(ie,is),numdat,egrid_dat, f0_values_dat,yp,sigma)  / &
                                (f0_values(ie,is) * spec(is)%temp)

          end do

          if (minval(f0_values(:,is)) .GT. epsilon(0.0)) then 
             acceptable_spline = .true.
          else
             sigma = sigma * sigma_fac
          end if
       end do

       f0prim(:,is) = - spec(is)%fprim

       ! Diagnostic output
       do ie = 1,negrid
          write(gen_f0_output_file,*) ie, egrid(ie,is), f0_values(ie,is), df0dE(ie,is), f0prim(ie,is)
       end do

    else if (mode .EQ. 2) then
       ! Read both f0 and df0/dE
       do ie = 1,numdat
          read(f0in_unit,*) egrid_dat(ie), f0_values_dat(ie), df0dE_dat(ie)
       end do
       close(f0in_unit)

       do while (.NOT. acceptable_spline)
          ! Generate spline parameters for f0
          call fitp_curv1(numdat,egrid_dat,f0_values_dat,0.0,0.0,3,yp,temp,sigma,ierr)

          if (ierr .NE. 0) then
             write(*,*) "fitp_curv1 returned error code ", ierr
             stop 1
          end if

          do ie = 1,negrid

             ! Interpolate to get F0 at grid points
             f0_values(ie,is) = fitp_curv2(egrid(ie,is),numdat,egrid_dat, &
                                f0_values_dat,yp,sigma)
          end do

          ! Generate spline parameters for df0/dE
          call fitp_curv1(numdat,egrid_dat,df0dE_dat,0.0,0.0,3,yp,temp,sigma,ierr)
          if (ierr .NE. 0) then
             write(*,*) "fitp_curv1 returned error code ", ierr
             stop 1
          end if

          do ie = 1,negrid
             ! Interpolate to get f0 at grid points
             df0dE(ie,is)     = fitp_curv2(egrid(ie,is),numdat,egrid_dat, &
                                df0dE_dat,yp,sigma)  / spec(is)%temp
          end do
 
          if (minval(f0_values(:,is)) .GT. epsilon(0.0)) then 
             acceptable_spline = .true.
          else
             sigma = sigma * sigma_fac
          end if

       end do

       f0prim(:,is) = - spec(is)%fprim

       ! Diagnostic output
       do ie = 1,negrid
          write(gen_f0_output_file,*) ie, egrid(ie,is), f0_values(ie,is), df0dE(ie,is), f0prim(ie,is)
       end do

    else if (mode .EQ. 3) then
       ! Read both f0 and df0/drho
       do ie = 1,numdat
          read(f0in_unit,*) egrid_dat(ie), f0_values_dat(ie), df0drho_dat(ie)
       end do
       close(f0in_unit)

       do while (.NOT. acceptable_spline)
          ! Generate spline parameters for f0
          call fitp_curv1(numdat,egrid_dat,f0_values_dat,0.0,0.0,3,yp,temp,sigma,ierr)

          if (ierr .NE. 0) then
             write(*,*) "fitp_curv1 returned error code ", ierr
             stop 1
          end if

          do ie = 1,negrid
             ! Interpolate to get F0 at grid points
             f0_values(ie,is) = fitp_curv2(egrid(ie,is),numdat,egrid_dat, &
                                 f0_values_dat,yp,sigma)
             ! Calculate d/dE lnF0 to get generalised temperature
             df0dE(ie,is) = fitp_curvd(egrid(ie,is),numdat,egrid_dat, &
                             f0_values_dat,yp,sigma) / (f0_values(ie,is) * spec(is)%temp)
          end do

          ! Generate spline parameters for df0/drho
          call fitp_curv1(numdat,egrid_dat,df0drho_dat,0.0,0.0,3,yp,temp,sigma,ierr)
          if (ierr .NE. 0) then
             write(*,*) "fitp_curv1 returned error code ", ierr
             stop 1
          end if

          do ie = 1,negrid
             ! Interpolate to get f0 at grid points
             f0prim(ie,is)    = fitp_curv2(egrid(ie,is),numdat,egrid_dat, &
                                df0dE_dat,yp,sigma)
          end do
 
          if (minval(f0_values(:,is)) .GT. epsilon(0.0)) then 
             acceptable_spline = .true.
          else
             sigma = sigma * sigma_fac
          end if

       end do

       ! Diagnostic output
       do ie = 1,negrid
          write(gen_f0_output_file,*) ie, egrid(ie,is), f0_values(ie,is), df0dE(ie,is), f0prim(ie,is)
       end do    
    else if (mode .EQ. 4) then
       ! Read f0, df0/dE, df0/drho
       do ie = 1,numdat
          read(f0in_unit,*) egrid_dat(ie), f0_values_dat(ie), df0dE_dat(ie), df0drho_dat(ie)
       end do
       close(f0in_unit)

       do while (.NOT. acceptable_spline)
          ! Generate spline parameters for f0
          call fitp_curv1(numdat,egrid_dat,f0_values_dat,0.0,0.0,3,yp,temp,sigma,ierr)

          if (ierr .NE. 0) then
             write(*,*) "fitp_curv1 returned error code ", ierr
             stop 1
          end if

          do ie = 1,negrid

             ! Interpolate to get F0 at grid points
             f0_values(ie,is) = fitp_curv2(egrid(ie,is),numdat,egrid_dat, &
                                f0_values_dat,yp,sigma)
          end do

          ! Generate spline parameters for df0/dE
          call fitp_curv1(numdat,egrid_dat,df0dE_dat,0.0,0.0,3,yp,temp,sigma,ierr)
          if (ierr .NE. 0) then
             write(*,*) "fitp_curv1 returned error code ", ierr
             stop 1
          end if

          do ie = 1,negrid
             ! Interpolate to get f0 at grid points
             df0dE(ie,is)     = fitp_curv2(egrid(ie,is),numdat,egrid_dat, &
                                df0dE_dat,yp,sigma)  / spec(is)%temp
          end do
 
          ! Generate spline parameters for df0/drho
          call fitp_curv1(numdat,egrid_dat,df0drho_dat,0.0,0.0,3,yp,temp,sigma,ierr)
          if (ierr .NE. 0) then
             write(*,*) "fitp_curv1 returned error code ", ierr
             stop 1
          end if

          do ie = 1,negrid
             ! Interpolate to get f0 at grid points
             f0prim(ie,is)     = fitp_curv2(egrid(ie,is),numdat,egrid_dat, &
                                 df0drho_dat,yp,sigma) 
          end do
 
          if (minval(f0_values(:,is)) .GT. epsilon(0.0)) then 
             acceptable_spline = .true.
          else
             sigma = sigma * sigma_fac
          end if

       end do

       ! Diagnostic output
       do ie = 1,negrid
          write(gen_f0_output_file,*) ie, egrid(ie,is), f0_values(ie,is), df0dE(ie,is), f0prim(ie,is)
       end do

       f0prim(:,is) = - spec(is)%fprim

    else
       write(*,*) "ERROR. First line in f0 input file should be mode" 
       stop 1
    end if

    ! Now calculate moments of F0 from trapezoidal rule.
    ! Should probably do this by quadrature, 
    ! but le_grids is not defined yet

    moment0 = 0.0
!    do ie = 1,negrid-1
!       moment0  = moment0 + 0.5*(grid(ie+1,is)-egrid(ie,is)) * &
!                           ( sqrt(egrid(ie,is))*f0_values(ie,is) + &
!                             sqrt(egrid(ie+1,is))*f0_values(ie+1,is) )
!    end do
    do ie = 1,numdat-1
       moment0  = moment0 + 0.5*(egrid_dat(ie+1)-egrid_dat(ie)) * &
                           ( sqrt(egrid_dat(ie))*f0_values_dat(ie) + &
                             sqrt(egrid_dat(ie+1))*f0_values_dat(ie+1) )
    end do
    moment0 = moment0*4.0*pi

!    ! Input parameter rescale_f0 determines the priority between the input species
!    ! density or the input F0.
!    if (rescale_f0) then
!       write(*,*) "rescale_f0 = T: Rescaling magnitude of F0 to agree with input dens"
!       write(*,*) "rescale factor is ", spec(is)%dens / moment0
!       f0_values = f0_values * spec(is)%dens / moment0
!    else
!       write(*,*) "rescale_f0 = F: Overriding species dens to agree with input F0"
!       write(*,*) "New dens = ", moment0
!       spec(is)%dens = moment0
!    end if

    gtempoz(:,is) = - 1.0/ (df0dE(ie,is)*spec(is)%z)
    zogtemp(:,is) = - df0dE(ie,is)* spec(is)%z 
    
    !if (.NOT. genquad_flag)
    weights(:,is) = weights(:,is) * f0_values(:,is)

    deallocate(egrid_dat,f0_values_dat,df0dE_dat,df0drho_dat,yp,temp)

  end subroutine calculate_f0_arrays_external

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Analytic Slowing-down
!! 
!! This subsection contains routines for handling 
!! the analytic approximation to the fast-particle 
!! slowing-down distribution 
!!
!! f0(v) = A/(vc^3 + v^3)                                           (if v < vinj)
!! 
!! vc    = critical speed, input
!! A     = normalization chosen so that f0 integrates to unity
!! vinj  = injection speed
!!
  subroutine calculate_f0_arrays_sdanalytic(is)
    use constants, only: pi
    use species, only: spec, nspec
    use mp, only: mp_abort, proc0
    use species, only: ion_species, electron_species, nonmaxw_species
    implicit none
    integer,intent(in):: is
    real:: v
    integer:: ie, i, electron_spec

    do ie = 1,negrid 
       v = sqrt(egrid(ie,is))
       call eval_f0_sdanalytic(is,v,f0_values(ie,is),df0dE(ie,is),f0prim(ie,is))
    end do

    gtempoz(:,is) = -1.0/(df0dE(:,is) * spec(is)%z)
    zogtemp(:,is) = - df0dE(:,is) * spec(is)%z 

    !if (.NOT. genquad_flag) 
    weights(:,is) = weights(:,is) * f0_values(:,is)

  end subroutine calculate_f0_arrays_sdanalytic

  subroutine eval_f0_sdanalytic(is,v,f0,df0dE_out,f0prim)
    use species, only: ion_species, electron_species, nonmaxw_species, spec, nspec
    use species, only: calculate_slowingdown_parameters
    use constants, only: pi
    use mp, only: mp_abort, proc0
    use spfunc, only: erf => erf_ext
    implicit none
    integer,intent(in):: is
    real,intent(in):: v
    real,intent(inout):: f0,df0dE_out, f0prim
    real:: A, df0dv,vcprim, vcrit
    integer:: ie, i, electron_spec

    if (vcrit .LE. 0.0) then
       call calculate_slowingdown_parameters(is,vcrit_opt,vcrit,vcprim)
    else
       !> If vcrit is specified, there is no way of knowing its radial derivative, so set dvc/drho = 0
       vcprim = 0.0
    end if

    A = (3.0/(4.0*pi))/log( (vcrit**3 + 1.0)/(vcrit**3))

    f0 = A/(vcrit**3 + v**3)
    df0dv = -A*3.0*v**2/(vcrit**3 + v**3)**2
    df0dE_out =  (0.5/v)*df0dv / ( f0 * spec(is)%temp )
    f0prim = -spec(is)%fprim - (-(vcrit**3/(vcrit**3 + v**3)) + (1.0/( log(1.0 + vcrit**(-3)) * (1.0 + vcrit**(3)))))*vcprim

  end subroutine eval_f0_sdanalytic

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Kappa distribution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 
  subroutine eval_f0_kappa(is,v,f0,df0dE_out,f0prim)
    use species, only: ion_species, electron_species, nonmaxw_species, spec, nspec
    use constants, only: pi
    use mp, only: mp_abort, proc0
    use spfunc, only: lgamma
    implicit none
    integer,intent(in):: is
    real,intent(in):: v
    real,intent(inout):: f0,df0dE_out, f0prim
    real:: gammakmh,gammak,A,df0dv
   
    gammakmh = exp(lgamma(kappa-0.5))  
    gammak = exp(lgamma(kappa))

    A = 1.0/(pi**1.5*gammakmh/gammak)

    f0 = A/(1.0 + v**2/kappa)**(kappa+1.0)

    !if (.NOT. genquad_flag) weights(:,is) = weights(:,is) * f0_values(:,is)
    weights(:,is) = weights(:,is) * f0_values(:,is)
 
    df0dv = - 2.0*v*A*(kappa+1.0)/( kappa* (1.0 + v**2/kappa))
    df0dE_out = (0.5/v)*df0dv / (f0*spec(is)%temp)
 
    f0prim = -spec(is)%fprim
     
  end subroutine eval_f0_kappa

  subroutine calculate_f0_arrays_kappa(is)
    use constants, only: pi
    use species, only: spec, nspec
    use mp, only: mp_abort, proc0
    implicit none
    integer,intent(in):: is
    real:: v, df0dv
    integer:: ie, i, electron_spec


    do ie = 1,negrid 
       v = sqrt(egrid(ie,is))
       call eval_f0_kappa(is,v,f0_values(ie,is),df0dE(ie,is),f0prim(ie,is))
       gtempoz(ie,is) = -1.0 / (df0dE(ie,is)* spec(is)%z)
       zogtemp(ie,is) = - df0dE(ie,is) *spec(is)%z 
    end do


    !if (.NOT. genquad_flag) weights(:,is) = weights(:,is) * f0_values(:,is)
    weights(:,is) = weights(:,is) * f0_values(:,is)

  end subroutine calculate_f0_arrays_kappa

!  function general_f0_unit_test_init_general_f0()
!    use unit_tests
!    use species
!    logical :: general_f0_unit_test_init_general_f0
!
!    call init_general_f0
!
!    call announce_check('alpha_f0_switch')
!    general_f0_unit_test_init_general_f0 = &
!      agrees_with(alpha_f0_switch, alpha_f0_split)
!    call process_check(general_f0_unit_test_init_general_f0, 'alpha_f0_switch')
!
!    call announce_check('spec 1 type')
!    general_f0_unit_test_init_general_f0 = &
!      general_f0_unit_test_init_general_f0 .and. &
!      agrees_with(spec(1)%type, ion_species)
!    call process_check(general_f0_unit_test_init_general_f0, 'spec 1 type')
!    call announce_check('spec 2 type')
!    general_f0_unit_test_init_general_f0 = &
!      general_f0_unit_test_init_general_f0 .and. &
!      agrees_with(spec(2)%type, electron_species)
!    call process_check(general_f0_unit_test_init_general_f0, 'spec 2 type')
!    call announce_check('spec 3 type')
!    general_f0_unit_test_init_general_f0 = &
!      general_f0_unit_test_init_general_f0 .and. &
!      agrees_with(spec(3)%type, nonmaxw_species)
!    call process_check(general_f0_unit_test_init_general_f0, 'spec 3 type')
!    
!  end function general_f0_unit_test_init_general_f0
!
!  function general_f0_unit_test_calculate_f0_arrays(epoints, wgts, vcut_in, rslts, err)
!    use unit_tests
!    real, dimension(:,:), intent(inout) :: epoints
!    real, dimension(:,:), intent(inout) :: wgts
!    real, dimension(:,:,:), intent(in) :: rslts
!    real, intent(in) :: err
!    real, intent(in) :: vcut_in
!    logical :: general_f0_unit_test_calculate_f0_arrays
!    logical :: tr ! test results
!    tr = .true.
!    call calculate_f0_arrays(epoints, wgts, vcut_in,genquad_flag)
!
!    call announce_check('ion energy grid')
!    tr = tr .and. agrees_with(egrid(:,1), epoints(:,1), err)
!    call process_check(tr, 'ion energy grid')
!
!    call announce_check('electron energy grid')
!    tr = tr .and. agrees_with(egrid(:,2), epoints(:,1), err)
!    call process_check(tr, 'electron energy grid')
!
!    call announce_check('ion f0')
!    tr = tr .and. agrees_with(f0_values(:,1), rslts(:,1,1), err)
!    call process_check(tr, 'ion f0')
!    call announce_check('electron f0')
!    tr = tr .and. agrees_with(f0_values(:,2), rslts(:,2,1), err)
!    call process_check(tr, 'electron f0')
!    call announce_check('alpha f0')
!    tr = tr .and. agrees_with(f0_values(:,3), rslts(:,3,1), err)
!    call process_check(tr, 'alpha f0')
!
!    call announce_check('ion gentemp')
!    tr = tr .and. agrees_with(generalised_temperature(:,1), rslts(:,1,2), err)
!    call process_check(tr,' ion gentemp')
!    call announce_check('electron gentemp')
!    tr = tr .and. agrees_with(generalised_temperature(:,2), rslts(:,2,2), err)
!    call process_check(tr,' electron gentemp')
!    call announce_check('alpha gentemp')
!    tr = tr .and. agrees_with(generalised_temperature(:,3), rslts(:,3,2), err)
!    call process_check(tr,' alpha gentemp')
!
!    !write (*,*) 'f0prim(:,3),', f0prim(:,3)
!    call announce_check('ion f0prim')
!    !tr = tr .and. agrees_with(f0prim(:,1), rslts(:,1,3), err)
!    call process_check(tr,' ion f0prim')
!    call announce_check('electron f0prim')
!    !tr = tr .and. agrees_with(f0prim(:,2), rslts(:,2,3), err)
!    call process_check(tr,' electron f0prim')
!    call announce_check('alpha f0prim')
!    tr = tr .and. agrees_with(f0prim(:,3), rslts(:,3,3), err*10.0)
!    call process_check(tr,' alpha f0prim')
!    !write (*,*) 'f0prim(:,3),', f0prim(:,3)
!
!
!    general_f0_unit_test_calculate_f0_arrays = tr
!
!  end function general_f0_unit_test_calculate_f0_arrays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Semianalytic
!! 
!! This subsection contains routines for handling 
!! distributions given by the semianalytical form,
!! ammended such that all alphas, including ash,
!! are combined in one species.

!  subroutine calculate_f0_arrays_sdsemianalytic(is)
!    use semianalytic_falpha, only: semianalytic_falpha_parameters_type
!    use semianalytic_falpha, only: calculate_arrays
!    use constants, only: pi
!    use species, only: spec 
!    use species, only: nspec
!    use species, only: electron_species
!    use species, only: ion_species
!    use species, only: has_electron_species
!    use mp, only: mp_abort, proc0
!    type(semianalytic_falpha_parameters_type) :: parameters
!    integer, intent(in) :: is
!    integer :: ie
!    integer :: electron_spec
!    integer :: i
!    real :: shift, scal ! For modifying the energy grid for alphas
!
!    !> Determine which species is the electrons
!    !main_ion_species = -1
!    do i = 1,nspec
!      if (spec(i)%type .eq. electron_species) electron_spec = i
!      if (main_ion_species < 1 .and. spec(i)%type .eq. ion_species) &
!        main_ion_species = i
!    end do
!
!    write(*,*) 'main_ion_species', main_ion_species, 'electron_spec', electron_spec
!
!    !if (main_ion_species < 1)  call mp_abort(&
!      !'main_ion_species < 1: please set main_ion_species in general_f0_parameters') 
!
!    parameters%source_prim = spec(is)%sprim
!    
!    parameters%alpha_density = spec(is)%dens      
!
!    parameters%alpha_is = is
!
!    parameters%alpha_ion_collision_rate = spec(is)%gamma_ai
!    parameters%alpha_electron_collision_rate = spec(is)%gamma_ae
!
!    parameters%alpha_vth       = spec(is)%stm
!
!    parameters%source          = spec(is)%source
!    parameters%ash_fraction    = spec(is)%ash_fraction
!    parameters%source_prim     = spec(is)%sprim
!    !write (*,*) 'Source prim is', parameters%source_prim
!
!    parameters%alpha_injection_energy = spec(is)%temp
!
!    parameters%alpha_charge    = spec(is)%z
!
!    parameters%alpha_mass      = spec(is)%mass
!
!    parameters%negrid = negrid
!    if (has_electron_species(spec)) then 
!      parameters%electron_fprim   = spec(electron_spec)%fprim
!      parameters%electron_tprim   = spec(electron_spec)%tprim
!      parameters%electron_mass   = spec(electron_spec)%mass
!      parameters%electron_charge = spec(electron_spec)%z
!      parameters%electron_vth    = spec(electron_species)%stm
!      parameters%electron_temp   = spec(electron_species)%temp
!      if (main_ion_species>0) then 
!        parameters%ash_temp        = spec(main_ion_species)%temp
!        parameters%ion_mass        = spec(main_ion_species)%mass
!        parameters%ion_tprim        = spec(main_ion_species)%tprim
!        parameters%ion_fprim        = spec(main_ion_species)%fprim
!        parameters%ion_charge      = spec(main_ion_species)%z
!        parameters%ion_vth         = spec(main_ion_species)%stm
!        parameters%ion_temp        = spec(main_ion_species)%temp
!      else
!        ! Assume main ions are deuterium with ti = te
!        parameters%ion_mass       = spec(electron_spec)%mass * 2.0 * 1836.0
!        parameters%ion_tprim      = 0.0
!        parameters%ion_fprim      = spec(electron_spec)%fprim * 0.0
!        parameters%ion_charge     = -spec(electron_spec)%z
!        parameters%ion_temp       = spec(electron_spec)%temp
!        parameters%ion_vth        = sqrt(parameters%ion_temp/parameters%ion_mass)
!      end if
!
!
!    else
!      if (main_ion_species>0) then 
!
!        parameters%electron_fprim   =  0.0
!        parameters%electron_tprim   =  0.0
!        ! Assume main ions are deuterium with ti = te
!        parameters%electron_mass   = spec(main_ion_species)%mass/1836.0/2.0
!        parameters%electron_charge = -spec(main_ion_species)%z
!        parameters%electron_temp   = spec(main_ion_species)%temp
!        parameters%electron_vth    = sqrt(parameters%electron_temp/parameters%electron_mass)
!      else
!        write (*,*) 'You have no electrons and no ions: not sure how to set'
!        write (*,*) 'alpha parameters!'
!        call mp_abort('')
!      end if
!    end if
!
!    !if (parameters%source .eq. 0.0) then 
!      !write (*,*) 'You have source = 0.0 for alphas!!!'
!      !call mp_abort(' ')
!    !end if
!
!    ! This bit of code shifts the grid so that vcut for alphas is 1
!    ! and the lowest point is energy min
!    !scal =  (energy_min - 1.0) / (egrid(1,1) - vcut)
!    !shift = energy_min - egrid(1,1) * scal
!    !egrid(:,is) = egrid_maxwell(:) *(1.0 - energy_min)/vcut + energy_min
!    !egrid(:,is) = egrid_maxwell(:) *scal + shift
!
!    !write (*,*) 'parameters', parameters
!    !write(*,*) 'calling calculate_arrays'
!    !write (*,*) 'egrid', egrid(:,is)
!    call calculate_arrays(parameters,&
!                          egrid, &
!                          f0_values, &
!                          generalised_temperature, &
!                          f0prim)
!    spec(is)%source = parameters%source
!    !write (*,*) 'f0prim', ',is', f0prim(:,is), is
!    !write (*,*) 'f0prim(:,3),', f0prim(:,3)
!    !write (*,*) 'f0_values(:,3),', f0_values(:,3)
!    !weights(:,is) =  weights_maxwell(:)*scal * f0_values(:,is)
!    if (.NOT. genquad_flag) weights(:,is) = weights(:,is) * f0_values(:,is)
!    gtempoz(:,is) = generalised_temperature(:,is) / spec(is)%z
!    zogtemp(:,is) = spec(is)%z / generalised_temperature(:,is)
!
!    if (proc0) write(*,*) "source = ", parameters%source
!    if (proc0) write(*,*) "ash_fraction = ", parameters%ash_fraction
!!    do ie = 1,negrid
!!       if (proc0) write(*,*) sqrt(egrid(ie,is)), f0_values(ie,is), generalised_temperature(ie,is)
!!    end do
!    
!  end subroutine calculate_f0_arrays_sdsemianalytic
!
!


end module general_f0
