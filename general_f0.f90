

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


  !> Initialises the module, chiefly reading the parameters.
  !! NB does not allocate arrays, as negrid must be provided
  !! by the egrid module.
  public :: init_general_f0

  !> Deallocate arrays and close output files. 
  public :: finish_general_f0

  !> Grid of f0 as function of energy and species
  public :: f0_values
  real, dimension (:,:), allocatable :: f0_values

  !> Grid of generalised temperature 
  !! = -1/T^*_sN d(F_s/d energy_s)_N
  !! = - T_r/T^*_s d(F_s/d energy_s) T*_s/F_s  
  !! (where T*_s is just temperature for Maxwellian species)
  !!  as function of energy and species.
  !! For Maxwellian species this grid is just equal to T_r/T_s
  !! For alphas, T^*_s = E_alpha, the injection energy
  !! and this grid is calculated in this module

  public :: generalised_temperature
  real, dimension (:,:), allocatable :: generalised_temperature

  ! These arrays below replace the species quantites of the same name

  !> Equal to sqrt(generalised_temperature/mass) as a 
  !! function of energy and species
  public :: stm
  real, dimension (:,:), allocatable :: stm
  
  !> Generalized dF0/drho  for Maxwellian speices
  public :: f0prim
  real,dimension(:,:), allocatable:: f0prim

  ! get rid of this!
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
  !> Equal to abs(sqrt(generalised_temperature/mass)/Z) as a 
  !! function of energy and species
  public :: smz
  real, dimension (:,:), allocatable :: smz

  !> Arrays are initially only calculated from proc0 as 
  !! calculate_f0_arrays is called from setvgrid. Later
  !! this function is called to put them on all procs
  public :: broadcast_arrays

  private

  integer :: alpha_f0_switch, &
             beam_f0_switch
  
  real :: alpha_Einj

  logical :: print_egrid
 
  !> Flag that controls whether or not an externally-supplied f0 
  !! is rescaled to fit the given species parameters.
  !! f0_rescale = T -- f0 is rescaled to fit spec(is)%dens
  !!            = F -- f0 is taken literally, spec(is)%dens is changed                     
  !! This option does not rescale according to spec(is)%temp!
  logical :: rescale_f0

  integer, parameter :: alpha_f0_maxwellian = 1, &
                        alpha_f0_analytic = 2, &
                        alpha_f0_external =3

  integer, parameter :: beam_f0_maxwellian = 1, &
                        beam_f0_analytic = 2, &
                        beam_f0_external =3

  integer :: gen_f0_output_file
  
  integer :: negrid

  real, dimension(:,:), allocatable :: egrid

  real, dimension(:), allocatable :: egrid_maxwell

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Module admin
!!
!! This subsection contains routines for initializing/finalizing 
!! the module and reading parameters. 

  subroutine init_general_f0
    use file_utils, only: open_output_file
    use mp, only: proc0
    logical, save :: initialized = .false.
    !write (*,*) "initializing parameter_scan"

    if (initialized) return
    initialized = .true.


    call read_parameters


    if (proc0) call open_output_file(gen_f0_output_file, ".general_f0")
      

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
    type (text_option), dimension (3), parameter :: alpha_f0_opts = &
         (/ text_option('maxwellian', alpha_f0_maxwellian), &
            text_option('analytic', alpha_f0_analytic), &
            text_option('external', alpha_f0_external) /)
    character(20) :: alpha_f0
    type (text_option), dimension (3), parameter :: beam_f0_opts = &
         (/ text_option('maxwellian', beam_f0_maxwellian), &
            text_option('analytic', beam_f0_analytic), &
            text_option('external', beam_f0_external) /)
    character(20) :: beam_f0
    namelist /general_f0_parameters/ &
            alpha_f0, &
            beam_f0, &
            alpha_Einj, rescale_f0, print_egrid

    integer :: ierr, in_file
    logical :: exist

    if (proc0) then

       alpha_f0 = 'maxwellian'
       beam_f0 = 'maxwellian'

       in_file = input_unit_exist ("general_f0_parameters", exist)
       if (exist) read (unit=in_file, nml=general_f0_parameters)

       ierr = error_unit()
       call get_option_value &
            (alpha_f0, alpha_f0_opts, alpha_f0_switch, &
            ierr, "alpha_f0 in general_f0_parameters")
       call get_option_value &
            (beam_f0, beam_f0_opts, beam_f0_switch, &
            ierr, "beam_f0 in general_f0_parameters")

    end if

    call broadcast (alpha_f0_switch)
    call broadcast (beam_f0_switch)

  end subroutine read_parameters

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! General functions
!  Functions which apply to all species

  !> NB this only gets called by proc0
  !! Arrays are broadcast later
  subroutine calculate_f0_arrays(epoints)
    use species, only: nspec, spec
    use species, only: ion_species, electron_species, alpha_species
    use species, only: beam_species
    real, dimension(:,:), intent(inout) :: epoints
    integer :: is
    negrid = size(epoints)
    call allocate_arrays
    egrid = epoints
    egrid_maxwell(:) = epoints(:,1)
    do is = 1,nspec
      select case (spec(is)%type)
      case (ion_species)
        call calculate_f0_arrays_maxwellian(is)
      case (electron_species)
        call calculate_f0_arrays_maxwellian(is)
      case (alpha_species)
        select case (alpha_f0_switch)
        case (alpha_f0_maxwellian)
          call calculate_f0_arrays_maxwellian(is)
        case (alpha_f0_external)
          call calculate_f0_arrays_external(is)
        end select
      case (beam_species)
        select case (beam_f0_switch)
        case (beam_f0_maxwellian)
          call calculate_f0_arrays_maxwellian(is)
        end select
      end select
    end do
    epoints = egrid
    
  end subroutine calculate_f0_arrays

  subroutine allocate_arrays
    use species, only: nspec

    allocate(egrid(negrid,nspec))
    allocate(egrid_maxwell(negrid))
    allocate(f0_values(negrid,nspec))
    allocate(generalised_temperature(negrid,nspec))
    allocate(stm(negrid,nspec))
    allocate(zstm(negrid,nspec))
    allocate(zogtemp(negrid,nspec))
    allocate(gtempoz(negrid,nspec))
    allocate(smz(negrid,nspec))
    allocate(F0prim(negrid,nspec))

  end subroutine allocate_arrays

  subroutine broadcast_arrays
    use mp, only: proc0, broadcast
    
    call broadcast(negrid)
    if (.not. proc0) call allocate_arrays
    call broadcast(egrid)
    call broadcast(egrid_maxwell)
    call broadcast(f0_values)
    call broadcast(generalised_temperature)
    call broadcast(stm)
    call broadcast(zstm)
    call broadcast(zogtemp)
    call broadcast(gtempoz)
    call broadcast(smz)
    call broadcast(F0prim)
  end subroutine broadcast_arrays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Maxwellian
!! 
!! This subsection contains routines for handling 
!! Maxwellian distributions. 


  subroutine calculate_f0_arrays_maxwellian(is)
    use species, only: spec
    integer, intent(in) :: is
    integer :: ie
    egrid(:,is) = egrid_maxwell(:)
    f0_values(:, is) = exp(-egrid(:,is))
    do ie = 1,negrid
       generalised_temperature(:,is) = spec(is)%temp
       if (print_egrid) write(*,*) ie, egrid(ie,is), f0_values(ie,is), & 
                                   generalised_temperature(ie,is)
    end do
    gtempoz(:,is) = generalised_temperature(:,is) / spec(is)%z
    zogtemp(:,is) = spec(is)%z / generalised_temperature(:,is)
    
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
!! of gs2. A cubic-spline is used to interpolate between the two.
!! This is a logarithmic interpolation on the splines to ensure positive-
!! definiteness of F0. The user can either specify df0/dE or it 
!! can be estimated internally.
  subroutine calculate_f0_arrays_external(is)
    !
    ! The egrid is already given by get_legendre_grids. The input values
    ! are cubic-spline interpolated to the energy grid of gs2.
    !
    ! The input file can take two forms, as controlled by the
    ! control paramter num_cols, the first integer 
    ! read from the file
    !  - Single-Column Mode: Only F0(E) is given. To calculate
    !      generalized temperature, a cubic spline is used to 
    !      interpolate and the slope is taken from that.
    !  - Two-Column Mode: First column is F0(E), the second is
    !      dF0/dE.
    ! 
    use le_grids, only: w
    use mp, only: broadcast
    use species, only: spec, nspec
    use file_utils, only: run_name
    use splines, only: fitp_curvd, fitp_curv1, fitp_curv2
    implicit none
    integer, intent(in) :: is
    integer:: f0in_unit = 21, num_cols, ie, ierr, numdat, js
    real:: df0dE, n0_alpha, test
    real:: pick_spec(nspec)
    real, dimension(:), allocatable:: f0_values_dat, df0dE_dat, egrid_dat, &
                                      f0_values_dat_log, df0dE_dat_log, yp, temp
    
    ! Open file and read column option
    open(unit=f0in_unit,file=trim(run_name)//'.f0in',status='old',action='read')
    read(f0in_unit,*) num_cols
    read(f0in_unit,*) numdat

    allocate(f0_values_dat(numdat))
    allocate(df0dE_dat(numdat))
    allocate(egrid_dat(numdat))
    allocate(f0_values_dat_log(numdat))
    allocate(df0dE_dat_log(numdat))
    allocate(yp(numdat))
    allocate(temp(numdat))

    if (num_cols .EQ. 2) then
       ! Read f0 values
       do ie = 1,numdat
          read(f0in_unit,*) egrid_dat(ie), f0_values_dat(ie)
          ! Interpolate the *logarithm* of f0 to ensure positive-definiteness
          f0_values_dat_log(ie) = log(f0_values_dat(ie))
       end do
       close(f0in_unit)

       ! Perform cubic spline to get F0 and its slope at grid points

       ! Generate spline parameters
       call fitp_curv1(numdat,egrid_dat,f0_values_dat_log,0.0,0.0,3,yp,temp,1.0,ierr)
       if (ierr .NE. 0) then
          write(*,*) "fitp_curv1 returned error code ", ierr
          stop 1
       end if

       do ie = 1,negrid
          ! Interpolate to get f0 at grid points
          f0_values(ie,is) = fitp_curv2(egrid(ie,is),numdat,egrid_dat, &
                             f0_values_dat_log,yp,1.0)

          ! Recover F0 from its logarithm
          f0_values(ie,is) = exp(f0_values(ie,is))

          ! Calculate d/dE lnF0 to get generalised temperature
          df0dE = fitp_curvd(egrid(ie,is),numdat,egrid_dat, &
                             f0_values_dat_log,yp,1.0)
          generalised_temperature(ie,is) = -alpha_Einj/df0dE

          ! Diagnostic output
          if (print_egrid) write(*,*) ie, egrid(ie,is), f0_values(ie,is), df0dE, & 
                                   generalised_temperature(ie,is)
       end do

    else if (num_cols .EQ. 3) then
       ! Read both f0 and df0/dE
       do ie = 1,numdat
          read(f0in_unit,*) egrid_dat(ie), f0_values_dat(ie), df0dE_dat(ie)
          f0_values_dat_log(ie) = log(f0_values_dat(ie))
          df0dE_dat_log(ie) = log(abs(df0dE_dat(ie)))
       end do
       close(f0in_unit)

       ! Generate spline parameters for f0
       call fitp_curv1(numdat,egrid_dat,f0_values_dat_log,0.0,0.0,3,yp,temp,1.0,ierr)

       if (ierr .NE. 0) then
          write(*,*) "fitp_curv1 returned error code ", ierr
          stop 1
       end if

       do ie = 1,negrid

          ! Interpolate to get F0 at grid points
          f0_values(ie,is) = fitp_curv2(egrid(ie,is),numdat,egrid_dat, &
                             f0_values_dat_log,yp,1.0)
       end do

       ! Recover F0 from its logarithm
       f0_values = exp(f0_values)

       ! Generate spline parameters for df0/dE
       call fitp_curv1(numdat,egrid_dat,df0dE_dat_log,0.0,0.0,3,yp,temp,1.0,ierr)
       if (ierr .NE. 0) then
          write(*,*) "fitp_curv1 returned error code ", ierr
          stop 1
       end if

       do ie = 1,negrid
          ! Interpolate to get f0 at grid points
          df0dE            = fitp_curv2(egrid(ie,is),numdat,egrid_dat, &
                             df0dE_dat_log,yp,1.0)
 
          ! Recover df0/dE from its logarithm (maintaining whatever sign it had before)
          df0dE = sign(1.0,df0dE_dat(ie))* exp(df0dE)

          generalised_temperature(ie,is) = -alpha_Einj*f0_values(ie,is)/df0dE

          ! Diagnostic output
          if (print_egrid) write(*,*) ie, egrid(ie,is), f0_values(ie,is), df0dE, & 
                                   generalised_temperature(ie,is)

       end do
    else
       write(*,*) "ERROR. First line in f0 input file should be num_cols: " 
       write(*,*) " num_cols = 1 if only f0 is to be input, "
       write(*,*) " num_cols=2 if f0 and df0/dE are input."
       stop 1
    end if

    ! This is wrong. Need to calculate 0th moment properly. For now, just trust that
    ! f0 as input agrees with spec(is)%dens
    ! Calculate 0th moment of f0 as input, and rescale either f0 or n0 consistently

!    f0_values = f0_values/sqrt(2.0)

!    n0_alpha = 0.0
!    do ie = 1,negrid
!       n0_alpha = n0_alpha + w(ie,is)*f0_values(ie,is)*4.0*egrid(ie,is)
!    end do
!    write(*,*) "n0_alpha = ", n0_alpha

!    if (rescale_f0) then
!       ! Calculate 0th moment of f0 as input, and rescale according to n0
!       f0_values(ie,is) = f0_values(ie,is) * spec(is)%dens / n0_alpha
!    else
!       spec(is)%dens = n0_alpha
!    end if
!    call broadcast(spec(is)%dens)

    gtempoz(:,is) = generalised_temperature(:,is) / spec(is)%z
    zogtemp(:,is) = spec(is)%z / generalised_temperature(:,is)
    
    f0prim(:,is) = -( spec(is)%fprim + (egrid(:,is) - 1.5)*spec(is)%tprim)


    deallocate(egrid_dat,f0_values_dat,df0dE_dat,f0_values_dat_log,df0dE_dat_log,yp,temp)

  end subroutine calculate_f0_arrays_external

end module general_f0
