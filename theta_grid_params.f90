module theta_grid_params
  implicit none

  private

  public :: init_theta_grid_params
  public :: finish_theta_grid_params
  public :: wnml_theta_grid_params, write_trinity_parameters
  public :: rhoc, rmaj, r_geo, eps, epsl, qinp, shat, alpmhd
  public :: pk, shift, akappa, akappri, tri, tripri, asym, asympri
  public :: btor_slab, betaprim, ntheta, nperiod
  public :: set_overrides

  real :: rhoc, rmaj, r_geo, eps, epsl
  real :: qinp, shat, alpmhd, pk, shift, akappa, akappri, tri, tripri
  real :: asym, asympri, btor_slab, betaprim

  integer :: ntheta, nperiod

  logical :: initialized = .false.
  real :: kp = -1.
  logical :: exist

contains
  subroutine init_theta_grid_params
    use unit_tests, only: debug_message
    implicit none
    integer, parameter :: verb=3
!    logical, save :: initialized = .false.

    
    call debug_message(verb, "theta_grid_params::init_theta_grid_params start")
    if (initialized) return
    initialized = .true.

    call debug_message(verb, &
      "theta_grid_params::init_theta_grid_params call read_parameters")
    call read_parameters

    call debug_message(verb, "theta_grid_params::init_theta_grid_params end")
  end subroutine init_theta_grid_params

  subroutine finish_theta_grid_params
    implicit none
    initialized = .false.
  end subroutine finish_theta_grid_params

  subroutine read_parameters
    use file_utils, only: input_unit, input_unit_exist
    use unit_tests, only: debug_message
    implicit none
    integer, parameter :: verb=3
    integer :: in_file
    character(4) :: ntheta_char

!CMR,2/2/2011: add btor_slab
! btor_slab = btor/bpol defines direction of a flow relative to B in slab 
! geometry, where flow is by definition in the toroidal direction.

    namelist /theta_grid_parameters/ rhoc, rmaj, r_geo, eps, epsl, &
         qinp, shat, alpmhd, pk, shift, akappa, akappri, tri, tripri, &
         ntheta, nperiod, kp, asym, asympri, btor_slab
    
       call debug_message(verb, "theta_grid_params::read_parameters start")


    rhoc = 0.5
    rmaj = 3.0
    r_geo = 3.0
    eps = 0.3
    epsl = 0.3
    qinp = 1.5
    shat = 0.75
    pk = 0.3
    shift = 0.0
    akappa = 1.0
    akappri = 0.0
    tri = 0.0
    tripri = 0.0
    asym = 0.0
    asympri = 0.0
    btor_slab = 0.0
    ntheta = 24
    nperiod = 2
    in_file = input_unit_exist("theta_grid_parameters", exist)
    if (exist) read (unit=in_file, nml=theta_grid_parameters)

    if (kp > 0.) pk = 2.*kp
    
  end subroutine read_parameters

  subroutine wnml_theta_grid_params(unit)
    implicit none
    integer, intent(in) :: unit
    if (.not.exist) return
    write (unit, *)
    write (unit, fmt="(' &',a)") "theta_grid_parameters"
    write (unit, fmt="(' ntheta =  ',i4)") ntheta
    write (unit, fmt="(' nperiod = ',i4)") nperiod
    write (unit, fmt="(' rhoc =    ',f7.4)") rhoc
    write (unit, fmt="(' Rmaj =    ',f7.4)") rmaj
    write (unit, fmt="(' R_geo =   ',f7.4)") r_geo
    write (unit, fmt="(' eps =     ',f7.4)") eps
    write (unit, fmt="(' epsl =    ',f7.4)") epsl
    write (unit, fmt="(' qinp =    ',f7.4)") qinp
    write (unit, fmt="(' shat =    ',f7.4)") shat
    write (unit, fmt="(' alpmhd =  ',f7.4)") alpmhd
    write (unit, fmt="(' pk =      ',f7.4)") pk
    write (unit, fmt="(' kp =      ',f7.4)") kp
    write (unit, fmt="(' shift =   ',f7.4)") shift
    write (unit, fmt="(' akappa =  ',f7.4)") akappa
    write (unit, fmt="(' akappri = ',f7.4)") akappri
    write (unit, fmt="(' tri =     ',f7.4)") tri
    write (unit, fmt="(' tripri =  ',f7.4)") tripri
    write (unit, fmt="(' asym =     ',f7.4)") asym
    write (unit, fmt="(' asympri =  ',f7.4)") asympri
    write (unit, fmt="(' btor_slab =',f7.4)") btor_slab
    write (unit, fmt="(' /')")
  end subroutine wnml_theta_grid_params

  subroutine write_trinity_parameters(trinpars_unit)
    integer, intent(in) :: trinpars_unit
    write (trinpars_unit, "(A22)") '&theta_grid_parameters'
    write (trinpars_unit, *) ' rhoc = ', rhoc
    write (trinpars_unit, *) ' qinp = ', qinp
    write (trinpars_unit, *) ' shat = ', shat
    write (trinpars_unit, *) ' rmaj = ', rmaj
    write (trinpars_unit, *) ' r_geo = ', r_geo
    write (trinpars_unit, *) ' akappa = ', akappa
    write (trinpars_unit, *) ' akappri = ', akappri
    write (trinpars_unit, *) ' tri = ', tri
    write (trinpars_unit, *) ' tripri = ', tripri
    write (trinpars_unit, *) ' shift = ', shift
    write (trinpars_unit, *) ' betaprim = ', betaprim
    write (trinpars_unit, "(A1)") '/'
  end subroutine write_trinity_parameters


  subroutine set_overrides(mgeo_ov)
    use overrides, only: miller_geometry_overrides_type
    type(miller_geometry_overrides_type), intent(in) :: mgeo_ov
          !write (*,*) 'Calling tgpso'
    if (mgeo_ov%override_rhoc) rhoc = mgeo_ov%rhoc
    if (mgeo_ov%override_qinp) qinp = mgeo_ov%qinp
    if (mgeo_ov%override_shat) shat = mgeo_ov%shat
    if (mgeo_ov%override_rgeo_lcfs) r_geo = mgeo_ov%rgeo_lcfs
    if (mgeo_ov%override_rgeo_local) rmaj = mgeo_ov%rgeo_local
    if (mgeo_ov%override_akappa) akappa = mgeo_ov%akappa
    if (mgeo_ov%override_akappri) akappri = mgeo_ov%akappri
    if (mgeo_ov%override_tri) tri = mgeo_ov%tri
    if (mgeo_ov%override_tripri) tripri = mgeo_ov%tripri
    if (mgeo_ov%override_shift) shift = mgeo_ov%shift
    if (mgeo_ov%override_betaprim) betaprim = mgeo_ov%betaprim
  end subroutine set_overrides


end module theta_grid_params
