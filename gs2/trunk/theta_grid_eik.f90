module theta_grid_eik
  implicit none

  public :: init_theta_grid_eik
  public :: eik_get_sizes
  public :: eik_get_grids

  private

contains

  subroutine init_theta_grid_eik
    use geometry, only: init_theta
    use geometry, only: eikcoefs
    use theta_grid_params, only: init_theta_grid_params, ntheta
    implicit none
    logical, save :: initialized = .false.

    if (initialized) return
    initialized = .true.

    call init_theta_grid_params

    call read_parameters
    call init_theta (ntheta)
    call eikcoefs
  end subroutine init_theta_grid_eik

  subroutine eik_get_sizes (nthetaout, nperiodout, nbsetout)
    use geometry, only: nperiod
    use theta_grid_params, only: ntheta
    implicit none
    integer, intent (out) :: nthetaout, nperiodout, nbsetout

    nthetaout = ntheta
    nperiodout = nperiod
    nbsetout = ntheta/2+1 ! upper bound
  end subroutine eik_get_sizes

  subroutine eik_get_grids (nperiod, ntheta, ntgrid, nbset, theta, bset, bmag,&
            gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, gds2, gds21, gds22,&
            grho, shat, drhodpsi)
    use theta_grid_gridgen, only: theta_grid_gridgen_init, gridgen_get_grids
    use geometry, only: theta_out => theta
    use geometry, only: gradpar_out => gradpar
    use geometry, only: bmag_out => bmag
    use geometry, only: cvdrift_out => cvdrift
    use geometry, only: cvdrift0_out => cvdrift0
    use geometry, only: gbdrift_out => gbdrift
    use geometry, only: gbdrift0_out => gbdrift0
    use geometry, only: gds2_out => gds2
    use geometry, only: gds21_out => gds21
    use geometry, only: gds22_out => gds22
    use geometry, only: grho_out => grho
    use geometry, only: s_hat_input, drhodpsin
    implicit none
    integer, intent (in) :: nperiod
    integer, intent (in out) :: ntheta, ntgrid, nbset
    real, dimension (-ntgrid:ntgrid), intent (out) :: theta
    real, dimension (nbset), intent (out) :: bset
    real, dimension (-ntgrid:ntgrid), intent (out) :: &
         bmag, gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, &
         gds2, gds21, gds22, grho
    real, intent (out) :: shat, drhodpsi

    theta(-ntgrid:ntgrid) = theta_out(-ntgrid:ntgrid)
    gradpar(-ntgrid:ntgrid) = gradpar_out(-ntgrid:ntgrid)
    bmag(-ntgrid:ntgrid) = bmag_out(-ntgrid:ntgrid)
    cvdrift(-ntgrid:ntgrid) = cvdrift_out(-ntgrid:ntgrid)
    cvdrift0(-ntgrid:ntgrid) = cvdrift0_out(-ntgrid:ntgrid)
    gbdrift(-ntgrid:ntgrid) = gbdrift_out(-ntgrid:ntgrid)
    gbdrift0(-ntgrid:ntgrid) = gbdrift0_out(-ntgrid:ntgrid)
    gds2(-ntgrid:ntgrid) = gds2_out(-ntgrid:ntgrid)
    gds21(-ntgrid:ntgrid) = gds21_out(-ntgrid:ntgrid)
    gds22(-ntgrid:ntgrid) = gds22_out(-ntgrid:ntgrid)
    grho(-ntgrid:ntgrid) = grho_out(-ntgrid:ntgrid)

    call theta_grid_gridgen_init
    call gridgen_get_grids (nperiod, ntheta, ntgrid, nbset, &
         theta, bset, bmag, &
         gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, gds2, gds21, gds22, &
         grho)
    shat = s_hat_input
    drhodpsi = drhodpsin
  end subroutine eik_get_grids

  subroutine read_parameters
    use file_utils, only: input_unit, error_unit
    use geometry, only: nperiod
    use geometry, only: rhoc
    use geometry, only: itor, iflux, irho
    use geometry, only: ppl_eq, gen_eq, vmom_eq, efit_eq, eqfile, local_eq
    use geometry, only: equal_arc
    use geometry, only: bishop
    use geometry, only: s_hat_input
    use geometry, only: alpha_input, invLp_input, beta_prime_input, dp_mult
    use geometry, only: rmaj, r_geo
    use geometry, only: shift, qinp, akappa, akappri, tri, tripri
    use geometry, only: delrho, rmin, rmax
    use geometry, only: ismooth, ak0, k1, k2
    use geometry, only: isym, in_nt, writelots
    use theta_grid_params, only: ntheta, nperiod_in => nperiod
    use theta_grid_params, only: rhoc_in => rhoc
    use theta_grid_params, only: rmaj_in => rmaj, r_geo_in => r_geo
    use theta_grid_params, only: eps_in => eps, epsl_in => epsl
    use theta_grid_params, only: qinp_in => qinp, shat_in => shat
    use theta_grid_params, only: alpmhd_in => alpmhd
    use theta_grid_params, only: shift_in => shift
    use theta_grid_params, only: akappa_in => akappa, akappri_in => akappri
    use theta_grid_params, only: tri_in => tri, tripri_in => tripri
    implicit none

    namelist /theta_grid_eik_knobs/ itor, iflux, irho, &
         ppl_eq, gen_eq, vmom_eq, efit_eq, eqfile, &
         equal_arc, bishop, local_eq, &
         s_hat_input, alpha_input, invLp_input, beta_prime_input, dp_mult, &
         delrho, rmin, rmax, ismooth, ak0, k1, k2, isym, writelots

    nperiod = nperiod_in
    rhoc = rhoc_in
    s_hat_input = shat_in
    alpha_input = alpmhd_in
    rmaj = rmaj_in
    r_geo = r_geo_in
    shift = shift_in
    qinp = qinp_in
    akappa = akappa_in
    akappri = akappri_in
    tri = tri_in
    tripri = tripri_in

    itor = 1
    iflux = 0
    irho = 2
    equal_arc = .true.
    bishop = 5
    dp_mult = 1.0
    delrho = 1e-3
    rmin = 1e-3
    rmax = 1.0
    ismooth = 0
    isym = 0
    in_nt = .false.
    writelots = .false.
    local_eq = .true.

    read (unit=input_unit("theta_grid_eik_knobs"), nml=theta_grid_eik_knobs)
  end subroutine read_parameters
end module theta_grid_eik
