module gs2_io

  implicit none

  private

  public :: init_gs2_io, nc_eigenfunc, nc_final_fields, nc_final_epar
  public :: nc_final_moments, nc_final_an, nc_finish
  public :: nc_qflux, nc_vflux, nc_pflux, nc_loop, nc_loop_moments
  public :: nc_loop_movie, nc_write_fields

  logical, parameter :: serial_io = .true.
  integer :: ncid

  integer :: naky_dim, nakx_dim, nttot_dim, negrid_dim, nlambda_dim, nspec_dim
  integer :: nsign_dim, time_dim, char10_dim, char200_dim, ri_dim, nlines_dim, nheat_dim
  
  integer, dimension (5) :: field_dim, final_mom_dim, heatk_dim
  integer, dimension (4) :: omega_dim, fluxk_dim, final_field_dim, loop_mom_dim
  integer, dimension (3) :: mode_dim, phase_dim, loop_phi_dim, heat_dim
  integer, dimension (2) :: kx_dim, ky_dim, om_dim, flux_dim, nin_dim, fmode_dim

  ! added by EAB 03/05/04 for movies
  logical :: my_make_movie
  integer :: ncid_movie
  integer :: nx_dim, ny_dim, nth_dim, time_movie_dim
  integer, dimension (4) :: xmode_dim
  integer :: nx_id, ny_id, nth_id, x_id, y_id, th_id, time_movie_id
  integer :: phi_by_xmode_id, apar_by_xmode_id, bpar_by_xmode_id
  integer :: code_id_movie

  integer :: nakx_id, naky_id, nttot_id, akx_id, aky_id, theta_id, nspec_id
  integer :: time_id, phi2_id, apar2_id, bpar2_id, theta0_id, nproc_id, nmesh_id
  integer :: phi2_by_mode_id, apar2_by_mode_id, bpar2_by_mode_id, anorm_id
  integer :: current_scan_parameter_value_id
  integer :: phtot_id, dmix_id, kperpnorm_id
  integer :: phi2_by_kx_id, apar2_by_kx_id, bpar2_by_kx_id
  integer :: phi2_by_ky_id, apar2_by_ky_id, bpar2_by_ky_id
  integer :: phi0_id, apar0_id, bpar0_id, sourcefac_id
  integer :: omega_id, omegaavg_id, phase_id, phiavg_id
  integer :: es_heat_flux_id, es_mom_flux_id, es_part_flux_id
  integer :: es_heat_par_id, es_heat_perp_id
  integer :: apar_heat_flux_id, apar_mom_flux_id, apar_part_flux_id
  integer :: apar_heat_par_id, apar_heat_perp_id
  integer :: bpar_heat_flux_id, bpar_mom_flux_id, bpar_part_flux_id
  integer :: bpar_heat_par_id, bpar_heat_perp_id
  integer :: hflux_tot_id, zflux_tot_id, vflux_tot_id
  integer :: es_heat_by_k_id, es_mom_by_k_id, es_part_by_k_id
  integer :: apar_heat_by_k_id, apar_mom_by_k_id, apar_part_by_k_id
  integer :: bpar_heat_by_k_id, bpar_mom_by_k_id, bpar_part_by_k_id
  integer :: phi_t_id, apar_t_id, bpar_t_id
  integer :: phi_norm_id, apar_norm_id, bpar_norm_id
  integer :: phi_id, apar_id, bpar_id, epar_id
  integer :: antot_id, antota_id, antotp_id
  integer :: ntot_id, density_id, upar_id, tpar_id, tperp_id
  integer :: qparflux_id, pperpj1_id, qpperpj1_id
  integer :: hrateavg_id, hrate_by_k_id
  integer :: ntot2_id, ntot2_by_mode_id
  integer :: phi00_id, ntot00_id, density00_id, upar00_id, tpar00_id, tperp00_id
  integer :: input_id
  integer :: charge_id, mass_id, dens_id, temp_id, tprim_id, fprim_id
  integer :: uprim_id, uprim2_id, vnewk_id, spec_type_id
  integer :: bmag_id, gradpar_id, gbdrift_id, gbdrift0_id
  integer :: cvdrift_id, cvdrift0_id, gds2_id, gds21_id, gds22_id
  integer :: grho_id, jacob_id, shat_id, eps_id, drhodpsi_id, q_id
  integer :: code_id, datestamp_id, timestamp_id, timezone_id

  real :: zero
  
  include 'netcdf.inc'
  
contains

  subroutine init_gs2_io (write_nl_flux, write_omega, &
      write_phiavg, write_hrate, make_movie, nmovie_tot, write_fields)
!David has made some changes to this subroutine (may 2005) now should be able to do movies for 
!linear box runs as well as nonlinear box runs.

    use mp, only: proc0, barrier
    use file_utils, only: run_name
    use netcdf_mod
    use gs2_transforms, only: init_transforms
    use kt_grids, only: naky, ntheta0,nx,ny
    use gs2_layouts, only: yxf_lo
    use theta_grid, only: ntgrid
    use le_grids, only: nlambda, negrid
    use species, only: nspec

    logical :: write_nl_flux, write_omega, write_phiavg, write_hrate, make_movie, write_fields
    logical :: accelerated
    character (300) :: filename, filename_movie
    integer :: ierr         ! 0 if initialization is successful
    integer :: nmovie_tot
    integer :: status


    zero = epsilon(0.0)
    ierr = 0

    call netcdf_init(serial_io)

    status = nf_noerr

    filename = trim(trim(run_name)//'.out.nc')
    if(serial_io) then
       ! only proc0 opens the file:
       if (proc0) &
            status = nf_create(trim(filename), 0, ncid) 
    else
       ! all processors open the file
       call barrier
       status = nf_create(trim(filename), 0, ncid) ! overwrites old
       call barrier
    endif

    if(status /= nf_noerr) then
       write(*,*) 'Unable to create netcdf file ', filename
       write(*,*) nf_strerror(status)
       ierr = 1
    endif

    ! added by EAB 03/2004 for making movie netcdf files
    my_make_movie = make_movie
    if(my_make_movie) then
    !perhaps put a call to init_y_transform
!    call init_y_transform_layouts  
    call init_transforms(ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx, ny, accelerated)
       filename_movie = trim(trim(run_name)//'.movie.nc')
       if(serial_io) then
          ! only proc0 opens the file:
          if (proc0) &
               status = nf_create(trim(filename_movie), 0, ncid_movie) 
       else
          ! all processors open the file
          call barrier
          status = nf_create(trim(filename_movie), 0, ncid_movie) ! overwrites old
          call barrier
       endif
       
       if(status /= nf_noerr) then
          write(*,*) 'Unable to create netcdf file for movies', filename_movie
          write(*,*) nf_strerror(status)
          ierr = 1
       endif
    endif

    if (proc0) then
       call define_dims (nmovie_tot)
       call define_vars (write_nl_flux, write_omega, write_phiavg, write_hrate, write_fields)
       call nc_grids
       call nc_species
       call nc_geo
    end if

  end subroutine init_gs2_io

  subroutine define_dims (nmovie_tot)

    use file_utils, only: num_input_lines
    use kt_grids, only: naky, ntheta0,nx,ny
    use gs2_layouts, only: yxf_lo    
    use theta_grid, only: ntgrid
    use le_grids, only: nlambda, negrid
    use species, only: nspec
    use netcdf_mod

    integer, intent (in) :: nmovie_tot
    integer :: status


    status = netcdf_def_dim(ncid, 'ky', naky, naky_dim)
    status = netcdf_def_dim(ncid, 'kx', ntheta0, nakx_dim)
    status = netcdf_def_dim(ncid, 'theta', 2*ntgrid+1, nttot_dim)
    status = netcdf_def_dim(ncid, 'lambda', nlambda, nlambda_dim)
    status = netcdf_def_dim(ncid, 'egrid', negrid, negrid_dim)
    status = netcdf_def_dim(ncid, 'species', nspec, nspec_dim)
    status = netcdf_def_dim(ncid, 'sign', 2, nsign_dim)
    status = netcdf_def_dim(ncid, 't', NF_UNLIMITED, time_dim)
    status = netcdf_def_dim(ncid, 'char10', 10, char10_dim)
    status = netcdf_def_dim(ncid, 'char200', 200, char200_dim)
    status = netcdf_def_dim(ncid, 'nlines', num_input_lines, nlines_dim)
    status = netcdf_def_dim(ncid, 'ri', 2, ri_dim)
    status = netcdf_def_dim(ncid, 'nheat', 3, nheat_dim)

    ! added by EAB 03/05/04 for GS2 movies
    if(my_make_movie) then



       
       nx=yxf_lo%nx
       ny=yxf_lo%ny
       status = netcdf_def_dim(ncid_movie, 'x', nx, nx_dim)
       status = netcdf_def_dim(ncid_movie, 'y', ny, ny_dim)
       status = netcdf_def_dim(ncid_movie, 'theta', 2*ntgrid+1, nth_dim)
!
! can only have one NF_UNLIMITED, so:
!
!       status = netcdf_def_dim(ncid_movie, 't', NF_UNLIMITED, time_movie_dim)
       status = netcdf_def_dim(ncid_movie, 'tm', nmovie_tot, time_movie_dim)
    endif

  end subroutine define_dims

  subroutine nc_grids

    use theta_grid, only: ntgrid, theta, eps, ntheta
    use kt_grids, only: naky, ntheta0, theta0, akx, aky, nx, ny
    use gs2_layouts, only: yxf_lo
    use species, only: nspec
    use le_grids, only: negrid, nlambda
    use nonlinear_terms, only: nonlin
    use netcdf_mod
    use constants

    integer :: status
    real :: nmesh

    real, dimension(:), allocatable :: x, y
    integer :: ik, it


    status = netcdf_put_var (ncid, nttot_id, 2*ntgrid+1)
    status = netcdf_put_var (ncid, naky_id, naky)
    status = netcdf_put_var (ncid, nakx_id, ntheta0)
    status = netcdf_put_var (ncid, nspec_id, nspec)

    status = netcdf_put_var (ncid, akx_id, akx)
    status = netcdf_put_var (ncid, aky_id, aky)
    status = netcdf_put_var (ncid, theta_id, theta)
    status = netcdf_put_var (ncid, theta0_id, theta0)

    if (nonlin) then
       nmesh = (2*ntgrid+1)*2*nlambda*negrid*nx*ny*nspec
    else
       nmesh = (2*ntgrid+1)*2*nlambda*negrid*ntheta0*naky*nspec
    end if

    status = netcdf_put_var (ncid, nmesh_id, nmesh)

    ! added by EAB 03/05/04 for GS2 movies
    if(my_make_movie) then

       allocate(x(yxf_lo%nx))
       allocate(y(yxf_lo%ny))
       do it = 1, yxf_lo%nx 
          x(it) = 2.0*pi/akx(2)*(-0.5+ real(it-1)/real(yxf_lo%nx))
       end do
       do ik = 1, yxf_lo%ny
          y(ik) = 2.0*pi/aky(2)*(-0.5+ real(ik-1)/real(yxf_lo%ny))
       end do
       
       status = netcdf_put_var (ncid_movie, nx_id, yxf_lo%nx)
       status = netcdf_put_var (ncid_movie, ny_id, yxf_lo%ny)
       status = netcdf_put_var (ncid_movie, nth_id, 2*ntgrid+1)
       status = netcdf_put_var (ncid_movie, x_id, x)
       status = netcdf_put_var (ncid_movie, y_id, y)
       status = netcdf_put_var (ncid_movie, th_id, theta)
       deallocate(x)
       deallocate(y)
    endif

  end subroutine nc_grids

  subroutine nc_finish

    use mp, only: proc0
    integer :: status

    if (proc0) then
       call save_input
       status = nf_close (ncid)
       if(my_make_movie) then
          status = nf_close (ncid_movie)
       endif
    end if

  end subroutine nc_finish

  subroutine save_input
    
    use file_utils, only: num_input_lines, get_input_unit
    use netcdf_mod

    character(200) line
    integer, dimension (2) :: nin_start, nin_count

    integer :: status, n, unit, iostat

    nin_start(1) = 1
    nin_start(2) = 1

    nin_count(2) = 1

    call get_input_unit (unit)
    rewind (unit=unit)
    do n = 1, num_input_lines
       read (unit=unit, fmt="(a)") line
       nin_count(1) = len(trim(line))
       status = nf_put_vara_text (ncid, input_id, nin_start, nin_count, line)
       nin_start(2) = nin_start(2) + 1
    end do

  end subroutine save_input

  subroutine define_vars (write_nl_flux, write_omega, write_phiavg, write_hrate, write_fields)

    use mp, only: nproc
    use species, only: nspec
    use kt_grids, only: naky, ntheta0, theta0
    use run_parameters, only: fphi, fapar, fbpar
    use parameter_scan_arrays, only: write_scan_parameter
    use netcdf_mod
    logical :: write_nl_flux, write_omega, write_phiavg, write_hrate, write_fields

    character (5) :: ci
    character (20) :: datestamp, timestamp, timezone
    logical :: d_fields_per = .false.
    
    integer :: status

    fmode_dim(1) = nakx_dim
    fmode_dim(2) = naky_dim

    mode_dim (1) = nakx_dim
    mode_dim (2) = naky_dim
    mode_dim (3) = time_dim

    kx_dim (1) = nakx_dim
    kx_dim (2) = time_dim
    
    ky_dim (1) = naky_dim
    ky_dim (2) = time_dim
    
    om_dim (1) = ri_dim
    om_dim (2) = time_dim

    omega_dim (1) = ri_dim
    omega_dim (2) = nakx_dim
    omega_dim (3) = naky_dim
    omega_dim (4) = time_dim

    phase_dim (1) = ri_dim
    phase_dim (2) = nakx_dim
    phase_dim (3) = naky_dim
    
    nin_dim(1) = char200_dim
    nin_dim(2) = nlines_dim
    
    flux_dim (1) = nspec_dim
    flux_dim (2) = time_dim

    fluxk_dim (1) = nakx_dim
    fluxk_dim (2) = naky_dim
    fluxk_dim (3) = nspec_dim
    fluxk_dim (4) = time_dim

    heat_dim (1) = nspec_dim
    heat_dim (2) = nheat_dim
    heat_dim (3) = time_dim

    heatk_dim (1) = nakx_dim
    heatk_dim (2) = naky_dim
    heatk_dim (3) = nspec_dim
    heatk_dim (4) = nheat_dim
    heatk_dim (5) = time_dim

    field_dim (1) = ri_dim
    field_dim (2) = nttot_dim
    field_dim (3) = nakx_dim
    field_dim (4) = naky_dim
    field_dim (5) = time_dim
    
    final_field_dim (1) = ri_dim
    final_field_dim (2) = nttot_dim
    final_field_dim (3) = nakx_dim
    final_field_dim (4) = naky_dim

    final_mom_dim (1) = ri_dim
    final_mom_dim (2) = nttot_dim
    final_mom_dim (3) = nakx_dim
    final_mom_dim (4) = naky_dim
    final_mom_dim (5) = nspec_dim

    loop_mom_dim (1) = ri_dim
    loop_mom_dim (2) = nakx_dim
    loop_mom_dim (3) = nspec_dim
    loop_mom_dim (4) = time_dim

    loop_phi_dim (1) = ri_dim
    loop_phi_dim (2) = nakx_dim
    loop_phi_dim (3) = time_dim

    status = netcdf_put_att (ncid, nf_global, 'title', 'GS2 Simulation Data')
    status = netcdf_put_att (ncid, nf_global, 'Conventions', &
         'http://gs2.sourceforge.net')

    datestamp(:) = ' '
    timestamp(:) = ' '
    timezone(:) = ' '
    call date_and_time (datestamp, timestamp, timezone)
    
    status = netcdf_def_var (ncid, 'code_info', nf_char, 1, char10_dim, code_id)
    status = netcdf_put_att (ncid, code_id, 'long_name', 'GS2')

    ci = 'c1'
    status = netcdf_put_att (ncid, code_id, trim(ci), 'Date: '//trim(datestamp))

    ci = 'c2'
    status = netcdf_put_att (ncid, code_id, trim(ci), &
         'Time: '//trim(timestamp)//' '//trim(timezone))

    ci = 'c3'
    status = netcdf_put_att (ncid, code_id, trim(ci), &
         'netCDF version '//trim(NF_INQ_LIBVERS()))

    ci = 'c4'
    status = netcdf_put_att (ncid, code_id, trim(ci), &
         'Units are determined with respect to reference temperature (T_ref),')

    ci = 'c5'
    status = netcdf_put_att (ncid, code_id, trim(ci), &
         'reference charge (q_ref), reference mass (mass_ref),')

    ci = 'c6'
    status = netcdf_put_att (ncid, code_id, trim(ci), &
         'reference field (B_ref), and reference length (a_ref)')

    ci = 'c7'
    status = netcdf_put_att (ncid, code_id, trim(ci), &
         'from which one may construct rho_ref and vt_ref/a,')

    ci = 'c8'
    status = netcdf_put_att (ncid, code_id, trim(ci), &
         'which are the basic units of perpendicular length and time.')

    ci = 'c9'
    status = netcdf_put_att (ncid, code_id, trim(ci), &
         'Macroscopic lengths are normalized to the minor radius.')

    ci = 'c10'
    status = netcdf_put_att (ncid, code_id, trim(ci), &
         'The difference between rho (normalized minor radius) and rho (gyroradius)')

    ci = 'c11'
    status = netcdf_put_att (ncid, code_id, trim(ci), &
         'should be clear from the context in which they appear below.')

    status = netcdf_def_var (ncid, 'nproc', nf_int, 0, 0, nproc_id)
    status = netcdf_put_att (ncid, nproc_id, 'long_name', 'Number of processors')

    status = netcdf_def_var (ncid, 'nmesh', nf_double, 0, 0, nmesh_id)
    status = netcdf_put_att (ncid, nmesh_id, 'long_name', 'Number of meshpoints')

    status = netcdf_def_var (ncid, 'nkx', nf_int, 0, 0, nakx_id)
    status = netcdf_def_var (ncid, 'nky', nf_int, 0, 0, naky_id)
    status = netcdf_def_var (ncid, 'ntheta_tot', nf_int, 0, 0, nttot_id)
    status = netcdf_def_var (ncid, 'nspecies', nf_int, 0, 0, nspec_id)

    status = netcdf_def_var (ncid, 't', nf_double, 1, time_dim, time_id)
    status = netcdf_put_att (ncid, time_id, 'long_name', 'Time')
    status = netcdf_put_att (ncid, time_id, 'units', 'L/vt')

    status = netcdf_def_var (ncid, 'charge', nf_int, 1, nspec_dim, charge_id)
    status = netcdf_put_att (ncid, charge_id, 'long_name', 'Charge')
    status = netcdf_put_att (ncid, charge_id, 'units', 'q')

    status = netcdf_def_var (ncid, 'mass', nf_double, 1, nspec_dim, mass_id)
    status = netcdf_put_att (ncid, mass_id, 'long_name', 'Atomic mass')
    status = netcdf_put_att (ncid, mass_id, 'units', 'm')

    status = netcdf_def_var (ncid, 'dens', nf_double, 1, nspec_dim, dens_id)
    status = netcdf_put_att (ncid, dens_id, 'long_name', 'Density')
    status = netcdf_put_att (ncid, dens_id, 'units', 'n_e')

    status = netcdf_def_var (ncid, 'temp', nf_double, 1, nspec_dim, temp_id)
    status = netcdf_put_att (ncid, temp_id, 'long_name', 'Temperature')
    status = netcdf_put_att (ncid, temp_id, 'units', 'T')

    status = netcdf_def_var (ncid, 'tprim', nf_double, 1, nspec_dim, tprim_id)
    status = netcdf_put_att (ncid, tprim_id, 'long_name', '-1/rho dT/drho')

    status = netcdf_def_var (ncid, 'fprim', nf_double, 1, nspec_dim, fprim_id) 
    status = netcdf_put_att (ncid, fprim_id, 'long_name', '-1/rho dn/drho')

    status = netcdf_def_var (ncid, 'uprim', nf_double, 1, nspec_dim, uprim_id)
    status = netcdf_put_att (ncid, uprim_id, 'long_name', '-1/v_t du_par/drho')

    status = netcdf_def_var (ncid, 'uprim2', nf_double, 1, nspec_dim, uprim2_id)

    status = netcdf_def_var (ncid, 'vnewk', nf_double, 1, nspec_dim, vnewk_id)
    status = netcdf_put_att (ncid, vnewk_id, 'long_name', 'Collisionality')
    status = netcdf_put_att (ncid, vnewk_id, 'units', 'v_t/L')
    
    status = netcdf_def_var (ncid, 'type_of_species', nf_int, 1, nspec_dim, spec_type_id)

    status = netcdf_def_var (ncid, 'theta0', nf_double, 2, fmode_dim, theta0_id)
    status = netcdf_put_att (ncid, theta0_id, 'long_name', 'Theta_0')

    status = netcdf_def_var (ncid, 'kx', nf_double, 1, nakx_dim, akx_id)
    status = netcdf_put_att (ncid, akx_id, 'long_name', 'kx rho')

    status = netcdf_def_var (ncid, 'ky', nf_double, 1, naky_dim, aky_id)
    status = netcdf_put_att (ncid, aky_id, 'long_name', 'ky rho')

    status = netcdf_def_var (ncid, 'theta', nf_double, 1, nttot_dim, theta_id)

    status = netcdf_def_var (ncid, 'bmag', nf_double, 1, nttot_dim, bmag_id)
    status = netcdf_put_att (ncid, bmag_id, 'long_name', '|B|(theta)')
    status = netcdf_put_att (ncid, bmag_id, 'units', 'B_0')

    status = netcdf_def_var (ncid, 'gradpar', nf_double, 1, nttot_dim, gradpar_id)
    status = netcdf_def_var (ncid, 'gbdrift', nf_double, 1, nttot_dim, gbdrift_id) 
    status = netcdf_def_var (ncid, 'gbdrift0', nf_double, 1, nttot_dim, gbdrift0_id)
    status = netcdf_def_var (ncid, 'cvdrift', nf_double, 1, nttot_dim, cvdrift_id)
    status = netcdf_def_var (ncid, 'cvdrift0', nf_double, 1, nttot_dim, cvdrift0_id)
    status = netcdf_def_var (ncid, 'gds2', nf_double, 1, nttot_dim, gds2_id)
    status = netcdf_def_var (ncid, 'gds21', nf_double, 1, nttot_dim, gds21_id)
    status = netcdf_def_var (ncid, 'gds22', nf_double, 1, nttot_dim, gds22_id)
    status = netcdf_def_var (ncid, 'grho', nf_double, 1, nttot_dim, grho_id)
    status = netcdf_def_var (ncid, 'jacob', nf_double, 1, nttot_dim, jacob_id)

    status = netcdf_def_var (ncid, 'q', nf_double, 0, 0, q_id)
    status = netcdf_def_var (ncid, 'eps', nf_double, 0, 0, eps_id)
    status = netcdf_def_var (ncid, 'shat', nf_double, 0, 0, shat_id)
    status = netcdf_put_att (ncid, shat_id, 'long_name', '(rho/q) dq/drho')

    status = netcdf_def_var (ncid, 'drhodpsi', nf_double, 0, 0, drhodpsi_id)
    status = netcdf_put_att (ncid, drhodpsi_id, 'long_name', 'drho/dPsi')

    if (write_scan_parameter) then
       status = netcdf_def_var (ncid, 'scan_parameter_value', &
                                nf_double, 1, time_dim, &
                                current_scan_parameter_value_id)
       status = netcdf_put_att (ncid, &
                         current_scan_parameter_value_id, &
                        'long_name',&
                      'The current value of the scan parameter')
     end if

    if (fphi > zero) then
       status = netcdf_def_var (ncid, 'phi2',         nf_double, 1, time_dim, phi2_id)
       status = netcdf_put_att (ncid, phi2_id, 'long_name', '|Potential**2|')
       status = netcdf_put_att (ncid, phi2_id, 'units', '(T/q rho/L)**2')

       status = netcdf_def_var (ncid, 'phi2_by_mode', nf_double, 3, mode_dim, phi2_by_mode_id)
       if (ntheta0 > 1) &
       status = netcdf_def_var (ncid, 'phi2_by_kx',   nf_double, 2, kx_dim, phi2_by_kx_id)
       if (naky > 1) &
       status = netcdf_def_var (ncid, 'phi2_by_ky',   nf_double, 2, ky_dim, phi2_by_ky_id)
       status = netcdf_def_var (ncid, 'phi0',         nf_double, 4, omega_dim, phi0_id)
       if (write_phiavg) status = netcdf_def_var (ncid, 'phiavg',       nf_double, 4, omega_dim, phiavg_id)
       if (write_nl_flux) then
          status = netcdf_def_var (ncid, 'es_heat_par',  nf_double, 2, flux_dim, es_heat_par_id)
          status = netcdf_def_var (ncid, 'es_heat_perp', nf_double, 2, flux_dim, es_heat_perp_id)
          status = netcdf_def_var (ncid, 'es_heat_flux', nf_double, 2, flux_dim, es_heat_flux_id)
          status = netcdf_def_var (ncid, 'es_mom_flux',  nf_double, 2, flux_dim, es_mom_flux_id)
          status = netcdf_def_var (ncid, 'es_part_flux', nf_double, 2, flux_dim, es_part_flux_id)
          status = netcdf_def_var (ncid, 'es_heat_by_k', nf_double, 4, fluxk_dim, es_heat_by_k_id)
          status = netcdf_def_var (ncid, 'es_mom_by_k',  nf_double, 4, fluxk_dim, es_mom_by_k_id)
          status = netcdf_def_var (ncid, 'es_part_by_k', nf_double, 4, fluxk_dim, es_part_by_k_id)
       end if
       if (write_fields) status = netcdf_def_var (ncid, 'phi_t', nf_double, 5, field_dim, phi_t_id)  !MR
       status = netcdf_def_var (ncid, 'phi_norm',     nf_double, 4, final_field_dim, phi_norm_id)
       status = netcdf_def_var (ncid, 'phi',          nf_double, 4, final_field_dim, phi_id)
       status = netcdf_put_att (ncid, phi_id, 'long_name', 'Electrostatic Potential')
       status = netcdf_put_att (ncid, phi_id, 'idl_name', '!7U!6')
       status = netcdf_put_att (ncid, phi_id, 'units', 'T/q rho/L')
       status = netcdf_def_var (ncid, 'antot',        nf_double, 4, final_field_dim, antot_id)
!       if (d_fields_per) status = netcdf_def_var (ncid, 'phi_t', nf_double, 5, field_dim, phi_t_id)
    end if

    if (fapar > zero) then
       status = netcdf_def_var (ncid, 'apar2',          nf_double, 1, time_dim, apar2_id)
       status = netcdf_def_var (ncid, 'apar2_by_mode',  nf_double, 3, mode_dim, apar2_by_mode_id)
       if (ntheta0 > 1) &
       status = netcdf_def_var (ncid, 'apar2_by_kx',    nf_double, 2, kx_dim, apar2_by_kx_id)
       if (naky > 1) &
       status = netcdf_def_var (ncid, 'apar2_by_ky',    nf_double, 2, ky_dim, apar2_by_ky_id)
       status = netcdf_def_var (ncid, 'apar0',          nf_double, 4, omega_dim, apar0_id)
       if (write_nl_flux) then
          status = netcdf_def_var (ncid, 'apar_heat_flux', nf_double, 2, flux_dim, apar_heat_flux_id)
          status = netcdf_def_var (ncid, 'apar_heat_par',  nf_double, 2, flux_dim, apar_heat_par_id)
          status = netcdf_def_var (ncid, 'apar_heat_perp', nf_double, 2, flux_dim, apar_heat_perp_id)
          status = netcdf_def_var (ncid, 'apar_mom_flux',  nf_double, 2, flux_dim, apar_mom_flux_id)
          status = netcdf_def_var (ncid, 'apar_part_flux', nf_double, 2, flux_dim, apar_part_flux_id)
          status = netcdf_def_var (ncid, 'apar_heat_by_k', nf_double, 4, fluxk_dim, apar_heat_by_k_id)
          status = netcdf_def_var (ncid, 'apar_mom_by_k',  nf_double, 4, fluxk_dim, apar_mom_by_k_id)
          status = netcdf_def_var (ncid, 'apar_part_by_k', nf_double, 4, fluxk_dim, apar_part_by_k_id)
       end if
       status = netcdf_def_var (ncid, 'apar_norm',      nf_double, 4, final_field_dim, apar_norm_id)
       status = netcdf_def_var (ncid, 'apar',           nf_double, 4, final_field_dim, apar_id)
       status = netcdf_def_var (ncid, 'antota',         nf_double, 4, final_field_dim, antota_id)
!       if (d_fields_per) status = netcdf_def_var (ncid, 'apar_t',  nf_double, 5, field_dim, apar_t_id)
       status = netcdf_put_att (ncid, apar2_by_mode_id, 'long_name', 'Apar squared')
       status = netcdf_put_att (ncid, apar_id, 'long_name', 'Parallel Magnetic Potential')      
       status = netcdf_put_att (ncid, apar_id, 'idl_name', '!6A!9!D#!N!6')      
       status = netcdf_put_att (ncid, apar2_id, 'long_name', 'Total A_par squared')
    end if

    if (fbpar > zero) then
       status = netcdf_def_var (ncid, 'bpar2',          nf_double, 1, time_dim, bpar2_id)
       status = netcdf_def_var (ncid, 'bpar2_by_mode',  nf_double, 3, mode_dim, bpar2_by_mode_id)
       if (ntheta0 > 1) &
       status = netcdf_def_var (ncid, 'bpar2_by_kx',    nf_double, 2, kx_dim, bpar2_by_kx_id)
       if (naky > 1) &
       status = netcdf_def_var (ncid, 'bpar2_by_ky',    nf_double, 2, ky_dim, bpar2_by_ky_id)
       status = netcdf_def_var (ncid, 'bpar0',          nf_double, 4, omega_dim, bpar0_id)
       if (write_nl_flux) then
          status = netcdf_def_var (ncid, 'bpar_heat_flux', nf_double, 2, flux_dim, bpar_heat_flux_id)
          status = netcdf_def_var (ncid, 'bpar_heat_par',  nf_double, 2, flux_dim, bpar_heat_par_id)
          status = netcdf_def_var (ncid, 'bpar_heat_perp', nf_double, 2, flux_dim, bpar_heat_perp_id)
          status = netcdf_def_var (ncid, 'bpar_mom_flux',  nf_double, 2, flux_dim, bpar_mom_flux_id)
          status = netcdf_def_var (ncid, 'bpar_part_flux', nf_double, 2, flux_dim, bpar_part_flux_id)
          status = netcdf_def_var (ncid, 'bpar_heat_by_k', nf_double, 4, fluxk_dim, bpar_heat_by_k_id)
          status = netcdf_def_var (ncid, 'bpar_mom_by_k',  nf_double, 4, fluxk_dim, bpar_mom_by_k_id)
          status = netcdf_def_var (ncid, 'bpar_part_by_k', nf_double, 4, fluxk_dim, bpar_part_by_k_id)
       end if
       status = netcdf_def_var (ncid, 'bpar_norm',      nf_double, 4, final_field_dim, bpar_norm_id)
       status = netcdf_def_var (ncid, 'bpar',           nf_double, 4, final_field_dim, bpar_id)
       status = netcdf_def_var (ncid, 'antotp',          nf_double, 4, final_field_dim, antotp_id)
!       if (d_fields_per) status = netcdf_def_var (ncid, 'bpar_t', nf_double, 5, field_dim, bpar_t_id)
       status = netcdf_put_att (ncid, bpar2_by_mode_id, 'long_name', 'A_perp squared')
       status = netcdf_put_att (ncid, bpar_id, 'long_name', 'delta B Parallel')      
       status = netcdf_put_att (ncid, bpar_id, 'idl_name', '!6B!9!D#!N!6')          
       status = netcdf_put_att (ncid, bpar2_id, 'long_name', 'Total A_perp squared')
    end if

    status = netcdf_def_var (ncid, 'anorm', nf_double, 1, time_dim, anorm_id)
    status = netcdf_put_att (ncid, anorm_id, 'long_name', 'Normalizing factor')

    status = netcdf_def_var (ncid, 'phase', nf_double, 3, phase_dim, phase_id)
    status = netcdf_put_att (ncid, phase_id, 'long_name', 'Normalizing phase')

!    status = netcdf_def_var (ncid, 'phtot', nf_double, 3, mode_dim, phtot_id)
!    status = netcdf_def_var (ncid, 'dmix',  nf_double, 3, mode_dim, dmix_id)
!    status = netcdf_def_var (ncid, 'kperpnorm', nf_double, 3, mode_dim, kperpnorm_id)

    if (write_omega) then
       status = netcdf_def_var (ncid, 'omega',   nf_double, 4, omega_dim, omega_id)
       status = netcdf_def_var (ncid, 'omegaavg',   nf_double, 4, omega_dim, omegaavg_id)
    end if

    if (write_nl_flux) then
       status = netcdf_def_var (ncid, 'hflux_tot', nf_double, 1, time_dim, hflux_tot_id)
       status = netcdf_def_var (ncid, 'vflux_tot', nf_double, 1, time_dim, vflux_tot_id)
       status = netcdf_def_var (ncid, 'zflux_tot', nf_double, 1, time_dim, zflux_tot_id)
    end if

    status = netcdf_def_var (ncid, 'epar',  nf_double, 4, final_field_dim, epar_id)

!    <phi>

    status = netcdf_def_var (ncid, 'ntot2',         nf_double, 2, flux_dim,  ntot2_id)
    status = netcdf_def_var (ncid, 'ntot2_by_mode', nf_double, 4, fluxk_dim, ntot2_by_mode_id)

    status = netcdf_def_var (ncid, 'ntot',    nf_double, 5, final_mom_dim, ntot_id)
    status = netcdf_def_var (ncid, 'density', nf_double, 5, final_mom_dim, density_id)
    status = netcdf_def_var (ncid, 'upar',    nf_double, 5, final_mom_dim, upar_id)
    status = netcdf_def_var (ncid, 'tpar',    nf_double, 5, final_mom_dim, tpar_id)
    status = netcdf_def_var (ncid, 'tperp',   nf_double, 5, final_mom_dim, tperp_id)
    status = netcdf_def_var (ncid, 'qparflux',   nf_double, 5, final_mom_dim, qparflux_id)
    status = netcdf_def_var (ncid, 'pperpj1',   nf_double, 5, final_mom_dim, pperpj1_id)
    status = netcdf_def_var (ncid, 'qpperpj1',   nf_double, 5, final_mom_dim, qpperpj1_id)

    status = netcdf_def_var (ncid, 'phi00',     nf_double, 3, loop_phi_dim, phi00_id)
    status = netcdf_def_var (ncid, 'ntot00',    nf_double, 4, loop_mom_dim, ntot00_id)
    status = netcdf_def_var (ncid, 'density00', nf_double, 4, loop_mom_dim, density00_id)
    status = netcdf_def_var (ncid, 'upar00',    nf_double, 4, loop_mom_dim, upar00_id)
    status = netcdf_def_var (ncid, 'tpar00',    nf_double, 4, loop_mom_dim, tpar00_id)
    status = netcdf_def_var (ncid, 'tperp00',   nf_double, 4, loop_mom_dim, tperp00_id)
    
    if (write_hrate) then
       status = netcdf_def_var (ncid, 'hrate_tot',  nf_double, 3, heat_dim, hrateavg_id)
       status = netcdf_def_var (ncid, 'hrate_by_k', nf_double, 5, heatk_dim, hrate_by_k_id)
    end if

!    status = netcdf_put_att (ncid, phtot_id, 'long_name', 'Field amplitude')

    status = netcdf_def_var (ncid, 'input_file', nf_char, 2, nin_dim, input_id)
    status = netcdf_put_att (ncid, input_id, 'long_name', 'Input file')

    status = nf_enddef (ncid)  ! out of definition mode

    status = netcdf_put_var (ncid, nproc_id, nproc)

        ! added by EAB 03/05/04 for GS2 movies
    if(my_make_movie) then
       xmode_dim (1) = nx_dim
       xmode_dim (2) = ny_dim
       xmode_dim (3) = nth_dim
       xmode_dim (4) = time_movie_dim
       status = netcdf_put_att (ncid_movie, nf_global, 'title', 'GS2 Simulation x,y,theta Data for Movies')
       status = netcdf_put_att (ncid_movie, nf_global, 'Conventions', &
            'http://gs2.sourceforge.net')
       !    status = netcdf_def_var (ncid_movie, 'code_info_movie', nf_char, 1, char10_dim, code_id_movie)
       !    status = netcdf_put_att (ncid_movie, code_id_movie, 'long_name', 'GS2')
       !    ci = 'c1'
       !    status = netcdf_put_att (ncid_movie, code_id_movie, trim(ci), 'Date: '//trim(datestamp))
       !    ci = 'c2'
       !    status = netcdf_put_att (ncid_movie, code_id_movie, trim(ci), &
       !         'Time: '//trim(timestamp)//' '//trim(timezone))
       !    ci = 'c3'
       !    status = netcdf_put_att (ncid_movie, code_id_movie, trim(ci), &
       !         'netCDF version '//trim(NF_INQ_LIBVERS()))
       
       status = netcdf_def_var (ncid_movie, 'nx', nf_int, 0, 0, nx_id)
       status = netcdf_def_var (ncid_movie, 'ny', nf_int, 0, 0, ny_id)
       status = netcdf_def_var (ncid_movie, 'ntheta', nf_int, 0, 0, nth_id)    
       status = netcdf_def_var (ncid_movie, 'tm', nf_double, 1, time_movie_dim, time_movie_id)
       status = netcdf_put_att (ncid_movie, time_movie_id, 'long_name', 'Time')
       status = netcdf_put_att (ncid_movie, time_movie_id, 'units', 'L/vt')
       status = netcdf_def_var (ncid_movie, 'x', nf_double, 1, nx_dim, x_id)
       status = netcdf_put_att (ncid_movie, x_id, 'long_name', 'x / rho')
       status = netcdf_def_var (ncid_movie, 'y', nf_double, 1, ny_dim, y_id)
       status = netcdf_put_att (ncid_movie, y_id, 'long_name', 'y / rho')
       status = netcdf_def_var (ncid_movie, 'theta', nf_double, 1, nth_dim, th_id)
       if(fphi > zero) then
          status = netcdf_def_var (ncid_movie, 'phi_by_xmode', nf_double, 4, xmode_dim, phi_by_xmode_id)
       end if
       if(fapar > zero) then
          status = netcdf_def_var (ncid_movie, 'apar_by_xmode', nf_double, 4, xmode_dim, apar_by_xmode_id)
       end if
       if(fbpar > zero) then
          status = netcdf_def_var (ncid_movie, 'bpar_by_xmode', nf_double, 4, xmode_dim, bpar_by_xmode_id)
       end if
       status = nf_enddef (ncid_movie)  ! out of definition mode
    endif

  end subroutine define_vars

  subroutine nc_eigenfunc (phase)

    use netcdf_mod, only: netcdf_put_var
    use convert, only: c2r
    use run_parameters, only: fphi, fapar, fbpar
    use fields_arrays, only: phi, apar, bpar
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    complex, dimension(:,:), intent (in) :: phase

    complex, dimension(-ntgrid:ntgrid, ntheta0, naky) :: tmp
    real, dimension(2, ntheta0, naky) :: ri2
    real, dimension (2, 2*ntgrid+1, ntheta0, naky) :: ri3
    integer :: status, ig

    call c2r (phase, ri2)
    status = netcdf_put_var (ncid, phase_id, ri2)

    if (fphi > zero) then
       do ig = -ntgrid, ntgrid
          tmp(ig,:,:) = phi(ig,:,:)/phase(:,:)
       end do
       call c2r (tmp, ri3)
       status = netcdf_put_var(ncid, phi_norm_id, ri3)
    end if

    if (fapar > zero) then
       do ig = -ntgrid, ntgrid
          tmp(ig,:,:) = apar(ig,:,:)/phase(:,:)
       end do
       call c2r (tmp, ri3)
       status = netcdf_put_var(ncid, apar_norm_id, ri3)
    end if

    if (fbpar > zero) then
       do ig = -ntgrid, ntgrid
          tmp(ig,:,:) = bpar(ig,:,:)/phase(:,:)
       end do
       call c2r (tmp, ri3)
       status = netcdf_put_var(ncid, bpar_norm_id, ri3)
    end if

  end subroutine nc_eigenfunc

  subroutine nc_final_fields

    use netcdf_mod, only: netcdf_put_var
    use convert, only: c2r
    use run_parameters, only: fphi, fapar, fbpar
    use fields_arrays, only: phi, apar, bpar
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    real, dimension (2, 2*ntgrid+1, ntheta0, naky) :: ri3
    integer :: status

    if (fphi > zero) then
       call c2r (phi, ri3)
       status = netcdf_put_var(ncid, phi_id, ri3)
    end if

    if (fapar > zero) then
       call c2r (apar, ri3)
       status = netcdf_put_var(ncid, apar_id, ri3)
    end if

    if (fbpar > zero) then
       call c2r (bpar, ri3)
       status = netcdf_put_var(ncid, bpar_id, ri3)
    end if

  end subroutine nc_final_fields

!MR begin
  subroutine nc_write_fields (nout, phinew, aparnew, bparnew)
    use netcdf_mod, only: netcdf_put_var
    use convert, only: c2r
    use run_parameters, only: fphi, fapar, fbpar
!    use fields_arrays, only: phi, apar, bpar
!    use fields, only: phinew, aparnew, bparnew
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    
    real, dimension (:,:,:), intent (in) :: phinew, aparnew, bparnew
    integer, intent (in) :: nout

    real, dimension (2, 2*ntgrid+1, ntheta0, naky) :: ri3
    integer, dimension (5) :: start5, count5
    integer :: status

    start5(1) = 1
    start5(2) = 1
    start5(3) = 1
    start5(4) = 1
    start5(5) = nout
    
    count5(1) = 2
    count5(2) = ntheta0
    count5(3) = naky
    count5(4) = nspec
    count5(5) = 1

    if (fphi > zero) then
       call c2r (phinew, ri3)
       status = netcdf_put_vara(ncid, phi_t_id, start5, count5, ri3)
    end if

    if (fapar > zero) then
       call c2r (aparnew, ri3)
       status = netcdf_put_vara(ncid, apar_t_id, start5, count5, ri3)
    end if

    if (fbpar > zero) then
       call c2r (bparnew, ri3)
       status = netcdf_put_vara(ncid, bpar_t_id, start5, count5, ri3)
    end if
  end subroutine nc_write_fields
!MR end

  subroutine nc_final_epar (epar)

    use netcdf_mod, only: netcdf_put_var
    use convert, only: c2r
    use run_parameters, only: fphi, fapar, fbpar
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    complex, dimension (:,:,:), intent (in) :: epar
    real, dimension (2, 2*ntgrid+1, ntheta0, naky) :: ri3
    integer :: status

    call c2r (epar, ri3)
    status = netcdf_put_var(ncid, epar_id, ri3)

  end subroutine nc_final_epar

  subroutine nc_final_moments (ntot, density, upar, tpar, tperp, qparflux, pperpj1, qpperpj1)

    use netcdf_mod, only: netcdf_put_var
    use convert, only: c2r

    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use species, only: nspec

    complex, dimension (:,:,:,:), intent (in) :: ntot, density, upar, tpar, tperp
    real, dimension (2, 2*ntgrid+1, ntheta0, naky, nspec) :: ri4
    integer :: status

    call c2r (ntot, ri4)
    status = netcdf_put_var(ncid, ntot_id, ri4)
    
    call c2r (density, ri4)
    status = netcdf_put_var(ncid, density_id, ri4)
    
    call c2r (upar, ri4)
    status = netcdf_put_var(ncid, upar_id, ri4)
    
    call c2r (tpar, ri4)
    status = netcdf_put_var(ncid, tpar_id, ri4)

    call c2r (tperp, ri4)
    status = netcdf_put_var(ncid, tperp_id, ri4)

    call c2r (qparflux, ri4)
    status = netcdf_put_var(ncid, qparflux_id, ri4)

    call c2r (pperpj1, ri4)
    status = netcdf_put_var(ncid, pperpj1_id, ri4)

    call c2r (qpperpj1, ri4)
    status = netcdf_put_var(ncid, qpperpj1_id, ri4)

  end subroutine nc_final_moments

  subroutine nc_loop_moments (nout, ntot2, ntot2_by_mode, &
       phi00, ntot00, density00, upar00, tpar00, tperp00)

    use netcdf_mod, only: netcdf_put_vara
    use convert, only: c2r

    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use species, only: nspec

    integer, intent (in) :: nout
    real, dimension (:), intent (in) :: ntot2
    real, dimension (:,:,:), intent (in) :: ntot2_by_mode
    complex, dimension (:), intent (in) :: phi00
    complex, dimension (:,:), intent (in) :: ntot00, density00, upar00, tpar00, tperp00
    real, dimension (2, ntheta0, nspec) :: ri2
    real, dimension (2, ntheta0) :: ri1
    integer, dimension (2) :: start, count
    integer, dimension (3) :: start3, count3
    integer, dimension (4) :: start00, count00, start4, count4
    integer :: status

    start00(1) = 1
    start00(2) = 1
    start00(3) = 1
    start00(4) = nout
    
    count00(1) = 2
    count00(2) = ntheta0
    count00(3) = nspec
    count00(4) = 1

    start3(1) = 1
    start3(2) = 1
    start3(3) = nout
    
    count3(1) = 2
    count3(2) = ntheta0
    count3(3) = 1

    start(1) = 1
    start(2) = nout
    
    count(1) = nspec
    count(2) = 1

    start4(1) = 1
    start4(2) = 1
    start4(3) = 1
    start4(4) = nout
    
    count4(1) = ntheta0
    count4(2) = naky
    count4(3) = nspec
    count4(4) = 1

    status = netcdf_put_vara(ncid, ntot2_id, start, count, ntot2)
    status = netcdf_put_vara(ncid, ntot2_by_mode_id, start4, count4, ntot2_by_mode)

    call c2r (phi00, ri1)
    status = netcdf_put_vara(ncid, phi00_id, start3, count3, ri1)
    
    call c2r (ntot00, ri2)
    status = netcdf_put_vara(ncid, ntot00_id, start00, count00, ri2)
    
    call c2r (density00, ri2)
    status = netcdf_put_vara(ncid, density00_id, start00, count00, ri2)
    
    call c2r (upar00, ri2)
    status = netcdf_put_vara(ncid, upar00_id, start00, count00, ri2)
    
    call c2r (tpar00, ri2)
    status = netcdf_put_vara(ncid, tpar00_id, start00, count00, ri2)

    call c2r (tperp00, ri2)
    status = netcdf_put_vara(ncid, tperp00_id, start00, count00, ri2)

  end subroutine nc_loop_moments

  subroutine nc_final_an (antot, antota, antotp)

    use netcdf_mod, only: netcdf_put_var
    use convert, only: c2r

    use run_parameters, only: fphi, fapar, fbpar
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    complex, dimension (:,:,:) :: antot, antota, antotp
    real, dimension (2, 2*ntgrid+1, ntheta0, naky) :: ri3
    integer :: status

    if (fphi > zero) then
       call c2r (antot, ri3)
       status = netcdf_put_var(ncid, antot_id, ri3)
    end if

    if (fapar > zero) then
       call c2r (antota, ri3)
       status = netcdf_put_var(ncid, antota_id, ri3)
    end if

    if (fbpar > zero) then
       call c2r (antotp, ri3)
       status = netcdf_put_var(ncid, antotp_id, ri3)
    end if

  end subroutine nc_final_an

  subroutine nc_qflux (nout, qheat, qmheat, qbheat, &
       heat_par,  mheat_par,  bheat_par, &
       heat_perp, mheat_perp, bheat_perp, &
       heat_fluxes, mheat_fluxes, bheat_fluxes, hflux_tot, anorm)

    use species, only: nspec
    use kt_grids, only: naky, ntheta0
    use run_parameters, only: fphi, fapar, fbpar
    use netcdf_mod, only: netcdf_put_vara, netcdf_put_var1

    integer, intent (in) :: nout
    real, dimension (:,:,:), intent (in) :: qheat, qmheat, qbheat
    real, dimension (:), intent (in) :: heat_par, mheat_par, bheat_par
    real, dimension (:), intent (in) :: heat_perp, mheat_perp, bheat_perp
    real, dimension (:), intent (in) :: heat_fluxes, mheat_fluxes, bheat_fluxes
    real, intent (in) :: hflux_tot, anorm
    integer, dimension (2) :: start, count
    integer, dimension (4) :: start4, count4
    integer :: status

    start4(1) = 1
    start4(2) = 1
    start4(3) = 1
    start4(4) = nout
    
    count4(1) = ntheta0
    count4(2) = naky
    count4(3) = nspec
    count4(4) = 1

    start(1) = 1
    start(2) = nout
    
    count(1) = nspec
    count(2) = 1

    status = netcdf_put_var1(ncid, anorm_id, nout, anorm)

    if (fphi > zero) then
       status = netcdf_put_vara(ncid, es_heat_flux_id, start, count, heat_fluxes)
       status = netcdf_put_vara(ncid, es_heat_par_id,  start, count, heat_par)
       status = netcdf_put_vara(ncid, es_heat_perp_id, start, count, heat_perp)
       status = netcdf_put_vara(ncid, es_heat_by_k_id, start4, count4, qheat)
    end if

    if (fapar > zero) then
       status = netcdf_put_vara(ncid, apar_heat_flux_id, start, count, mheat_fluxes)
       status = netcdf_put_vara(ncid, apar_heat_par_id,  start, count, mheat_par)
       status = netcdf_put_vara(ncid, apar_heat_perp_id, start, count, mheat_perp)
       status = netcdf_put_vara(ncid, apar_heat_by_k_id, start4, count4, qmheat)
    end if

    if (fbpar > zero) then
       status = netcdf_put_vara(ncid, bpar_heat_flux_id, start, count, bheat_fluxes)
       status = netcdf_put_vara(ncid, bpar_heat_par_id,  start, count, bheat_par)
       status = netcdf_put_vara(ncid, bpar_heat_perp_id, start, count, bheat_perp)
       status = netcdf_put_vara(ncid, bpar_heat_by_k_id, start4, count4, qbheat)
    end if
    
    status = netcdf_put_var1(ncid, hflux_tot_id, nout, hflux_tot)

  end subroutine nc_qflux

  subroutine nc_pflux (nout, pflux, pmflux, pbflux, &
       part_fluxes, mpart_fluxes, bpart_fluxes, zflux_tot)

    use species, only: nspec
    use kt_grids, only: naky, ntheta0
    use run_parameters, only: fphi, fapar, fbpar
    use netcdf_mod, only: netcdf_put_vara, netcdf_put_var1

    integer, intent (in) :: nout
    real, dimension (:,:,:), intent (in) :: pflux, pmflux, pbflux
    real, dimension(:), intent (in) :: part_fluxes, mpart_fluxes, bpart_fluxes
    real, intent (in) :: zflux_tot
    integer, dimension (2) :: start, count
    integer, dimension (4) :: start4, count4
    integer :: status

    start4(1) = 1
    start4(2) = 1
    start4(3) = 1
    start4(4) = nout
    
    count4(1) = ntheta0
    count4(2) = naky
    count4(3) = nspec
    count4(4) = 1

    start(1) = 1
    start(2) = nout
    
    count(1) = nspec
    count(2) = 1

    if (fphi > zero) then
       status = netcdf_put_vara(ncid, es_part_flux_id, start, count, part_fluxes)
       status = netcdf_put_vara(ncid, es_part_by_k_id, start4, count4, pflux)
    end if

    if (fapar > zero) then
       status = netcdf_put_vara(ncid, apar_part_flux_id, start, count, mpart_fluxes)
       status = netcdf_put_vara(ncid, apar_part_by_k_id, start4, count4, pmflux)
    end if

    if (fbpar > zero) then
       status = netcdf_put_vara(ncid, bpar_part_flux_id, start, count, bpart_fluxes)
       status = netcdf_put_vara(ncid, bpar_part_by_k_id, start4, count4, pbflux)
    end if
    
    status = netcdf_put_var1(ncid, zflux_tot_id, nout, zflux_tot)

  end subroutine nc_pflux

  subroutine nc_vflux (nout, vflux, vmflux, vbflux, &
       mom_fluxes, mmom_fluxes, bmom_fluxes, vflux_tot)

    use species, only: nspec
    use kt_grids, only: naky, ntheta0
    use run_parameters, only: fphi, fapar, fbpar
    use netcdf_mod, only: netcdf_put_vara, netcdf_put_var1

    integer, intent (in) :: nout
    real, dimension (:,:,:), intent (in) :: vflux, vmflux, vbflux
    real, dimension(:), intent (in) :: mom_fluxes, mmom_fluxes, bmom_fluxes
    real, intent (in) :: vflux_tot
    integer, dimension (2) :: start, count
    integer, dimension (4) :: start4, count4
    integer :: status

    start4(1) = 1
    start4(2) = 1
    start4(3) = 1
    start4(4) = nout
    
    count4(1) = ntheta0
    count4(2) = naky
    count4(3) = nspec
    count4(4) = 1

    start(1) = 1
    start(2) = nout
    
    count(1) = nspec
    count(2) = 1

    if (fphi > zero) then
       status = netcdf_put_vara(ncid, es_mom_flux_id, start, count, mom_fluxes)
       status = netcdf_put_vara(ncid, es_mom_by_k_id, start4, count4, vflux)
    end if

    if (fapar > zero) then
       status = netcdf_put_vara(ncid, apar_mom_flux_id, start, count, mmom_fluxes)
       status = netcdf_put_vara(ncid, apar_mom_by_k_id, start4, count4, vmflux)
    end if

    if (fbpar > zero) then
       status = netcdf_put_vara(ncid, bpar_mom_flux_id, start, count, bmom_fluxes)
       status = netcdf_put_vara(ncid, bpar_mom_by_k_id, start4, count4, vbflux)
    end if
    
    status = netcdf_put_var1(ncid, vflux_tot_id, nout, vflux_tot)

  end subroutine nc_vflux

  subroutine nc_loop (nout, time, fluxfac, &
       phi0,   phi2,   phi2_by_mode, &! phiavg, &
       apar0,  apar2,  apar2_by_mode, &
       bpar0, bpar2, bpar2_by_mode, &
       hrateavg, rate_by_k, &
       omega, omegaavg, woutunits, phitot, write_omega, write_hrate)

    use run_parameters, only: fphi, fapar, fbpar
    use netcdf_mod, only: netcdf_put_var1, netcdf_put_vara
    use kt_grids, only: naky, ntheta0
    use species, only: nspec
    use convert, only: c2r
    use parameter_scan_arrays, only: write_scan_parameter,&
                                     current_scan_parameter_value

    integer, intent (in) :: nout
    real, intent (in) :: time, phi2, apar2, bpar2
    real, dimension (:), intent (in) :: fluxfac, woutunits
    complex, dimension(:,:), intent (in) :: phi0, apar0, bpar0, omega, omegaavg !, phiavg
    real, dimension(:,:), intent (in) :: phi2_by_mode, apar2_by_mode, bpar2_by_mode, phitot
    real, dimension (:,:), intent (in) :: hrateavg
    real, dimension (:,:,:,:), intent (in) :: rate_by_k
    logical :: write_omega, write_hrate
    real, dimension (ntheta0) :: field2_by_kx
    real, dimension (naky) :: field2_by_ky
    real, dimension (2, ntheta0, naky) :: ri2
    complex, dimension (ntheta0, naky) :: tmp
    integer, dimension (5) :: start5, count5
    integer, dimension (4) :: start0, count0, start4, count4
    integer, dimension (3) :: start, count, starth, counth
    integer, dimension (2) :: startx, countx, starty, county, starts, counts
    integer :: status, it, ik

    status = netcdf_put_var1(ncid, time_id, nout, time)

    start(1) = 1
    start(2) = 1
    start(3) = nout

    count(1) = ntheta0
    count(2) = naky
    count(3) = 1

    start0(1) = 1
    start0(2) = 1
    start0(3) = 1
    start0(4) = nout

    count0(1) = 2
    count0(2) = ntheta0
    count0(3) = naky
    count0(4) = 1

    start4(1) = 1
    start4(2) = 1
    start4(3) = 1
    start4(4) = nout

    count4(1) = ntheta0
    count4(2) = naky
    count4(3) = nspec
    count4(4) = 1

    start5(1) = 1
    start5(2) = 1
    start5(3) = 1
    start5(4) = 1
    start5(5) = nout

    count5(1) = ntheta0
    count5(2) = naky
    count5(3) = nspec
    count5(4) = 3
    count5(5) = 1

    starty(1) = 1
    starty(2) = nout

    county(1) = naky
    county(2) = 1

    startx(1) = 1
    startx(2) = nout

    countx(1) = ntheta0
    countx(2) = 1

    starts(1) = 1
    starts(2) = nout

    counts(1) = nspec
    counts(2) = 1

    starth(1) = 1 
    starth(2) = 1
    starth(3) = nout

    counth(1) = nspec
    counth(2) = 3
    counth(3) = 1

    if (fphi > zero) then

       if (ntheta0 > 1) then
          do it = 1, ntheta0
             field2_by_kx(it) = sum(phi2_by_mode(it,:)*fluxfac(:))
          end do
          status = netcdf_put_vara(ncid, phi2_by_kx_id, startx, countx, field2_by_kx)
       end if

       if (naky > 1) then
          do ik = 1, naky
             field2_by_ky(ik) = sum(phi2_by_mode(:,ik)*fluxfac(ik))
          end do
          status = netcdf_put_vara(ncid, phi2_by_ky_id, starty, county, field2_by_ky)
       end if

       call c2r (phi0, ri2) 
       status = netcdf_put_vara(ncid, phi0_id, start0, count0, ri2)
!       call c2r (phiavg, ri2) 
!       status = netcdf_put_vara(ncid, phiavg_id, start0, count0, ri2)
       status = netcdf_put_vara(ncid, phi2_by_mode_id, start, count, phi2_by_mode)
       status = netcdf_put_var1(ncid, phi2_id, nout, phi2)
    end if

    if (fapar > zero) then

       if (ntheta0 > 1) then
          do it = 1, ntheta0
             field2_by_kx(it) = sum(apar2_by_mode(it,:)*fluxfac(:))
          end do
          status = netcdf_put_vara(ncid, apar2_by_kx_id, startx, countx, field2_by_kx)
       end if

       if (naky > 1) then
          do ik = 1, naky
             field2_by_ky(ik) = sum(apar2_by_mode(:,ik)*fluxfac(ik))
          end do
          status = netcdf_put_vara(ncid, apar2_by_ky_id, starty, county, field2_by_ky)
       end if

       call c2r (apar0, ri2) 
       status = netcdf_put_vara(ncid, apar0_id, start0, count0, ri2)
       status = netcdf_put_vara(ncid, apar2_by_mode_id, start, count, apar2_by_mode)
       status = netcdf_put_var1(ncid, apar2_id, nout, apar2)
    end if

    if (fbpar > zero) then

       if (ntheta0 > 1) then
          do it = 1, ntheta0
             field2_by_kx(it) = sum(bpar2_by_mode(it,:)*fluxfac(:))
          end do
          status = netcdf_put_vara(ncid, bpar2_by_kx_id, startx, countx, field2_by_kx)
       end if

       if (naky > 1) then
          do ik = 1, naky
             field2_by_ky(ik) = sum(bpar2_by_mode(:,ik)*fluxfac(ik))
          end do
          status = netcdf_put_vara(ncid, bpar2_by_ky_id, starty, county, field2_by_ky)
       end if

       call c2r (bpar0, ri2) 
       status = netcdf_put_vara(ncid, bpar0_id, start0, count0, ri2)
       status = netcdf_put_vara(ncid, bpar2_by_mode_id, start, count, bpar2_by_mode)
       status = netcdf_put_var1(ncid, bpar2_id, nout, bpar2)
    end if
        
    if (write_hrate) then
       status = netcdf_put_vara (ncid, hrateavg_id, starth, counth, hrateavg)
       status = netcdf_put_vara (ncid, hrate_by_k_id, start5, count5, rate_by_k)
    end if

    if (write_omega) then
       do it = 1, ntheta0
          tmp(it, :) = omega(it, :) * woutunits
       end do
       
       call c2r (tmp, ri2) 
       status = netcdf_put_vara(ncid, omega_id, start0, count0, ri2)
       
       do it = 1, ntheta0
          tmp(it, :) = omegaavg(it, :) * woutunits
       end do
       
       call c2r (tmp, ri2) 
       status = netcdf_put_vara(ncid, omegaavg_id, start0, count0, ri2)
    end if

    if (write_scan_parameter) then
       status = netcdf_put_var1(ncid, &
                                current_scan_parameter_value_id, &
                                nout, &
                                current_scan_parameter_value)
    end if
!    status = netcdf_put_vara(ncid, phtot_id, start, count, phitot)
    
    if (mod(nout, 10) == 0) status = nf_sync (ncid)

  end subroutine nc_loop

  ! added by EAB on 03/05/04
  subroutine nc_loop_movie (nout_movie, time, yxphi, yxapar, yxbpar)
    use gs2_layouts, only: yxf_lo
    use theta_grid, only: ntgrid
    use run_parameters, only: fphi, fapar, fbpar
    use netcdf_mod, only: netcdf_put_var1, netcdf_put_vara
    integer, intent (in) :: nout_movie
    real, intent (in) :: time
    real, dimension (:,:,:), intent (in):: yxphi, yxapar, yxbpar 
    integer :: status
    integer, dimension (4) :: start, count

      
    status = netcdf_put_var1(ncid_movie, time_movie_id, nout_movie, time)
    start(1) = 1
    start(2) = 1
    start(3) = 1
    start(4) = nout_movie
    count(1) = yxf_lo%nx
    count(2) = yxf_lo%ny
    count(3) = 2*ntgrid+1
    count(4) = 1


    if ( fphi > zero) then
       status = netcdf_put_vara(ncid_movie, phi_by_xmode_id, start, count, yxphi)
    end if
    if ( fapar > zero) then
       status = netcdf_put_vara(ncid_movie, apar_by_xmode_id, start, count, yxapar)
    end if

    if ( fbpar > zero) then
       status = netcdf_put_vara(ncid_movie, bpar_by_xmode_id, start, count, yxbpar)
    end if

    if (mod(nout_movie, 10) == 0) status = nf_sync (ncid_movie)


  end subroutine nc_loop_movie

  subroutine nc_species
    
    use run_parameters, only: beta
    use netcdf_mod, only: netcdf_put_var
    use species, only: spec
    real :: betatot
    integer :: status

    status = netcdf_put_var (ncid, charge_id, spec%z)
    status = netcdf_put_var (ncid, mass_id, spec%mass)
    status = netcdf_put_var (ncid, dens_id, spec%dens)
    status = netcdf_put_var (ncid, temp_id, spec%temp)
    status = netcdf_put_var (ncid, tprim_id, spec%tprim)
    status = netcdf_put_var (ncid, fprim_id, spec%fprim)
    status = netcdf_put_var (ncid, uprim_id, spec%uprim)
    status = netcdf_put_var (ncid, uprim2_id, spec%uprim2)
    status = netcdf_put_var (ncid, vnewk_id, spec%vnewk)
    status = netcdf_put_var (ncid, spec_type_id, spec%type)

!    betatot = beta * ...
!    status = netcdf_put_var (ncid, betatot_id, betatot)

  end subroutine nc_species

  subroutine nc_geo

    use netcdf_mod, only: netcdf_put_var
    use theta_grid, only: bmag, gradpar, gbdrift, gbdrift0, &
         cvdrift, cvdrift0, gds2, gds21, gds22, grho, jacob, &
         shat, drhodpsi, eps
    integer :: status

    status = netcdf_put_var (ncid, bmag_id, bmag)
    status = netcdf_put_var (ncid, gradpar_id, gradpar)
    status = netcdf_put_var (ncid, gbdrift_id, gbdrift)
    status = netcdf_put_var (ncid, gbdrift0_id, gbdrift0)
    status = netcdf_put_var (ncid, cvdrift_id, cvdrift)
    status = netcdf_put_var (ncid, cvdrift0_id, cvdrift0)
    status = netcdf_put_var (ncid, gds2_id, gds2)
    status = netcdf_put_var (ncid, gds21_id, gds21)
    status = netcdf_put_var (ncid, gds22_id, gds22)
    status = netcdf_put_var (ncid, grho_id, grho)
    status = netcdf_put_var (ncid, jacob_id, jacob)

    status = netcdf_put_var (ncid, shat_id, shat)
    status = netcdf_put_var (ncid, eps_id, eps)
!    status = netcdf_put_var (ncid, q_id, eps)  ! find the q variable: qsf, q_val, or qinp?
!    status = netcdf_put_var (ncid, q_id, eps)  ! write both beta_prime variables with comment
    status = netcdf_put_var (ncid, drhodpsi_id, drhodpsi)   

  end subroutine nc_geo

end module gs2_io

