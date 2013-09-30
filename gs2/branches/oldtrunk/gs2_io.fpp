# include "define.inc"

module gs2_io

# ifdef NETCDF
  use netcdf, only: NF90_NOERR
  use netcdf_utils, only: netcdf_error, kind_nf
# endif

  implicit none

  private

  public :: init_gs2_io, nc_eigenfunc, nc_final_fields, nc_final_epar
  public :: nc_final_moments, nc_final_an, nc_finish
  public :: nc_qflux, nc_vflux, nc_pflux, nc_loop, nc_loop_moments
  public :: nc_loop_stress, nc_loop_vres
  public :: nc_loop_movie
  public :: nc_loop_fullmom

  logical, parameter :: serial_io = .true.
# ifdef NETCDF
!  integer :: ncid
  logical :: proc_write
  integer (kind_nf) :: ncid

!  integer :: naky_dim, nakx_dim, nttot_dim, negrid_dim, nlambda_dim, nspec_dim, ncoord_dim, ncoordt_dim
!  integer :: nsign_dim, time_dim, char10_dim, char200_dim, ri_dim, nlines_dim, nheat_dim
  integer (kind_nf) :: naky_dim, nakx_dim, nttot_dim, negrid_dim, nlambda_dim, nspec_dim, ncoord_dim, ncoordt_dim
  integer (kind_nf) :: nsign_dim, time_dim, char10_dim, char200_dim, ri_dim, nlines_dim, nheat_dim
  
  integer, dimension (6) :: lp_dim
  integer, dimension (5) :: field_dim, final_mom_dim, heatk_dim, err_dim
  integer, dimension (4) :: omega_dim, fluxk_dim, final_field_dim, loop_mom_dim
  integer, dimension (3) :: fluxx_dim
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
  integer :: phi2_by_mode_id, apar2_by_mode_id, bpar2_by_mode_id, errest_by_mode_id, lpcoef_by_mode_id
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
  integer :: apar_heat_by_x_id
  integer :: bpar_heat_by_k_id, bpar_mom_by_k_id, bpar_part_by_k_id
  integer :: phi_t_id, apar_t_id, bpar_t_id
  integer :: phi_norm_id, apar_norm_id, bpar_norm_id
  integer :: phi_id, apar_id, bpar_id, epar_id
  integer :: antot_id, antota_id, antotp_id
  integer :: ntot_id, density_id, upar_id, tpar_id, tperp_id
  integer :: rstress_id, ustress_id
  integer :: ntot2_id, ntot2_by_mode_id, ntot20_id, ntot20_by_mode_id
  integer :: phi00_id, ntot00_id, density00_id, upar00_id, tpar00_id, tperp00_id
  integer :: qflux_neo_by_k_id, pflux_neo_by_k_id, input_id
  integer :: charge_id, mass_id, dens_id, temp_id, tprim_id, fprim_id
  integer :: uprim_id, uprim2_id, vnewk_id, spec_type_id
  integer :: bmag_id, gradpar_id, gbdrift_id, gbdrift0_id
  integer :: cvdrift_id, cvdrift0_id, gds2_id, gds21_id, gds22_id
  integer :: grho_id, jacob_id, shat_id, eps_id, drhodpsi_id, q_id
  integer :: code_id, datestamp_id, timestamp_id, timezone_id
  integer :: h_energy_id, h_energy_dot_id, h_antenna_id
  integer :: h_eapar_id, h_ebpar_id
  integer :: h_delfs2_id, h_hs2_id, h_phis2_id
  integer :: h_hypervisc_id, h_hyperres_id, h_collisions_id
  integer :: h_gradients_id, h_hypercoll_id, h_heating_id
  integer :: h_imp_colls_id
  integer :: hk_energy_id, hk_energy_dot_id, hk_antenna_id
  integer :: hk_eapar_id, hk_ebpar_id
  integer :: hk_delfs2_id, hk_hs2_id, hk_phis2_id
  integer :: hk_hypervisc_id, hk_hyperres_id, hk_collisions_id
  integer :: hk_gradients_id, hk_hypercoll_id, hk_heating_id, hk_imp_colls_id
  integer, dimension (5) :: mom_dim
  integer :: ntot0_id, density0_id, upar0_id, tpar0_id, tperp0_id
  logical :: write_apar_t, write_phi_t, write_bpar_t ! Should the fields be written out every nwrite?

# endif
  real :: zero
  
!  include 'netcdf.inc'
  
contains

  subroutine init_gs2_io (write_nl_flux, write_omega, write_stress, &
      write_phiavg, write_hrate, write_final_antot, write_eigenfunc, &
      make_movie, nmovie_tot, write_verr, &
      write_full_moments_notgc, write_phi_over_time, write_apar_over_time, &
      write_bpar_over_time) 

    !write_fields_over_time added by EGH 08/2009
!David has made some changes to this subroutine (may 2005) now should 
!be able to do movies for linear box runs as well as nonlinear box runs.

    use mp, only: proc0, barrier
    use file_utils, only: run_name
    use gs2_transforms, only: init_transforms
    use kt_grids, only: naky, ntheta0,nx,ny
    use gs2_layouts, only: yxf_lo
    use theta_grid, only: ntgrid
    use le_grids, only: nlambda, negrid
    use species, only: nspec
!    use netcdf_mod
# ifdef NETCDF
    use netcdf, only: NF90_CLOBBER, nf90_create
    use netcdf_utils, only: get_netcdf_code_precision, netcdf_real
# endif
!    logical :: write_nl_flux, write_omega, write_stress, write_phiavg, write_hrate, make_movie
!    logical :: write_final_antot, write_eigenfunc, write_verr
    logical, intent(in) :: write_nl_flux, write_omega, write_stress
    logical, intent(in) :: write_phiavg, write_hrate, make_movie
    logical, intent(in) :: write_final_antot, write_eigenfunc, write_verr
    logical, intent(in) :: write_full_moments_notgc, write_phi_over_time
    logical, intent(in) :: write_apar_over_time,  write_bpar_over_time

    integer :: nmovie_tot
# ifdef NETCDF
    logical :: accelerated
    character (300) :: filename, filename_movie
!    integer :: ierr         ! 0 if initialization is successful
    integer :: status

    !EGH
    write_phi_t = write_phi_over_time
    write_apar_t = write_apar_over_time
    write_bpar_t = write_bpar_over_time


    zero = epsilon(0.0)
!    ierr = 0

    call netcdf_init (serial_io)
    if (netcdf_real == 0) netcdf_real = get_netcdf_code_precision()

    status = NF90_NOERR

    filename = trim(trim(run_name)//'.out.nc')
    if (serial_io) then
       ! only proc0 opens the file:
       if (proc0) then
          status = nf90_create (trim(filename), NF90_CLOBBER, ncid) 
          if (status /= NF90_NOERR) call netcdf_error (status, file=filename)
       end if
    else
       ! all processors open the file
       call barrier
       ! TT: do we want NF90_SHARE?
       status = nf90_create (trim(filename), NF90_CLOBBER, ncid) ! overwrites old
       if (proc0 .and. status /= NF90_NOERR) &
            call netcdf_error (status, file=filename)
       call barrier
    endif

!    if(status /= nf_noerr) then
!       write(*,*) 'Unable to create netcdf file ', filename
!       write(*,*) nf_strerror(status)
!       ierr = 1
!    endif

    ! added by EAB 03/2004 for making movie netcdf files
    my_make_movie = make_movie
    if (my_make_movie) then
    !perhaps put a call to init_y_transform
!    call init_y_transform_layouts  
       call init_transforms(ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx, ny, accelerated)
       filename_movie = trim(trim(run_name)//'.movie.nc')
       if (serial_io) then
          ! only proc0 opens the file:
          if (proc0) then 
             status = nf90_create (trim(filename_movie), NF90_CLOBBER, ncid_movie) 
             if (status /= NF90_NOERR) &
                  call netcdf_error (status, file=filename_movie)
          end if
       else
          ! all processors open the file
          call barrier
          status = nf90_create (trim(filename_movie), NF90_CLOBBER, ncid_movie) ! overwrites old
          if (proc0 .and. status /= NF90_NOERR) &
               call netcdf_error (status, file=filename_movie)
          call barrier
       endif
       
!       if(status /= nf_noerr) then
!          write(*,*) 'Unable to create netcdf file for movies', filename_movie
!          write(*,*) nf_strerror(status)
!          ierr = 1
!       endif
    endif
# endif
    if (proc0) then
       call define_dims (nmovie_tot)
       call define_vars (write_nl_flux, write_omega, write_stress, write_phiavg, &
            write_hrate, write_final_antot, write_eigenfunc, write_verr, &
            write_full_moments_notgc)
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
 !   use netcdf_mod
# ifdef NETCDF
    use netcdf, only: NF90_UNLIMITED
    use netcdf, only: nf90_def_dim
# endif

    integer, intent (in) :: nmovie_tot
# ifdef NETCDF
    integer :: status


    !<wkdoc> Associate the grid variables, e.g. ky, kx, with their size, e.g. naky, ntheta0 (= nakx) , and a variable which is later used to store these sizes in the NetCDF file, e.g. naky_dim, nakx_dim </wkdoc>
    status = nf90_def_dim (ncid, 'ky', naky, naky_dim)
    if (status /= NF90_NOERR) call netcdf_error (status, dim='ky')
    status = nf90_def_dim (ncid, 'kx', ntheta0, nakx_dim)
    if (status /= NF90_NOERR) call netcdf_error (status, dim='kx')
    status = nf90_def_dim (ncid, 'theta', 2*ntgrid+1, nttot_dim)
    if (status /= NF90_NOERR) call netcdf_error (status, dim='theta')
    status = nf90_def_dim (ncid, 'lambda', nlambda, nlambda_dim)
    if (status /= NF90_NOERR) call netcdf_error (status, dim='lambda')
    status = nf90_def_dim (ncid, 'egrid', negrid, negrid_dim)
    if (status /= NF90_NOERR) call netcdf_error (status, dim='egrid')
    status = nf90_def_dim (ncid, 'species', nspec, nspec_dim)
    if (status /= NF90_NOERR) call netcdf_error (status, dim='species')
    status = nf90_def_dim (ncid, 'sign', 2, nsign_dim)
    if (status /= NF90_NOERR) call netcdf_error (status, dim='sign')
    status = nf90_def_dim (ncid, 't', NF90_UNLIMITED, time_dim)
    if (status /= NF90_NOERR) call netcdf_error (status, dim='t')
    status = nf90_def_dim (ncid, 'char10', 10, char10_dim)
    if (status /= NF90_NOERR) call netcdf_error (status, dim='char10')
    status = nf90_def_dim (ncid, 'char200', 200, char200_dim)
    if (status /= NF90_NOERR) call netcdf_error (status, dim='char200')
    status = nf90_def_dim (ncid, 'nlines', num_input_lines, nlines_dim)
    if (status /= NF90_NOERR) call netcdf_error (status, dim='nlines')
    status = nf90_def_dim (ncid, 'ri', 2, ri_dim)
    if (status /= NF90_NOERR) call netcdf_error (status, dim='ri')
    status = nf90_def_dim (ncid, 'nheat', 7, nheat_dim)
    if (status /= NF90_NOERR) call netcdf_error (status, dim='nheat')

    ! added by EAB 03/05/04 for GS2 movies
    ! <wkdoc> If make_movie = .true., define some extra dimensions for movie output; x with dimension yxf_lo%nx; y with dimension yxf_lo%ny, and tm, (i.e. time) with dimension nmovie_tot = nstep/nmovie </wkdoc>
    if (my_make_movie) then
       nx=yxf_lo%nx
       ny=yxf_lo%ny
       status = nf90_def_dim (ncid_movie, 'x', nx, nx_dim)
       if (status /= NF90_NOERR) call netcdf_error (status, dim='x')
       status = nf90_def_dim (ncid_movie, 'y', ny, ny_dim)
       if (status /= NF90_NOERR) call netcdf_error (status, dim='y')
       status = nf90_def_dim (ncid_movie, 'theta', 2*ntgrid+1, nth_dim)
       if (status /= NF90_NOERR) call netcdf_error (status, dim='theta')
!
! can only have one NF_UNLIMITED, so:
!
!       status = nf90_def_dim(ncid_movie, 't', NF90_UNLIMITED, time_movie_dim)
       status = nf90_def_dim (ncid_movie, 'tm', nmovie_tot, time_movie_dim)
       if (status /= NF90_NOERR) call netcdf_error (status, dim='tm')
    endif
# endif
  end subroutine define_dims

  subroutine nc_grids

    use theta_grid, only: ntgrid, theta, eps, ntheta
    use kt_grids, only: naky, ntheta0, theta0, akx_out, aky_out, nx, ny
    use gs2_layouts, only: yxf_lo
    use species, only: nspec
    use le_grids, only: negrid, nlambda
    use nonlinear_terms, only: nonlin
# ifdef NETCDF
!    use netcdf_mod
    use netcdf, only: nf90_put_var
    use constants, only: pi
    
    integer :: status
    real :: nmesh

    real, dimension(:), allocatable :: x, y
    integer :: ik, it

    !<wkdoc> Store the size of the grid dimensions (as defined in def_dims), in the NetCDF file </wkdoc>
    status = nf90_put_var (ncid, nttot_id, 2*ntgrid+1)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, nttot_id)
    status = nf90_put_var (ncid, naky_id, naky)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, naky_id)
    status = nf90_put_var (ncid, nakx_id, ntheta0)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, nakx_id)
    status = nf90_put_var (ncid, nspec_id, nspec)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, nspec_id)

    status = nf90_put_var (ncid, akx_id, akx_out)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, akx_id)
    status = nf90_put_var (ncid, aky_id, aky_out)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, aky_id)
    status = nf90_put_var (ncid, theta_id, theta)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, theta_id)
    status = nf90_put_var (ncid, theta0_id, theta0)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, theta0_id)

    if (nonlin) then
       nmesh = (2*ntgrid+1)*2*nlambda*negrid*nx*ny*nspec
    else
       nmesh = (2*ntgrid+1)*2*nlambda*negrid*ntheta0*naky*nspec
    end if

    status = nf90_put_var (ncid, nmesh_id, nmesh)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, nmesh_id)

    ! added by EAB 03/05/04 for GS2 movies
    if(my_make_movie) then

       allocate(x(yxf_lo%nx))
       allocate(y(yxf_lo%ny))
       do it = 1, yxf_lo%nx 
          x(it) = 2.0*pi/akx_out(2)*(-0.5+ real(it-1)/real(yxf_lo%nx))
       end do
       do ik = 1, yxf_lo%ny
          y(ik) = 2.0*pi/aky_out(2)*(-0.5+ real(ik-1)/real(yxf_lo%ny))
       end do
       
       status = nf90_put_var (ncid_movie, nx_id, yxf_lo%nx)
       if (status /= NF90_NOERR) &
            call netcdf_error (status, ncid_movie, nx_id)
       status = nf90_put_var (ncid_movie, ny_id, yxf_lo%ny)
       if (status /= NF90_NOERR) &
            call netcdf_error (status, ncid_movie, ny_id)
       status = nf90_put_var (ncid_movie, nth_id, 2*ntgrid+1)
       if (status /= NF90_NOERR) &
            call netcdf_error (status, ncid_movie, nth_id)
       status = nf90_put_var (ncid_movie, x_id, x)
       if (status /= NF90_NOERR) &
            call netcdf_error (status, ncid_movie, x_id)
       status = nf90_put_var (ncid_movie, y_id, y)
       if (status /= NF90_NOERR) &
            call netcdf_error (status, ncid_movie, y_id)
       status = nf90_put_var (ncid_movie, th_id, theta)
       if (status /= NF90_NOERR) &
            call netcdf_error (status, ncid_movie, th_id)
       deallocate(x)


       deallocate(y)
    endif
# endif
  end subroutine nc_grids

  subroutine netcdf_init (serial_io2)
    use mp, only: proc0
    logical, intent(in) :: serial_io2
# ifdef NETCDF
    proc_write = proc0 .or. (.not.serial_io2)
# endif
  end subroutine netcdf_init

  subroutine nc_finish
    use mp, only: proc0
# ifdef NETCDF
    use netcdf, only: nf90_close
    use netcdf_utils, only: netcdf_error

    integer :: status

    if (proc0) then
       call save_input
       status = nf90_close (ncid)
       if (status /= NF90_NOERR) call netcdf_error (status)
       if(my_make_movie) then
          status = nf90_close (ncid_movie)
          if (status /= NF90_NOERR) call netcdf_error (status)
       endif
    end if
# endif
  end subroutine nc_finish

  subroutine save_input
    !<wkdoc> Save the input file in the NetCDF file </wkdoc>
# ifdef NETCDF    
    use file_utils, only: num_input_lines, get_input_unit
!    use netcdf_mod
    use netcdf, only: nf90_put_var

    character(200) line
    integer, dimension (2) :: nin_start, nin_count

    integer :: status, n, unit

    nin_start(1) = 1
    nin_start(2) = 1

    nin_count(2) = 1

    call get_input_unit (unit)
    rewind (unit=unit)
    do n = 1, num_input_lines
       read (unit=unit, fmt="(a)") line
       nin_count(1) = len(trim(line))
!       status = nf_put_vara_text (ncid, input_id, nin_start, nin_count, line)
       status = nf90_put_var (ncid, input_id, line, start=nin_start, count=nin_count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, input_id)
       nin_start(2) = nin_start(2) + 1
    end do
# endif
  end subroutine save_input

  subroutine define_vars (write_nl_flux, write_omega, write_stress, write_phiavg, &
       write_hrate, write_final_antot, write_eigenfunc, write_verr, &
       write_full_moments_notgc)

    use mp, only: nproc
    use species, only: nspec
    use kt_grids, only: naky, ntheta0, theta0
    use run_parameters, only: fphi, fapar, fbpar
# ifdef NETCDF
    use netcdf, only: NF90_CHAR, NF90_INT, NF90_GLOBAL
    use netcdf, only: nf90_def_var, nf90_put_att, nf90_enddef, nf90_put_var
    use netcdf, only: nf90_inq_libvers
    use netcdf_utils, only: netcdf_real
# endif
!    use netcdf_mod
!    logical :: write_nl_flux, write_omega, write_stress, write_phiavg, write_hrate
!    logical :: write_final_antot, write_eigenfunc, write_verr
    logical, intent(in) :: write_nl_flux, write_omega, write_stress
    logical, intent(in) :: write_phiavg, write_hrate, write_final_antot
    logical, intent(in) :: write_eigenfunc, write_verr
    logical, intent(in) :: write_full_moments_notgc
# ifdef NETCDF
    character (5) :: ci
    character (20) :: datestamp, timestamp, timezone
    !logical :: d_fields_per = .false. - unnecessary - now set in input file
    logical :: d_neo = .false.
    
    integer :: status

    fmode_dim(1) = nakx_dim
    fmode_dim(2) = naky_dim

    mode_dim (1) = nakx_dim
    mode_dim (2) = naky_dim
    mode_dim (3) = time_dim

    err_dim (1) = nttot_dim
    err_dim (2) = nakx_dim
    err_dim (3) = naky_dim
    err_dim (4) = ncoordt_dim
    err_dim (5) = time_dim

    lp_dim (1) = nttot_dim
    lp_dim (2) = nakx_dim
    lp_dim (3) = naky_dim
    lp_dim (4) = nspec_dim
    lp_dim (5) = ncoord_dim
    lp_dim (6) = time_dim

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

    fluxx_dim (1) = nakx_dim
    fluxx_dim (2) = nspec_dim
    fluxx_dim (3) = time_dim

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

    mom_dim(1) = ri_dim
    mom_dim(2) = nakx_dim
    mom_dim(3) = naky_dim
    mom_dim(4) = nspec_dim
    mom_dim(5) = time_dim


    !<wkdoc> Write some useful general information into the NetCDF file </wkdoc>
    status = nf90_put_att (ncid, NF90_GLOBAL, 'title', 'GS2 Simulation Data')
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, NF90_GLOBAL, att='title')
    status = nf90_put_att (ncid, NF90_GLOBAL, 'Conventions', &
         'http://gs2.sourceforge.net')
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, NF90_GLOBAL, att='Conventions')

    datestamp(:) = ' '
    timestamp(:) = ' '
    timezone(:) = ' '
    call date_and_time (datestamp, timestamp, timezone)

    
    status = nf90_def_var (ncid, 'code_info', NF90_CHAR, char10_dim, code_id)
    if (status /= NF90_NOERR) call netcdf_error (status, var='code_info')
    status = nf90_put_att (ncid, code_id, 'long_name', 'GS2')
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, code_id, att='long_name')

    ci = 'c1'
    status = nf90_put_att (ncid, code_id, trim(ci), 'Date: '//trim(datestamp))
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, code_id, att=ci)

    ci = 'c2'
    status = nf90_put_att (ncid, code_id, trim(ci), &
         'Time: '//trim(timestamp)//' '//trim(timezone))
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, code_id, att=ci)

    ci = 'c3'
    status = nf90_put_att (ncid, code_id, trim(ci), &
         'netCDF version '//trim(nf90_inq_libvers()))
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, code_id, att=ci)

    ci = 'c4'
    status = nf90_put_att (ncid, code_id, trim(ci), &
         'Units are determined with respect to reference temperature (T_ref),')
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, code_id, att=ci)

    ci = 'c5'
    status = nf90_put_att (ncid, code_id, trim(ci), &
         'reference charge (q_ref), reference mass (mass_ref),')
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, code_id, att=ci)

    ci = 'c6'
    status = nf90_put_att (ncid, code_id, trim(ci), &
         'reference field (B_ref), and reference length (a_ref)')
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, code_id, att=ci)

    ci = 'c7'
    status = nf90_put_att (ncid, code_id, trim(ci), &
         'from which one may construct rho_ref and vt_ref/a,')
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, code_id, att=ci)

    ci = 'c8'
    status = nf90_put_att (ncid, code_id, trim(ci), &
         'which are the basic units of perpendicular length and time.')
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, code_id, att=ci)

    ci = 'c9'
    status = nf90_put_att (ncid, code_id, trim(ci), &
         'Macroscopic lengths are normalized to the minor radius.')
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, code_id, att=ci)

    ci = 'c10'
    status = nf90_put_att (ncid, code_id, trim(ci), &
         'The difference between rho (normalized minor radius) and rho (gyroradius)')
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, code_id, att=ci)

    ci = 'c11'
    status = nf90_put_att (ncid, code_id, trim(ci), &
         'should be clear from the context in which they appear below.')
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, code_id, att=ci)

    !<wkdoc> Write lots of input variables into the NetCDF file </wkdoc>
    status = nf90_def_var (ncid, 'nproc', NF90_INT, nproc_id)
    if (status /= NF90_NOERR) call netcdf_error (status, var='nproc')
    status = nf90_put_att (ncid, nproc_id, 'long_name', 'Number of processors')
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, nproc_id, att='long_name')

    status = nf90_def_var (ncid, 'nmesh', netcdf_real, nmesh_id)
    if (status /= NF90_NOERR) call netcdf_error (status, var='nmesh')
    status = nf90_put_att (ncid, nmesh_id, 'long_name', 'Number of meshpoints')
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, nmesh_id, att='long_name')

    status = nf90_def_var (ncid, 'nkx', NF90_INT, nakx_id)
    if (status /= NF90_NOERR) call netcdf_error (status, var='nkx')
    status = nf90_def_var (ncid, 'nky', NF90_INT, naky_id)
    if (status /= NF90_NOERR) call netcdf_error (status, var='nky')
    status = nf90_def_var (ncid, 'ntheta_tot', NF90_INT, nttot_id)
    if (status /= NF90_NOERR) call netcdf_error (status, var='ntheta_tot')
    status = nf90_def_var (ncid, 'nspecies', NF90_INT, nspec_id)
    if (status /= NF90_NOERR) call netcdf_error (status, var='nspecies')

    status = nf90_def_var (ncid, 't', netcdf_real, time_dim, time_id)
    if (status /= NF90_NOERR) call netcdf_error (status, var='t')
    status = nf90_put_att (ncid, time_id, 'long_name', 'Time')
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, time_id, att='long_name')
    status = nf90_put_att (ncid, time_id, 'units', 'L/vt')
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, time_id, att='units')

    status = nf90_def_var (ncid, 'charge', NF90_INT, nspec_dim, charge_id)
    if (status /= NF90_NOERR) call netcdf_error (status, var='charge')
    status = nf90_put_att (ncid, charge_id, 'long_name', 'Charge')
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, charge_id, att='long_name')
    status = nf90_put_att (ncid, charge_id, 'units', 'q')
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, charge_id, att='units')

    status = nf90_def_var (ncid, 'mass', netcdf_real, nspec_dim, mass_id)
    if (status /= NF90_NOERR) call netcdf_error (status, var='mass')
    status = nf90_put_att (ncid, mass_id, 'long_name', 'Atomic mass')
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, mass_id, att='long_name')
    status = nf90_put_att (ncid, mass_id, 'units', 'm')
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, mass_id, att='units')

    status = nf90_def_var (ncid, 'dens', netcdf_real, nspec_dim, dens_id)
    if (status /= NF90_NOERR) call netcdf_error (status, var='dens')
    status = nf90_put_att (ncid, dens_id, 'long_name', 'Density')
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, dens_id, att='long_name')
    status = nf90_put_att (ncid, dens_id, 'units', 'n_e')
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, dens_id, att='units')

    status = nf90_def_var (ncid, 'temp', netcdf_real, nspec_dim, temp_id)
    if (status /= NF90_NOERR) call netcdf_error (status, var='temp')
    status = nf90_put_att (ncid, temp_id, 'long_name', 'Temperature')
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, temp_id, att='long_name')
    status = nf90_put_att (ncid, temp_id, 'units', 'T')
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, temp_id, att='units')

    status = nf90_def_var (ncid, 'tprim', netcdf_real, nspec_dim, tprim_id)
    if (status /= NF90_NOERR) call netcdf_error (status, var='tprim')
    status = nf90_put_att (ncid, tprim_id, 'long_name', '-1/rho dT/drho')
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, tprim_id, att='long_name')

    status = nf90_def_var (ncid, 'fprim', netcdf_real, nspec_dim, fprim_id) 
    if (status /= NF90_NOERR) call netcdf_error (status, var='fprim')
    status = nf90_put_att (ncid, fprim_id, 'long_name', '-1/rho dn/drho')
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, fprim_id, att='long_name')

    status = nf90_def_var (ncid, 'uprim', netcdf_real, nspec_dim, uprim_id)
    if (status /= NF90_NOERR) call netcdf_error (status, var='uprim')
    status = nf90_put_att (ncid, uprim_id, 'long_name', '-1/v_t du_par/drho')
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, uprim_id, att='long_name')

    status = nf90_def_var (ncid, 'uprim2', netcdf_real, nspec_dim, uprim2_id)
    if (status /= NF90_NOERR) call netcdf_error (status, var='uprim2')

    status = nf90_def_var (ncid, 'vnewk', netcdf_real, nspec_dim, vnewk_id)
    if (status /= NF90_NOERR) call netcdf_error (status, var='vnewk')
    status = nf90_put_att (ncid, vnewk_id, 'long_name', 'Collisionality')
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, vnewk_id, att='long_name')
    status = nf90_put_att (ncid, vnewk_id, 'units', 'v_t/L')
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, vnewk_id, att='units')
    
    status = nf90_def_var (ncid, 'type_of_species', NF90_INT, nspec_dim, spec_type_id)
    if (status /= NF90_NOERR) call netcdf_error (status, var='type_of_species')

    status = nf90_def_var (ncid, 'theta0', netcdf_real, fmode_dim, theta0_id)
    if (status /= NF90_NOERR) call netcdf_error (status, var='theta0')
    status = nf90_put_att (ncid, theta0_id, 'long_name', 'Theta_0')
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, theta0_id, att='long_name')

    status = nf90_def_var (ncid, 'kx', netcdf_real, nakx_dim, akx_id)
    if (status /= NF90_NOERR) call netcdf_error (status, var='kx')
    status = nf90_put_att (ncid, akx_id, 'long_name', 'kx rho')
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, akx_id, att='long_name')

    status = nf90_def_var (ncid, 'ky', netcdf_real, naky_dim, aky_id)
    if (status /= NF90_NOERR) call netcdf_error (status, var='ky')
    status = nf90_put_att (ncid, aky_id, 'long_name', 'ky rho')
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, aky_id, att='long_name')

    status = nf90_def_var (ncid, 'theta', netcdf_real, nttot_dim, theta_id)
    if (status /= NF90_NOERR) call netcdf_error (status, var='theta')

    status = nf90_def_var (ncid, 'bmag', netcdf_real, nttot_dim, bmag_id)
    if (status /= NF90_NOERR) call netcdf_error (status, var='bmag')
    status = nf90_put_att (ncid, bmag_id, 'long_name', '|B|(theta)')
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, bmag_id, att='long_name')
    status = nf90_put_att (ncid, bmag_id, 'units', 'B_0')
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, bmag_id, att='units')

    status = nf90_def_var (ncid, 'gradpar', netcdf_real, nttot_dim, gradpar_id)
    if (status /= NF90_NOERR) call netcdf_error (status, var='gradpar')
    status = nf90_def_var (ncid, 'gbdrift', netcdf_real, nttot_dim, gbdrift_id) 
    if (status /= NF90_NOERR) call netcdf_error (status, var='gbdrift')
    status = nf90_def_var (ncid, 'gbdrift0', netcdf_real, nttot_dim, gbdrift0_id)
    if (status /= NF90_NOERR) call netcdf_error (status, var='gbdrift0')
    status = nf90_def_var (ncid, 'cvdrift', netcdf_real, nttot_dim, cvdrift_id)
    if (status /= NF90_NOERR) call netcdf_error (status, var='cvdrift')
    status = nf90_def_var (ncid, 'cvdrift0', netcdf_real, nttot_dim, cvdrift0_id)
    if (status /= NF90_NOERR) call netcdf_error (status, var='cvdrift0')

    status = nf90_def_var (ncid, 'gds2', netcdf_real, nttot_dim, gds2_id)
    if (status /= NF90_NOERR) call netcdf_error (status, var='gds2')
    status = nf90_def_var (ncid, 'gds21', netcdf_real, nttot_dim, gds21_id)
    if (status /= NF90_NOERR) call netcdf_error (status, var='gds21')
    status = nf90_def_var (ncid, 'gds22', netcdf_real, nttot_dim, gds22_id)
    if (status /= NF90_NOERR) call netcdf_error (status, var='gds22')
    status = nf90_def_var (ncid, 'grho', netcdf_real, nttot_dim, grho_id)
    if (status /= NF90_NOERR) call netcdf_error (status, var='grho')
    status = nf90_def_var (ncid, 'jacob', netcdf_real, nttot_dim, jacob_id)
    if (status /= NF90_NOERR) call netcdf_error (status, var='jacob')

    status = nf90_def_var (ncid, 'q', netcdf_real, q_id)
    if (status /= NF90_NOERR) call netcdf_error (status, var='q')
    status = nf90_def_var (ncid, 'eps', netcdf_real, eps_id)
    if (status /= NF90_NOERR) call netcdf_error (status, var='eps')
    status = nf90_def_var (ncid, 'shat', netcdf_real, shat_id)
    if (status /= NF90_NOERR) call netcdf_error (status, var='shat')
    status = nf90_put_att (ncid, shat_id, 'long_name', '(rho/q) dq/drho')
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, shat_id, att='long_name')

    status = nf90_def_var (ncid, 'drhodpsi', netcdf_real, drhodpsi_id)
    if (status /= NF90_NOERR) call netcdf_error (status, var='drhodpsi')
    status = nf90_put_att (ncid, drhodpsi_id, 'long_name', 'drho/dPsi')
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, drhodpsi_id, att='long_name')

    if (write_verr) then
!       status = nf90_def_var (ncid, 'errest_by_mode', netcdf_real, err_dim, errest_by_mode_id)
!       if (status /= NF90_NOERR) call netcdf_error (status, var='errest_by_mode')
!       status = nf90_def_var (ncid, 'lpcoef_by_mode', netcdf_real, lp_dim, lpcoef_by_mode_id)
!       if (status /= NF90_NOERR) call netcdf_error (status, var='lpcoef_by_mode')
    end if

    if (fphi > zero) then
       status = nf90_def_var (ncid, 'phi2', netcdf_real, time_dim, phi2_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='phi2')
       status = nf90_put_att (ncid, phi2_id, 'long_name', '|Potential**2|')
       if (status /= NF90_NOERR) &
            call netcdf_error (status, ncid, phi2_id, att='long_name')
       status = nf90_put_att (ncid, phi2_id, 'units', '(T/q rho/L)**2')
       if (status /= NF90_NOERR) &
            call netcdf_error (status, ncid, phi2_id, att='units')

       status = nf90_def_var &
            (ncid, 'phi2_by_mode', netcdf_real, mode_dim, phi2_by_mode_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='phi2_by_mode')
       if (ntheta0 > 1) then
          status = nf90_def_var &
               (ncid, 'phi2_by_kx', netcdf_real, kx_dim, phi2_by_kx_id)
          if (status /= NF90_NOERR) &
               call netcdf_error (status, var='phi2_by_kx')
       end if

       if (naky > 1) then
          status = nf90_def_var &
               (ncid, 'phi2_by_ky', netcdf_real, ky_dim, phi2_by_ky_id)
          if (status /= NF90_NOERR) &
               call netcdf_error (status, var='phi2_by_ky')
       end if

       status = nf90_def_var (ncid, 'phi0', netcdf_real, omega_dim, phi0_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='phi0')
       if (write_phiavg) then
          status = nf90_def_var (ncid, 'phiavg', netcdf_real, omega_dim, phiavg_id)
          if (status /= NF90_NOERR) call netcdf_error (status, var='phiavg')
       endif
       if (write_nl_flux) then
          status = nf90_def_var (ncid, 'es_heat_par',  netcdf_real, flux_dim, es_heat_par_id)
          if (status /= NF90_NOERR) call netcdf_error (status, var='es_heat_par')
          status = nf90_def_var (ncid, 'es_heat_perp', netcdf_real, flux_dim, es_heat_perp_id)
          if (status /= NF90_NOERR) call netcdf_error (status, var='es_heat_perp')
          status = nf90_def_var (ncid, 'es_heat_flux', netcdf_real, flux_dim, es_heat_flux_id)
          if (status /= NF90_NOERR) call netcdf_error (status, var='es_heat_flux')
          status = nf90_def_var (ncid, 'es_mom_flux',  netcdf_real, flux_dim, es_mom_flux_id)
          if (status /= NF90_NOERR) call netcdf_error (status, var='es_mom_flux')
          status = nf90_def_var (ncid, 'es_part_flux', netcdf_real, flux_dim, es_part_flux_id)
          if (status /= NF90_NOERR) call netcdf_error (status, var='es_part_flux')
          status = nf90_def_var (ncid, 'es_heat_by_k', netcdf_real, fluxk_dim, es_heat_by_k_id)
          if (status /= NF90_NOERR) call netcdf_error (status, var='es_heat_by_k')
          status = nf90_def_var (ncid, 'es_mom_by_k',  netcdf_real, fluxk_dim, es_mom_by_k_id)
          if (status /= NF90_NOERR) call netcdf_error (status, var='es_mom_by_k')
          status = nf90_def_var (ncid, 'es_part_by_k', netcdf_real, fluxk_dim, es_part_by_k_id)
          if (status /= NF90_NOERR) call netcdf_error (status, var='es_part_by_k')
       end if
       if (write_eigenfunc) then
          status = nf90_def_var (ncid, 'phi_norm',  netcdf_real, final_field_dim, phi_norm_id)
          if (status /= NF90_NOERR) call netcdf_error (status, var='phi_norm')
       endif
       status = nf90_def_var (ncid, 'phi', netcdf_real, final_field_dim, phi_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='phi')
       status = nf90_put_att &
            (ncid, phi_id, 'long_name', 'Electrostatic Potential')
       if (status /= NF90_NOERR) &
            call netcdf_error (status, ncid, phi_id, att='long_name')
       status = nf90_put_att (ncid, phi_id, 'idl_name', '!7U!6')
       if (status /= NF90_NOERR) &
            call netcdf_error (status, ncid, phi_id, att='idl_name')
       status = nf90_put_att (ncid, phi_id, 'units', 'T/q rho/L')
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, phi_id, att='units')
       if (write_final_antot) then
          status = nf90_def_var (ncid, 'antot', netcdf_real, final_field_dim, antot_id)
          if (status /= NF90_NOERR) call netcdf_error (status, var='antot')
       endif
       if (write_phi_t) then
          status = nf90_def_var &
               (ncid, 'phi_t', netcdf_real, field_dim, phi_t_id)
          if (status /= NF90_NOERR) call netcdf_error (status, var='phi_t')
          status = nf90_put_att (ncid, phi_t_id, 'long_name', 'Electrostatic Potential over time')
          if (status /= NF90_NOERR) call netcdf_error (status, ncid, phi_t_id, att='long_name')
       end if
    end if

    if (fapar > zero) then
       status = nf90_def_var (ncid, 'apar2', netcdf_real, time_dim, apar2_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='apar2')
       status = nf90_def_var &
            (ncid, 'apar2_by_mode', netcdf_real, mode_dim, apar2_by_mode_id)
       if (status /= NF90_NOERR) &
            call netcdf_error (status, var='apar2_by_mode')
       if (ntheta0 > 1) then
          status = nf90_def_var &
               (ncid, 'apar2_by_kx', netcdf_real, kx_dim, apar2_by_kx_id)
          if (status /= NF90_NOERR) &
               call netcdf_error (status, var='apar2_by_kx')
       end if
       if (naky > 1) then
          status = nf90_def_var &
               (ncid, 'apar2_by_ky', netcdf_real, ky_dim, apar2_by_ky_id)
          if (status /= NF90_NOERR) &
               call netcdf_error (status, var='apar2_by_ky')
       end if

       status = nf90_def_var (ncid, 'apar0', netcdf_real, omega_dim, apar0_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='apar0')
       if (write_nl_flux) then
          status = nf90_def_var &
               (ncid,'apar_heat_flux',netcdf_real, flux_dim, apar_heat_flux_id)
          if (status /= NF90_NOERR) &
               call netcdf_error (status, var='apar_heat_flux')
          status = nf90_def_var &
               (ncid, 'apar_heat_par', netcdf_real, flux_dim, apar_heat_par_id)
          if (status /= NF90_NOERR) call netcdf_error (status, var='apar_heat_par')
          status = nf90_def_var (ncid, 'apar_heat_perp', netcdf_real, flux_dim, apar_heat_perp_id)
          if (status /= NF90_NOERR) call netcdf_error (status, var='apar_heat_perp')
          status = nf90_def_var (ncid, 'apar_mom_flux',  netcdf_real, flux_dim, apar_mom_flux_id)
          if (status /= NF90_NOERR) call netcdf_error (status, var='apar_mom_flux')
          status = nf90_def_var (ncid, 'apar_part_flux', netcdf_real, flux_dim, apar_part_flux_id)
          if (status /= NF90_NOERR) call netcdf_error (status, var='apar_part_flux')
          status = nf90_def_var (ncid, 'apar_heat_by_k', netcdf_real, fluxk_dim, apar_heat_by_k_id)
          if (status /= NF90_NOERR) call netcdf_error (status, var='apar_heat_by_k')
          status = nf90_def_var (ncid, 'apar_heat_by_x', netcdf_real, fluxx_dim, apar_heat_by_x_id)
          if (status /= NF90_NOERR) call netcdf_error (status, var='apar_heat_by_x')
          status = nf90_def_var (ncid, 'apar_mom_by_k',  netcdf_real, fluxk_dim, apar_mom_by_k_id)
          if (status /= NF90_NOERR) call netcdf_error (status, var='apar_mom_by_k')
          status = nf90_def_var (ncid, 'apar_part_by_k', netcdf_real, fluxk_dim, apar_part_by_k_id)
          if (status /= NF90_NOERR) call netcdf_error (status, var='apar_part_by_k')
       end if
       if (write_eigenfunc) then
          status = nf90_def_var (ncid, 'apar_norm', netcdf_real, final_field_dim, apar_norm_id)
          if (status /= NF90_NOERR) call netcdf_error (status, var='apar_norm')
       endif
       status = nf90_def_var (ncid, 'apar', netcdf_real, final_field_dim, apar_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='apar')
       if (write_final_antot) then
          status = nf90_def_var (ncid, 'antota', netcdf_real, final_field_dim, antota_id)
          if (status /= NF90_NOERR) call netcdf_error (status, var='antota')
       endif
       if (write_apar_t) then
          status = nf90_def_var (ncid, 'apar_t', netcdf_real, field_dim, apar_t_id)
          if (status /= NF90_NOERR) call netcdf_error (status, var='apar_t')
          status = nf90_put_att (ncid, apar_t_id, 'long_name', 'Parallel Magnetic Potential over Time')
          if (status /= NF90_NOERR) call netcdf_error (status, ncid, apar_t_id, att='long_name')
       end if
       status = nf90_put_att (ncid, apar2_by_mode_id, 'long_name', 'Apar squared')
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, apar2_by_mode_id, att='long_name')
       status = nf90_put_att (ncid, apar_id, 'long_name', 'Parallel Magnetic Potential')
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, apar_id, att='long_name')
       status = nf90_put_att (ncid, apar_id, 'idl_name', '!6A!9!D#!N!6')
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, apar_id, att='idl_name')
       status = nf90_put_att (ncid, apar2_id, 'long_name', 'Total A_par squared')
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, apar2_id, att='long_name')
    end if

    if (fbpar > zero) then
       status = nf90_def_var (ncid, 'bpar2', netcdf_real, time_dim, bpar2_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='bpar2')
       status = nf90_def_var (ncid, 'bpar2_by_mode', netcdf_real, mode_dim, bpar2_by_mode_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='bpar2_by_mode')
       if (ntheta0 > 1) then
          status = nf90_def_var (ncid, 'bpar2_by_kx', netcdf_real, kx_dim, bpar2_by_kx_id)
          if (status /= NF90_NOERR) call netcdf_error (status, var='bpar2_by_kx')
       end if
       if (naky > 1) then
          status = nf90_def_var (ncid, 'bpar2_by_ky', netcdf_real, ky_dim, bpar2_by_ky_id)
          if (status /= NF90_NOERR) call netcdf_error (status, var='bpar2_by_ky')
       end if
       status = nf90_def_var (ncid, 'bpar0', netcdf_real, omega_dim, bpar0_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='bpar0')
       if (write_nl_flux) then
          status = nf90_def_var (ncid, 'bpar_heat_flux', netcdf_real, flux_dim, bpar_heat_flux_id)
          if (status /= NF90_NOERR) call netcdf_error (status, var='bpar_heat_flux')
          status = nf90_def_var (ncid, 'bpar_heat_par', netcdf_real, flux_dim, bpar_heat_par_id)
          if (status /= NF90_NOERR) call netcdf_error (status, var='bpar_heat_par')
          status = nf90_def_var (ncid, 'bpar_heat_perp', netcdf_real, flux_dim, bpar_heat_perp_id)
          if (status /= NF90_NOERR) call netcdf_error (status, var='bpar_heat_perp')
          status = nf90_def_var (ncid, 'bpar_mom_flux', netcdf_real, flux_dim, bpar_mom_flux_id)
          if (status /= NF90_NOERR) call netcdf_error (status, var='bpar_mom_flux')
          status = nf90_def_var (ncid, 'bpar_part_flux', netcdf_real, flux_dim, bpar_part_flux_id)
          if (status /= NF90_NOERR) call netcdf_error (status, var='bpar_part_flux')
          status = nf90_def_var (ncid, 'bpar_heat_by_k', netcdf_real, fluxk_dim, bpar_heat_by_k_id)
          if (status /= NF90_NOERR) call netcdf_error (status, var='bpar_heat_by_k')
          status = nf90_def_var (ncid, 'bpar_mom_by_k', netcdf_real, fluxk_dim, bpar_mom_by_k_id)
          if (status /= NF90_NOERR) call netcdf_error (status, var='bpar_mom_by_k')
          status = nf90_def_var (ncid, 'bpar_part_by_k', netcdf_real, fluxk_dim, bpar_part_by_k_id)
          if (status /= NF90_NOERR) call netcdf_error (status, var='bpar_part_by_k')
       end if
       if (write_eigenfunc) then
          status = nf90_def_var (ncid, 'bpar_norm', netcdf_real, final_field_dim, bpar_norm_id)
          if (status /= NF90_NOERR) call netcdf_error (status, var='bpar_norm')
       endif
       status = nf90_def_var (ncid, 'bpar', netcdf_real, final_field_dim, bpar_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='bpar')
       if (write_final_antot) then
          status = nf90_def_var (ncid, 'antotp', netcdf_real, final_field_dim, antotp_id)
          if (status /= NF90_NOERR) call netcdf_error (status, var='antotp')
       endif
       if (write_bpar_t) then
          status = nf90_def_var (ncid, 'bpar_t', netcdf_real, field_dim, bpar_t_id)
          if (status /= NF90_NOERR) call netcdf_error (status, var='bpar_t')
          status = nf90_put_att (ncid, bpar_t_id, 'long_name', 'delta B Parallel over time')
          if (status /= NF90_NOERR) call netcdf_error (status, ncid, bpar_t_id, att='long_name')
       end if

       status = nf90_put_att (ncid, bpar2_by_mode_id, 'long_name', 'A_perp squared')
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, bpar2_by_mode_id, att='long_name')
       status = nf90_put_att (ncid, bpar_id, 'long_name', 'delta B Parallel')
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, bpar_id, att='long_name')
       status = nf90_put_att (ncid, bpar_id, 'idl_name', '!6B!9!D#!N!6')
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, bpar_id, att='idl_name')
       status = nf90_put_att (ncid, bpar2_id, 'long_name', 'Total A_perp squared')
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, bpar2_id, att='long_name')
    end if

    status = nf90_def_var (ncid, 'phase', netcdf_real, phase_dim, phase_id)
    if (status /= NF90_NOERR) call netcdf_error (status, var='phase')
    status = nf90_put_att (ncid, phase_id, 'long_name', 'Normalizing phase')
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, phase_id, att='long_name')

!    status = nf90_def_var (ncid, 'phtot', netcdf_real, mode_dim, phtot_id)
!    if (status /= NF90_NOERR) call netcdf_error (status, var='phtot')
!    status = nf90_def_var (ncid, 'dmix', netcdf_real, mode_dim, dmix_id)
!    if (status /= NF90_NOERR) call netcdf_error (status, var='dmix')
!    status = nf90_def_var (ncid, 'kperpnorm', netcdf_real, mode_dim, kperpnorm_id)
!    if (status /= NF90_NOERR) call netcdf_error (status, var='kperpnorm')

    if (write_omega) then
       status = nf90_def_var (ncid, 'omega', netcdf_real, omega_dim, omega_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='omega')
       status = nf90_def_var (ncid, 'omegaavg', netcdf_real, omega_dim, omegaavg_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='omegaavg')
    end if

    if (write_nl_flux) then
       status = nf90_def_var (ncid, 'hflux_tot', netcdf_real, time_dim, hflux_tot_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='kflux_tot')
       status = nf90_def_var (ncid, 'vflux_tot', netcdf_real, time_dim, vflux_tot_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='vflux_tot')
       status = nf90_def_var (ncid, 'zflux_tot', netcdf_real, time_dim, zflux_tot_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='zflux_tot')
    end if

    if (d_neo) then
       status = nf90_def_var (ncid, 'qflux_neo_by_k', netcdf_real, fluxk_dim, qflux_neo_by_k_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='qflux_neo_by_k')
       status = nf90_def_var (ncid, 'pflux_neo_by_k', netcdf_real, fluxk_dim, pflux_neo_by_k_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='pflux_neo_by_k')
       status = nf90_def_var (ncid, 'sourcefac', netcdf_real, om_dim, sourcefac_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='sourcefac')
    end if

    status = nf90_def_var (ncid, 'epar', netcdf_real, final_field_dim, epar_id)
    if (status /= NF90_NOERR) call netcdf_error (status, var='epar')

!    <phi>

    status = nf90_def_var (ncid, 'ntot2', netcdf_real, flux_dim,  ntot2_id)
    if (status /= NF90_NOERR) call netcdf_error (status, var='ntot2')
    status = nf90_def_var (ncid, 'ntot2_by_mode', netcdf_real, fluxk_dim, ntot2_by_mode_id)
    if (status /= NF90_NOERR) call netcdf_error (status, var='ntot2_by_mode')

    status = nf90_def_var (ncid, 'ntot20', netcdf_real, flux_dim,  ntot20_id)
    if (status /= NF90_NOERR) call netcdf_error (status, var='ntot20')
    status = nf90_put_att (ncid, ntot20_id, 'long_name', 'Density**2 at theta=0')
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, ntot20_id, att='long_name')
    status = nf90_def_var (ncid, 'ntot20_by_mode', netcdf_real, fluxk_dim, ntot20_by_mode_id)
    if (status /= NF90_NOERR) call netcdf_error (status, var='ntot20_by_mode')

    status = nf90_def_var (ncid, 'ntot', netcdf_real, final_mom_dim, ntot_id)
    if (status /= NF90_NOERR) call netcdf_error (status, var='ntot')
    status = nf90_def_var (ncid, 'density', netcdf_real, final_mom_dim, density_id)
    if (status /= NF90_NOERR) call netcdf_error (status, var='density')
    status = nf90_def_var (ncid, 'upar', netcdf_real, final_mom_dim, upar_id)
    if (status /= NF90_NOERR) call netcdf_error (status, var='upar')
    status = nf90_def_var (ncid, 'tpar', netcdf_real, final_mom_dim, tpar_id)
    if (status /= NF90_NOERR) call netcdf_error (status, var='tpar')
    status = nf90_def_var (ncid, 'tperp', netcdf_real, final_mom_dim, tperp_id)
    if (status /= NF90_NOERR) call netcdf_error (status, var='tperp')

    status = nf90_def_var (ncid, 'phi00', netcdf_real, loop_phi_dim, phi00_id)
    if (status /= NF90_NOERR) call netcdf_error (status, var='phi00')
    status = nf90_def_var (ncid, 'ntot00', netcdf_real, loop_mom_dim, ntot00_id)
    if (status /= NF90_NOERR) call netcdf_error (status, var='ntot00')
    status = nf90_def_var (ncid, 'density00', netcdf_real, loop_mom_dim, density00_id)
    if (status /= NF90_NOERR) call netcdf_error (status, var='density00')
    status = nf90_def_var (ncid, 'upar00', netcdf_real, loop_mom_dim, upar00_id)
    if (status /= NF90_NOERR) call netcdf_error (status, var='upar00')
    status = nf90_def_var (ncid, 'tpar00', netcdf_real, loop_mom_dim, tpar00_id)
    if (status /= NF90_NOERR) call netcdf_error (status, var='tpar00')
    status = nf90_def_var (ncid, 'tperp00', netcdf_real, loop_mom_dim, tperp00_id)
    if (status /= NF90_NOERR) call netcdf_error (status, var='tperp00')
    
    if (write_stress) then
       status = nf90_def_var (ncid, 'rstress', netcdf_real, loop_mom_dim, rstress_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='rstress')
       status = nf90_def_var (ncid, 'ustress', netcdf_real, loop_mom_dim, ustress_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='ustress')
    end if

    if (write_hrate) then
       status = nf90_def_var (ncid, 'h_energy', netcdf_real, time_dim, h_energy_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='h_energy')
       status = nf90_def_var (ncid, 'h_energy_dot', netcdf_real, time_dim, h_energy_dot_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='h_energy_dot')
       status = nf90_def_var (ncid, 'h_antenna', netcdf_real, time_dim, h_antenna_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='h_antenna')
       status = nf90_def_var (ncid, 'h_eapar', netcdf_real, time_dim, h_eapar_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='h_eapar')
       status = nf90_def_var (ncid, 'h_ebpar', netcdf_real, time_dim, h_ebpar_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='h_ebpar')

       status = nf90_def_var (ncid, 'h_delfs2', netcdf_real, flux_dim, h_delfs2_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='h_delfs2')
       status = nf90_def_var (ncid, 'h_hs2', netcdf_real, flux_dim, h_hs2_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='h_hs2')
       status = nf90_def_var (ncid, 'h_phis2', netcdf_real, flux_dim, h_phis2_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='h_phis2')
       status = nf90_def_var (ncid, 'h_hypervisc', netcdf_real, flux_dim, h_hypervisc_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='h_hypervisc')
       status = nf90_def_var (ncid, 'h_hyperres', netcdf_real, flux_dim, h_hyperres_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='h_hyperres')
       status = nf90_def_var (ncid, 'h_collisions', netcdf_real, flux_dim, h_collisions_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='h_collisions')
       status = nf90_def_var (ncid, 'h_gradients', netcdf_real, flux_dim, h_gradients_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='h_gradients')
       status = nf90_def_var (ncid, 'h_hypercoll', netcdf_real, flux_dim, h_hypercoll_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='h_hypercoll')
       status = nf90_def_var (ncid, 'h_imp_colls', netcdf_real, flux_dim, h_imp_colls_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='h_imp_colls')
       status = nf90_def_var (ncid, 'h_heating', netcdf_real, flux_dim, h_heating_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='h_heating')
       status = nf90_put_att (ncid, h_heating_id, 'long_name', 'Total heating by species')
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, h_heating_id, att='long_name')

       status = nf90_def_var (ncid, 'hk_energy', netcdf_real, mode_dim, hk_energy_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='h_heating')
       status = nf90_def_var (ncid, 'hk_energy_dot', netcdf_real, mode_dim, hk_energy_dot_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='h_heating')
       status = nf90_def_var (ncid, 'hk_antenna', netcdf_real, mode_dim, hk_antenna_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='h_heating')
       status = nf90_def_var (ncid, 'hk_eapar', netcdf_real, mode_dim, hk_eapar_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='h_heating')
       status = nf90_def_var (ncid, 'hk_ebpar', netcdf_real, mode_dim, hk_ebpar_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='h_heating')

       status = nf90_def_var (ncid, 'hk_delfs2', netcdf_real, fluxk_dim, hk_delfs2_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='hk_delfs2')
       status = nf90_def_var (ncid, 'hk_hs2', netcdf_real, fluxk_dim, hk_hs2_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='hk_hs2')
       status = nf90_def_var (ncid, 'hk_phis2', netcdf_real, fluxk_dim, hk_phis2_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='hk_phis2')
       status = nf90_def_var (ncid, 'hk_hypervisc', netcdf_real, fluxk_dim, hk_hypervisc_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='hk_hypervisc')
       status = nf90_def_var (ncid, 'hk_hyperres', netcdf_real, fluxk_dim, hk_hyperres_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='hk_hyperres')
       status = nf90_def_var (ncid, 'hk_collisions', netcdf_real, fluxk_dim, hk_collisions_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='hk_collisions')
       status = nf90_def_var (ncid, 'hk_gradients', netcdf_real, fluxk_dim, hk_gradients_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='hk_gradients')
       status = nf90_def_var (ncid, 'hk_imp_colls', netcdf_real, fluxk_dim, hk_imp_colls_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='hk_imp_colls')
       status = nf90_def_var (ncid, 'hk_hypercoll', netcdf_real, fluxk_dim, hk_hypercoll_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='hk_hypercoll')
       status = nf90_def_var (ncid, 'hk_heating', netcdf_real, fluxk_dim, hk_heating_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='hk_heating')
       status = nf90_put_att (ncid, hk_heating_id, 'long_name', 'Total heating by species and mode')
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, hk_heating_id, att='long_name')
    end if

    ! RN> not guiding center moments
    if (write_full_moments_notgc) then
       status = nf90_def_var (ncid, 'ntot0', netcdf_real, mom_dim, ntot0_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='ntot0')
       status = nf90_def_var (ncid, 'density0', netcdf_real, mom_dim, density0_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='density0')
       status = nf90_def_var (ncid, 'upar0', netcdf_real, mom_dim, upar0_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='upar0')
       status = nf90_def_var (ncid, 'tpar0', netcdf_real, mom_dim, tpar0_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='tpar0')
       status = nf90_def_var (ncid, 'tperp0', netcdf_real, mom_dim, tperp0_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='tperp0')
    endif
    ! <RN

!    status = nf90_put_att (ncid, phtot_id, 'long_name', 'Field amplitude')

    status = nf90_def_var (ncid, 'input_file', NF90_CHAR, nin_dim, input_id)
    if (status /= NF90_NOERR) call netcdf_error (status, var='input_file')
    status = nf90_put_att (ncid, input_id, 'long_name', 'Input file')
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, input_id, att='long_name')

    status = nf90_enddef (ncid)  ! out of definition mode
    if (status /= NF90_NOERR) call netcdf_error (status)

    status = nf90_put_var (ncid, nproc_id, nproc)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, nproc_id)

        ! added by EAB 03/05/04 for GS2 movies
    if(my_make_movie) then
       xmode_dim (1) = nx_dim
       xmode_dim (2) = ny_dim
       xmode_dim (3) = nth_dim
       xmode_dim (4) = time_movie_dim
       status = nf90_put_att (ncid_movie, NF90_GLOBAL, 'title', 'GS2 Simulation x,y,theta Data for Movies')
       if (status /= NF90_NOERR) call netcdf_error (status, ncid_movie, NF90_GLOBAL, att='title')
       status = nf90_put_att (ncid_movie, NF90_GLOBAL, 'Conventions', &
            'http://gs2.sourceforge.net')
       if (status /= NF90_NOERR) call netcdf_error (status, ncid_movie, NF90_GLOBAL, att='Conventions')
       !    status = nf90_def_var (ncid_movie, 'code_info_movie', NF90_CHAR, char10_dim, code_id_movie)
       !    if (status /= NF90_NOERR) call netcdf_error (status, var='code_info')
       !    status = nf90_put_att (ncid_movie, code_id_movie, 'long_name', 'GS2')
       !    if (status /= NF90_NOERR) call netcdf_error (status, ncid_movie, code_id_movie, att='long_name')
       !    ci = 'c1'
       !    status = nf90_put_att (ncid_movie, code_id_movie, trim(ci), 'Date: '//trim(datestamp))
       !    if (status /= NF90_NOERR) call netcdf_error (status, ncid_movie, code_id_movie, att=ci)
       !    ci = 'c2'
       !    status = nf90_put_att (ncid_movie, code_id_movie, trim(ci), &
       !         'Time: '//trim(timestamp)//' '//trim(timezone))
       !    if (status /= NF90_NOERR) call netcdf_error (status, ncid_movie, code_id_movie, att=ci)
       !    ci = 'c3'
       !    status = nf90_put_att (ncid_movie, code_id_movie, trim(ci), &
       !         'netCDF version '//trim(nf90_inq_libvers()))
       !    if (status /= NF90_NOERR) call netcdf_error (status, ncid_movie, code_id_movie, att=ci)
       
       status = nf90_def_var (ncid_movie, 'nx', NF90_INT, nx_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='nx')
       status = nf90_def_var (ncid_movie, 'ny', NF90_INT, ny_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='ny')
       status = nf90_def_var (ncid_movie, 'ntheta', NF90_INT, nth_id)    
       if (status /= NF90_NOERR) call netcdf_error (status, var='ntheta')
       status = nf90_def_var (ncid_movie, 'tm', netcdf_real, time_movie_dim, time_movie_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='tm')
       status = nf90_put_att (ncid_movie, time_movie_id, 'long_name', 'Time')
       if (status /= NF90_NOERR) call netcdf_error (status, ncid_movie, time_movie_id, att='long_name')
       status = nf90_put_att (ncid_movie, time_movie_id, 'units', 'L/vt')
       if (status /= NF90_NOERR) call netcdf_error (status, ncid_movie, time_movie_id, att='units')
       status = nf90_def_var (ncid_movie, 'x', netcdf_real, nx_dim, x_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='x')
       status = nf90_put_att (ncid_movie, x_id, 'long_name', 'x / rho')
       if (status /= NF90_NOERR) call netcdf_error (status, ncid_movie, x_id, att='long_name')
       status = nf90_def_var (ncid_movie, 'y', netcdf_real, ny_dim, y_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='y')
       status = nf90_put_att (ncid_movie, y_id, 'long_name', 'y / rho')
       if (status /= NF90_NOERR) call netcdf_error (status, ncid_movie, y_id, att='long_name')
       status = nf90_def_var (ncid_movie, 'theta', netcdf_real, nth_dim, th_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='theta')
       if(fphi > zero) then
          status = nf90_def_var (ncid_movie, 'phi_by_xmode', netcdf_real, xmode_dim, phi_by_xmode_id)
          if (status /= NF90_NOERR) call netcdf_error (status, var='phi_by_xmode')
       end if
! TEMPORARY BD
       status = nf90_def_var (ncid_movie, 'density_by_xmode', netcdf_real, xmode_dim, apar_by_xmode_id)
       if (status /= NF90_NOERR) call netcdf_error (status, var='density_by_xmode')
!       if(fapar > zero) then
!          status = nf90_def_var (ncid_movie, 'apar_by_xmode', netcdf_real, xmode_dim, apar_by_xmode_id)
!          if (status /= NF90_NOERR) call netcdf_error (status, var='apar_by_xmode')
!       end if
       if(fbpar > zero) then
          status = nf90_def_var (ncid_movie, 'bpar_by_xmode', netcdf_real, xmode_dim, bpar_by_xmode_id)
          if (status /= NF90_NOERR) call netcdf_error (status, var='bpar_by_xmode')
       end if

       status = nf90_enddef (ncid_movie)  ! out of definition mode
       if (status /= NF90_NOERR) call netcdf_error (status)
    endif

# endif
  end subroutine define_vars

  subroutine nc_eigenfunc (phase)

!    use nf90_mod, only: nf90_put_var
    use convert, only: c2r
    use run_parameters, only: fphi, fapar, fbpar
    use fields_arrays, only: phi, apar, bpar
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
# ifdef NETCDF
    use netcdf, only: nf90_put_var
# endif
    complex, dimension(:,:), intent (in) :: phase
# ifdef NETCDF
    complex, dimension(-ntgrid:ntgrid, ntheta0, naky) :: tmp
    real, dimension(2, ntheta0, naky) :: ri2
    real, dimension (2, 2*ntgrid+1, ntheta0, naky) :: ri3
    integer :: status, ig

    call c2r (phase, ri2)
    status = nf90_put_var (ncid, phase_id, ri2)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, phase_id)

    if (fphi > zero) then
       do ig = -ntgrid, ntgrid
          tmp(ig,:,:) = phi(ig,:,:)/phase(:,:)
       end do
       call c2r (tmp, ri3)
       status = nf90_put_var(ncid, phi_norm_id, ri3)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, phi_norm_id)
    end if

    if (fapar > zero) then
       do ig = -ntgrid, ntgrid
          tmp(ig,:,:) = apar(ig,:,:)/phase(:,:)
       end do
       call c2r (tmp, ri3)
       status = nf90_put_var(ncid, apar_norm_id, ri3)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, apar_norm_id)
    end if

    if (fbpar > zero) then
       do ig = -ntgrid, ntgrid
          tmp(ig,:,:) = bpar(ig,:,:)/phase(:,:)
       end do
       call c2r (tmp, ri3)
       status = nf90_put_var(ncid, bpar_norm_id, ri3)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, bpar_norm_id)
    end if
# endif
  end subroutine nc_eigenfunc

  subroutine nc_final_fields

!    use nf90_mod, only: nf90_put_var
    use convert, only: c2r
    use run_parameters, only: fphi, fapar, fbpar
    use fields_arrays, only: phi, apar, bpar
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
# ifdef NETCDF
    use netcdf, only: nf90_put_var

    real, dimension (2, 2*ntgrid+1, ntheta0, naky) :: ri3
    integer :: status

    if (fphi > zero) then
       call c2r (phi, ri3)
       status = nf90_put_var (ncid, phi_id, ri3)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, phi_id)
    end if

    if (fapar > zero) then
       call c2r (apar, ri3)
       status = nf90_put_var (ncid, apar_id, ri3)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, apar_id)
    end if

    if (fbpar > zero) then
       call c2r (bpar, ri3)
       status = nf90_put_var (ncid, bpar_id, ri3)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, bpar_id)
    end if
# endif
  end subroutine nc_final_fields

  subroutine nc_final_epar (epar)

!    use nf90_mod, only: nf90_put_var
    use convert, only: c2r
    use run_parameters, only: fphi, fapar, fbpar
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
# ifdef NETCDF
    use netcdf, only: nf90_put_var
# endif
    complex, dimension (:,:,:), intent (in) :: epar
# ifdef NETCDF
    real, dimension (2, 2*ntgrid+1, ntheta0, naky) :: ri3
    integer :: status

    call c2r (epar, ri3)
    status = nf90_put_var (ncid, epar_id, ri3)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, epar_id)
# endif
  end subroutine nc_final_epar

  subroutine nc_final_moments (ntot, density, upar, tpar, tperp)

!    use nf90_mod, only: nf90_put_var
    use convert, only: c2r
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use species, only: nspec
# ifdef NETCDF
    use netcdf, only: nf90_put_var
# endif
    complex, dimension (:,:,:,:), intent (in) :: ntot, density, upar, tpar, tperp
# ifdef NETCDF
    real, dimension (2, 2*ntgrid+1, ntheta0, naky, nspec) :: ri4
    integer :: status

    call c2r (ntot, ri4)
    status = nf90_put_var (ncid, ntot_id, ri4)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, ntot_id)

    call c2r (density, ri4)
    status = nf90_put_var (ncid, density_id, ri4)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, density_id)

    call c2r (upar, ri4)
    status = nf90_put_var (ncid, upar_id, ri4)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, upar_id)

    call c2r (tpar, ri4)
    status = nf90_put_var (ncid, tpar_id, ri4)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, tpar_id)

    call c2r (tperp, ri4)
    status = nf90_put_var (ncid, tperp_id, ri4)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, tperp_id)
# endif
  end subroutine nc_final_moments

  subroutine nc_loop_stress (nout, rstress, ustress)

!    use nf90_mod, only: nf90_put_vara
    use convert, only: c2r
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use species, only: nspec
# ifdef NETCDF
    use netcdf, only: nf90_put_var
# endif
    integer, intent (in) :: nout
    complex, dimension (:,:), intent (in) :: rstress, ustress
# ifdef NETCDF
    real, dimension (2, ntheta0, nspec) :: ri2
    integer, dimension (4) :: start, count
    integer :: status

    start(1) = 1
    start(2) = 1
    start(3) = 1
    start(4) = nout
    
    count(1) = 2
    count(2) = ntheta0
    count(3) = nspec
    count(4) = 1

    call c2r (rstress, ri2)
    status = nf90_put_var (ncid, rstress_id, ri2, start=start, count=count)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, nttot_id)

    call c2r (ustress, ri2)
    status = nf90_put_var (ncid, ustress_id, ri2, start=start, count=count)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, ustress_id)
# endif
  end subroutine nc_loop_stress

  subroutine nc_loop_vres (nout, errest_by_mode, lpcoef_by_mode)

!    use nf90_mod, only: nf90_put_vara
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use species, only: nspec

    integer, intent (in) :: nout
    real, dimension(:,:,:,:), intent (in) :: errest_by_mode
    real, dimension(:,:,:,:,:), intent (in) :: lpcoef_by_mode
# ifdef NETCDF
    integer, dimension (6) :: start2, count2
    integer, dimension (5) :: start1, count1
!    integer :: status

    start1(1) = 1
    start1(2) = 1
    start1(3) = 1
    start1(4) = 1
    start1(5) = nout

    start2(1) = 1
    start2(2) = 1
    start2(3) = 1
    start2(4) = 1
    start2(5) = 1
    start2(6) = nout
  
    count1(1) = 2*ntgrid+1
    count1(2) = ntheta0
    count1(3) = naky
    count1(4) = 4
    count1(5) = 1

    count2(1) = 2*ntgrid+1
    count2(2) = ntheta0
    count2(3) = naky
    count2(4) = nspec
    count2(5) = 2
    count2(6) = 1

!    status = nf90_put_var (ncid, errest_by_mode_id, start=start1, count=count1, errest_by_mode)
!    if (status /= NF90_NOERR) call netcdf_error (status, ncid, errest_by_mode_id)
!    status = nf90_put_var (ncid, lpcoef_by_mode_id, start=start2, count=count2, lpcoef_by_mode)
!    if (status /= NF90_NOERR) call netcdf_error (status, ncid, lpcoef_by_mode_id)
# endif
  end subroutine nc_loop_vres

  subroutine nc_loop_moments (nout, ntot2, ntot2_by_mode, ntot20, ntot20_by_mode, &
       phi00, ntot00, density00, upar00, tpar00, tperp00)

!    use nf90_mod, only: nf90_put_vara
    use convert, only: c2r
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use species, only: nspec
# ifdef NETCDF
    use netcdf, only: nf90_put_var
# endif
    integer, intent (in) :: nout
    real, dimension (:), intent (in) :: ntot2, ntot20
    real, dimension (:,:,:), intent (in) :: ntot2_by_mode, ntot20_by_mode
    complex, dimension (:), intent (in) :: phi00
    complex, dimension (:,:), intent (in) :: ntot00, density00, upar00, tpar00, tperp00
# ifdef NETCDF
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

    status = nf90_put_var (ncid, ntot2_id, ntot2, start=start, count=count)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, ntot2_id)
    status = nf90_put_var (ncid, ntot2_by_mode_id, ntot2_by_mode, start=start4, count=count4)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, ntot2_by_mode_id)

    status = nf90_put_var (ncid, ntot20_id, ntot20, start=start, count=count)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, ntot20_id)
    status = nf90_put_var (ncid, ntot20_by_mode_id, ntot20_by_mode, start=start4, count=count4)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, ntot20_by_mode_id)

    call c2r (phi00, ri1)
    status = nf90_put_var (ncid, phi00_id, ri1, start=start3, count=count3)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, phi00_id)
    
    call c2r (ntot00, ri2)
    status = nf90_put_var (ncid, ntot00_id, ri2, start=start00, count=count00)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, ntot00_id)
    
    call c2r (density00, ri2)
    status = nf90_put_var (ncid, density00_id, ri2, start=start00, count=count00)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, density00_id)
    
    call c2r (upar00, ri2)
    status = nf90_put_var (ncid, upar00_id, ri2, start=start00, count=count00)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, upar00_id)
    
    call c2r (tpar00, ri2)
    status = nf90_put_var (ncid, tpar00_id, ri2, start=start00, count=count00)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, tpar00_id)

    call c2r (tperp00, ri2)
    status = nf90_put_var (ncid, tperp00_id, ri2, start=start00, count=count00)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, tperp00_id)
# endif
  end subroutine nc_loop_moments

  subroutine nc_final_an (antot, antota, antotp)

!    use nf90_mod, only: nf90_put_var
    use convert, only: c2r
    use run_parameters, only: fphi, fapar, fbpar
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
# ifdef NETCDF
    use netcdf, only: nf90_put_var
# endif
    complex, dimension (:,:,:) :: antot, antota, antotp ! intent?
# ifdef NETCDF
    real, dimension (2, 2*ntgrid+1, ntheta0, naky) :: ri3
    integer :: status

    if (fphi > zero) then
       call c2r (antot, ri3)
       status = nf90_put_var(ncid, antot_id, ri3)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, antot_id)
    end if

    if (fapar > zero) then
       call c2r (antota, ri3)
       status = nf90_put_var(ncid, antota_id, ri3)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, antota_id)
    end if

    if (fbpar > zero) then
       call c2r (antotp, ri3)
       status = nf90_put_var(ncid, antotp_id, ri3)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, antotp_id)
    end if
# endif
  end subroutine nc_final_an

  subroutine nc_qflux (nout, qheat, qmheat, qbheat, &
       heat_par,  mheat_par,  bheat_par, &
       heat_perp, mheat_perp, bheat_perp, &
       heat_fluxes, mheat_fluxes, bheat_fluxes, x_qmflux, hflux_tot)

    use species, only: nspec
    use kt_grids, only: naky, ntheta0
    use run_parameters, only: fphi, fapar, fbpar
!    use nf90_mod, only: nf90_put_vara, nf90_put_var1
# ifdef NETCDF
    use netcdf, only: nf90_put_var
# endif
    integer, intent (in) :: nout
    real, dimension (:,:,:), intent (in) :: qheat, qmheat, qbheat
    real, dimension (:), intent (in) :: heat_par, mheat_par, bheat_par
    real, dimension (:), intent (in) :: heat_perp, mheat_perp, bheat_perp
    real, dimension (:), intent (in) :: heat_fluxes, mheat_fluxes, bheat_fluxes
    real, dimension (:,:), intent (in) :: x_qmflux
    real, intent (in) :: hflux_tot
# ifdef NETCDF
    integer, dimension (2) :: start, count
    integer, dimension (3) :: start3, count3
    integer, dimension (4) :: start4, count4
    integer :: status

    start3(1) = 1
    start3(2) = 1
    start3(3) = nout

    count3(1) = ntheta0
    count3(2) = nspec
    count3(3) = 1    

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
       status = nf90_put_var (ncid, es_heat_flux_id, heat_fluxes, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, es_heat_flux_id)
       status = nf90_put_var (ncid, es_heat_par_id, heat_par, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, es_heat_par_id)
       status = nf90_put_var (ncid, es_heat_perp_id, heat_perp, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, es_heat_perp_id)
       status = nf90_put_var (ncid, es_heat_by_k_id, qheat, start=start4, count=count4)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, es_heat_by_k_id)
    end if

    if (fapar > zero) then
       status = nf90_put_var (ncid, apar_heat_flux_id, mheat_fluxes, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, apar_heat_flux_id)
       status = nf90_put_var (ncid, apar_heat_par_id, mheat_par, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, apar_heat_par_id)
       status = nf90_put_var (ncid, apar_heat_perp_id, mheat_perp, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, apar_heat_perp_id)
       status = nf90_put_var (ncid, apar_heat_by_k_id, qmheat, start=start4, count=count4)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, apar_heat_by_k_id)
       status = nf90_put_var (ncid, apar_heat_by_x_id, qmheat, start=start3, count=count3)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, apar_heat_by_x_id)
    end if

    if (fbpar > zero) then
       status = nf90_put_var (ncid, bpar_heat_flux_id, bheat_fluxes, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, bpar_heat_flux_id)
       status = nf90_put_var (ncid, bpar_heat_par_id, bheat_par, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, bpar_heat_par_id)
       status = nf90_put_var (ncid, bpar_heat_perp_id, bheat_perp, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, bpar_heat_perp_id)
       status = nf90_put_var (ncid, bpar_heat_by_k_id, qbheat, start=start4, count=count4)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, bpar_heat_by_k_id)
    end if
    
    status = nf90_put_var (ncid, hflux_tot_id, hflux_tot, start=(/ nout /))
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, hflux_tot_id)
# endif
  end subroutine nc_qflux

  subroutine nc_pflux (nout, pflux, pmflux, pbflux, &
       part_fluxes, mpart_fluxes, bpart_fluxes, zflux_tot)

    use species, only: nspec
    use kt_grids, only: naky, ntheta0
    use run_parameters, only: fphi, fapar, fbpar
!    use nf90_mod, only: nf90_put_vara, nf90_put_var1
# ifdef NETCDF
    use netcdf, only: nf90_put_var
# endif
    integer, intent (in) :: nout
    real, dimension (:,:,:), intent (in) :: pflux, pmflux, pbflux
    real, dimension(:), intent (in) :: part_fluxes, mpart_fluxes, bpart_fluxes
    real, intent (in) :: zflux_tot
# ifdef NETCDF
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
       status = nf90_put_var (ncid, es_part_flux_id, part_fluxes, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, es_part_flux_id)
       status = nf90_put_var (ncid, es_part_by_k_id, pflux, start=start4, count=count4)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, es_part_by_k_id)
    end if

    if (fapar > zero) then
       status = nf90_put_var (ncid, apar_part_flux_id, mpart_fluxes, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, apar_part_flux_id)
       status = nf90_put_var (ncid, apar_part_by_k_id, pmflux, start=start4, count=count4)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, apar_part_by_k_id)
    end if

    if (fbpar > zero) then
       status = nf90_put_var (ncid, bpar_part_flux_id, bpart_fluxes, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, bpar_part_flux_id)
       status = nf90_put_var (ncid, bpar_part_by_k_id, pbflux, start=start4, count=count4)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, bpar_part_by_k_id)
    end if
    
    status = nf90_put_var (ncid, zflux_tot_id, zflux_tot, start=(/ nout /))
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, zflux_tot_id)
# endif
  end subroutine nc_pflux

  subroutine nc_vflux (nout, vflux, vmflux, vbflux, &
       mom_fluxes, mmom_fluxes, bmom_fluxes, vflux_tot)

    use species, only: nspec
    use kt_grids, only: naky, ntheta0
    use run_parameters, only: fphi, fapar, fbpar
!    use nf90_mod, only: nf90_put_vara, nf90_put_var1
# ifdef NETCDF
    use netcdf, only: nf90_put_var
# endif
    integer, intent (in) :: nout
    real, dimension (:,:,:), intent (in) :: vflux, vmflux, vbflux
    real, dimension(:), intent (in) :: mom_fluxes, mmom_fluxes, bmom_fluxes
    real, intent (in) :: vflux_tot
# ifdef NETCDF
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
       status = nf90_put_var (ncid, es_mom_flux_id, mom_fluxes, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, es_mom_flux_id)
       status = nf90_put_var (ncid, es_mom_by_k_id, vflux, start=start4, count=count4)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, es_mom_by_k_id)
    end if

    if (fapar > zero) then
       status = nf90_put_var (ncid, apar_mom_flux_id, mmom_fluxes, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, apar_mom_flux_id)
       status = nf90_put_var (ncid, apar_mom_by_k_id, vmflux, start=start4, count=count4)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, apar_mom_by_k_id)
    end if

    if (fbpar > zero) then
       status = nf90_put_var (ncid, bpar_mom_flux_id, bmom_fluxes, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, bpar_mom_flux_id)
       status = nf90_put_var (ncid, bpar_mom_by_k_id, vbflux, start=start4, count=count4)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, bpar_mom_by_k_id)
    end if
    
    status = nf90_put_var (ncid, vflux_tot_id, vflux_tot, start=(/ nout /))
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, vflux_tot_id)
# endif
  end subroutine nc_vflux

  subroutine nc_loop (nout, time, fluxfac, &
       phi0,   phi2,   phi2_by_mode, &! phiavg, &
       apar0,  apar2,  apar2_by_mode, &
       bpar0, bpar2, bpar2_by_mode, &
       h, hk, omega, omegaavg, woutunits, phitot, write_omega, write_hrate)

    use gs2_heating, only: heating_diagnostics, hk_repack
    use run_parameters, only: fphi, fapar, fbpar
    use kt_grids, only: naky, ntheta0
    use theta_grid, only: ntgrid
    use species, only: nspec
    use convert, only: c2r
    use fields_arrays, only: phi, apar, bpar

!    use nf90_mod, only: nf90_put_var1, nf90_put_vara
# ifdef NETCDF
    use netcdf, only: nf90_put_var, nf90_sync
# endif

    integer, intent (in) :: nout
    real, intent (in) :: time, phi2, apar2, bpar2
    real, dimension (:), intent (in) :: fluxfac, woutunits
    complex, dimension(:,:), intent (in) :: phi0, apar0, bpar0, omega, omegaavg !, phiavg
    real, dimension(:,:), intent (in) :: phi2_by_mode, apar2_by_mode, bpar2_by_mode, phitot
    type(heating_diagnostics), intent (in) :: h
    type(heating_diagnostics), dimension(:,:), intent (in) :: hk
    logical :: write_omega, write_hrate
# ifdef NETCDF
    real, dimension (ntheta0) :: field2_by_kx
    real, dimension (naky) :: field2_by_ky
    real, dimension (2, ntheta0, naky) :: ri2
    real, dimension (ntheta0, naky, nspec) :: tmps
    complex, dimension (ntheta0, naky) :: tmp
    integer, dimension (4) :: start0, count0, start4, count4
    integer, dimension (3) :: start, count, starth, counth
    integer, dimension (2) :: startx, countx, starty, county, starts, counts
    integer :: status, it, ik


    !Added by EGH ; for writing phi_t, the whole potential vs time
    real, dimension (2, 2*ntgrid+1, ntheta0, naky) :: ri3
    real, dimension (2, 2*ntgrid+1, ntheta0, naky, 1) :: ri_phi_t
    integer, dimension (5) :: start5, count5

    start5(1) = 1
    start5(2) = 1
    start5(3) = 1
    start5(4) = 1
    start5(5) = nout
  
    count5(1) = 2
    count5(2) = 2*ntgrid+1
    count5(3) = ntheta0
    count5(4) = naky
    count5(5) = 1
  
    !End EGH



    status = nf90_put_var (ncid, time_id, time, start=(/ nout /))

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
    counth(2) = 7
    counth(3) = 1

    if (fphi > zero) then


	!<wkdoc> Write fields at the current timestep, if [[write_phi_over_time]], [[write_apar_over_time]], [[write_bpar_over_time]] are set in the input file</wkdoc>
	if(write_phi_t) then
          call c2r (phi, ri3)
          !ri_phi_t(:,:,:,:,1) = ri3(:,:,:,:)
	  status = nf90_put_var (ncid, phi_t_id, ri3, start=start5, count=count5)
          if (status /= NF90_NOERR) call netcdf_error (status, ncid, phi_id)
	end if

	if(write_apar_t) then
          call c2r (apar, ri3)
          !ri_apar_t(:,:,:,:,1) = ri3(:,:,:,:)
	  status = nf90_put_var (ncid, apar_t_id, ri3, start=start5, count=count5)
          if (status /= NF90_NOERR) call netcdf_error (status, ncid, apar_id)
	end if

	if(write_bpar_t) then
          call c2r (bpar, ri3)
          !ri_bpar_t(:,:,:,:,1) = ri3(:,:,:,:)
	  status = nf90_put_var (ncid, bpar_t_id, ri3, start=start5, count=count5)
          if (status /= NF90_NOERR) call netcdf_error (status, ncid, bpar_id)
	end if


       if (ntheta0 > 1) then
          do it = 1, ntheta0
             field2_by_kx(it) = sum(phi2_by_mode(it,:)*fluxfac(:))
          end do
          status = nf90_put_var (ncid, phi2_by_kx_id, field2_by_kx, start=startx, count=countx)
          if (status /= NF90_NOERR) call netcdf_error (status, ncid, phi2_by_kx_id)
       end if

       if (naky > 1) then
          do ik = 1, naky
             field2_by_ky(ik) = sum(phi2_by_mode(:,ik)*fluxfac(ik))
          end do
          status = nf90_put_var (ncid, phi2_by_ky_id, field2_by_ky, start=starty, count=county)
          if (status /= NF90_NOERR) call netcdf_error (status, ncid, phi2_by_ky_id)
       end if

       call c2r (phi0, ri2) 
       status = nf90_put_var (ncid, phi0_id, ri2, start=start0, count=count0)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, phi0_id)
!       call c2r (phiavg, ri2) 
!       status = nf90_put_var (ncid, phiavg_id, ri2, start=start0, count=count0)
!       if (status /= NF90_NOERR) call netcdf_error (status, ncid, phiavg_id)
       status = nf90_put_var (ncid, phi2_by_mode_id, phi2_by_mode, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, phi2_by_mode_id)
       status = nf90_put_var (ncid, phi2_id, phi2, start=(/nout/))
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, phi2_id)
    end if

    if (fapar > zero) then

       if (ntheta0 > 1) then
          do it = 1, ntheta0
             field2_by_kx(it) = sum(apar2_by_mode(it,:)*fluxfac(:))
          end do
          status = nf90_put_var (ncid, apar2_by_kx_id, field2_by_kx, start=startx, count=countx)
          if (status /= NF90_NOERR) call netcdf_error (status, ncid, apar2_by_kx_id)
       end if

       if (naky > 1) then
          do ik = 1, naky
             field2_by_ky(ik) = sum(apar2_by_mode(:,ik)*fluxfac(ik))
          end do
          status = nf90_put_var (ncid, apar2_by_ky_id, field2_by_ky, start=starty, count=county)
          if (status /= NF90_NOERR) call netcdf_error (status, ncid, apar2_by_ky_id)
       end if

       call c2r (apar0, ri2) 
       status = nf90_put_var (ncid, apar0_id, ri2, start=start0, count=count0)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, apar0_id)
       status = nf90_put_var (ncid, apar2_by_mode_id, apar2_by_mode, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, apar2_by_mode_id)
       status = nf90_put_var (ncid, apar2_id, apar2, start=(/nout/))
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, apar2_id)
    end if

    if (fbpar > zero) then

       if (ntheta0 > 1) then
          do it = 1, ntheta0
             field2_by_kx(it) = sum(bpar2_by_mode(it,:)*fluxfac(:))
          end do
          status = nf90_put_var (ncid, bpar2_by_kx_id, field2_by_kx, start=startx, count=countx)
          if (status /= NF90_NOERR) call netcdf_error (status, ncid, bpar2_by_kx_id)
       end if

       if (naky > 1) then
          do ik = 1, naky
             field2_by_ky(ik) = sum(bpar2_by_mode(:,ik)*fluxfac(ik))
          end do
          status = nf90_put_var (ncid, bpar2_by_ky_id, field2_by_ky, start=starty, count=county)
          if (status /= NF90_NOERR) call netcdf_error (status, ncid, bpar2_by_ky_id)
       end if

       call c2r (bpar0, ri2) 
       status = nf90_put_var (ncid, bpar0_id, ri2, start=start0, count=count0)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, bpar0_id)
       status = nf90_put_var (ncid, bpar2_by_mode_id, bpar2_by_mode, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, bpar2_by_mode_id)
       status = nf90_put_var (ncid, bpar2_id, bpar2, start=(/nout/))
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, bpar2_id)
    end if

    if (write_hrate) then
       status = nf90_put_var (ncid, h_energy_id, h%energy, start=(/nout/))
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, h_energy_id)
       status = nf90_put_var (ncid, h_energy_dot_id, h%energy_dot, start=(/nout/))
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, h_energy_dot_id)
       status = nf90_put_var (ncid, h_antenna_id, h%antenna, start=(/nout/))
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, h_antenna_id)
       status = nf90_put_var (ncid, h_eapar_id, h%eapar, start=(/nout/))
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, h_eapar_id)
       status = nf90_put_var (ncid, h_ebpar_id, h%ebpar, start=(/nout/))
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, h_ebpar_id)

       status = nf90_put_var (ncid, h_delfs2_id, h%delfs2, start=starts, count=counts)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, h_delfs2_id)
       status = nf90_put_var (ncid, h_hs2_id, h%hs2, start=starts, count=counts)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, h_hs2_id)
       status = nf90_put_var (ncid, h_phis2_id, h%phis2, start=starts, count=counts)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, h_phis2_id)
       status = nf90_put_var (ncid, h_hypervisc_id, h%hypervisc, start=starts, count=counts)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, h_hypervisc_id)
       status = nf90_put_var (ncid, h_hyperres_id, h%hyperres, start=starts, count=counts)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, h_hyperres_id)
       status = nf90_put_var (ncid, h_collisions_id, h%collisions, start=starts, count=counts)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, h_collisions_id)
       status = nf90_put_var (ncid, h_gradients_id, h%gradients, start=starts, count=counts)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, h_gradients_id)
       status = nf90_put_var (ncid, h_imp_colls_id, h%imp_colls, start=starts, count=counts)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, h_imp_colls_id)
       status = nf90_put_var (ncid, h_hypercoll_id, h%hypercoll, start=starts, count=counts)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, h_hypercoll_id)
       status = nf90_put_var (ncid, h_heating_id, h%heating, start=starts, count=counts)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, h_heating_id)

       status = nf90_put_var (ncid, hk_energy_id, hk%energy, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, hk_energy_id)
       status = nf90_put_var (ncid, hk_energy_dot_id, hk%energy_dot, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, hk_energy_dot_id)
       status = nf90_put_var (ncid, hk_antenna_id, hk%antenna, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, hk_antenna_id)
       status = nf90_put_var (ncid, hk_eapar_id, hk%eapar, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, hk_eapar_id)
       status = nf90_put_var (ncid, hk_ebpar_id, hk%ebpar, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, hk_ebpar_id)

       call hk_repack (hk, 1, tmps)
       status = nf90_put_var (ncid, hk_hypervisc_id, tmps, start=start4, count=count4)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, hk_hypervisc_id)

       call hk_repack (hk, 2, tmps)
       status = nf90_put_var (ncid, hk_hyperres_id, tmps, start=start4, count=count4)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, hk_hyperres_id)

       call hk_repack (hk, 3, tmps)
       status = nf90_put_var (ncid, hk_collisions_id, tmps, start=start4, count=count4)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, hk_collisions_id)

       call hk_repack (hk, 4, tmps)
       status = nf90_put_var (ncid, hk_gradients_id, tmps, start=start4, count=count4)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, hk_gradients_id)

       call hk_repack (hk, 5, tmps)
       status = nf90_put_var (ncid, hk_heating_id, tmps, start=start4, count=count4)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, hk_heating_id)

       call hk_repack (hk, 6, tmps)
       status = nf90_put_var (ncid, hk_hypercoll_id, tmps, start=start4, count=count4)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, hk_hypercoll_id)

       call hk_repack (hk, 7, tmps)
       status = nf90_put_var (ncid, hk_delfs2_id, tmps, start=start4, count=count4)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, hk_delfs2_id)

       call hk_repack (hk, 8, tmps)
       status = nf90_put_var (ncid, hk_hs2_id, tmps, start=start4, count=count4)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, hk_hs2_id)

       call hk_repack (hk, 9, tmps)
       status = nf90_put_var (ncid, hk_phis2_id, tmps, start=start4, count=count4)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, hk_phis2_id)

       call hk_repack (hk, 10, tmps)
       status = nf90_put_var (ncid, hk_imp_colls_id, tmps, start=start4, count=count4)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, hk_imp_colls_id)
    end if

    if (write_omega) then
       do it = 1, ntheta0
          tmp(it, :) = omega(it, :) * woutunits
       end do
       
       call c2r (tmp, ri2)
       status = nf90_put_var (ncid, omega_id, ri2, start=start0, count=count0)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, omega_id)
       
       do it = 1, ntheta0
          tmp(it, :) = omegaavg(it, :) * woutunits
       end do

       call c2r (tmp, ri2)
       status = nf90_put_var (ncid, omegaavg_id, ri2, start=start0, count=count0)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid, omegaavg_id)
    end if

!    status = nf90_put_var (ncid, phtot_id, phitot, start=start, count=count)
!    if (status /= NF90_NOERR) call netcdf_error (status, ncid, phtot_id)
    
    if (mod(nout, 10) == 0) then
       status = nf90_sync (ncid)
       if (status /= NF90_NOERR) call netcdf_error (status)
    end if
# endif
  end subroutine nc_loop

  ! RN> output full not guiding center moments 
  subroutine nc_loop_fullmom (nout, time, &
       & ntot0, density0, upar0, tpar0, tperp0)
    use kt_grids, only: naky, nakx => ntheta0
    use species, only: nspec
    use convert, only: c2r
# ifdef NETCDF
    use netcdf, only: nf90_put_var, nf90_sync
    use netcdf_utils, only: netcdf_error
# endif
    implicit none

    integer, intent (in) :: nout
    real, intent (in) :: time
    complex, intent(in) :: ntot0(:,:,:), density0(:,:,:), upar0(:,:,:)
    complex, intent(in) :: tpar0(:,:,:), tperp0(:,:,:)
# ifdef NETCDF
    real, dimension (2, nakx, naky, nspec) :: ri3
    integer, dimension (5) :: startmom, countmom
    integer :: status, it, ik, is

    status = nf90_put_var (ncid, time_id, time, start=(/nout/))
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, time_id)

    startmom(1) = 1;    countmom(1)=2
    startmom(2) = 1;    countmom(2)=nakx
    startmom(3) = 1;    countmom(3)=naky
    startmom(4) = 1;    countmom(4)=nspec
    startmom(5) = nout; countmom(5)=1

    call c2r (ntot0, ri3)
    status = nf90_put_var (ncid, ntot0_id, ri3, start=startmom, count=countmom)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, ntot0_id)
    call c2r (density0, ri3)
    status = nf90_put_var (ncid, density0_id, ri3, start=startmom, count=countmom)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, density0_id)
    call c2r (upar0, ri3)
    status = nf90_put_var (ncid, upar0_id, ri3, start=startmom, count=countmom)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, upar0_id)
    call c2r (tpar0, ri3)
    status = nf90_put_var (ncid, tpar0_id, ri3, start=startmom, count=countmom)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, tpar0_id)
    call c2r (tperp0, ri3)
    status = nf90_put_var (ncid, tperp0_id, ri3, start=startmom, count=countmom)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, tperp0_id)

    if (mod(nout, 10) == 0) then
       status = nf90_sync (ncid)
       if (status /= NF90_NOERR) call netcdf_error (status)
    end if
# endif
  end subroutine nc_loop_fullmom
  ! <RN

  ! added by EAB on 03/05/04
  subroutine nc_loop_movie (nout_movie, time, yxphi, yxapar, yxbpar)
    use gs2_layouts, only: yxf_lo
    use theta_grid, only: ntgrid
    use run_parameters, only: fphi, fapar, fbpar
!    use nf90_mod, only: nf90_put_var1, nf90_put_vara
# ifdef NETCDF
    use netcdf, only: nf90_put_var, nf90_sync
# endif
    integer, intent (in) :: nout_movie
    real, intent (in) :: time
    real, dimension (:,:,:), intent (in):: yxphi, yxapar, yxbpar 
# ifdef NETCDF
    integer :: status
    integer, dimension (4) :: start, count

    status = nf90_put_var (ncid_movie, time_movie_id, time, start=(/nout_movie/))
    if (status /= NF90_NOERR) call netcdf_error (status, ncid_movie, time_movie_id)
    start(1) = 1
    start(2) = 1
    start(3) = 1
    start(4) = nout_movie
    count(1) = yxf_lo%nx
    count(2) = yxf_lo%ny
    count(3) = 2*ntgrid+1
    count(4) = 1


    if ( fphi > zero) then
       status = nf90_put_var (ncid_movie, phi_by_xmode_id, yxphi, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid_movie, phi_by_xmode_id)
    end if
!    if ( fapar > zero) then
       status = nf90_put_var (ncid_movie, apar_by_xmode_id, yxapar, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid_movie, apar_by_xmode_id)
!    end if

    if ( fbpar > zero) then
       status = nf90_put_var (ncid_movie, bpar_by_xmode_id, yxbpar, start=start, count=count)
       if (status /= NF90_NOERR) call netcdf_error (status, ncid_movie, bpar_by_xmode_id)
    end if

    if (mod(nout_movie, 10) == 0) then
       status = nf90_sync (ncid_movie)
       if (status /= NF90_NOERR) call netcdf_error (status)
    end if
# endif
  end subroutine nc_loop_movie

  subroutine nc_species

    use run_parameters, only: beta
    use species, only: spec
!    use nf90_mod, only: nf90_put_var
# ifdef NETCDF
    use netcdf, only: nf90_put_var

!    real :: betatot
    integer :: status

    status = nf90_put_var (ncid, charge_id, spec%z)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, charge_id)
    status = nf90_put_var (ncid, mass_id, spec%mass)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, mass_id)
    status = nf90_put_var (ncid, dens_id, spec%dens)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, dens_id)
    status = nf90_put_var (ncid, temp_id, spec%temp)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, temp_id)
    status = nf90_put_var (ncid, tprim_id, spec%tprim)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, tprim_id)
    status = nf90_put_var (ncid, fprim_id, spec%fprim)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, fprim_id)
    status = nf90_put_var (ncid, uprim_id, spec%uprim)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, uprim_id)
    status = nf90_put_var (ncid, uprim2_id, spec%uprim2)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, uprim2_id)
    status = nf90_put_var (ncid, vnewk_id, spec%vnewk)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, vnewk_id)
    status = nf90_put_var (ncid, spec_type_id, spec%type)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, spec_type_id)

!    betatot = beta * ...
!    status = nf90_put_var (ncid, betatot_id, betatot)
!    if (status /= NF90_NOERR) call netcdf_error (status, ncid_movie, betatot_id)
# endif
  end subroutine nc_species

  subroutine nc_geo

    use theta_grid, only: bmag, gradpar, gbdrift, gbdrift0, &
         cvdrift, cvdrift0, gds2, gds21, gds22, grho, jacob, &
         shat, drhodpsi, eps
!    use nf90_mod, only: nf90_put_var
# ifdef NETCDF
    use netcdf, only: nf90_put_var

    integer :: status

    status = nf90_put_var (ncid, bmag_id, bmag)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, bmag_id)
    status = nf90_put_var (ncid, gradpar_id, gradpar)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, gradpar_id)
    status = nf90_put_var (ncid, gbdrift_id, gbdrift)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, gbdrift_id)
    status = nf90_put_var (ncid, gbdrift0_id, gbdrift0)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, gbdrift0_id)
    status = nf90_put_var (ncid, cvdrift_id, cvdrift)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, cvdrift_id)
    status = nf90_put_var (ncid, cvdrift0_id, cvdrift0)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, cvdrift0_id)
    status = nf90_put_var (ncid, gds2_id, gds2)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, gds2_id)
    status = nf90_put_var (ncid, gds21_id, gds21)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, gds21_id)
    status = nf90_put_var (ncid, gds22_id, gds22)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, gds22_id)
    status = nf90_put_var (ncid, grho_id, grho)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, grho_id)
    status = nf90_put_var (ncid, jacob_id, jacob)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, jacob_id)

    status = nf90_put_var (ncid, shat_id, shat)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, shat_id)
    status = nf90_put_var (ncid, eps_id, eps)
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, eps_id)
!    status = nf90_put_var (ncid, q_id, eps)  ! find the q variable: qsf, q_val, or qinp?
!    if (status /= NF90_NOERR) call netcdf_error (status, ncid, q_id)
!    status = nf90_put_var (ncid, q_id, eps)  ! write both beta_prime variables with comment
!    if (status /= NF90_NOERR) call netcdf_error (status, ncid, q_id)
    status = nf90_put_var (ncid, drhodpsi_id, drhodpsi)   
    if (status /= NF90_NOERR) call netcdf_error (status, ncid, drhodpsi_id)
# endif
  end subroutine nc_geo

end module gs2_io