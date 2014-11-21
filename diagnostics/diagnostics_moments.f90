!> A module for writing out moments of the distribution function,
!! for example density, temperature, parallel flow
module diagnostics_moments
  implicit none
  private
  public :: write_moments, write_full_moments_notgc
contains
  subroutine write_moments(gnostics)
    use gs2_io, only: nc_write_moments
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use species, only: nspec
    use dist_fn, only: getmoms
    use fields_arrays, only: phinew, bparnew
    use diagnostics_config, only: diagnostics_type
    implicit none
    type(diagnostics_type), intent(in) :: gnostics
    complex, dimension (-ntgrid:ntgrid,ntheta0,naky,nspec) :: ntot, density, &
         upar, tpar, tperp, qparflux, pperpj1, qpperpj1

    call getmoms (phinew, bparnew, ntot, density, upar, tpar, tperp, qparflux, pperpj1, qpperpj1)

    call write_standard_moment_properties(gnostics, &
         'ntot',  'The total perturbed species density', 'n_r', ntot, gnostics%distributed)
    call write_standard_moment_properties(gnostics, &
         'density',  'The non-adiabatic part of the perturbed species density', 'n_r', density, gnostics%distributed)
    call write_standard_moment_properties(gnostics, &
         'upar',  'The perturbed parallel flow', 'v_ths', upar, gnostics%distributed)
    call write_standard_moment_properties(gnostics, &
         'tpar',  'The perturbed parallel temperature', 'T_r', tpar, gnostics%distributed)
    call write_standard_moment_properties(gnostics, &
         'tperp',  'The perturbed perpendicular temperature', 'T_r', tperp, gnostics%distributed)
    call write_standard_moment_properties(gnostics, &
         'qparflux',  'The parallel heat flux', 'n_r T_r v_thr', qparflux, gnostics%distributed)
    call write_standard_moment_properties(gnostics, &
         'pperpj1',  'A modified perpendicular pressure which give the particle flux from bpar',&
         'n_r T_r / e', pperpj1, gnostics%distributed)
    call write_standard_moment_properties(gnostics, &
         'qpperpj1',  'A modified perpendicular pressure * energy which give the heat flux from bpar',&
         'n_r T_r^2 / e', qpperpj1, gnostics%distributed)
    !call write_standard_moment_properties(gnostics, &
    !'upar',  'The non-adiabatic part of the perturbed species density', 'n_r', ntot, gnostics%distributed)
  end subroutine write_moments

  subroutine write_full_moments_notgc(gnostics)
    use gs2_io, only: nc_write_moments
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use species, only: nspec
    use dist_fn, only: getmoms_notgc
    use fields_arrays, only: phinew, bparnew
    use diagnostics_config, only: diagnostics_type
    implicit none
    type(diagnostics_type), intent(in) :: gnostics
    complex, dimension (-ntgrid:ntgrid,ntheta0,naky,nspec) :: ntot, density, &
         upar, tpar, tperp !, qparflux, pperpj1, qpperpj1

    call getmoms_notgc(phinew, bparnew, density,upar,tpar,tperp,ntot)

    call write_standard_moment_properties(gnostics, &
         'ntot_notgc',  'The total perturbed species density &
         & in non-guiding centre coordinates', 'n_r', ntot, gnostics%distributed)
    call write_standard_moment_properties(gnostics, &
         'density_notgc',  'The non-adiabatic part of the perturbed species density &
         & in non-guiding centre coordinates', 'n_r', density, gnostics%distributed)
    call write_standard_moment_properties(gnostics, &
         'upar_notgc',  'The perturbed parallel flow &
         & in non-guiding centre coordinates', 'v_ths', upar, gnostics%distributed)
    call write_standard_moment_properties(gnostics, &
         'tpar_notgc',  'The perturbed parallel temperature &
         & in non-guiding centre coordinates', 'T_r', tpar, gnostics%distributed)
    call write_standard_moment_properties(gnostics, &
         'tperp_notgc',  'The perturbed perpendicular temperature &
         & in non-guiding centre coordinates', 'T_r', tperp, gnostics%distributed)
  end subroutine write_full_moments_notgc

  subroutine write_standard_moment_properties(gnostics, moment_name, moment_description, &
       moment_units, moment_value, distributed)
    use diagnostics_create_and_write, only: create_and_write_variable
    use diagnostics_create_and_write, only: create_and_write_distributed_fieldlike_variable
    use volume_averages, only: average_theta, average_kx, average_ky
    use fields_parallelization, only: field_k_local
    use mp, only: sum_allreduce
    use theta_grid, only: ntgrid
    use kt_grids, only: ntheta0, naky
    use species, only: nspec
    use diagnostics_config, only: diagnostics_type
    type(diagnostics_type), intent(in) :: gnostics
    character(*), intent(in) :: moment_name, moment_description, moment_units
    complex, dimension(-ntgrid:,:,:,:), intent(in) :: moment_value
    logical, intent(in) :: distributed
    real, dimension(ntheta0, naky, nspec) :: moment2_by_mode
    real, dimension(naky,nspec) :: moment2_by_ky
    real, dimension(ntheta0,nspec) :: moment2_by_kx
    complex, dimension(ntheta0, naky, nspec) :: moment_by_mode
    complex, dimension(ntheta0, naky, nspec) :: moment_igomega_by_mode
    complex, dimension(ntheta0, nspec) :: moment_flx_surfavg
    logical :: write_moment_by_time
    integer :: it !, ik

    call average_theta(moment_value, moment_value, moment2_by_mode, distributed)
    !call create_and_write_moment(gnostics%sfile, moment_name, moment_description, moment_units, moment_value)
    !call create_and_write_moment_by_mode(gnostics, moment_name, moment_description, moment_units, &
    !moment_value, moment2_by_mode, distributed)

    call average_kx(moment2_by_mode, moment2_by_ky, distributed)
    call create_and_write_variable(gnostics, gnostics%rtype, moment_name//"2_by_ky", "Yst", &
         moment_description//" squared and averaged over theta and kx", &
         "("//moment_units//")^2", moment2_by_ky)

    call average_ky(moment2_by_mode, moment2_by_kx, distributed)
    call create_and_write_variable(gnostics, gnostics%rtype, moment_name//"2_by_kx", "Xst", &
         moment_description//" squared and averaged over theta and ky", &
         "("//moment_units//")^2", moment2_by_kx)

    call create_and_write_variable(gnostics, gnostics%rtype, moment_name//"2", "st", &
         moment_description//" squared and averaged over theta, kx and ky", &
         "("//moment_units//")^2", sum(moment2_by_kx, 1))

    call average_theta(moment_value, moment_by_mode, gnostics%distributed)
    moment_igomega_by_mode(:,:,:) = moment_value(gnostics%igomega,:,:,:)

    ! moment_by_mode could be distributed so we have to be careful here
    moment_flx_surfavg(:,:) = 0.0
    do it = 1,ntheta0
       if ((.not. gnostics%distributed).or.field_k_local(it, 1))&
            moment_flx_surfavg(it, :) = moment_by_mode(it, 1, :)
    end do
    if (gnostics%distributed) call sum_allreduce(moment_flx_surfavg)
    
    call create_and_write_variable(gnostics, gnostics%rtype, moment_name//"_flxsurf_avg", "rXst", &
         moment_description//" flux surface averaged: &
         & averaged over theta, at ky=0 (actally ik==1),  as a function of kx" , moment_units, moment_flx_surfavg)
    
    call create_and_write_distributed_fieldlike_variable( &
         gnostics, gnostics%rtype, moment_name, "rzXYs", moment_description, moment_units, moment_value)
    call create_and_write_distributed_fieldlike_variable( &
         gnostics, gnostics%rtype, moment_name//"_igomega_by_mode", "rXYst", &
         moment_description//" at ig=igomega, as a function of kx and ky" , moment_units, moment_igomega_by_mode)
    call create_and_write_distributed_fieldlike_variable( &
         gnostics, gnostics%rtype, moment_name//"_by_mode", "rXYst", &
         moment_description//"  averaged over theta, as a function of kx and ky" , &
         moment_units, moment_by_mode)
    call create_and_write_distributed_fieldlike_variable( &
         gnostics, gnostics%rtype, moment_name//"2_by_mode", "XYst", &
         moment_description//" squared and averaged over theta, as a function of kx and ky" , &
         moment_units, moment2_by_mode)

    write_moment_by_time = .false.
    if (moment_name .eq. 'ntot'  .and. gnostics%write_ntot_over_time ) write_moment_by_time = .true.
    if (moment_name .eq. 'density' .and. gnostics%write_density_over_time) write_moment_by_time = .true.
    if (moment_name .eq. 'upar' .and. gnostics%write_upar_over_time) write_moment_by_time = .true.
    if (moment_name .eq. 'tperp' .and. gnostics%write_tperp_over_time) write_moment_by_time = .true.

    if (write_moment_by_time) & 
         call create_and_write_distributed_fieldlike_variable( &
         gnostics, gnostics%rtype, moment_name//"_t", "rzXYst", &
         moment_description//": the whole moment, as a function of time" , &
         moment_units, moment_value)

  end subroutine write_standard_moment_properties

!Is this dead?
  !subroutine create_and_write_moment_by_mode(gnostics, moment_name, moment_description, moment_units, &
  !val, moment2_by_mode, distributed)
  !use fields_parallelization, only: field_k_local
  !use volume_averages, only: average_theta
  !type(diagnostics_type), intent(in) :: gnostics
  !character(*), intent(in) :: moment_name
  !character(*), intent(in) :: moment_description
  !character(*), intent(in) :: moment_units
  !character(len=len(moment_name)+8) :: moment_by_mode_name
  !character(len=len(moment_name)+12) :: moment_flx_surfavg_name
  !character(len=len(moment_name)+16) :: moment_igomega_by_mode_name
  !character(len=len(moment_name)+9) :: moment2_by_mode_name
  !character(len=len(moment_name)+2) :: moment_t_name
  !complex, intent(in), dimension(-ntgrid:ntgrid, ntheta0, naky, nspec) :: val
  !complex, dimension(ntheta0, naky, nspec) :: moment_by_mode
  !complex, dimension(ntheta0, naky, nspec) :: moment_igomega_by_mode
  !complex, dimension(ntheta0, nspec) :: moment_flx_surfavg
  !! It is unnecessary to pass moment2_by_mode in to this routine... 
  !! all it does is save one volume average among thousands. Should
  !! be changed at some point?
  !real, dimension(ntheta0, naky, nspec), intent(in) :: moment2_by_mode
  !logical, intent(in) :: distributed
  !complex, dimension(-ntgrid:ntgrid, 1, 1) :: dummyc
  !real, dimension(1,1) :: dummyr
  !complex, dimension(1,1) :: dummyc1
  !logical :: write_moment_by_time
  !integer :: it,ik

  !moment_flx_surfavg_name =  moment_name//"_flxsurf_avg" 
  !moment_by_mode_name =  moment_name//"_by_mode" 
  !moment_igomega_by_mode_name =  moment_name//"_igomega_by_mode" 
  !moment2_by_mode_name =  moment_name//"2_by_mode" 
  !moment_t_name = moment_name//"_t"


  !write_moment_by_time = .false.
  !if (moment_name .eq. 'ntot'  .and. gnostics%write_ntot_over_time ) write_moment_by_time = .true.
  !if (moment_name .eq. 'density' .and. gnostics%write_density_over_time) write_moment_by_time = .true.
  !if (moment_name .eq. 'upar' .and. gnostics%write_upar_over_time) write_moment_by_time = .true.
  !if (moment_name .eq. 'tperp' .and. gnostics%write_tperp_over_time) write_moment_by_time = .true.

  !if (gnostics%create) then 
  !call create_variable(gnostics%sfile, gnostics%rtype, moment_name, "rzXYs", moment_description, moment_units)
  !call create_variable(gnostics%sfile, gnostics%rtype, moment_flx_surfavg_name, "rXst", &
  !moment_description//"flux surface averaged: &
  !& averaged over theta, at ky=0 (actally ik==1),  as a function of kx" , moment_units)
  !call create_variable(gnostics%sfile, gnostics%rtype, moment_by_mode_name, "rXYst", &
  !moment_description//" averaged over theta, as a function of kx and ky" , moment_units)
  !call create_variable(gnostics%sfile, gnostics%rtype, moment_igomega_by_mode_name, "rXYst", &
  !moment_description//" at ig=igomega, as a function of kx and ky" , moment_units)
  !call create_variable(gnostics%sfile, gnostics%rtype, moment2_by_mode_name, "XYst", &
  !moment_description//" squared and averaged over theta, as a function of kx and ky" , "("//moment_units//")^2")
  !if (write_moment_by_time) then 
  !call create_variable(gnostics%sfile, gnostics%rtype, moment_t_name, "rzXYst", &
  !moment_description//" as a function of time", moment_units)
  !end if


  !end if

  !if (gnostics%create .or. .not. gnostics%wryte) return

  !call average_theta(val, moment_by_mode, gnostics%distributed)
  !moment_igomega_by_mode(:,:,:) = val(gnostics%igomega,:,:,:)
  !moment_flx_surfavg(:,:) = moment_by_mode(:,1,:)

  !if (.not. distributed) then

  !! The easy case where all fields and all moments are on all processors:
  !! proc0 just writes the whole thing. NB gnostics%wryte is only true for
  !! proc0 in this instance.
  !call write_variable(gnostics%sfile, moment_name, val)
  !call write_variable(gnostics%sfile, moment_flx_surfavg_name, moment_flx_surfavg)
  !call write_variable(gnostics%sfile, moment_by_mode_name, moment_by_mode)
  !call write_variable(gnostics%sfile, moment_igomega_by_mode_name, moment_igomega_by_mode)
  !call write_variable(gnostics%sfile, moment2_by_mode_name, moment2_by_mode)
  !if (write_moment_by_time) call write_variable(gnostics%sfile, moment_t_name, val)

  !else

  !call set_count(gnostics%sfile, moment_name, "X", 1)
  !call set_count(gnostics%sfile, moment_name, "Y", 1)
  !call set_count(gnostics%sfile, moment_flx_surfavg_name, "X", 1)
  !call set_count(gnostics%sfile, moment_by_mode_name, "X", 1)
  !call set_count(gnostics%sfile, moment_by_mode_name, "Y", 1)
  !call set_count(gnostics%sfile, moment_igomega_by_mode_name, "X", 1)
  !call set_count(gnostics%sfile, moment_igomega_by_mode_name, "Y", 1)
  !call set_count(gnostics%sfile, moment2_by_mode_name, "X", 1)
  !call set_count(gnostics%sfile, moment2_by_mode_name, "Y", 1)
  !if (write_moment_by_time) then
  !call set_count(gnostics%sfile, moment_t_name, "X", 1)
  !call set_count(gnostics%sfile, moment_t_name, "Y", 1)
  !end if

  !! For some reason every process has to make at least
  !! one write to a variable with an infinite dimension.
  !! Here we make some dummy writes to satisfy that
  !if (write_moment_by_time) call write_variable_with_offset(gnostics%sfile, moment_t_name, dummyc)
  !call write_variable_with_offset(gnostics%sfile, moment_flx_surfavg_name, dummyc1)
  !call write_variable_with_offset(gnostics%sfile, moment_by_mode_name, dummyc1)
  !call write_variable_with_offset(gnostics%sfile, moment_igomega_by_mode_name, dummyc1)
  !call write_variable_with_offset(gnostics%sfile, moment2_by_mode_name, dummyr)

  !do ik = 1,naky
  !do it = 1,ntheta0
  !if (field_k_local(it,ik)) then

  !call set_start(gnostics%sfile, moment_name, "X", it)
  !call set_start(gnostics%sfile, moment_name, "Y", ik)
  !call write_variable_with_offset(gnostics%sfile, moment_name, val)

  !if (ik.eq.1) then
  !call set_start(gnostics%sfile, moment_flx_surfavg_name, "X", it)
  !call write_variable_with_offset(gnostics%sfile, moment_flx_surfavg_name, moment_flx_surfavg)
  !end if

  !call set_start(gnostics%sfile, moment_by_mode_name, "X", it)
  !call set_start(gnostics%sfile, moment_by_mode_name, "Y", ik)
  !call write_variable_with_offset(gnostics%sfile, moment_by_mode_name, moment_by_mode)

  !call set_start(gnostics%sfile, moment_igomega_by_mode_name, "X", it)
  !call set_start(gnostics%sfile, moment_igomega_by_mode_name, "Y", ik)
  !call write_variable_with_offset(gnostics%sfile, moment_igomega_by_mode_name, moment_igomega_by_mode)

  !call set_start(gnostics%sfile, moment2_by_mode_name, "X", it)
  !call set_start(gnostics%sfile, moment2_by_mode_name, "Y", ik)
  !call write_variable_with_offset(gnostics%sfile, moment2_by_mode_name, moment2_by_mode)

  !if (write_moment_by_time) then
  !call set_start(gnostics%sfile, moment_t_name, "X", it)
  !call set_start(gnostics%sfile, moment_t_name, "Y", ik)
  !call write_variable_with_offset(gnostics%sfile, moment_t_name, val)
  !end if
  !end if
  !end do
  !end do
  !end if
  !end subroutine create_and_write_moment_by_mode

end module diagnostics_moments
