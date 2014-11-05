module diagnostics_moments
   use simpledataio
   use simpledataio_write
  use diagnostics_config
   use theta_grid, only: ntgrid
   use kt_grids, only: naky, ntheta0
   use species, only: nspec
  public :: write_moments
contains
  subroutine write_moments(gnostics)
    use gs2_io, only: nc_write_moments
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use species, only: nspec
    use dist_fn, only: getmoms
    use mp, only: proc0
    use fields_arrays, only: phinew, bparnew
    !use diagnostics_write_fields, only: write_standard_field_properties
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

  subroutine write_standard_moment_properties(gnostics, moment_name, moment_description, &
    moment_units, moment_value, distributed)
    use diagnostics_create_and_write, only: create_and_write_variable
    use volume_averages
    type(diagnostics_type), intent(in) :: gnostics
    character(*), intent(in) :: moment_name, moment_description, moment_units
    complex, dimension(:,:,:,:), intent(in) :: moment_value
    logical, intent(in) :: distributed
    real, dimension(ntheta0, naky, nspec) :: moment2_by_mode
    real, dimension(naky,nspec) :: moment2_by_ky
    real, dimension(ntheta0,nspec) :: moment2_by_kx

    call average_theta(moment_value, moment_value, moment2_by_mode, distributed)
    !call create_and_write_moment(gnostics%sfile, moment_name, moment_description, moment_units, moment_value)
    call create_and_write_moment_by_mode(gnostics, moment_name, moment_description, moment_units, &
      moment_value, moment2_by_mode, distributed)

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

  end subroutine write_standard_moment_properties

  subroutine create_and_write_moment_by_mode(gnostics, moment_name, moment_description, moment_units, &
      val, moment2_by_mode, distributed)
   use fields_parallelization, only: field_k_local
   type(diagnostics_type), intent(in) :: gnostics
   character(*), intent(in) :: moment_name
   character(*), intent(in) :: moment_description
   character(*), intent(in) :: moment_units
   character(len=len(moment_name)+9) :: moment2_by_mode_name
   character(len=len(moment_name)+2) :: moment_t_name
   complex, intent(in), dimension(-ntgrid:ntgrid, ntheta0, naky, nspec) :: val
   real, dimension(ntheta0, naky, nspec), intent(in) :: moment2_by_mode
   logical, intent(in) :: distributed
   complex, dimension(-ntgrid:ntgrid, 1, 1) :: dummyc
   real, dimension(1,1) :: dummyr
   logical :: write_moment_by_time
   integer :: it,ik

   moment2_by_mode_name =  moment_name//"2_by_mode" 
   moment_t_name = moment_name//"_t"


   write_moment_by_time = .false.
   if (moment_name .eq. 'phi'  .and. gnostics%write_phi_over_time ) write_moment_by_time = .true.
   if (moment_name .eq. 'apar' .and. gnostics%write_apar_over_time) write_moment_by_time = .true.
   if (moment_name .eq. 'bpar' .and. gnostics%write_bpar_over_time) write_moment_by_time = .true.
   
   if (gnostics%create) then 
     call create_variable(gnostics%sfile, gnostics%rtype, moment_name, "rzXYs", moment_description, moment_units)
     call create_variable(gnostics%sfile, gnostics%rtype, moment2_by_mode_name, "XYst", &
       moment_description//" squared and averaged over theta, as a function of kx and ky" , "("//moment_units//")^2")
     if (write_moment_by_time) then 
       call create_variable(gnostics%sfile, gnostics%rtype, moment_t_name, "rzXYst", &
         moment_description//" as a function of time", moment_units)
     end if


   end if
   
   if (gnostics%create .or. .not. gnostics%wryte) return
   
   if (.not. distributed) then

     call write_variable(gnostics%sfile, moment_name, val)
     call write_variable(gnostics%sfile, moment2_by_mode_name, moment2_by_mode)
     if (write_moment_by_time) call write_variable(gnostics%sfile, moment_t_name, val)

   else

     call set_count(gnostics%sfile, moment_name, "X", 1)
     call set_count(gnostics%sfile, moment_name, "Y", 1)
     call set_count(gnostics%sfile, moment2_by_mode_name, "X", 1)
     call set_count(gnostics%sfile, moment2_by_mode_name, "Y", 1)
     if (write_moment_by_time) then
       call set_count(gnostics%sfile, moment_t_name, "X", 1)
       call set_count(gnostics%sfile, moment_t_name, "Y", 1)
     end if

     ! For some reason every process has to make at least
     ! one write to a variable with an infinite dimension.
     ! Here we make some dummy writes to satisfy that
     if (write_moment_by_time) call write_variable_with_offset(gnostics%sfile, moment_t_name, dummyc)
     call write_variable_with_offset(gnostics%sfile, moment2_by_mode_name, dummyr)

     do ik = 1,naky
       do it = 1,ntheta0
         if (field_k_local(it,ik)) then

           call set_start(gnostics%sfile, moment_name, "X", it)
           call set_start(gnostics%sfile, moment_name, "Y", ik)
           call write_variable_with_offset(gnostics%sfile, moment_name, val)

           call set_start(gnostics%sfile, moment2_by_mode_name, "X", it)
           call set_start(gnostics%sfile, moment2_by_mode_name, "Y", ik)
           call write_variable_with_offset(gnostics%sfile, moment2_by_mode_name, moment2_by_mode)

           if (write_moment_by_time) then
             call set_start(gnostics%sfile, moment_t_name, "X", it)
             call set_start(gnostics%sfile, moment_t_name, "Y", ik)
             call write_variable_with_offset(gnostics%sfile, moment_t_name, val)
           end if
         end if
       end do
     end do
   end if
  end subroutine create_and_write_moment_by_mode

end module diagnostics_moments
