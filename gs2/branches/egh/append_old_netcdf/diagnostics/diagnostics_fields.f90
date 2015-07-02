!> This module contains functions for writing quantities
!! that have the dimensions of the fields, including
!! volume averaged quantities, etc.
module diagnostics_fields
  implicit none
  
  private
  
  public :: write_fields, dump_fields_periodically, write_movie

contains

  subroutine write_fields(gnostics)
    use fields_arrays, only: phinew, aparnew, bparnew
    use run_parameters, only: fphi, fapar, fbpar
    use diagnostics_config, only: diagnostics_type
    implicit none
    type(diagnostics_type), intent(inout) :: gnostics

    if (gnostics%replay) return 

    if (fphi >epsilon(0.0)) call write_standard_field_properties(gnostics, &
      'phi',  'The electrostatic potential', 'T_r/e', phinew)
    if (fapar>epsilon(0.0)) call write_standard_field_properties(gnostics, &
      'apar', 'The parallel magnetic potential', '...', aparnew)
    if (fbpar>epsilon(0.0)) call write_standard_field_properties(gnostics, &
      'bpar', 'The parallel magnetic potential', '...', bparnew)
  end subroutine write_fields

  subroutine write_standard_field_properties(gnostics, field_name, field_description, &
       field_units, field_value)
    use diagnostics_create_and_write, only: create_and_write_variable
    use diagnostics_create_and_write, only: create_and_write_variable_noread
    use diagnostics_create_and_write, only: create_and_write_distributed_fieldlike_variable
    use diagnostics_dimensions, only: dim_string
    use diagnostics_config, only: diagnostics_type
    use volume_averages, only: average_theta, average_kx, average_ky
    use kt_grids, only: ntheta0, naky
    use theta_grid, only: ntgrid
    use mp,only: broadcast
    implicit none
    type(diagnostics_type), intent(inout) :: gnostics
    character(*), intent(in) :: field_name, field_description, field_units
    complex, dimension(-ntgrid:,:,:), intent(inout) :: field_value
    logical :: distributed
    !Should these be allocatable to avoid them hanging around?
    complex, dimension(ntheta0, naky) :: field_igomega_by_mode
    real, dimension(ntheta0, naky) :: field2_by_mode
    real, dimension(naky) :: field2_by_ky
    real, dimension(ntheta0) :: field2_by_kx
    logical :: write_field_by_time    

    distributed = gnostics%distributed

    call average_theta(field_value, field_value, field2_by_mode, distributed)

    call average_kx(field2_by_mode, field2_by_ky, distributed)
    call create_and_write_variable(gnostics, gnostics%rtype, field_name//"2_by_ky", &
         dim_string([gnostics%dims%ky,gnostics%dims%time]), &
         field_description//" squared and averaged over theta and kx, as a function of time", &
         "("//field_units//")^2", field2_by_ky)

    call average_ky(field2_by_mode, field2_by_kx, distributed)
    call create_and_write_variable(gnostics, gnostics%rtype, field_name//"2_by_kx", &
         dim_string([gnostics%dims%kx,gnostics%dims%time]), &
         field_description//" squared and averaged over theta and ky, as a function of time", &
         "("//field_units//")^2", field2_by_kx)

    call create_and_write_variable_noread(gnostics, gnostics%rtype, field_name//"2", &
         dim_string([gnostics%dims%time]), &
         field_description//" squared and averaged over theta, kx and ky, as a function of time", &
         "("//field_units//")^2", sum(field2_by_kx))

    if (field_name .eq. 'phi') gnostics%current_results%phi2 = sum(field2_by_kx)
    if (field_name .eq. 'apar') gnostics%current_results%apar2 = sum(field2_by_kx)
    if (field_name .eq. 'bpar') gnostics%current_results%bpar2 = sum(field2_by_kx)

    ! Below we deal with the larger arrays which are functions of ky and kx and
    ! may be distributed

    write_field_by_time = .false.

    if (field_name .eq. 'phi'  .and. gnostics%write_phi_over_time ) write_field_by_time = .true.
    if (field_name .eq. 'apar' .and. gnostics%write_apar_over_time) write_field_by_time = .true.
    if (field_name .eq. 'bpar' .and. gnostics%write_bpar_over_time) write_field_by_time = .true.
 
    field_igomega_by_mode(:,:) = field_value(gnostics%igomega, :, :)

    !write (*,*) 'dim string is ', dim_string([gnostics%dims%ri,gnostics%dims%theta,gnostics%dims%kx,gnostics%dims%ky])
    call create_and_write_distributed_fieldlike_variable( &
         gnostics, gnostics%rtype, field_name, &
         dim_string([gnostics%dims%ri,gnostics%dims%theta,gnostics%dims%kx,gnostics%dims%ky]), &
         field_description, field_units, field_value)
    call create_and_write_distributed_fieldlike_variable( &
         gnostics, gnostics%rtype, field_name//"_igomega_by_mode", &
         dim_string([gnostics%dims%ri,gnostics%dims%kx,gnostics%dims%ky,gnostics%dims%time]), &
         field_description//" at ig=igomega, as a function of kx and ky" , field_units, field_igomega_by_mode)
    call create_and_write_distributed_fieldlike_variable( &
         gnostics, gnostics%rtype, field_name//"2_by_mode", &
         dim_string([gnostics%dims%kx,gnostics%dims%ky,gnostics%dims%time]), &
         field_description//" squared and averaged over theta, as a function of kx and ky" , &
         field_units, field2_by_mode)

    if (write_field_by_time) then 
       call create_and_write_distributed_fieldlike_variable( &
            gnostics, gnostics%rtype, field_name//"_t", &
            dim_string([gnostics%dims%ri,gnostics%dims%theta,&
            gnostics%dims%kx,gnostics%dims%ky,gnostics%dims%time]), &
            field_description//": the whole field as a function time " , &
            field_units, field_value)
    endif
  end subroutine write_standard_field_properties

  !> Dump the fields to a time-labelled ascii file:
  !! suggestion... use the netcdf file with 
  !! write_phi_over_time etc on instead :-)
  subroutine dump_fields_periodically(gnostics)
    use file_utils, only: get_unused_unit
    use kt_grids, only: naky, ntheta0, aky, theta0, akx
    use theta_grid, only: theta, ntgrid
    use fields_arrays, only: phinew, aparnew, bparnew
    use mp, only: proc0
    use diagnostics_config, only: diagnostics_type
    implicit none
    type(diagnostics_type), intent(in) :: gnostics
    real :: t
    character(200) :: filename
    integer :: ik, it, ig, unit
    integer :: ntg_out 

    !<DD>Adding guard to only let proc0 do this writing
    if(.not.proc0) return

    ! EGH I'm pretty sure that there is no difference
    ! between ntg_out and ntgrid in the old diagnostics 
    ! module... anyone disagree?
    !<DD>I agree... we should probably remove it. 
    ntg_out = ntgrid
    
    t = gnostics%user_time
    
    !<DD>Should this only occur for write_ascii=T?
    ! EGH Yes...
    call get_unused_unit (unit)
    write (filename, "('dump.fields.t=',e13.6)") t
    open (unit=unit, file=filename, status="unknown")
    do ik = 1, naky
       do it = 1, ntheta0
          do ig = -ntg_out, ntg_out
             write (unit=unit, fmt="(20(1x,e12.5))") &
                  theta(ig), aky(ik), theta0(it,ik), &
                  phinew(ig,it,ik), aparnew(ig,it,ik), &
                  bparnew(ig,it,ik), t, akx(it)
          end do
          write (unit, "()")
       end do
    end do
    close (unit=unit)
  end subroutine dump_fields_periodically

  subroutine write_movie(gnostics)
    use gs2_transforms, only: init_transforms
    use gs2_layouts, only: yxf_lo
    use theta_grid, only: ntgrid
    use run_parameters, only: fphi, fapar, fbpar
    use gs2_transforms, only: transform2
    use gs2_io, only: nc_loop_movie
    use fields_arrays, only: phinew, aparnew, bparnew
    use mp, only: proc0, mp_abort
    use diagnostics_create_and_write, only: create_and_write_variable
    use diagnostics_dimensions, only: dim_string
    use le_grids, only: nlambda, negrid
    use species, only: nspec
    use kt_grids, only: nx, ny, naky, ntheta0
    use file_utils, only: error_unit
    use diagnostics_config, only: diagnostics_type
    implicit none
    type(diagnostics_type), intent(in) :: gnostics
    real :: t
    integer :: nnx, nny
    real, dimension(:,:,:), allocatable :: yxphi, yxapar, yxbpar
    logical :: accelerated ! dummy variable here
    logical, save :: warning = .false.

    call init_transforms (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx, ny, accelerated)

    t = gnostics%user_time
    
    ! EAB 09/17/03 -- modify dump_fields_periodically to print out inverse fft of fields in x,y
    nnx = yxf_lo%nx
    nny = yxf_lo%ny
    
    if (proc0 .and. .not. warning) then
       write (error_unit(),*) "WARNING: write_movie writes arrays of size nx*ny which contain a lot of &
            & redundancy, leading to unnecessarily large output files... consider doing &
            & the Fourier transforms of the fields in post-processing instead"
       warning = .true.
    endif
    
    if (fphi > epsilon(0.0)) then
       allocate (yxphi(nnx,nny,-ntgrid:ntgrid))
       call transform2 (phinew, yxphi, nny, nnx)
       call create_and_write_variable(gnostics, gnostics%rtype, "phi_by_xmode", &
            dim_string([gnostics%dims%xx,gnostics%dims%yy,gnostics%dims%theta,gnostics%dims%time]), &
            "The electric potential in real space", "Tr/e", yxphi)
       deallocate (yxphi)
    end if

    if (fapar > epsilon(0.0)) then
       allocate (yxapar(nnx,nny,-ntgrid:ntgrid))
       call transform2 (aparnew, yxapar, nny, nnx)
       call create_and_write_variable(gnostics, gnostics%rtype, "apar_by_xmode", &
            dim_string([gnostics%dims%xx,gnostics%dims%yy,gnostics%dims%theta,gnostics%dims%time]), &
            "The parallel vector potential in real space", "TBC", yxapar)
       deallocate (yxapar)
    end if

    if (fbpar > epsilon(0.0)) then 
       allocate (yxbpar(nnx,nny,-ntgrid:ntgrid))
       call transform2 (bparnew, yxbpar, nny, nnx)
       call create_and_write_variable(gnostics, gnostics%rtype, "bpar_by_xmode", &
            dim_string([gnostics%dims%xx,gnostics%dims%yy,gnostics%dims%theta,gnostics%dims%time]), &
            "The parallel magnetic field in real space", "TBC", yxbpar)
       deallocate (yxbpar)
    end if
  end subroutine write_movie
end module diagnostics_fields

