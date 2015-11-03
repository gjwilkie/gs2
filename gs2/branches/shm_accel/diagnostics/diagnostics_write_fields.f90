!> This module contains functions for writing quantities
!! that have the dimensions of the fields, including
!! volume averaged quantities, etc.
module diagnostics_write_fields
   use simpledataio
   use simpledataio_write
   use theta_grid, only: ntgrid
   use kt_grids, only: naky, ntheta0
   use diagnostics_config, only: diagnostics_type
   implicit none
   !integer, parameter :: gnostics%rtype = SDATIO_DOUBLE

contains

  subroutine write_fields(gnostics)
    use fields_arrays, only: phinew, aparnew, bparnew
    use run_parameters, only: fphi, fapar, fbpar
    type(diagnostics_type), intent(in) :: gnostics

    
    if (fphi >epsilon(0.0)) call write_standard_field_properties(gnostics, &
      'phi',  'The electrostatic potential', 'T_r/e', phinew, gnostics%distributed)
    if (fapar>epsilon(0.0)) call write_standard_field_properties(gnostics, &
      'apar', 'The parallel magnetic potential', '...', aparnew, gnostics%distributed)
    if (fbpar>epsilon(0.0)) call write_standard_field_properties(gnostics, &
      'bpar', 'The parallel magnetic potential', '...', bparnew, gnostics%distributed)
  end subroutine write_fields

  subroutine write_standard_field_properties(gnostics, field_name, field_description, &
    field_units, field_value, distributed)
    use diagnostics_create_and_write, only: create_and_write_variable
    use volume_averages
    type(diagnostics_type), intent(in) :: gnostics
    character(*), intent(in) :: field_name, field_description, field_units
    complex, dimension(:,:,:), intent(in) :: field_value
    logical, intent(in) :: distributed
    real, dimension(ntheta0, naky) :: field2_by_mode
    real, dimension(naky) :: field2_by_ky
    real, dimension(ntheta0) :: field2_by_kx

    call average_theta(field_value, field_value, field2_by_mode, distributed)
    !call create_and_write_field(gnostics%sfile, field_name, field_description, field_units, field_value)
    call create_and_write_field_by_mode(gnostics, field_name, field_description, field_units, &
      field_value, field2_by_mode, distributed)

    call average_kx(field2_by_mode, field2_by_ky, distributed)
    call create_and_write_variable(gnostics, gnostics%rtype, field_name//"2_by_ky", "Yt", &
      field_description//" squared and averaged over theta and kx, as a function of time", &
      "("//field_units//")^2", field2_by_ky)

    call average_ky(field2_by_mode, field2_by_kx, distributed)
    call create_and_write_variable(gnostics, gnostics%rtype, field_name//"2_by_kx", "Xt", &
      field_description//" squared and averaged over theta and ky, as a function of time", &
      "("//field_units//")^2", field2_by_kx)

    call create_and_write_variable(gnostics, gnostics%rtype, field_name//"2", "t", &
      field_description//" squared and averaged over theta, kx and ky, as a function of time", &
      "("//field_units//")^2", sum(field2_by_kx))

  end subroutine write_standard_field_properties

  subroutine create_and_write_field_by_mode(gnostics, field_name, field_description, field_units, &
      val, field2_by_mode, distributed)
   use fields_parallelization, only: field_k_local
   type(diagnostics_type), intent(in) :: gnostics
   character(*), intent(in) :: field_name
   character(*), intent(in) :: field_description
   character(*), intent(in) :: field_units
   character(len=len(field_name)+9) :: field2_by_mode_name
   character(len=len(field_name)+2) :: field_t_name
   complex, intent(in), dimension(-ntgrid:ntgrid, ntheta0, naky) :: val
   real, dimension(ntheta0, naky), intent(in) :: field2_by_mode
   logical, intent(in) :: distributed
   complex, dimension(-ntgrid:ntgrid, 1, 1) :: dummyc
   real, dimension(1,1) :: dummyr
   logical :: write_field_by_time
   integer :: it,ik

   field2_by_mode_name =  field_name//"2_by_mode" 
   field_t_name = field_name//"_t"


   write_field_by_time = .false.
   if (field_name .eq. 'phi'  .and. gnostics%write_phi_over_time ) write_field_by_time = .true.
   if (field_name .eq. 'apar' .and. gnostics%write_apar_over_time) write_field_by_time = .true.
   if (field_name .eq. 'bpar' .and. gnostics%write_bpar_over_time) write_field_by_time = .true.
   
   if (gnostics%create) then 
     call create_variable(gnostics%sfile, gnostics%rtype, field_name, "rzXY", field_description, field_units)
     call create_variable(gnostics%sfile, gnostics%rtype, field2_by_mode_name, "XYt", &
       field_description//" squared and averaged over theta, as a function of kx and ky" , "("//field_units//")^2")
     if (write_field_by_time) then 
       call create_variable(gnostics%sfile, gnostics%rtype, field_t_name, "rzXYt", &
         field_description//" as a function of time", field_units)
     end if


   end if
   
   if (gnostics%create .or. .not. gnostics%wryte) return
   
   if (.not. distributed) then

     call write_variable(gnostics%sfile, field_name, val)
     call write_variable(gnostics%sfile, field2_by_mode_name, field2_by_mode)
     if (write_field_by_time) call write_variable(gnostics%sfile, field_t_name, val)

   else

     call set_count(gnostics%sfile, field_name, "X", 1)
     call set_count(gnostics%sfile, field_name, "Y", 1)
     call set_count(gnostics%sfile, field2_by_mode_name, "X", 1)
     call set_count(gnostics%sfile, field2_by_mode_name, "Y", 1)
     if (write_field_by_time) then
       call set_count(gnostics%sfile, field_t_name, "X", 1)
       call set_count(gnostics%sfile, field_t_name, "Y", 1)
     end if

     ! For some reason every process has to make at least
     ! one write to a variable with an infinite dimension.
     ! Here we make some dummy writes to satisfy that
     if (write_field_by_time) call write_variable_with_offset(gnostics%sfile, field_t_name, dummyc)
     call write_variable_with_offset(gnostics%sfile, field2_by_mode_name, dummyr)

     do ik = 1,naky
       do it = 1,ntheta0
         if (field_k_local(it,ik)) then

           call set_start(gnostics%sfile, field_name, "X", it)
           call set_start(gnostics%sfile, field_name, "Y", ik)
           call write_variable_with_offset(gnostics%sfile, field_name, val)

           call set_start(gnostics%sfile, field2_by_mode_name, "X", it)
           call set_start(gnostics%sfile, field2_by_mode_name, "Y", ik)
           call write_variable_with_offset(gnostics%sfile, field2_by_mode_name, field2_by_mode)

           if (write_field_by_time) then
             call set_start(gnostics%sfile, field_t_name, "X", it)
             call set_start(gnostics%sfile, field_t_name, "Y", ik)
             call write_variable_with_offset(gnostics%sfile, field_t_name, val)
           end if
         end if
       end do
     end do
   end if
  end subroutine create_and_write_field_by_mode

end module diagnostics_write_fields

