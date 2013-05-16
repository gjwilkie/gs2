
program test_analytical_falpha
  use analytical_falpha
  real :: energy, rslt

  type(analytical_falpha_parameters_type) :: parameters
  write (*,*) '*************************'
  write (*,*) 'Testing analytical_falpha'
  write (*,*) '*************************'

  call announce_test('is_converged')
  call process_test(unit_test_is_converged(), 'is_converged')
  call announce_test('chandrasekhar')
  call process_test(analytical_falpha_unit_test_chandrasekhar(), 'chandrasekhar')

  parameters%alpha_ion_collision_rate = 0.01
  parameters%alpha_electron_collision_rate = 0.006
  parameters%ion_vth = 1.0
  parameters%electron_vth = 42.0
  parameters%alpha_vth = (3.6e2 / 4.0)**0.5

  call announce_test('nu_parallel')
  energy = 0.8
  rslt = 6.9741586014337e-5
  call process_test(analytical_falpha_unit_test_nu_parallel(parameters, energy, &
    rslt), 'nu_parallel')




contains
  subroutine announce_test(test_name)
    character(*), intent(in) :: test_name
    write (*,*) '--> Testing ', test_name
  end subroutine announce_test
  subroutine process_test(rslt, test_name)
    logical, intent (in) :: rslt
    character(*), intent(in) :: test_name
    if (.not. rslt) then 
      write(*,*) '--> ', test_name, ' failed'
      stop 1
    end if

    write (*,*) '--> ', test_name, ' passed'
  end subroutine process_test


end program test_analytical_falpha
