!> A module which generates an analytical form
!! for the alpha particle equilibrium distribution
!! function by integrating an expression given in...
!!
!! This is free software released under the GPLv3
!! Written by
!!    Edmund Highcock (edmundhighcock@sourceforge.net)
!!    Ian Abel ()


module analytical_falpha

  !> This subroutine does the work of the module.
  !! It takes in a grid of energy values and the species index
  !! and then, using the species parameter arrays, calculates
  !! f0, generalised_temperature (a function of df0/denergy)
  !! and df0/drho for that species.
  public :: calculate_arrays

  public :: unit_test_is_converged
  public :: analytical_falpha_unit_test_chandrasekhar
  public :: analytical_falpha_unit_test_nu_parallel
  public :: analytical_falpha_unit_test_falpha_integrand
  public :: analytical_falpha_unit_test_falpha

  public :: analytical_falpha_parameters_type
  type analytical_falpha_parameters_type
    integer :: alpha_is
    real :: alpha_mass
    real :: alpha_vth
    real :: alpha_injection_energy
    real :: ion_temp
    real :: ion_vth
    real :: electron_temp
    real :: electron_vth
    real :: alpha_ion_collision_rate
    real :: alpha_electron_collision_rate
    integer :: negrid
  end type analytical_falpha_parameters_type

  private 


contains

  subroutine calculate_arrays(parameters, &
                              egrid, &
                              f0_values, &
                              generalised_temperature, &
                              f0prim)
    type(analytical_falpha_parameters_type), intent(in) :: parameters
    real, dimension(:,:), intent(in) :: egrid
    real, dimension(:,:), intent(out) :: f0_values
    real, dimension(:,:), intent(out) :: generalised_temperature
    real, dimension(:,:), intent(out) :: f0prim

    integer :: ie
    integer :: is
    integer :: resolution
    logical :: converged = .false.
    real, dimension(2) :: f0
    real, dimension(2) :: gentemp
    real, dimension(2) :: f0pr



    is = parameters%alpha_is


    do ie = 1,parameters%negrid
      ! Initialise arrays to test for convergence 
      f0 = -1.0
      gentemp = 0.0
      f0pr = 0.0
      resolution = 64
      do while (.not. converged)
        f0(1)      = f0(2)
        gentemp(1) = gentemp(2)
        f0pr(1)    = f0pr(2) 

        f0(2)      = falpha(parameters, egrid(ie, is), egrid(1, is), resolution)
        gentemp(2) = dfalpha_denergy(parameters, & 
                        egrid(ie, is), f0_values(ie, is), resolution)
        f0pr(2)    = falpha_prim(parameters,  &
                        egrid(ie, is), f0_values(ie, is), resolution)

        converged  = (is_converged(f0) .and.  &
                     is_converged(gentemp) .and. &
                     is_converged(f0pr))
        resolution = resolution * 2
      end do 
      f0_values(ie, is) = f0(2) 
      generalised_temperature(ie, is) = gentemp(2)
      f0prim(ie, is) = f0pr(2)
    end do

  end subroutine calculate_arrays

  function is_converged(list)
    real, dimension(2), intent(in) :: list
    logical :: is_converged

    if (abs(list(2) - list(1)/list(2)) .lt. 1.0e-6) then 
      is_converged = .true.
    else
      is_converged = .false.
    end if
  end function is_converged

  function unit_test_is_converged()
    real, dimension(2) :: test_array
    logical :: unit_test_is_converged

    test_array = (/1.0, 1.0/)
    unit_test_is_converged = is_converged(test_array)
    test_array = (/1.0e-5, 1.0003e-5/)
    unit_test_is_converged = (unit_test_is_converged .and. .not.  &
      is_converged(test_array))
  end function unit_test_is_converged


  !> Calculates the normalised alpha distribution function
  !! f_alpha/(n_alpha/v_inj^3)
  function falpha(parameters, energy, energy_0, resolution)
    implicit none
    type(analytical_falpha_parameters_type), intent(in) :: parameters
    integer, intent(in) :: resolution
    real, intent(in) :: energy
    real, intent(in) :: energy_0
    real :: falpha
    real :: integral
    real :: dv
    real :: v_2j, v_2jm1
    real :: energy_top
    integer :: j

    energy_top = min(energy, 1.0)
    integral = 0.0
    dv = (energy_top**0.5-energy_0**0.5) / real(resolution)

    ! Let's use the composite Simpson's rule (see Wikipedia!). 
    ! wrt the Wikipedia article
    ! resolution = n
    ! b = energy_top**0.5
    ! a = energy_0**0.5
    ! x_j = v_j = a + j*h = energy_0**).5 + j*dv
    ! dv = h = (b-a)/n = (energy_top**0.5-energy_0**0.5)/n
    integral = falpha_integrand(parameters, energy, energy_0)
    write (*,*) 'i1,', integral
    do j = 1, resolution/2-1
      v_2j = energy_0**0.5 + real(j*2)*dv 
      integral = integral + 2.0 * &
        falpha_integrand(parameters, energy, v_2j**2.0) 
      write (*,*) 'v_2j', v_2j, 'energy', v_2j**2.0, falpha_integrand(parameters, energy, v_2j**2.0)
    end do
    write (*,*) 'i2,', integral
    do j = 1, resolution/2
      v_2jm1 = energy_0**0.5 + real(j*2 -1)*dv 
      write (*,*) 'energy', v_2jm1**2.0
      integral = integral + 4.0 * &
        falpha_integrand(parameters, energy, v_2jm1**2.0) 
    end do
    write (*,*) 'i3,', integral, 'energy_top', energy_top
    integral = integral + falpha_integrand(parameters, energy, energy_top)
    write (*,*) 'i4,', integral, 'dv', dv
    integral = integral * dv/3.0
    write (*,*) 'i5,', integral

    falpha = integral




  end function falpha

  function analytical_falpha_unit_test_falpha(parameters, &
                                              energy, &
                                              energy_0, &
                                              resolution, &
                                              rslt, &
                                              err)
    use unit_tests
    type(analytical_falpha_parameters_type), intent(in) :: parameters
    integer, intent(in) :: resolution
    real, intent(in) :: energy
    real, intent(in) :: energy_0
    real, intent(in) :: rslt
    real, intent(in) :: err
    logical :: analytical_falpha_unit_test_falpha

    analytical_falpha_unit_test_falpha = &
      agrees_with(falpha(parameters, energy, energy_0, resolution), rslt, err)

  end function analytical_falpha_unit_test_falpha


  !> The integrand for the falpha integral, see eq
  function falpha_integrand(parameters, energy, energy_dummy_var)
    type(analytical_falpha_parameters_type), intent(in) :: parameters
    real, intent(in) :: energy
    real, intent(in) :: energy_dummy_var
    real :: falpha_integrand

    falpha_integrand = 2.0 / &
                       (nu_parallel(parameters, energy_dummy_var) * &
                         energy_dummy_var**2.0) * &
                       exp(parameters%alpha_mass * &
                         (energy_dummy_var - energy) / &
                         (2.0*parameters%ion_temp) &
                       )


  end function falpha_integrand

  !> Unit test used by test_analytical_falpha, testing the private
  !! function falpha_integrand
  function analytical_falpha_unit_test_falpha_integrand(parameters, &
                                                        energy, &
                                                        energy_dummy_var, &
                                                        rslt, &
                                                        err)
    use unit_tests
    type(analytical_falpha_parameters_type), intent(in) :: parameters
    real, intent(in) :: energy
    real, intent(in) :: energy_dummy_var
    real, intent(in) :: rslt
    real, intent(in) :: err
    logical :: analytical_falpha_unit_test_falpha_integrand
    
    analytical_falpha_unit_test_falpha_integrand = &
      agrees_with(falpha_integrand(parameters, energy, energy_dummy_var), &
      rslt, err)

  end function analytical_falpha_unit_test_falpha_integrand

  !> Normalised nu_parallel, see eq () of  
  function nu_parallel(parameters, energy)
    type(analytical_falpha_parameters_type), intent(in) :: parameters
    real, intent(in) :: energy
    real :: nu_parallel

    nu_parallel = 2.0/energy**1.5 * &
                  (parameters%alpha_ion_collision_rate * &
                    chandrasekhar(energy**0.5 * &
                      parameters%alpha_vth / parameters%ion_vth) - & 
                  parameters%alpha_electron_collision_rate * &
                    chandrasekhar(energy**0.5 * &
                      parameters%alpha_vth / parameters%electron_vth))
                  
                      


  end function nu_parallel

  !> Unit test used by test_analytical_falpha, testing the private
  !! function nu_parallel
  function analytical_falpha_unit_test_nu_parallel(parameters, energy, rslt, err)
    use unit_tests
    type(analytical_falpha_parameters_type), intent(in) :: parameters
    real, intent(in) :: rslt, energy, err
    logical :: analytical_falpha_unit_test_nu_parallel

    analytical_falpha_unit_test_nu_parallel = agrees_with(nu_parallel(parameters, energy), rslt, err)
  end function analytical_falpha_unit_test_nu_parallel



  function chandrasekhar(argument)
    real, intent(in) :: argument
    real ::  chandrasekhar
    chandrasekhar = (erf(argument) - &
        argument * 2.0 * exp(-argument**2.0) / 1.7724538509055159) / &
        (2.0 * argument**2.0)
  end function chandrasekhar

  function analytical_falpha_unit_test_chandrasekhar()
    logical :: analytical_falpha_unit_test_chandrasekhar
    real :: result1, result2

    result1 = 0.0373878
    result2 = 0.213797
    write (*,*) 'chandrasekhar(0.1): ', chandrasekhar(0.1), &
      ' should be ', result1
    write (*,*) 'chandrasekhar(1.0): ', chandrasekhar(1.0), &
      ' should be ', result2
    analytical_falpha_unit_test_chandrasekhar = &
      (abs(chandrasekhar(0.1)-result1)/result1 .lt. 1.0e-6 .and. &
      abs(chandrasekhar(1.0)-result2)/result2 .lt. 1.0e-5)
  end function analytical_falpha_unit_test_chandrasekhar


  !> Calculates the normalised d f_alpha / d energy
  !! which replaces temperature for Maxwellian species 
  !! in the GK equation and field equations
  !! = T_r/E_alpha * (d f_alpha / d energy) * E_alpha/f_alpha
  function dfalpha_denergy(parameters, energy, falph, resolution)
    type(analytical_falpha_parameters_type), intent(in) :: parameters
    integer, intent(in) :: resolution
    real, intent(in) :: energy
    real, intent(in) :: falph
    real :: dfalpha_denergy

  end function dfalpha_denergy

  !> Calculates the normalised gradient of f_alpha
  !! d(ln f_alpha) / d rho = d (ln f_alpha)/ d psi * d psi/d rho
  !! where rho is the GS2 flux label (usually r/a, see docs for irho)
  function falpha_prim(parameters, energy, falph, resolution)
    type(analytical_falpha_parameters_type), intent(in) :: parameters
    integer, intent(in) :: resolution
    real, intent(in) :: energy
    real, intent(in) :: falph
    real :: falpha_prim

  end function falpha_prim

end module analytical_falpha
