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

  !> Unit tests.. see general_f0 for a better example of 
  !! the latest unit testing practice
  public :: unit_test_is_converged
  public :: analytical_falpha_unit_test_chandrasekhar
  public :: analytical_falpha_unit_test_chandrasekhar_prime
  public :: analytical_falpha_unit_test_nu_parallel
  public :: analytical_falpha_unit_test_nu_parallel_prime
  public :: analytical_falpha_unit_test_falpha_integrand
  public :: analytical_falpha_unit_test_dfalpha_dti_integrand
  public :: analytical_falpha_unit_test_falpha
  public :: analytical_falpha_unit_test_dfalpha_dti
  public :: analytical_falpha_unit_test_calculate_arrays
  public :: analytical_falpha_unit_test_simpson

  public :: analytical_falpha_parameters_type
  type analytical_falpha_parameters_type
    integer :: alpha_is
    real :: energy_0
    real :: source
    real :: source_prim
    real :: alpha_mass
    real :: alpha_injection_energy
    real :: alpha_vth
    real :: alpha_charge
    real :: ion_mass
    real :: ion_temp
    real :: ion_vth
    real :: ion_charge
    real :: ion_tprim
    real :: ion_fprim
    real :: electron_mass
    real :: electron_temp
    real :: electron_vth
    real :: electron_charge
    real :: electron_tprim
    real :: electron_fprim
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
      resolution = 8
      converged = .false.
      do while (.not. converged)
        f0(1)      = f0(2)
        gentemp(1) = gentemp(2)
        f0pr(1)    = f0pr(2) 

        f0(2)      = falpha(parameters, egrid(ie, is), parameters%energy_0, resolution)
!        gentemp(2) = dfalpha_denergy(parameters, & 
!                        egrid(ie, is), f0(2), resolution)
!
! Changed definition here to agree with the definition of generalised temperature
! 1/T* = - (Tref/Eref) * d/dE ln F0
! (GW)
        gentemp(2) = - parameters%alpha_injection_energy / dfalpha_denergy(parameters, & 
                        egrid(ie, is), f0(2), resolution)
        f0pr(2)    = falpha_prim(parameters,  &
                        egrid(ie, is), f0(2), resolution)

        !write (*,*) 'egrid ', egrid(ie, is), ' f0 ', f0(2), ' gentemp ', gentemp(2), egrid(ie,is), parameters%energy_0
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

  function analytical_falpha_unit_test_calculate_arrays(parameters, &
                              egrid, &
                              f0_values, &
                              generalised_temperature, &
                              f0prim, &
                              f0_rslt, &
                              gentemp_rslt, &
                              f0prim_rslt, &
                              err)
    use unit_tests
    type(analytical_falpha_parameters_type), intent(in) :: parameters
    real, dimension(:,:), intent(in) :: egrid
    real, dimension(:,:), intent(out) :: f0_values
    real, dimension(:,:), intent(out) :: generalised_temperature
    real, dimension(:,:), intent(out) :: f0prim
    real, dimension(:,:), intent(in) :: f0_rslt
    real, dimension(:,:), intent(in) :: gentemp_rslt
    real, dimension(:,:), intent(in) :: f0prim_rslt
    real, intent(in) :: err
    logical :: analytical_falpha_unit_test_calculate_arrays

    integer :: i

    call calculate_arrays(parameters, egrid, f0_values, &
        generalised_temperature, f0prim)
    analytical_falpha_unit_test_calculate_arrays = .true.

    call announce_check('f0_values')
    analytical_falpha_unit_test_calculate_arrays = &
      analytical_falpha_unit_test_calculate_arrays .and. &
       agrees_with(f0_values(:, parameters%alpha_is), f0_rslt(:, parameters%alpha_is), err)
    call process_check(analytical_falpha_unit_test_calculate_arrays, 'f0 values')
    call announce_check('gentemp')
    analytical_falpha_unit_test_calculate_arrays = &
      analytical_falpha_unit_test_calculate_arrays .and. &
       agrees_with(generalised_temperature(:, parameters%alpha_is), gentemp_rslt(:, parameters%alpha_is), err)
    call process_check(analytical_falpha_unit_test_calculate_arrays, 'gentemp')
    do i = 1,parameters%negrid
      analytical_falpha_unit_test_calculate_arrays = &
        analytical_falpha_unit_test_calculate_arrays .and. &
        agrees_with(generalised_temperature(i, parameters%alpha_is), &
                    gentemp_rslt(i, parameters%alpha_is), err) .and. &
        agrees_with(f0prim(i, parameters%alpha_is), &
                    f0prim_rslt(i, parameters%alpha_is), err) 
    end do

  end function analytical_falpha_unit_test_calculate_arrays

  function is_converged(list)
    real, dimension(2), intent(in) :: list
    logical :: is_converged
    if (list(2) .eq. 0.0) then 
      is_converged = abs(list(1)) .lt. 1.0e-9
    else
      is_converged = (abs((list(1)-list(2))/list(2)) .lt. 1.0e-9)
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
    test_array = (/-1.0, 0.0/)
    unit_test_is_converged = (unit_test_is_converged .and. .not.  &
      is_converged(test_array))
  end function unit_test_is_converged


  !> Calculates the integral of func using the composite Simpson rule
  function simpson(parameters, func, lower_lim, upper_lim, energy, resolution)
    implicit none
    type(analytical_falpha_parameters_type), intent(in) :: parameters
    integer, intent(in) :: resolution
    real, external :: func
    real, intent(in) :: lower_lim
    real, intent(in) :: upper_lim
    real, intent(in) :: energy
    real :: simpson
    real :: integral
    real :: dv
    real :: v_2j, v_2jm1
    real :: energy_top
    integer :: j

    integral = 0.0
    dv = (upper_lim-lower_lim) / real(resolution)

    ! Let's use the composite Simpson's rule (see Wikipedia!). 
    ! wrt the Wikipedia article
    ! resolution = n
    ! b = upper_lim
    ! a = lower_lim
    ! x_j = v_j = a + j*h = lower_lim + j*dv
    ! dv = h = (b-a)/n = (upper_lim-lower_lim)/n
    ! NB the integral is wrt v, not energy, so we have to put v^2 in everywhere
    integral = func(parameters, energy, lower_lim**2.0)
    do j = 1, resolution/2-1
      v_2j = lower_lim + real(j*2)*dv 
      integral = integral + 2.0 * &
        func(parameters, energy, v_2j**2.0) 
    end do
    do j = 1, resolution/2
      v_2jm1 = lower_lim + real(j*2 -1)*dv 
      integral = integral + 4.0 * &
        func(parameters, energy, v_2jm1**2.0) 

    end do
    integral = integral + func(parameters, energy, upper_lim**2.0)
    integral = integral * dv/3.0

    simpson = integral





  end function simpson

  function simpson_test_func(parameters, dummy, energy)
    type(analytical_falpha_parameters_type), intent(in) :: parameters
    real, intent(in) :: dummy
    real, intent(in) :: energy
    real ::  simpson_test_func
    simpson_test_func = energy**1.5 + energy**2.0 ! v^3 + v^4
  end function simpson_test_func

  function analytical_falpha_unit_test_simpson(eps)
    use unit_tests
    real, intent(in) :: eps
    logical :: analytical_falpha_unit_test_simpson
    type(analytical_falpha_parameters_type) :: parameters
    real :: dummy

    analytical_falpha_unit_test_simpson = agrees_with(&
      simpson(parameters, simpson_test_func, 0.0, 5.3, dummy, 1024),&
      5.3**4.0/4.0 + 5.3**5.0/5.0,  eps)
  end function analytical_falpha_unit_test_simpson
      


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

    ! Let's use the composite Simpson's rule (see Wikipedia!). 
    ! wrt the Wikipedia article
    ! resolution = n
    ! b = energy_top**0.5
    ! a = energy_0**0.5
    ! x_j = v_j = a + j*h = energy_0**).5 + j*dv
    ! dv = h = (b-a)/n = (energy_top**0.5-energy_0**0.5)/n
    ! NB the integral is wrt v, not energy

    integral = simpson(parameters, &
                       falpha_integrand, &
                       energy_0**0.5,&
                       energy_top**0.5,&
                       energy,&
                       resolution)

    falpha = integral * parameters%source / 4.0 / 3.14159265358979




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
                       exp(parameters%alpha_injection_energy * &
                         (energy_dummy_var - energy) / &
                         parameters%ion_temp &
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

  function dfalpha_dti(parameters, falph, energy, energy_0, resolution)
    implicit none
    type(analytical_falpha_parameters_type), intent(in) :: parameters
    integer, intent(in) :: resolution
    real, intent(in) :: energy
    real, intent(in) :: energy_0
    real, intent(in) :: falph
    real :: dfalpha_dti
    real :: integral
    real :: dv
    real :: v_2j, v_2jm1
    real :: energy_top
    integer :: j

    energy_top = min(energy, 1.0)

    ! Let's use the composite Simpson's rule (see Wikipedia!). 
    ! wrt the Wikipedia article
    ! resolution = n
    ! b = energy_top**0.5
    ! a = energy_0**0.5
    ! x_j = v_j = a + j*h = energy_0**).5 + j*dv
    ! dv = h = (b-a)/n = (energy_top**0.5-energy_0**0.5)/n
    ! NB the integral is wrt v, not energy

    integral = simpson(parameters, &
                       dfalpha_dti_integrand, &
                       energy_0**0.5,&
                       energy_top**0.5,&
                       energy,&
                       resolution)

    dfalpha_dti = - parameters%ion_tprim  &
      * parameters%alpha_injection_energy  &
      / parameters%ion_temp * (&
      energy - &
      integral * parameters%source / 2.0 / 3.14159265358979 / falph)
    




  end function dfalpha_dti

  function analytical_falpha_unit_test_dfalpha_dti(parameters, &
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
    logical :: analytical_falpha_unit_test_dfalpha_dti

    analytical_falpha_unit_test_dfalpha_dti = &
      agrees_with(&
      dfalpha_dti(parameters, &
                  falpha(parameters, energy, energy_0, resolution), &
                  energy, energy_0, resolution), rslt, err)

  end function analytical_falpha_unit_test_dfalpha_dti
  !> The integrand for the dfalpha/d Ti integral, see eq
  function dfalpha_dti_integrand(parameters, energy, energy_dummy_var)
    type(analytical_falpha_parameters_type), intent(in) :: parameters
    real, intent(in) :: energy
    real, intent(in) :: energy_dummy_var
    real :: dfalpha_dti_integrand

    dfalpha_dti_integrand = 1.0 / &
                       (nu_parallel(parameters, energy_dummy_var) * &
                         energy_dummy_var) * &
                       exp(parameters%alpha_injection_energy * &
                         (energy_dummy_var - energy) / &
                         parameters%ion_temp &
                       )


  end function dfalpha_dti_integrand

  !> Unit test used by test_analytical_falpha, testing the private
  !! function dfalpha_dti_integrand
  function analytical_falpha_unit_test_dfalpha_dti_integrand(parameters, &
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
    logical :: analytical_falpha_unit_test_dfalpha_dti_integrand
    
    analytical_falpha_unit_test_dfalpha_dti_integrand = &
      agrees_with(dfalpha_dti_integrand(parameters, energy, energy_dummy_var), &
      rslt, err)

  end function analytical_falpha_unit_test_dfalpha_dti_integrand


  !> Returns gamma_aiN * (Za/Zi)^2 * (mi/ma)^2 * (vthi/vtha)^3
  !! where gamma_ai is the alpha-ion collisionality parameter
  !! from eq ()
  function gamma_fac_ions(parameters)
    type(analytical_falpha_parameters_type), intent(in) :: parameters
    real gamma_fac_ions
    gamma_fac_ions = parameters%alpha_ion_collision_rate * &
                    (parameters%ion_vth / parameters%alpha_vth)**3.0 * &
                    (parameters%ion_mass / parameters%alpha_mass * &
                     parameters%alpha_charge / parameters%ion_charge)**2.0
  end function gamma_fac_ions

  function gamma_fac_electrons(parameters)
    type(analytical_falpha_parameters_type), intent(in) :: parameters
    real gamma_fac_electrons
    gamma_fac_electrons = parameters%alpha_electron_collision_rate * &
                    (parameters%electron_vth / parameters%alpha_vth)**3.0 * &
                    (parameters%electron_mass / parameters%alpha_mass * &
                     parameters%alpha_charge / parameters%electron_charge)**2.0  
  end function gamma_fac_electrons

  !> Calculates d nu_parallel_N/d rho, see eq () of  
  function nu_parallel_prime(parameters, energy)
    type(analytical_falpha_parameters_type), intent(in) :: parameters
    real, intent(in) :: energy
    real :: nu_parallel_prime

    nu_parallel_prime = 2.0/energy**1.5 * &
                  (0.5 * gamma_fac_ions(parameters) * &
                    chandrasekhar_prime(energy**0.5 * &
                      parameters%alpha_vth / parameters%ion_vth) * &
                    energy**0.5 * parameters%alpha_vth / parameters%ion_vth * &
                    parameters%ion_tprim & 
                  + 0.5 * gamma_fac_electrons(parameters) * &
                    chandrasekhar_prime(energy**0.5 * &
                      parameters%alpha_vth / parameters%electron_vth) * &
                    energy**0.5 * parameters%alpha_vth / parameters%electron_vth * &
                    parameters%electron_tprim &
                  - gamma_fac_ions(parameters) * &
                    parameters%ion_fprim * &
                    chandrasekhar(energy**0.5 * &
                      parameters%alpha_vth / parameters%ion_vth) & 
                  - gamma_fac_electrons(parameters) * &
                    parameters%electron_fprim * &
                    chandrasekhar(energy**0.5 * &
                      parameters%alpha_vth / parameters%electron_vth))
                  
                      


  end function nu_parallel_prime

  !> Unit test used by test_analytical_falpha, testing the private
  !! function nu_parallel
  function analytical_falpha_unit_test_nu_parallel_prime(parameters, energy, rslt, err)
    use unit_tests
    type(analytical_falpha_parameters_type), intent(in) :: parameters
    real, intent(in) :: rslt, energy, err
    logical :: analytical_falpha_unit_test_nu_parallel_prime

    analytical_falpha_unit_test_nu_parallel_prime = agrees_with(nu_parallel_prime(parameters, energy), rslt, err)
  end function analytical_falpha_unit_test_nu_parallel_prime
  !> Normalised nu_parallel, see eq () of  
  function nu_parallel(parameters, energy)
    type(analytical_falpha_parameters_type), intent(in) :: parameters
    real, intent(in) :: energy
    real :: nu_parallel

    nu_parallel = 2.0/energy**1.5 * &
                  (gamma_fac_ions(parameters) * &
                    chandrasekhar(energy**0.5 * &
                      parameters%alpha_vth / parameters%ion_vth) + & 
                     gamma_fac_electrons(parameters) * &
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

  function chandrasekhar_prime(argument)
    use constants, only: pi
    real, intent(in) :: argument
    real ::  chandrasekhar_prime
    chandrasekhar_prime = -2.0 / argument * chandrasekhar(argument) &
      + 2.0/pi**0.5*exp(-argument**2.0)
  end function chandrasekhar_prime

  function analytical_falpha_unit_test_chandrasekhar_prime(arg, rslt, eps)
    use unit_tests
    logical :: analytical_falpha_unit_test_chandrasekhar_prime
    real, intent(in) :: arg, rslt, eps

    analytical_falpha_unit_test_chandrasekhar_prime = &
      agrees_with(chandrasekhar_prime(arg), rslt, eps)
      
  end function analytical_falpha_unit_test_chandrasekhar_prime


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
    dfalpha_denergy = 0.0
    if (energy .gt. 1.0) then 
      dfalpha_denergy = - parameters%alpha_injection_energy / parameters%ion_temp
    else if (falph .eq. 0.0) then 
      dfalpha_denergy = 0.0
    else
      dfalpha_denergy = parameters%source / falph / (&
        4.0 * 3.14159265358979 * nu_parallel(parameters, energy) *  &
        energy ** (5.0/2.0) ) - &
        parameters%alpha_injection_energy / parameters%ion_temp
    end if
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
    falpha_prim = 0.0
    falpha_prim = - parameters%source_prim

  end function falpha_prim

end module analytical_falpha
