!> A module which generates an analytical form
!! for the alpha particle equilibrium distribution
!! function by integrating an expression given in...
!!
!! This is free software released under the GPLv3
!! Written by
!!    George Wilkie
!!    Edmund Highcock (edmundhighcock@sourceforge.net)
!!    Ian Abel ()


module semianalytic_falpha

  !> This subroutine does the work of the module.
  !! It takes in a grid of energy values and the species index
  !! and then, using the species parameter arrays, calculates
  !! f0, generalised_temperature (a function of df0/denergy)
  !! and df0/drho for that species.
  public :: calculate_arrays

  public :: semianalytic_falpha_parameters_type
  type semianalytic_falpha_parameters_type
    integer :: alpha_is
    real :: source
    real :: source_prim
    real :: alpha_mass
    real :: alpha_injection_energy
    real :: alpha_vth
    real :: alpha_charge
    real :: alpha_density
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
  end type semianalytic_falpha_parameters_type

  private 

  integer :: good_resolution

contains

  subroutine calculate_arrays(parameters, &
                              egrid, &
                              f0_values, &
                              generalised_temperature, &
                              f0prim)
    use mp, only: nproc, iproc, sum_allreduce, mp_abort
    use unit_tests, only: print_with_stars
    type(semianalytic_falpha_parameters_type), intent(inout) :: parameters
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
    real :: T_ash, df0dE

    is = parameters%alpha_is

    f0_values(:,is) = 0.0
    generalised_temperature(:,is) = 0.0
    f0prim(:,is) = 0.0
    do ie = 1,parameters%negrid
      if (.not. mod(ie, nproc) .eq. iproc) cycle

      !write (*,*) 'iproc ', iproc, ' calculating ', ie, ' energy ',  egrid(ie, is)
      ! Initialise arrays to test for convergence 
      f0 = -1.0
      gentemp = 0.0
      f0pr = 0.0
      resolution = 64
      converged = .false.
      do while (.not. converged)
        !write (*,*) 'resolution,', resolution
        f0(1)      = f0(2)
        f0(2)      = falpha(parameters, egrid(ie, is), resolution)
          converged  = is_converged(f0)
        !end if
        resolution = resolution * 2
        !write (*,*) ' f0', f0, ' f0pr', f0pr, ' gentemp', gentemp
      end do 
      resolution = 64
      converged = .false.
      do while (.not. converged)
        !write (*,*) 'resolution,', resolution
        f0(1)      = f0(2)
        gentemp(1) = gentemp(2)
        f0pr(1)    = f0pr(2) 

        f0(2)      = falpha(parameters, egrid(ie, is), resolution)
        gentemp(2) = - parameters%alpha_injection_energy *f0(2) / dfalpha_denergy(parameters,egrid(ie, is),  resolution)
        f0pr(2)    = falpha_prim(parameters,  &
                        egrid(ie, is), f0(2), resolution)

          converged  = (is_converged(f0) .and.  &
                       is_converged(gentemp) .and. &
                       is_converged(f0pr))
        resolution = resolution * 2
      end do 
      f0_values(ie, is) = f0(2)  * falpha_exponential(parameters , egrid(ie,is))
      generalised_temperature(ie, is) = gentemp(2)
      f0prim(ie, is) = f0pr(2)
    end do

    ! At this point we have H_alpha(v). We still need to add the ash.
    ! To do this, we now need to integrate the newly-found function
    ! Two options:
    !   - We can "is-converged" our way to a good Simpons grid
    !   - Import le_grids weights (which have already been calcualted!)
    ! The second option preempts the use of the generalized quadrature scheme, since that would result in a circular dependence 
    ! between the integration weights and F0. Use simpson instead.
    ! Note we lose the ability to say that the numerical integral of F0 is absolutely unity. This may cause problems.
    ! (GW)

    ! Take ash isothermal with main ions for now (GW)
    T_ash = parameters%ion_temp

    call add_ash(parameters,egrid(:,is),f0_values(:,is),generalised_temperature(:,is))

    call sum_allreduce(f0_values(:,is))
    call sum_allreduce(generalised_temperature(:,is))
    call sum_allreduce(f0prim(:,is))

  end subroutine calculate_arrays
  
  subroutine add_ash(parameters,egrid,f0,gentemp)
    use constants, only: pi
    implicit none
    type(semianalytic_falpha_parameters_type), intent(in) :: parameters
    real,dimension(:),intent(inout):: f0,gentemp
    real,dimension(:),intent(in):: egrid
    real:: f0moment, n_ash, T_ash,n_tot,Ealpha, energy_top, dF0dE
    real,dimension(2):: integral
    integer:: resolution, ie
    logical:: converged

    T_ash = parameters%ion_temp
    Ealpha = parameters%alpha_injection_energy
    n_tot = parameters%alpha_density
    energy_top = 1.2*Ealpha

    ! Find the 0th moment of H_alpha(v): the input distribution of fast alphas
    integral = -1.0
    resolution = 64
    converged = .false.
    do while (.not. converged)
       integral(1)      = integral(2)
       integral(2)      = simpson(parameters,falpha,0.0,energy_top**0.5,0.0,resolution)
       converged  = is_converged(integral)

       resolution = resolution * 2
    end do 
    ! Define n_ash
    n_ash = n_tot - f0moment

    do ie = 1,size(egrid)
       ! Unwrap gentemp, and 
       df0dE = -f0(ie) * Ealpha / gentemp(ie)
       df0dE = df0dE - n_ash*(Ealpha/T_ash)*exp(-egrid(ie)*Ealpha/T_ash)/(2.0*pi)**1.5
       gentemp(ie) = -f0(ie)*Ealpha/df0dE
      
       f0(ie) = f0(ie) + n_ash*exp(-egrid(ie)*Ealpha/T_ash)/(2.0*pi)**1.5
    end do

  end subroutine add_ash
 
  function is_converged(list)
    real, dimension(2), intent(in) :: list
    logical :: is_converged
    if (abs(list(2)) .lt. 1.0e-60) then 
      is_converged = abs(list(1)) .lt. 1.0e-60
    else
      is_converged = (abs((list(1)-list(2))/list(2)) .lt. 1.0e-7)
    end if

  end function is_converged

  !> Calculates the integral of func using the composite Simpson rule
  function simpson(parameters, func, lower_lim, upper_lim, energy, resolution)
    implicit none
    type(semianalytic_falpha_parameters_type), intent(in) :: parameters
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

  function falpha_integral_function(parameters, energy, energy_dummy_var)
    type(semianalytic_falpha_parameters_type), intent(in) :: parameters
    real, intent(in) :: energy
    real, intent(in) :: energy_dummy_var
    real :: falpha_integral_function
    falpha_integral_function = falpha(parameters, &
                                      energy_dummy_var, &
                                      good_resolution) * &
                               falpha_exponential(parameters, energy_dummy_var)*&
                               energy_dummy_var
   end function falpha_integral_function


  function falpha_exponential(parameters, energy)
    implicit none
    type(semianalytic_falpha_parameters_type), intent(in) :: parameters
    real, intent(in) :: energy
    real :: falpha_exponential
    falpha_exponential = exp(parameters%alpha_injection_energy * &
       (- energy) / parameters%ion_temp)
  end function falpha_exponential

  !> NB this does not calculate falpha, but falpha without the 
  !! exp(-E_alpha * energy/T_i) part, as this can be very small
  !! for some energies, and can mess up the calculation of the 
  !! gradients
  function falpha(parameters, energy, resolution)
    implicit none
    type(semianalytic_falpha_parameters_type), intent(in) :: parameters
    integer, intent(in) :: resolution
    real, intent(in) :: energy
    real :: falpha
    real :: integral
    real :: dv
    real :: v_2j, v_2jm1
    real :: energy_top
    integer :: j

!    energy_top = min(energy, 1.0)
    energy_top = energy

    integral = simpson(parameters, &
                       falpha_integrand, &
                       0.0,&
                       energy_top**0.5,&
                       0.0,& ! Get rid of exponential for better numerics
                       resolution)

    falpha = integral * parameters%source / 2.0 / 3.14159265358979 

  end function falpha


  !> The integrand for the falpha integral, see eq
  function falpha_integrand(parameters, energy, energy_dummy_var)
    type(semianalytic_falpha_parameters_type), intent(in) :: parameters
    real, intent(in) :: energy
    real, intent(in) :: energy_dummy_var
    real :: falpha_integrand
    real :: extra_fac, heaviside, arg

    if (energy_dummy_var .GT. 1.0) then
       heaviside = 0.0
    else
       heaviside = 1.0
    end if

    arg = sqrt(parameters%alpha_injection_energy*energy_dummy_var/parameters%ion_temp)
    extra_fac = 2.0*energy_dummy_var*chandrasekhar(arg) - 1.0 + heaviside

    falpha_integrand = 2.0 * extra_fac / &
                       (nu_parallel(parameters, energy_dummy_var) * &
                         energy_dummy_var**2.0) * &
                       exp(parameters%alpha_injection_energy * &
                         (energy_dummy_var - energy) / &
                         parameters%ion_temp &
                       )


  end function falpha_integrand

  !> Unit test used by test_semianalytic_falpha, testing the private
  !! function falpha_integrand
  function semianalytic_falpha_unit_test_falpha_integrand(parameters, &
                                                        energy, &
                                                        energy_dummy_var, &
                                                        rslt, &
                                                        err)
    use unit_tests
    type(semianalytic_falpha_parameters_type), intent(in) :: parameters
    real, intent(in) :: energy
    real, intent(in) :: energy_dummy_var
    real, intent(in) :: rslt
    real, intent(in) :: err
    logical :: semianalytic_falpha_unit_test_falpha_integrand
    
    semianalytic_falpha_unit_test_falpha_integrand = &
      agrees_with(falpha_integrand(parameters, energy, energy_dummy_var), &
      rslt, err)

  end function semianalytic_falpha_unit_test_falpha_integrand

  !> Calculates (1/Falpha * (dTi/drho) * (dFalpha/dTi)) 
  !! NB not dFalpha/dTi
  function dfalpha_dti(parameters, falph, energy, resolution)
    use constants, only: pi
    implicit none
    type(semianalytic_falpha_parameters_type), intent(in) :: parameters
    integer, intent(in) :: resolution
    real, intent(in) :: energy
    real, intent(in) :: falph
    real :: dfalpha_dti
    real :: integral
    real :: dv
    real :: v_2j, v_2jm1
    real :: energy_top
    integer :: j

!    energy_top = min(energy, 1.0)
    energy_top = energy

    integral = simpson(parameters, &
                       dfalpha_dti_integrand, &
                       0.0,&
                       energy_top**0.5,&
                       0.0,& ! Get rid of exp factor as it cancels with falpha
                       resolution) / falph !  (&
               !simpson(parameters, &
                       !falpha_integrand, &
                       !parameters%energy_0**0.5,&
                       !energy_top**0.5,&
                       !0.0,& ! Get rid of exp factor
                       !resolution) * parameters%source / 4.0 / pi)

    dfalpha_dti = - parameters%ion_tprim  &
      * parameters%alpha_injection_energy  &
      / parameters%ion_temp * (&
      energy - &
      integral * parameters%source / 2.0 / 3.14159265358979) ! / &
      !falpha(parameters, 0.0, parameters%energy_0, resolution))
    

  end function dfalpha_dti

  !> The integrand for the dfalpha/d Ti integral, see eq
  function dfalpha_dti_integrand(parameters, energy, energy_dummy_var)
    type(semianalytic_falpha_parameters_type), intent(in) :: parameters
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

  !> Unit test used by test_semianalytic_falpha, testing the private
  !! function dfalpha_dti_integrand
  function semianalytic_falpha_unit_test_dfalpha_dti_integrand(parameters, &
                                                        energy, &
                                                        energy_dummy_var, &
                                                        rslt, &
                                                        err)
    use unit_tests
    type(semianalytic_falpha_parameters_type), intent(in) :: parameters
    real, intent(in) :: energy
    real, intent(in) :: energy_dummy_var
    real, intent(in) :: rslt
    real, intent(in) :: err
    logical :: semianalytic_falpha_unit_test_dfalpha_dti_integrand
    
    semianalytic_falpha_unit_test_dfalpha_dti_integrand = &
      agrees_with(dfalpha_dti_integrand(parameters, energy, energy_dummy_var), &
      rslt, err)

  end function semianalytic_falpha_unit_test_dfalpha_dti_integrand

  !> Calculates 1/Falpha dFalpha/dnu_par d nu_par d rho
  !! where d nu_par / d rho goes under the integral for Falpha
  function dfalpha_dnupar(parameters, falph, energy, resolution)
    use constants, only: pi
    implicit none
    type(semianalytic_falpha_parameters_type), intent(in) :: parameters
    integer, intent(in) :: resolution
    real, intent(in) :: energy
    real, intent(in) :: falph
    real :: dfalpha_dnupar
    real :: integral
    real :: dv
    real :: v_2j, v_2jm1
    real :: energy_top
    integer :: j

    energy_top = energy

    ! Let's use the composite Simpson's rule (see Wikipedia!). 
    ! wrt the Wikipedia article
    ! resolution = n
    ! b = energy_top**0.5
    ! a = energy_0**0.5
    ! x_j = v_j = a + j*h = energy_0**).5 + j*dv
    ! dv = h = (b-a)/n = (energy_top**0.5-energy_0**0.5)/n
    ! NB the integral is wrt v, not energy

    integral = simpson(parameters, &
                       dfalpha_dnupar_integrand, &
                       0.0,&
                       energy_top**0.5,&
                       0.0,&
                       resolution)/  falph! 

    dfalpha_dnupar = - integral * parameters%source / 4.0 / pi

  end function dfalpha_dnupar

  !> The integrand for the dfalpha/d nu_parallel integral, see eq
  function dfalpha_dnupar_integrand(parameters, energy, energy_dummy_var)
    type(semianalytic_falpha_parameters_type), intent(in) :: parameters
    real, intent(in) :: energy
    real, intent(in) :: energy_dummy_var
    real :: dfalpha_dnupar_integrand

    !dfalpha_dnupar_integrand = 2.0 / &
                       !(nu_parallel(parameters, energy_dummy_var)**2.0 * &
                         !energy_dummy_var**2.0) * &
                       !exp(parameters%alpha_injection_energy * &
                         !(energy_dummy_var - energy) / &
                         !parameters%ion_temp &
                       !)
    dfalpha_dnupar_integrand = &
      falpha_integrand(parameters, energy, energy_dummy_var) / &
      nu_parallel(parameters, energy_dummy_var) &
        * nu_parallel_prime(parameters, energy_dummy_var)


  end function dfalpha_dnupar_integrand

  !> Returns gamma_aiN * (Za/Zi)^2 * (mi/ma)^2 * (vthi/vtha)^3
  !! where gamma_ai is the alpha-ion collisionality parameter
  !! from eq ()
  function gamma_fac_ions(parameters)
    type(semianalytic_falpha_parameters_type), intent(in) :: parameters
    real gamma_fac_ions
    gamma_fac_ions = parameters%alpha_ion_collision_rate * &
                    (parameters%ion_vth / parameters%alpha_vth)**3.0 * &
                    (parameters%ion_mass / parameters%alpha_mass * &
                     parameters%alpha_charge / parameters%ion_charge)**2.0
  end function gamma_fac_ions

  function gamma_fac_electrons(parameters)
    type(semianalytic_falpha_parameters_type), intent(in) :: parameters
    real gamma_fac_electrons
    gamma_fac_electrons = parameters%alpha_electron_collision_rate * &
                    (parameters%electron_vth / parameters%alpha_vth)**3.0 * &
                    (parameters%electron_mass / parameters%alpha_mass * &
                     parameters%alpha_charge / parameters%electron_charge)**2.0  
  end function gamma_fac_electrons

  !> Calculates d nu_parallel_N/d rho, see eq () of  
  function nu_parallel_prime(parameters, energy)
    type(semianalytic_falpha_parameters_type), intent(in) :: parameters
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

  !> Normalised nu_parallel, see eq () of  
  function nu_parallel(parameters, energy)
    type(semianalytic_falpha_parameters_type), intent(in) :: parameters
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

  function chandrasekhar(argument)
    real, intent(in) :: argument
    real ::  chandrasekhar
    chandrasekhar = (erf(argument) - &
        argument * 2.0 * exp(-argument**2.0) / 1.7724538509055159) / &
        (2.0 * argument**2.0)
  end function chandrasekhar

  function chandrasekhar_prime(argument)
    use constants, only: pi
    real, intent(in) :: argument
    real ::  chandrasekhar_prime
    chandrasekhar_prime = -2.0 / argument * chandrasekhar(argument) &
      + 2.0/pi**0.5*exp(-argument**2.0)
  end function chandrasekhar_prime

  !> Calculates the normalised d f_alpha / d energy
  !! which replaces temperature for Maxwellian species 
  !! in the GK equation and field equations
  !! = T_r/E_alpha * (d f_alpha / d energy) * E_alpha/f_alpha
  real function dfalpha_denergy(parameters, energy, resolution)
    type(semianalytic_falpha_parameters_type), intent(in) :: parameters
    integer, intent(in) :: resolution
    real, intent(in) :: energy
    real :: sum

    sum = - parameters%alpha_injection_energy / parameters%ion_temp
    sum = sum + (0.5/sqrt(energy))*parameters%source *falpha_integrand(parameters,energy,energy) / (2.0 * 3.14159265358979 ) 
 
    dfalpha_denergy = 1.0/sum
    return
  end function dfalpha_denergy

  !> Calculates the normalised gradient of f_alpha
  !! d(ln f_alpha) / d rho = d (ln f_alpha)/ d psi * d psi/d rho
  !! where rho is the GS2 flux label (usually r/a, see docs for irho)
  function falpha_prim(parameters, energy, falph, resolution)
    type(semianalytic_falpha_parameters_type), intent(in) :: parameters
    integer, intent(in) :: resolution
    real, intent(in) :: energy
    real, intent(in) :: falph
    real :: falpha_prim
    falpha_prim = 0.0
    falpha_prim = - parameters%source_prim &
      + dfalpha_dti(parameters, falph, energy, resolution) &
      + dfalpha_dnupar(parameters, falph, energy, resolution)  
  end function falpha_prim

end module semianalytic_falpha
