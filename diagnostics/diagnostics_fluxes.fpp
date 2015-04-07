!> Module which contains functions for calculating
!! and writing out the fluxes of heat and momentum etc.
module diagnostics_fluxes
 
  implicit none

  private

  !> Allocate arrays
  public :: init_diagnostics_fluxes

  !> Deallocate arrays
  public :: finish_diagnostics_fluxes

  !> Calculate and possibly write fluxes.  The fluxes are actually calculated
  !!  in the module dist_fn. They are returned as a function of
  !!  x, y and species. This function writes the whole array,
  !!  and also various averages of them.
  public :: calculate_fluxes
  public :: write_symmetry, write_parity

  real, dimension (:,:,:,:), allocatable ::  qheat, qmheat, qbheat
  real, dimension (:,:,:), allocatable ::  pflux,  vflux, vflux_par, vflux_perp
  real, dimension (:,:,:), allocatable ::  pflux_tormom
  real, dimension (:,:,:), allocatable :: vflux0, vflux1  ! low flow correction to turbulent momentum flux
  real, dimension (:,:,:), allocatable :: pmflux, vmflux
  real, dimension (:,:,:), allocatable :: pbflux, vbflux
  real, dimension (:,:,:), allocatable :: exchange
  real, dimension (:,:,:), allocatable :: exchange1

contains
  subroutine init_diagnostics_fluxes(gnostics)
    use kt_grids, only: naky, ntheta0
    use species, only: nspec
    use diagnostics_config, only: diagnostics_type
    implicit none
    type(diagnostics_type), intent(inout) :: gnostics
    
    gnostics%current_results%species_heat_flux_avg = 0.0
    gnostics%current_results%species_momentum_flux_avg = 0.0
    gnostics%current_results%species_particle_flux_avg = 0.0
  
    allocate (pflux (ntheta0,naky,nspec)) ; pflux = 0.
    allocate (pflux_tormom (ntheta0,naky,nspec)) ; pflux_tormom = 0. 
    allocate (qheat (ntheta0,naky,nspec,3)) ; qheat = 0.
    allocate (vflux (ntheta0,naky,nspec)) ; vflux = 0.
    
    allocate (exchange1 (ntheta0,naky,nspec)) ; exchange1 = 0.
    allocate (exchange (ntheta0,naky,nspec)) ; exchange = 0.

    allocate (vflux_par (ntheta0,naky,nspec)) ; vflux_par = 0.
    allocate (vflux_perp (ntheta0,naky,nspec)) ; vflux_perp = 0.
    
    allocate (vflux0 (ntheta0,naky,nspec)) ; vflux0 = 0.
    allocate (vflux1 (ntheta0,naky,nspec)) ; vflux1 = 0.
    
    allocate (pmflux(ntheta0,naky,nspec)) ; pmflux = 0.
    allocate (qmheat(ntheta0,naky,nspec,3)) ; qmheat = 0.
    allocate (vmflux(ntheta0,naky,nspec)) ; vmflux = 0.
    
    allocate (pbflux(ntheta0,naky,nspec)) ; pbflux = 0.
    allocate (qbheat(ntheta0,naky,nspec,3)) ; qbheat = 0.
    allocate (vbflux(ntheta0,naky,nspec)) ; vbflux = 0.   
  end subroutine init_diagnostics_fluxes
  
  subroutine finish_diagnostics_fluxes
    deallocate (pflux, pflux_tormom, qheat, vflux, vflux_par, vflux_perp, pmflux, qmheat, vmflux, &
         pbflux, qbheat, vbflux, vflux0, vflux1, exchange, exchange1)
  end subroutine finish_diagnostics_fluxes

  subroutine calculate_fluxes(gnostics)
    use dist_fn, only: flux, lf_flux, eexchange
    use dist_fn_arrays, only: g_adjust, gnew
    use species, only: nspec, spec
    use fields_arrays, only: phinew, bparnew, aparnew, phi
    use run_parameters, only: fphi, fapar, fbpar
    use nonlinear_terms, only: nonlinear_mode_switch, nonlinear_mode_none
    use diagnostics_create_and_write, only: create_and_write_variable
    use diagnostics_config, only: diagnostics_type
    use diagnostics_dimensions, only: dim_string
    implicit none
    type(diagnostics_type), intent(inout) :: gnostics
    integer :: is

    gnostics%current_results%total_heat_flux = 0.0
    gnostics%current_results%total_momentum_flux = 0.0
    gnostics%current_results%total_particle_flux = 0.0
    gnostics%current_results%species_heat_flux = 0.0
    gnostics%current_results%species_momentum_flux = 0.0
    gnostics%current_results%species_particle_flux = 0.0

    if (nonlinear_mode_switch .eq. nonlinear_mode_none) &
         call write_diffusive_estimates(gnostics)

    call g_adjust (gnew, phinew, bparnew, fphi, fbpar)
    call flux (phinew, aparnew, bparnew, &
         pflux,  qheat,  vflux, vflux_par, vflux_perp, &
         pmflux, qmheat, vmflux, pbflux, qbheat, vbflux, pflux_tormom)

#ifdef LOWFLOW
    !lowflow terms only implemented in electrostatic limit at present
    call lf_flux (phinew, vflux0, vflux1)
#endif

    call g_adjust (gnew, phinew, bparnew, -fphi, -fbpar)
    call eexchange (phinew, phi, exchange1, exchange)
    
    if (fphi > epsilon(0.0)) then
       do is = 1, nspec
          qheat(:,:,is,1) = qheat(:,:,is,1) * spec(is)%dens*spec(is)%temp
          qheat(:,:,is,2) = qheat(:,:,is,2) * spec(is)%dens*spec(is)%temp
          qheat(:,:,is,3) = qheat(:,:,is,3) * spec(is)%dens*spec(is)%temp
          
          pflux(:,:,is) = pflux(:,:,is) * spec(is)%dens
          pflux_tormom(:,:,is) = pflux_tormom(:,:,is) * spec(is)%dens  
          
          vflux(:,:,is) = vflux(:,:,is) * spec(is)%dens*sqrt(spec(is)%mass*spec(is)%temp)
          vflux_par(:,:,is) = vflux_par(:,:,is) * spec(is)%dens*sqrt(spec(is)%mass*spec(is)%temp)
          vflux_perp(:,:,is) = vflux_perp(:,:,is) * spec(is)%dens*spec(is)%mass*spec(is)%stm
          exchange(:,:,is) = exchange(:,:,is) * spec(is)%dens*spec(is)%z
          
#ifdef LOWFLOW
          vflux0(:,:,is) = vflux0(:,:,is) * spec(is)%dens*sqrt(spec(is)%mass*spec(is)%temp)
          vflux1(:,:,is) = vflux1(:,:,is) * spec(is)%dens*spec(is)%mass*spec(is)%temp/spec(is)%z
#endif
       end do
       call calculate_standard_flux_properties(gnostics, &
            'es_heat_flux',  'Turbulent flux of heat', 'Q_gB = ', qheat(:,:,:,1), gnostics%distributed)
       call calculate_standard_flux_properties(gnostics, &
            'es_heat_flux_par',  'Turbulent flux of parallel heat', 'Q_gB = ', qheat(:,:,:,2), gnostics%distributed)
       call calculate_standard_flux_properties(gnostics, &
            'es_heat_flux_perp',  'Turbulent flux of perpendicular heat', 'Q_gB = ', qheat(:,:,:,3), gnostics%distributed)
       call calculate_standard_flux_properties(gnostics, &
            'es_part_flux',  'Turbulent flux of particles', 'n_r? ', pflux, gnostics%distributed)
       call calculate_standard_flux_properties(gnostics, &
            'es_part_tormom_flux',  'Ask Jung-Pyo Lee...', '? ', pflux_tormom, gnostics%distributed)
       call calculate_standard_flux_properties(gnostics, &
            'es_mom_flux',  'Flux of toroidal angular momentum', 'Pi_gB =  ', vflux, gnostics%distributed)
       call calculate_standard_flux_properties(gnostics, &
            'es_mom_flux_par', 'Flux of the parallel component of the toroidal angular momentum', 'Pi_gB =  ', &
            vflux_par, gnostics%distributed)
       call calculate_standard_flux_properties(gnostics, &
            'es_energy_exchange', '??', 'Pi_gB =  ', &
            exchange, gnostics%distributed)
       call calculate_standard_flux_properties(gnostics, &
            'es_mom_flux_perp', 'Flux of the perpendicular component of the toroidal angular momentum', 'Pi_gB =  ', &
            vflux_perp, gnostics%distributed)
#ifdef LOWFLOW
       call calculate_standard_flux_properties(gnostics, &
            'es_mom0', 'Low-flow momentum flux 0 (Ask Michael Barnes)', 'Pi_gB =  ', &
            vflux0, gnostics%distributed)
       call calculate_standard_flux_properties(gnostics, &
            'es_mom1', 'Low-flow momentum flux 1 (Ask Michael Barnes)', 'Pi_gB =  ', &
            vflux0, gnostics%distributed)
#endif
    end if


    if (fapar > epsilon(0.0)) then
       do is = 1, nspec
          qmheat(:,:,is,1)=qmheat(:,:,is,1) * spec(is)%dens*spec(is)%temp
          qmheat(:,:,is,2)=qmheat(:,:,is,2) * spec(is)%dens*spec(is)%temp
          qmheat(:,:,is,3)=qmheat(:,:,is,3) * spec(is)%dens*spec(is)%temp
          pmflux(:,:,is)=pmflux(:,:,is) * spec(is)%dens
          vmflux(:,:,is)=vmflux(:,:,is) * spec(is)%dens*spec(is)%mass*spec(is)%stm
       end do
       call calculate_standard_flux_properties(gnostics, &
            'apar_heat_flux',  'Turbulent flux of heat resulting from &
            & perpendicular magnetic fluctuations', 'Q_gB = ', qmheat(:,:,:,1), gnostics%distributed)
       call calculate_standard_flux_properties(gnostics, &
            'apar_heat_flux_par',  'Turbulent flux of parallel heat resulting from &
            & perpendicular magnetic fluctuations', 'Q_gB = ', qmheat(:,:,:,2), gnostics%distributed)
       call calculate_standard_flux_properties(gnostics, &
            'apar_heat_flux_perp',  'Turbulent flux of perpendicular heat resulting from &
            & perpendicular magnetic fluctuations', 'Q_gB = ', qmheat(:,:,:,3), gnostics%distributed)
       call calculate_standard_flux_properties(gnostics, &
            'apar_part_flux',  'Turbulent flux of particles resulting from &
            & perpendicular magnetic fluctuations', 'TBC ', pmflux(:,:,:), gnostics%distributed)
       call calculate_standard_flux_properties(gnostics, &
            'apar_mom_flux',  'Turbulent flux of momentum resulting from &
            & perpendicular magnetic fluctuations', 'Pi_gB = ', vmflux(:,:,:), gnostics%distributed)
    end if
    if (fbpar > epsilon(0.0)) then
       do is = 1, nspec
          qbheat(:,:,is,1)=qbheat(:,:,is,1) * spec(is)%dens*spec(is)%temp
          qbheat(:,:,is,2)=qbheat(:,:,is,2) * spec(is)%dens*spec(is)%temp
          qbheat(:,:,is,3)=qbheat(:,:,is,3) * spec(is)%dens*spec(is)%temp
          pbflux(:,:,is)=pbflux(:,:,is) * spec(is)%dens
          vbflux(:,:,is)=vbflux(:,:,is) * spec(is)%dens*spec(is)%mass*spec(is)%stm
       end do
       call calculate_standard_flux_properties(gnostics, &
            'bpar_heat_flux',  'Turbulent flux of heat resulting from &
            & parallel magnetic fluctuations', 'Q_gB = ', qbheat(:,:,:,1), gnostics%distributed)
       call calculate_standard_flux_properties(gnostics, &
            'bpar_heat_flux_par',  'Turbulent flux of parallel heat resulting from &
            & parallel magnetic fluctuations', 'Q_gB = ', qbheat(:,:,:,2), gnostics%distributed)
       call calculate_standard_flux_properties(gnostics, &
            'bpar_heat_flux_perp',  'Turbulent flux of perpendicular heat resulting from &
            & parallel magnetic fluctuations', 'Q_gB = ', qbheat(:,:,:,3), gnostics%distributed)
       call calculate_standard_flux_properties(gnostics, &
            'bpar_part_flux',  'Turbulent flux of particles resulting from &
            & parallel magnetic fluctuations', 'TBC ', pbflux(:,:,:), gnostics%distributed)
       call calculate_standard_flux_properties(gnostics, &
            'bpar_mom_flux',  'Turbulent flux of momentum resulting from &
            & parallel magnetic fluctuations', 'Pi_gB = ', vbflux(:,:,:), gnostics%distributed)
    end if
    
    ! Update averages
    gnostics%current_results%species_heat_flux_avg = gnostics%current_results%species_heat_flux_avg + & 
         gnostics%current_results%species_heat_flux * (gnostics%user_time-gnostics%user_time_old)
    gnostics%current_results%species_particle_flux_avg = gnostics%current_results%species_particle_flux_avg + & 
         gnostics%current_results%species_particle_flux * (gnostics%user_time-gnostics%user_time_old)
    gnostics%current_results%species_momentum_flux_avg = gnostics%current_results%species_momentum_flux_avg + & 
         gnostics%current_results%species_momentum_flux * (gnostics%user_time-gnostics%user_time_old)
    ! Don't need to broadcast them any more as fluxes are now calculated on all
    ! processors
    
    ! Write totals
    if (gnostics%write_fluxes) then
       call create_and_write_variable(gnostics, gnostics%rtype, "heat_flux_tot", &
            dim_string([gnostics%dims%time]), &
            "Total heat flux, as a function of  time", &
            "Q_gB", gnostics%current_results%total_heat_flux)
       call create_and_write_variable(gnostics, gnostics%rtype, "hflux_tot", &
            dim_string([gnostics%dims%time]), &
            "Total heat flux, as a function of  time, same as heat_flux_tot, included &
            & for backwards compatiblity", &
            "Q_gB", gnostics%current_results%total_heat_flux)
       call create_and_write_variable(gnostics, gnostics%rtype, "mom_flux_tot", &
            dim_string([gnostics%dims%time]), &
            "Total momentum flux, as a function of  time", &
            "Pi_gB", gnostics%current_results%total_momentum_flux)
       call create_and_write_variable(gnostics, gnostics%rtype, "vflux_tot", &
            dim_string([gnostics%dims%time]), &
            "Total momentum flux, as a function of  time, same as mom_flux_tot, &
            & included for backwards compatiblity.", &
            "Pi_gB", gnostics%current_results%total_momentum_flux)
       call create_and_write_variable(gnostics, gnostics%rtype, "part_flux_tot", &
            dim_string([gnostics%dims%time]), &
            "Total particle flux, as a function of  time", &
            "TBC", gnostics%current_results%total_particle_flux)
       call create_and_write_variable(gnostics, gnostics%rtype, "zflux_tot", &
            dim_string([gnostics%dims%time]), &
            "Total particle flux, as a function of  time, same as part_flux_tot, &
            & included for backwards compatiblity.", &
            "TBC", gnostics%current_results%total_particle_flux)
    end if    
  end subroutine calculate_fluxes

  !> Calculate estimates of the heat and particles fluxes using
  !! gamma / k^2 estimate of the diffusivity
  subroutine write_diffusive_estimates(gnostics)
    use diagnostics_omega, only: omega_average
    use fields_parallelization, only: field_k_local
    use species, only: spec, nspec
    use kt_grids, only: kperp2, ntheta0, naky
    use theta_grid, only: grho
    use diagnostics_create_and_write, only: create_and_write_variable
    use diagnostics_dimensions, only: dim_string
    use diagnostics_config, only: diagnostics_type
    implicit none
    type(diagnostics_type), intent(inout) :: gnostics
    real, dimension(ntheta0, naky) :: diffusivity_by_k
    !real :: diffusivity
    real :: heat_flux_max
    real, dimension(nspec) :: heat_flux_max_by_spec
    real :: particle_flux_max
    real, dimension(nspec) :: particle_flux_max_by_spec
    !real :: grho
    real, dimension(ntheta0, naky, nspec) :: heat_flux
    real, dimension(ntheta0, naky, nspec) :: particle_flux
    real, dimension(naky, nspec) :: heat_flux_by_ky
    real, dimension(naky, nspec) :: particle_flux_by_ky
    real, dimension(ntheta0, naky) :: momentum_flux
    integer :: it, ik, is

    diffusivity_by_k = 0.0
    heat_flux = 0.0
    particle_flux = 0.0
    momentum_flux = 0.0
    heat_flux_max = 0.0

    do ik = 1,naky
       do it = 1,ntheta0
          if (.not. gnostics%distributed .or. field_k_local(it,ik)) then
             if (kperp2(gnostics%igomega,it,ik).eq.0.0) cycle
             diffusivity_by_k(it,ik) = &
                  max(aimag(omega_average(it,ik)),0.0)/kperp2(gnostics%igomega, it, ik)*2.0
             do is = 1,nspec
                ! Q = n chi grad T = n (gamma / k^2) dT / dr
                ! = dens  n_r (gamma_N v_thr / k_N**2 rho_r a) dT / drho drho/dr
                ! = dens  n_r (gamma_N v_thr rho_r **2 / k_N**2 a) T a / L_T drho/dr
                ! = dens  n_r (gamma_N v_thr rho_r **2 / k_N**2 a) temp T_r tprim drho/dr_N/a
                ! Q / (n_r  v_r T_r rho_r**2/a**2) 
                ! = dens (gamma_N / k_N**2) temp tprim grho
                !   
                heat_flux(it,ik,is) = diffusivity_by_k(it,ik) * &
                     spec(is)%dens * spec(is)%temp * spec(is)%tprim *  grho(gnostics%igomega)
                particle_flux(it,ik,is) = diffusivity_by_k(it,ik) * &
                     spec(is)%dens **2.0 * spec(is)%fprim * grho(gnostics%igomega)
             end do
          end if
       end do
    end do

    call calculate_standard_flux_properties(gnostics, &
         'heat_flux_diff',  'Diffusive estimate of turbulent flux of heat', &
         'Q_gB = ', heat_flux, gnostics%distributed)
    
    call calculate_standard_flux_properties(gnostics, &
         'part_flux_diff',  'Diffusive estimate of turbulent flux of particles', &
         'n_r? ', particle_flux, gnostics%distributed)
    
    heat_flux_by_ky = maxval(heat_flux, 1)
    particle_flux_by_ky = maxval(particle_flux, 1)
    
    heat_flux_max_by_spec = maxval(heat_flux_by_ky, 1)
    particle_flux_max_by_spec = maxval(particle_flux_by_ky, 1)
    
    heat_flux_max = sum(heat_flux_max_by_spec)
    particle_flux_max = sum(particle_flux_max_by_spec)
    
    if (gnostics%write_fluxes) then
       call create_and_write_variable(gnostics, gnostics%rtype, "es_heat_flux_diff_max", &
            dim_string([gnostics%dims%species,gnostics%dims%time]), &
            " A turbulent estimate of the heat flux, as a function of species and time", &
            "Q_gB", heat_flux_max_by_spec)
       call create_and_write_variable(gnostics, gnostics%rtype, "heat_flux_diff_max", &
            dim_string([gnostics%dims%time]), &
            " A turbulent estimate of the heat flux, as a function of  time", &
            "Q_gB", heat_flux_max)
       call create_and_write_variable(gnostics, gnostics%rtype, "es_particle_flux_diff_max", &
            dim_string([gnostics%dims%species,gnostics%dims%time]), &
            " A turbulent estimate of the particle flux, as a function of species and time", &
            "Q_gB", particle_flux_max_by_spec)
       call create_and_write_variable(gnostics, gnostics%rtype, "particle_flux_diff_max", &
            dim_string([gnostics%dims%time]), &
            " A turbulent estimate of the particle flux, as a function of  time", &
            "Q_gB", particle_flux_max)
    end if
  end subroutine write_diffusive_estimates

  !> Writes a range of different summed and averaged properties of the given
  !! flux... e.g. the flux summed over kx as a function of ky, species and time 
  subroutine calculate_standard_flux_properties(gnostics, flux_name, flux_description, &
       flux_units, flux_value, distributed)
    use kt_grids, only: ntheta0, naky
    use species, only: nspec
    use diagnostics_create_and_write, only: create_and_write_variable
    use diagnostics_create_and_write, only: create_and_write_distributed_fieldlike_variable
    use fields_parallelization, only: field_k_local
    use diagnostics_dimensions, only: dim_string
    use mp, only: broadcast, sum_allreduce
    use volume_averages, only: average_all, average_kx, average_ky
    use diagnostics_config, only: diagnostics_type
    implicit none
    type(diagnostics_type), intent(inout) :: gnostics
    character(*), intent(in) :: flux_name, flux_description, flux_units
    real, dimension(ntheta0,naky,nspec), intent(inout) :: flux_value
    real, dimension(ntheta0,naky,nspec) :: flux_value_temp
    logical, intent(in) :: distributed
    real, dimension(ntheta0, naky) :: total_flux_by_mode
    real, dimension(naky) :: total_flux_by_ky
    real, dimension(ntheta0) :: total_flux_by_kx
    real, dimension(nspec) :: flux_by_species
    integer :: is, ik, it

    call broadcast(flux_value)

    call average_all(flux_value, flux_by_species, distributed) 

    flux_value_temp = flux_value

    call broadcast(flux_value_temp)

    total_flux_by_mode = 0.
    do ik = 1,naky
      do it = 1,ntheta0
        if (.not. distributed .or. field_k_local(it, ik)) then
          do is = 1,nspec
            total_flux_by_mode(it, ik) = &
              total_flux_by_mode(it, ik) + flux_value(it, ik, is)
          end do
        end if
      end do
    end do

!Is this dead?
    !if (gnostics%write_fluxes_by_mode) then 
      !if (gnostics%write_fluxes) call create_and_write_flux_by_mode(gnostics, flux_name, flux_description, flux_units, &
        !flux_value, total_flux_by_mode, distributed)
    !end if

    call average_kx(total_flux_by_mode, total_flux_by_ky, distributed)
    call average_ky(total_flux_by_mode, total_flux_by_kx, distributed)
    if (gnostics%write_fluxes) then 
       call create_and_write_variable(gnostics, gnostics%rtype, flux_name, &
            dim_string([gnostics%dims%species,gnostics%dims%time]), &
            flux_description//" averaged over kx and ky, as a function of species and time", &
            flux_units, flux_by_species)
       call create_and_write_variable(gnostics, gnostics%rtype, "total_"//flux_name//"_by_ky", &
            dim_string([gnostics%dims%ky,gnostics%dims%time]), &
            flux_description//" summed over species and averaged over kx, as a function of ky and time", &
            flux_units, total_flux_by_ky)
       
       call create_and_write_variable(gnostics, gnostics%rtype, "total_"//flux_name//"_by_kx", &
            dim_string([gnostics%dims%kx,gnostics%dims%time]), &
            flux_description//" summed over species and averaged over ky, as a function of kx and time", &
            flux_units, total_flux_by_kx)
       
       call create_and_write_variable(gnostics, gnostics%rtype, "total_"//flux_name, &
            dim_string([gnostics%dims%time]), &
            flux_description//" summed over species and averaged over kx and ky, as a function of time", &
            flux_units, sum(total_flux_by_kx))
       
       if (gnostics%write_fluxes_by_mode) then 
          call create_and_write_distributed_fieldlike_variable( &
               gnostics, gnostics%rtype, flux_name//"_by_mode", &
               dim_string([gnostics%dims%kx,gnostics%dims%ky,gnostics%dims%species,gnostics%dims%time]), &
               flux_description//" as a function of species, kx and ky", &
               flux_units, flux_value)
          call create_and_write_distributed_fieldlike_variable( &
               gnostics, gnostics%rtype, "total_"//flux_name//"_by_mode", &
               dim_string([gnostics%dims%kx,gnostics%dims%ky,gnostics%dims%time]), &
               flux_description//" summed over species, as a function of kx and ky", &
               flux_units, total_flux_by_mode)
          
       end if
    end if
    
    if (flux_name .eq. 'es_heat_flux') gnostics%current_results%species_es_heat_flux = flux_by_species
    if (flux_name .eq. 'es_energy_exchange') gnostics%current_results%species_energy_exchange = flux_by_species
    if (flux_name .eq. 'apar_heat_flux') gnostics%current_results%species_apar_heat_flux = flux_by_species
    if (flux_name .eq. 'bpar_heat_flux') gnostics%current_results%species_bpar_heat_flux = flux_by_species
    
    if (flux_name .eq. 'es_heat_flux' &
         .or. flux_name .eq. 'apar_heat_flux' &
         .or. flux_name .eq. 'bpar_heat_flux') then 
       gnostics%current_results%total_heat_flux = gnostics%current_results%total_heat_flux + sum(flux_by_species)
       gnostics%current_results%species_heat_flux = gnostics%current_results%species_heat_flux + flux_by_species
    else if (flux_name .eq. 'es_mom_flux' &
         !.or.flux_name .eq. 'es_mom0' & ! Low flow fluxes, currently disabled
         .or.flux_name .eq. 'apar_mom_flux' &
         .or.flux_name .eq. 'bpar_mom_flux') then
       gnostics%current_results%total_momentum_flux = gnostics%current_results%total_momentum_flux + sum(flux_by_species)
       gnostics%current_results%species_momentum_flux = gnostics%current_results%species_momentum_flux + flux_by_species
    else if (flux_name .eq. 'es_part_flux' &
         .or.flux_name .eq. 'apar_part_flux' &
         .or.flux_name .eq. 'bpar_part_flux') then
       gnostics%current_results%total_particle_flux = gnostics%current_results%total_particle_flux + sum(flux_by_species)
       gnostics%current_results%species_particle_flux = gnostics%current_results%species_particle_flux + flux_by_species
    end if
  end subroutine calculate_standard_flux_properties


  !> 
  subroutine write_symmetry(gnostics)
    use dist_fn, only: flux_vs_theta_vs_vpa
    use dist_fn, only: pflux_vs_theta_vs_vpa
    use theta_grid, only: ntgrid
    use le_grids, only: nlambda, negrid
    use species, only: nspec
    use gs2_io, only: nc_loop_sym
    use fields_arrays, only: phinew
    use diagnostics_create_and_write, only: create_and_write_variable
    use diagnostics_dimensions, only: dim_string
    use diagnostics_config, only: diagnostics_type
    implicit none
    type (diagnostics_type), intent(in) :: gnostics
    real, dimension(:,:,:), allocatable :: vflx_sym
    real, dimension(:,:,:), allocatable :: pflux_sym

    allocate (pflux_sym(-ntgrid:ntgrid,nlambda*negrid,nspec))
    call pflux_vs_theta_vs_vpa (pflux_sym)
    call create_and_write_variable(gnostics, gnostics%rtype, "es_part_tormom_sym", &
         dim_string([gnostics%dims%theta,gnostics%dims%vpar,gnostics%dims%species,gnostics%dims%time]), &
         "Particle flux density as a function theta and vspace, used for looking at the effect of asymmetry, &
         & see Parra et al POP 18 062501 2011 and ask Jung-Pyo Lee", "See ref", pflux_sym)
    !if (proc0) call nc_loop_partsym_tormom (nout, pflux_sym)
    deallocate (pflux_sym)

    allocate (vflx_sym(-ntgrid:ntgrid,nlambda*negrid,nspec))
    call flux_vs_theta_vs_vpa (phinew,vflx_sym)
    call create_and_write_variable(gnostics, gnostics%rtype, "es_mom_sym", &
         dim_string([gnostics%dims%theta,gnostics%dims%vpar,gnostics%dims%species,gnostics%dims%time]), &
         "Momentum flux density as a function theta and vspace, used for looking at the effect of asymmetry, &
         & see Parra et al POP 18 062501 2011", "See ref", vflx_sym)
    !if (proc0) call nc_loop_sym (nout, vflx_sym)
    deallocate (vflx_sym)
  end subroutine write_symmetry

  subroutine write_parity(gnostics)
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use le_grids, only: negrid, nlambda, integrate_moment, integrate_kysum
    use species, only: nspec
    use gs2_layouts, only: init_parity_layouts, p_lo, g_lo, idx_local
    use gs2_layouts, only: ik_idx, it_idx, il_idx, ie_idx, is_idx, idx, proc_id
    use mp, only: send, receive, proc0
    use dist_fn_arrays, only: gnew, g_adjust, aj0
    use fields_arrays, only: phinew, bparnew
    use run_parameters, only: fphi, fbpar
    use constants, only: zi
    use volume_averages, only: average_theta, average_all
    use diagnostics_config, only: diagnostics_type
    implicit none
    type(diagnostics_type), intent(in) :: gnostics
    real :: t 
    integer :: iplo, iglo, sgn2, isgn, il, ie, ig, ik, it, is
    complex, dimension (:,:,:,:), allocatable :: gparity, gmx, gpx
    complex, dimension (:,:,:), allocatable :: g0, gm, gp
    complex, dimension (:,:,:), allocatable :: g_kint, gm_kint, gp_kint
    complex, dimension (:,:), allocatable :: g_avg, gnorm_avg, phim
    complex, dimension (:), allocatable :: g_all_tot, g_nokx_tot, g_nosig_tot, gtmp
    complex, dimension (:), allocatable :: gnorm_all_tot, gnorm_nokx_tot, gnorm_nosig_tot
    real, dimension (:,:,:), allocatable :: gmnorm, gmint, gpint, gmnormint, gmavg, gpavg, gmnormavg
    real, dimension (:), allocatable :: gmtot, gptot, gtot, gmnormtot
    real, dimension (:), allocatable :: gm_nokx_tot, gp_nokx_tot, gmnorm_nokx_tot
    real, dimension (:), allocatable :: gm_nosig_tot, gp_nosig_tot, gmnorm_nosig_tot
    logical :: first = .true.

    t = gnostics%user_time

    ! initialize layouts for parity diagnostic
    if (first) then
       call init_parity_layouts (naky, nlambda, negrid, nspec)
       first = .false.
    end if
    
    allocate (gparity(-ntgrid:ntgrid,ntheta0,2,p_lo%llim_proc:p_lo%ulim_alloc))
    allocate (g0(-ntgrid:ntgrid,2,p_lo%llim_proc:p_lo%ulim_alloc))
    allocate (gm(-ntgrid:ntgrid,2,p_lo%llim_proc:p_lo%ulim_alloc))
    allocate (gp(-ntgrid:ntgrid,2,p_lo%llim_proc:p_lo%ulim_alloc))
    allocate (gmnorm(-ntgrid:ntgrid,2,p_lo%llim_proc:p_lo%ulim_alloc))
    allocate (gmint(-ntgrid:ntgrid,naky,nspec))
    allocate (gpint(-ntgrid:ntgrid,naky,nspec))
    allocate (gmnormint(-ntgrid:ntgrid,naky,nspec))
    allocate (gmavg(ntheta0,naky,nspec))
    allocate (gpavg(ntheta0,naky,nspec))
    allocate (gmnormavg(ntheta0,naky,nspec))
    allocate (gmtot(nspec), gm_nokx_tot(nspec), gm_nosig_tot(nspec))
    allocate (gptot(nspec), gp_nokx_tot(nspec), gp_nosig_tot(nspec))
    allocate (gtot(nspec), gmnormtot(nspec), gmnorm_nokx_tot(nspec), gmnorm_nosig_tot(nspec))
    allocate (phim(-ntgrid:ntgrid,naky))
    allocate (g_kint(-ntgrid:ntgrid,2,nspec), gm_kint(-ntgrid:ntgrid,2,nspec), gp_kint(-ntgrid:ntgrid,2,nspec))
    allocate (g_avg(ntheta0,nspec), gnorm_avg(ntheta0,nspec))
    allocate (g_all_tot(nspec), g_nokx_tot(nspec), g_nosig_tot(nspec), gtmp(nspec))
    allocate (gnorm_all_tot(nspec), gnorm_nokx_tot(nspec), gnorm_nosig_tot(nspec))
    
    ! convert from g to h
    call g_adjust (gnew, phinew, bparnew, fphi, fbpar)
    
    ! below we define gparity = J0 * h, where delta f = h - (q phi / T) F0
    ! because we're ultimately interested in int J0 h * phi (i.e. particle flux)
    do isgn = 1, 2
       do ig = -ntgrid, ntgrid
          do iglo = g_lo%llim_world, g_lo%ulim_world
             ik = ik_idx(g_lo,iglo)
             ie = ie_idx(g_lo,iglo)
             il = il_idx(g_lo,iglo)
             is = is_idx(g_lo,iglo)
             it = it_idx(g_lo,iglo)
             ! count total index for given (ik,il,ie,is) combo (dependent on layout)
             iplo = idx(p_lo,ik,il,ie,is)
             ! check to see if value of g corresponding to iglo is on this processor
             if (idx_local(g_lo,iglo)) then
                if (idx_local(p_lo,iplo)) then
                   ! if g value corresponding to iplo should be on this processor, then get it
                   gparity(ig,it,isgn,iplo) = gnew(ig,isgn,iglo)*aj0(ig,iglo)
                else
                   ! otherwise, send g for this iplo value to correct processor
                   call send (gnew(ig,isgn,iglo)*aj0(ig,iglo),proc_id(p_lo,iplo))
                end if
             else if (idx_local(p_lo,iplo)) then
                ! if g for this iplo not on processor, receive it
                call receive (gparity(ig,it,isgn,iplo),proc_id(g_lo,iglo))
             end if
          end do
       end do
    end do
    
    ! convert from h back to g
    call g_adjust (gnew, phinew, bparnew, -fphi, -fbpar)
    
    ! first diagnostic is phase space average of 
    ! |J0*(h(z,vpa,kx) +/- h(-z,-vpa,-kx))|**2 / ( |J0*h(z,vpa,kx)|**2 + |J0*h(-z,-vpa,-kx)|**2 )
    do it = 1, ntheta0
       ! have to treat kx=0 specially
       if (it == 1) then
          do isgn = 1, 2
             sgn2 = mod(isgn,2)+1
             do ig = -ntgrid, ntgrid
                g0(ig,isgn,:) = gparity(-ig,1,sgn2,:)
             end do
          end do
       else
          do isgn = 1, 2
             sgn2 = mod(isgn,2)+1
             do ig = -ntgrid, ntgrid
                g0(ig,isgn,:) = gparity(-ig,ntheta0-it+2,sgn2,:)
             end do
          end do
       end if
       gm = gparity(:,it,:,:)-g0
       gp = gparity(:,it,:,:)+g0
       gmnorm = real(g0*conjg(g0))
       ! integrate out velocity dependence
       call integrate_moment (real(gm*conjg(gm)),gmint,1)
       call integrate_moment (real(gp*conjg(gp)),gpint,1)
       call integrate_moment (gmnorm,gmnormint,1)
       ! average along field line
       do is = 1, nspec
          do ik = 1, naky
             call average_theta (gmint(:,ik,is),gmavg(it,ik,is))
             call average_theta (gpint(:,ik,is),gpavg(it,ik,is))
             call average_theta (gmnormint(:,ik,is),gmnormavg(it,ik,is))
          end do
       end do

       ! phim(theta,kx) = phi(-theta,-kx)
       ! have to treat kx=0 specially
       if (it == 1) then
          do ig = -ntgrid, ntgrid
             phim(ig,:) = phinew(-ig,1,:)
          end do
       else
          do ig = -ntgrid, ntgrid
             phim(ig,:) = phinew(-ig,ntheta0-it+2,:)
          end do
       end if
       
       ! want int dtheta sum_{kx} int d3v | sum_{ky} [ J0*(h+ * conjg(phi-) + h- * conjg(phi+)) ky ] |**2
       ! J0*(h+ * conjg(phi-) + h- * conjg(phi+)) = h*conjg(phi) - h(-theta,-vpa,-kx)*conjg(phi(-theta,-kx))
       do iplo = p_lo%llim_proc, p_lo%ulim_proc
          ik = ik_idx(p_lo,iplo)
          do isgn = 1, 2
             gm(:,isgn,iplo) = gparity(:,it,isgn,iplo)*conjg(phinew(:,it,ik)) - g0(:,isgn,iplo)*conjg(phim(:,ik))
          end do
       end do
       do isgn = 1, 2
          do ig = -ntgrid, ntgrid
             call integrate_kysum (gm(ig,isgn,p_lo%llim_proc:p_lo%ulim_proc),ig,g_kint(ig,isgn,:),1)
          end do
       end do
       do is = 1, nspec
          call average_theta (g_kint(:,1,is)+g_kint(:,2,is),g_avg(it,is))
       end do
       
       ! get normalizing term for above diagnostic
       do iplo = p_lo%llim_proc, p_lo%ulim_proc
          ik = ik_idx(p_lo,iplo)
          do isgn = 1, 2
             gm(:,isgn,iplo) = gparity(:,it,isgn,iplo)*conjg(phinew(:,it,ik))
             gp(:,isgn,iplo) = g0(:,isgn,iplo)*conjg(phim(:,ik))
          end do
       end do
       do isgn = 1, 2
          do ig = -ntgrid, ntgrid
             call integrate_kysum (gm(ig,isgn,:),ig,gm_kint(ig,isgn,:),1)
             call integrate_kysum (gp(ig,isgn,:),ig,gp_kint(ig,isgn,:),1)
          end do
       end do
       do is = 1, nspec
          call average_theta (gm_kint(:,1,is)+gm_kint(:,2,is) &
               + gp_kint(:,1,is)+gp_kint(:,2,is),gnorm_avg(it,is))
       end do
    end do
    
    ! average over perp plane
    do is = 1, nspec
       call average_all (gmavg(:,:,is), gmtot(is), gnostics%distributed)
       call average_all (gpavg(:,:,is), gptot(is), gnostics%distributed)
       call average_all (gmnormavg(:,:,is), gmnormtot(is), gnostics%distributed)    
       g_all_tot(is) = sum(g_avg(:,is))
       gnorm_all_tot(is) = sum(gnorm_avg(:,is))
    end do
    
    allocate (gmx(-ntgrid:ntgrid,ntheta0,2,p_lo%llim_proc:p_lo%ulim_alloc))
    allocate (gpx(-ntgrid:ntgrid,ntheta0,2,p_lo%llim_proc:p_lo%ulim_alloc))
    
    ! now we want diagnostic of phase space average of
    ! |J0*(h(z,vpa) +/- h(-z,-vpa))|**2 / ( |J0*h(z,vpa)|**2 + |J0*h(-z,-vpa)|**2 )
    do it = 1, ntheta0
       do isgn = 1, 2
          sgn2 = mod(isgn,2)+1
          do ig = -ntgrid, ntgrid
             g0(ig,isgn,:) = gparity(-ig,it,sgn2,:)
          end do
       end do
       gm = gparity(:,it,:,:)-g0
       gp = gparity(:,it,:,:)+g0
       gmnorm = real(g0*conjg(g0))
       ! integrate out velocity dependence
       call integrate_moment (real(gm*conjg(gm)),gmint,1)
       call integrate_moment (real(gp*conjg(gp)),gpint,1)
       call integrate_moment (gmnorm,gmnormint,1)
       ! average along field line
       do is = 1, nspec
          do ik = 1, naky
             call average_theta (gmint(:,ik,is),gmavg(it,ik,is))
             call average_theta (gpint(:,ik,is),gpavg(it,ik,is))
             call average_theta (gmnormint(:,ik,is),gmnormavg(it,ik,is))
          end do
       end do
       
       ! phim(theta) = phi(-theta)
       ! have to treat kx=0 specially
       do ig = -ntgrid, ntgrid
          phim(ig,:) = phinew(-ig,it,:)
       end do
       
       ! want int dtheta int d3v | sum_{kx} sum_{ky} [ J0*(h+ * conjg(phi-) + h- * conjg(phi+)) ky ] |**2
       ! J0*(h+ * conjg(phi-) + h- * conjg(phi+)) = h*conjg(phi) - h(-theta,-vpa)*conjg(phi(-theta))
       do iplo = p_lo%llim_proc, p_lo%ulim_proc
          ik = ik_idx(p_lo,iplo)
          do isgn = 1, 2
             gmx(:,it,isgn,iplo) = g0(:,isgn,iplo)*conjg(phim(:,ik))
             gpx(:,it,isgn,iplo) = gparity(:,it,isgn,iplo)*conjg(phinew(:,it,ik)) - gmx(:,it,isgn,iplo) 
          end do
       end do
    end do
    
    ! sum over kx
    gp = sum(gpx,2)
    deallocate (gpx)
    do isgn = 1, 2
       do ig = -ntgrid, ntgrid
          call integrate_kysum (gp(ig,isgn,:),ig,g_kint(ig,isgn,:),1)
       end do
    end do
    do is = 1, nspec
       call average_theta (g_kint(:,1,is)+g_kint(:,2,is), g_nokx_tot(is))
    end do
    
    ! get normalizing terms for above diagnostic
    gm = sum(gmx,2)
    deallocate (gmx)
    do isgn = 1, 2
       do ig = -ntgrid, ntgrid
          call integrate_kysum (gm(ig,isgn,:),ig,gm_kint(ig,isgn,:),1)
          call integrate_kysum (gm(ig,isgn,:)+gp(ig,isgn,:),ig,gp_kint(ig,isgn,:),1)
       end do
    end do
    do is = 1, nspec
       call average_theta (gm_kint(:,1,is)+gm_kint(:,2,is) &
            + gp_kint(:,1,is)+gp_kint(:,2,is), gnorm_nokx_tot(is))
    end do
    
    ! average over perp plane
    do is = 1, nspec
       call average_all (gmavg(:,:,is), gm_nokx_tot(is), gnostics%distributed)
       call average_all (gpavg(:,:,is), gp_nokx_tot(is), gnostics%distributed)
       call average_all (gmnormavg(:,:,is), gmnorm_nokx_tot(is), gnostics%distributed)   
    end do
    
    ! final diagnostic is phase space average of 
    ! |J0*(h(z,kx) +/- h(-z,-kx))|**2 / ( |J0*h(z,kx)|**2 + |J0*h(-z,-kx)|**2 )
    do it = 1, ntheta0
       ! have to treat kx=0 specially
       if (it == 1) then
          do ig = -ntgrid, ntgrid
             g0(ig,:,:) = gparity(-ig,1,:,:)
          end do
       else
          do ig = -ntgrid, ntgrid
             g0(ig,:,:) = gparity(-ig,ntheta0-it+2,:,:)
          end do
       end if
       gm = gparity(:,it,:,:)-g0
       gp = gparity(:,it,:,:)+g0
       gmnorm = real(g0*conjg(g0))
       ! integrate out velocity dependence
       call integrate_moment (real(gm*conjg(gm)),gmint,1)
       call integrate_moment (real(gp*conjg(gp)),gpint,1)
       call integrate_moment (gmnorm,gmnormint,1)
       ! average along field line
       do is = 1, nspec
          do ik = 1, naky
             call average_theta (gmint(:,ik,is),gmavg(it,ik,is))
             call average_theta (gpint(:,ik,is),gpavg(it,ik,is))
             call average_theta (gmnormint(:,ik,is),gmnormavg(it,ik,is))
          end do
       end do
       
       ! phim(theta,kx) = phi(-theta,-kx)
       ! have to treat kx=0 specially
       if (it == 1) then
          do ig = -ntgrid, ntgrid
             phim(ig,:) = phinew(-ig,1,:)
          end do
       else
          do ig = -ntgrid, ntgrid
             phim(ig,:) = phinew(-ig,ntheta0-it+2,:)
          end do
       end if
       
       ! want int dtheta sum_{kx} int dE dmu | sum_{sigma} sum_{ky} [ J0*(h+ * conjg(phi-) + h- * conjg(phi+)) ky ] |**2
       ! J0*(h+ * conjg(phi-) + h- * conjg(phi+)) = h*conjg(phi) - h(-theta,-kx)*conjg(phi(-theta,-kx))
       do iplo = p_lo%llim_proc, p_lo%ulim_proc
          ik = ik_idx(p_lo,iplo)
          do isgn = 1, 2
             gm(:,isgn,iplo) = gparity(:,it,isgn,iplo)*conjg(phinew(:,it,ik)) - g0(:,isgn,iplo)*conjg(phim(:,ik))
          end do
       end do
       do ig = -ntgrid, ntgrid
          call integrate_kysum (gm(ig,1,:)+gm(ig,2,:),ig,g_kint(ig,1,:),1)
       end do
       do is = 1, nspec
          call average_theta (g_kint(:,1,is),g_avg(it,is))
       end do
       
       ! get normalizing term for above diagnostic
       do iplo = p_lo%llim_proc, p_lo%ulim_proc
          ik = ik_idx(p_lo,iplo)
          do isgn = 1, 2
             gm(:,isgn,iplo) = gparity(:,it,isgn,iplo)*conjg(phinew(:,it,ik))
             gp(:,isgn,iplo) = g0(:,isgn,iplo)*conjg(phim(:,ik))
          end do
       end do
       do ig = -ntgrid, ntgrid
          call integrate_kysum (gm(ig,1,:)+gm(ig,2,:),ig,gm_kint(ig,1,:),1)
          call integrate_kysum (gp(ig,1,:)+gp(ig,2,:),ig,gp_kint(ig,1,:),1)
       end do
       do is = 1, nspec
          call average_theta (gm_kint(:,1,is)+gp_kint(:,1,is),gnorm_avg(it,is))
       end do
    end do

    ! average over perp plane
    do is = 1, nspec
       call average_all (gmavg(:,:,is), gm_nosig_tot(is), gnostics%distributed)
       call average_all (gpavg(:,:,is), gp_nosig_tot(is), gnostics%distributed)
       call average_all (gmnormavg(:,:,is), gmnorm_nosig_tot(is), gnostics%distributed)    
       g_nosig_tot(is) = sum(g_avg(:,is))
       gnorm_nosig_tot(is) = sum(gnorm_avg(:,is))
    end do
    
    deallocate (gparity) ; allocate (gparity(-ntgrid:ntgrid,ntheta0,naky,nspec))
    ! obtain normalization factor = int over phase space of |g|**2
    call g_adjust (gnew, phinew, bparnew, fphi, fbpar)
    !<DD>Do all processors need to know the full result here? Only proc0 seems to do any writing.
    !If not then remove the last two arguments in the following call.
    call integrate_moment (spread(aj0,2,2)*spread(aj0,2,2)*gnew*conjg(gnew), gparity, .true., full_arr=.true.)
    call g_adjust (gnew, phinew, bparnew, -fphi, -fbpar)
    do is = 1, nspec
       do ik = 1, naky
          do it = 1, ntheta0
             call average_theta (real(gparity(:,it,ik,is)),gmavg(it,ik,is))
          end do
       end do
       call average_all (gmavg(:,:,is), gtot(is), gnostics%distributed)
    end do
    
    ! normalize g(theta,vpa,kx) - g(-theta,-vpa,-kx) terms
    where (gtot+gmnormtot > epsilon(0.0))
       gmtot = sqrt(gmtot/(gtot+gmnormtot)) ; gptot = sqrt(gptot/(gtot+gmnormtot))
    elsewhere
       gmtot = sqrt(gmtot) ; gptot = sqrt(gptot)
    end where
    
    where (real(gnorm_all_tot) > epsilon(0.0))
       gtmp = sqrt(real(g_all_tot)/real(gnorm_all_tot))
    elsewhere
       gtmp = sqrt(real(g_all_tot))
    end where
    where (aimag(gnorm_all_tot) > epsilon(0.0))
       g_all_tot = gtmp + zi*sqrt(aimag(g_all_tot)/aimag(gnorm_all_tot))
    elsewhere
       g_all_tot = gtmp + zi*sqrt(aimag(g_all_tot))
    end where
    
    ! normalize g(theta,vpa) +/- g(-theta,-vpa) terms
    where (gtot+gmnorm_nokx_tot > epsilon(0.0))
       gm_nokx_tot = sqrt(gm_nokx_tot/(gtot+gmnorm_nokx_tot)) ; gp_nokx_tot = sqrt(gp_nokx_tot/(gtot+gmnorm_nokx_tot))
    elsewhere
       gm_nokx_tot = sqrt(gm_nokx_tot) ; gp_nokx_tot = sqrt(gp_nokx_tot)
    end where
    
    where (real(gnorm_nokx_tot) > epsilon(0.0))
       gtmp = sqrt(real(g_nokx_tot)/real(gnorm_nokx_tot))
    elsewhere
       gtmp = sqrt(real(g_nokx_tot))
    end where
    where (aimag(gnorm_nokx_tot) > epsilon(0.0))
       g_nokx_tot = gtmp + zi*sqrt(aimag(g_nokx_tot)/aimag(gnorm_nokx_tot))
    elsewhere
       g_nokx_tot = gtmp + zi*sqrt(aimag(g_nokx_tot))
    end where
    
    ! normalize g(theta,kx) +/ g(-theta,-kx) terms
    where (gtot+gmnorm_nosig_tot > epsilon(0.0))
       gm_nosig_tot = sqrt(gm_nosig_tot/(gtot+gmnorm_nosig_tot)) ; gp_nosig_tot = sqrt(gp_nosig_tot/(gtot+gmnorm_nosig_tot))
    elsewhere
       gm_nosig_tot = sqrt(gm_nosig_tot) ; gp_nosig_tot = sqrt(gp_nosig_tot)
    end where
    
    where (real(gnorm_nosig_tot) > epsilon(0.0))
       gtmp = sqrt(real(g_nosig_tot)/real(gnorm_nosig_tot))
    elsewhere
       gtmp = sqrt(real(g_nosig_tot))
    end where
    where (aimag(gnorm_nosig_tot) > epsilon(0.0))
       g_nosig_tot = gtmp + zi*sqrt(aimag(g_nosig_tot)/aimag(gnorm_nosig_tot))
    elsewhere
       g_nosig_tot = gtmp + zi*sqrt(aimag(g_nosig_tot))
    end where
    
    if (proc0 .and. gnostics%write_ascii) write (gnostics%ascii_files%parity,"(19(1x,e12.5))") &
         t, gmtot, gptot, real(g_all_tot), aimag(g_all_tot), &
         real(gnorm_all_tot), aimag(gnorm_all_tot), gm_nokx_tot, gp_nokx_tot, real(g_nokx_tot), aimag(g_nokx_tot), &
         real(gnorm_nokx_tot), aimag(gnorm_nokx_tot), gm_nosig_tot, gp_nosig_tot, &
         real(g_nosig_tot), aimag(g_nosig_tot), real(gnorm_nosig_tot), aimag(gnorm_nosig_tot)
    
    deallocate (gparity, g0)
    deallocate (gm, gp, gmnorm, gmint, gpint, gmnormint)
    deallocate (gmavg, gpavg, gmnormavg)
    deallocate (gmtot, gm_nokx_tot, gm_nosig_tot)
    deallocate (gptot, gp_nokx_tot, gp_nosig_tot)
    deallocate (gtot, gmnormtot, gmnorm_nokx_tot, gmnorm_nosig_tot)
    deallocate (g_avg, gnorm_avg)
    deallocate (g_kint, gm_kint, gp_kint)
    deallocate (g_all_tot, g_nokx_tot, g_nosig_tot, gtmp)
    deallocate (gnorm_all_tot, gnorm_nokx_tot, gnorm_nosig_tot)
    deallocate (phim)
  end subroutine write_parity
end module diagnostics_fluxes
