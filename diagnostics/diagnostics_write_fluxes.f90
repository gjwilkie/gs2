!> Module which contains functions for calculating
!! and writing out the fluxes of heat and momentum etc.
module diagnostics_write_fluxes

  use kt_grids, only: naky, ntheta0
  use diagnostics_config, only: diagnostics_type
  use species, only: nspec
  use theta_grid, only: ntgrid
  use simpledataio
  use simpledataio_write
  

implicit none

!> Allocate arrays
public :: init_diagnostics_write_fluxes

!> Deallocate arrays
public :: finish_diagnostics_write_fluxes

!> Calculate and write fluxes.  The fluxes are actually calculated
!!  in the module dist_fn. They are returned as a function of
!!  x, y and species. This function writes the whole array,
!!  and also various averages of them.
public :: write_fluxes

private
  real, dimension (:,:,:,:), allocatable ::  qheat, qmheat, qbheat
  real, dimension (:,:,:), allocatable ::  pflux,  vflux, vflux_par, vflux_perp
  real, dimension (:,:,:), allocatable ::  pflux_tormom
  real, dimension (:,:,:), allocatable :: vflux0, vflux1  ! low flow correction to turbulent momentum flux
  real, dimension (:,:,:), allocatable :: pmflux, vmflux
  real, dimension (:,:,:), allocatable :: pbflux, vbflux
  real, dimension (:,:,:), allocatable :: exchange
  !integer, parameter :: gnostics%rtype = SDATIO_DOUBLE


contains
  subroutine init_diagnostics_write_fluxes
  
    allocate (pflux (ntheta0,naky,nspec)) ; pflux = 0.
    allocate (pflux_tormom (ntheta0,naky,nspec)) ; pflux_tormom = 0. 
    allocate (qheat (ntheta0,naky,nspec,3)) ; qheat = 0.
    allocate (vflux (ntheta0,naky,nspec)) ; vflux = 0.
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
  end subroutine init_diagnostics_write_fluxes

  subroutine finish_diagnostics_write_fluxes
    deallocate (pflux, qheat, vflux, vflux_par, vflux_perp, pmflux, qmheat, vmflux, &
         pbflux, qbheat, vbflux, vflux0, vflux1, exchange)
  end subroutine finish_diagnostics_write_fluxes

  subroutine write_fluxes(gnostics)
    use dist_fn, only: flux
    use dist_fn_arrays, only: g_adjust, gnew
    use species, only: nspec, spec
    use fields_arrays, only: phinew, bparnew, aparnew
    use run_parameters, only: fphi, fapar, fbpar
    use nonlinear_terms, only: nonlinear_mode_switch, nonlinear_mode_none
    type(diagnostics_type), intent(in) :: gnostics
    integer :: is
    !if (istep > 0) then

    if (nonlinear_mode_switch .eq. nonlinear_mode_none) &
      call write_diffusive_estimates(gnostics)
    call g_adjust (gnew, phinew, bparnew, fphi, fbpar)
    call flux (phinew, aparnew, bparnew, &
      pflux,  qheat,  vflux, vflux_par, vflux_perp, &
      pmflux, qmheat, vmflux, pbflux, qbheat, vbflux, pflux_tormom)
!#ifdef LOWFLOW
     ! lowflow terms only implemented in electrostatic limit at present
     !call lf_flux (phinew, vflux0, vflux1)
!#endif
    call g_adjust (gnew, phinew, bparnew, -fphi, -fbpar)
    !end if
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
      end do
      call write_standard_flux_properties(gnostics, &
        'heat_flux',  'Turbulent flux of heat', 'Q_gB = ', qheat(:,:,:,1), gnostics%distributed)
      call write_standard_flux_properties(gnostics, &
        'heat_flux_par',  'Turbulent flux of parallel heat', 'Q_gB = ', qheat(:,:,:,2), gnostics%distributed)
      call write_standard_flux_properties(gnostics, &
        'heat_flux_perp',  'Turbulent flux of perpendicular heat', 'Q_gB = ', qheat(:,:,:,3), gnostics%distributed)

      call write_standard_flux_properties(gnostics, &
        'part_flux',  'Turbulent flux of particles', 'n_r? ', pflux, gnostics%distributed)
      call write_standard_flux_properties(gnostics, &
        'part_tormom_flux',  'Ask Jung-Pyo Lee...', '? ', pflux_tormom, gnostics%distributed)

      call write_standard_flux_properties(gnostics, &
        'mom_flux',  'Flux of toroidal angular momentum', 'Pi_gB =  ', vflux, gnostics%distributed)
      call write_standard_flux_properties(gnostics, &
        'mom_flux_par',  &
        'Flux of the parallel component of the toroidal angular momentum', 'Pi_gB =  ', &
        vflux_par, gnostics%distributed)
      call write_standard_flux_properties(gnostics, &
        'mom_flux_perp', &
        'Flux of the perpendicular component of the toroidal angular momentum', 'Pi_gB =  ', &
        vflux_perp, gnostics%distributed)
    end if


          !call get_volume_average (vflux(:,:,is), mom_fluxes(is))

          !call get_volume_average (vflux_par(:,:,is), parmom_fluxes(is))

          !call get_volume_average (vflux_perp(:,:,is), perpmom_fluxes(is))

          !exchange(:,:,is) = exchange(:,:,is) * spec(is)%dens*spec(is)%z
          !call get_volume_average (exchange(:,:,is), energy_exchange(is))

!#ifdef LOWFLOW
          !vflux0(:,:,is) = vflux0(:,:,is) * spec(is)%dens*sqrt(spec(is)%mass*spec(is)%temp)
          !call get_volume_average (vflux0(:,:,is), lfmom_fluxes(is))

          !vflux1(:,:,is) = vflux1(:,:,is) * spec(is)%dens*spec(is)%mass*spec(is)%temp/spec(is)%z
          !call get_volume_average (vflux1(:,:,is), vflux1_avg(is))

!! TMP UNTIL VFLUX0 IS TESTED
!!                   mom_fluxes = mom_fluxes + lfmom_fluxes
!#endif
       !end do
    !end if
    if (fapar > epsilon(0.0)) then
       do is = 1, nspec
          !qmheat(:,:,is,1)=qmheat(:,:,is,1) * spec(is)%dens*spec(is)%temp
          !call get_volume_average (qmheat(:,:,is,1), mheat_fluxes(is))

          !call get_surf_average (qmheat(:,:,is,1), x_qmflux(:,is))

          !qmheat(:,:,is,2)=qmheat(:,:,is,2) * spec(is)%dens*spec(is)%temp
          !call get_volume_average (qmheat(:,:,is,2), mheat_par(is))

          !qmheat(:,:,is,3)=qmheat(:,:,is,3) * spec(is)%dens*spec(is)%temp
          !call get_volume_average (qmheat(:,:,is,3), mheat_perp(is))
          
          !pmflux(:,:,is)=pmflux(:,:,is) * spec(is)%dens
          !call get_volume_average (pmflux(:,:,is), mpart_fluxes(is))

          !vmflux(:,:,is)=vmflux(:,:,is) * spec(is)%dens*spec(is)%mass*spec(is)%stm
          !call get_volume_average (vmflux(:,:,is), mmom_fluxes(is))
       end do
    end if
    if (fbpar > epsilon(0.0)) then
       do is = 1, nspec
          !qbheat(:,:,is,1)=qbheat(:,:,is,1) * spec(is)%dens*spec(is)%temp
          !call get_volume_average (qbheat(:,:,is,1), bheat_fluxes(is))

          !qbheat(:,:,is,2)=qbheat(:,:,is,2) * spec(is)%dens*spec(is)%temp
          !call get_volume_average (qbheat(:,:,is,2), bheat_par(is))

          !qbheat(:,:,is,3)=qbheat(:,:,is,3) * spec(is)%dens*spec(is)%temp
          !call get_volume_average (qbheat(:,:,is,3), bheat_perp(is))
          
          !pbflux(:,:,is)=pbflux(:,:,is) * spec(is)%dens
          !call get_volume_average (pbflux(:,:,is), bpart_fluxes(is))

          !vbflux(:,:,is)=vbflux(:,:,is) * spec(is)%dens*spec(is)%mass*spec(is)%stm
          !call get_volume_average (vbflux(:,:,is), bmom_fluxes(is))
       end do
    end if
  end subroutine write_fluxes

  !> Calculate estimates of the heat and particles fluxes using
  !! gamma / k^2 estimate of the diffusivity
  subroutine write_diffusive_estimates(gnostics)
    use diagnostics_write_omega, only: omega_average
    use dist_fn_arrays, only: kperp2
    use fields_parallelization, only: field_k_local
    use species, only: spec
    use kt_grids, only: aky, akx
    !use geometry, only: surfarea, dvdrhon, grho
    use theta_grid, only: grho
    use diagnostics_create_and_write
    use volume_averages
    implicit none
    type(diagnostics_type), intent(in) :: gnostics
    real, dimension(ntheta0, naky) :: diffusivity_by_k
    real :: diffusivity
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
    real :: k2
    diffusivity_by_k = 0.0
    heat_flux = 0.0
    particle_flux = 0.0
    momentum_flux = 0.0
    heat_flux_max = 0.0
    !grho = 1.0
    !grho = surfarea/dvdrhon
    !write (*,*) 'grho', grho
    do ik = 1,naky
      do it = 1,ntheta0
        if (.not. gnostics%distributed .or. field_k_local(it,ik)) then
          !k2 = (aky(ik)**2.0 + akx(it)**2.0)**0.5
          !if (k2.eq.epsilon(0.0)) cycle
          !if (akx(it) .eq. 0.0) cycle
          !if (aky(ik) .eq. 0.0) cycle
          if (kperp2(gnostics%igomega,it,ik).eq.0.0) cycle
          diffusivity_by_k(it,ik) = &
            max(aimag(omega_average(it,ik)),0.0)/kperp2(gnostics%igomega, it, ik)*2.0
            !max(aimag(omega_average(it,ik)),0.0)/aky(ik)**2.0/2.0**0.5
            !max(aimag(omega_average(it,ik)),0.0)/akx(ik)**2.0/2.0**0.5
            !max(aimag(omega_average(it,ik)),0.0)/k2
            !aimag(omega_average(it,ik))
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
!
           !heat_flux(it,ik,is) = max(aimag(omega_average(it,ik)),0.0)

           !heat_flux_max = max(heat_flux(it,ik,is), heat_flux_max)
             
          end do
        end if
      end do
    end do
    !heat_flux = 0.0
    call write_standard_flux_properties(gnostics, &
      'heat_flux_diff',  'Diffusive estimate of turbulent flux of heat', &
      'Q_gB = ', heat_flux, gnostics%distributed)
    
    call write_standard_flux_properties(gnostics, &
      'part_flux_diff',  'Diffusive estimate of turbulent flux of particles', &
      'n_r? ', particle_flux, gnostics%distributed)


    !call average_kx(heat_flux, heat_flux_by_ky, gnostics%distributed)
    !call average_kx(particle_flux, particle_flux_by_ky, gnostics%distributed)
    heat_flux_by_ky = maxval(heat_flux, 1)
    particle_flux_by_ky = maxval(particle_flux, 1)

    heat_flux_max_by_spec = maxval(heat_flux_by_ky, 1)
    particle_flux_max_by_spec = maxval(particle_flux_by_ky, 1)

    heat_flux_max = sum(heat_flux_max_by_spec)
    particle_flux_max = sum(particle_flux_max_by_spec)
    
    call create_and_write_variable(gnostics, gnostics%rtype, "es_heat_flux_diff_max", "st", &
      " A turbulent estimate of the heat flux, as a function of species and  time", &
      "Q_gB", heat_flux_max_by_spec)
    call create_and_write_variable(gnostics, gnostics%rtype, "heat_flux_diff_max", "t", &
      " A turbulent estimate of the heat flux, as a function of  time", &
      "Q_gB", heat_flux_max)
    call create_and_write_variable(gnostics, gnostics%rtype, "es_particle_flux_diff_max", "st", &
      " A turbulent estimate of the particle flux, as a function of species and  time", &
      "Q_gB", particle_flux_max_by_spec)
    call create_and_write_variable(gnostics, gnostics%rtype, "particle_flux_diff_max", "t", &
      " A turbulent estimate of the particle flux, as a function of  time", &
      "Q_gB", particle_flux_max)

  end subroutine write_diffusive_estimates

  subroutine write_standard_flux_properties(gnostics, flux_name, flux_description, &
    flux_units, flux_value, distributed)
    use kt_grids, only: ntheta0, naky
    use species, only: nspec
    use diagnostics_create_and_write, only: create_and_write_variable
    use volume_averages
    use fields_parallelization, only: field_k_local
    type(diagnostics_type), intent(in) :: gnostics
    character(*), intent(in) :: flux_name, flux_description, flux_units
    real, dimension(ntheta0,naky,nspec), intent(in) :: flux_value
    logical, intent(in) :: distributed
    real, dimension(ntheta0, naky) :: total_flux_by_mode
    real, dimension(naky) :: total_flux_by_ky
    real, dimension(ntheta0) :: total_flux_by_kx
    real, dimension(nspec) :: flux_by_species
    integer :: is, ik, it
    !logical, external :: k_local

    !call average_theta(flux_value, flux_value, flux_by_mode, .true.)
    !!call create_and_write_flux(gnostics%sfile, flux_name, flux_description, flux_units, flux_value)
    call average_all(flux_value, flux_by_species, distributed) 
    call create_and_write_variable(gnostics, gnostics%rtype, "es_"//flux_name, "st", &
      flux_description//" averaged over kx and ky, as a function of species and time", &
      flux_units, flux_by_species)

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

    if (gnostics%write_fluxes_by_mode) then 
      call create_and_write_flux_by_mode(gnostics, flux_name, flux_description, flux_units, &
        flux_value, total_flux_by_mode, distributed)
    end if



    call average_kx(total_flux_by_mode, total_flux_by_ky, distributed)
    call create_and_write_variable(gnostics, gnostics%rtype, "total_"//flux_name//"_by_ky", "Yt", &
      flux_description//" summed over species and averaged over kx, as a function of ky and time", &
      flux_units, total_flux_by_ky)

    call average_ky(total_flux_by_mode, total_flux_by_kx, distributed)
    call create_and_write_variable(gnostics, gnostics%rtype, "total_"//flux_name//"_by_kx", "Xt", &
      flux_description//" summed over species and averaged over ky, as a function of kx and time", &
      flux_units, total_flux_by_kx)

    call create_and_write_variable(gnostics, gnostics%rtype, flux_name//"_tot", "t", &
      flux_description//" summed over species and averaged over kx and ky, as a function of time", &
      flux_units, sum(total_flux_by_kx))

  end subroutine write_standard_flux_properties

  subroutine create_and_write_flux_by_mode(gnostics, flux_name, flux_description, flux_units, &
      flux_value, total_flux_by_mode, distributed)
    use fields_parallelization, only: field_k_local
   type(diagnostics_type), intent(in) :: gnostics
   character(*), intent(in) :: flux_name
   character(*), intent(in) :: flux_description
   character(*), intent(in) :: flux_units
   character(len=len(flux_name)+14) :: total_flux_by_mode_name
   character(len=len(flux_name)+3+8) :: flux_by_mode_name
   real, intent(in), dimension(ntheta0, naky, nspec) :: flux_value
   real, dimension(ntheta0, naky), intent(in) :: total_flux_by_mode
   logical, intent(in) :: distributed
   !logical, external :: k_local

   real, dimension(1, 1, nspec) :: dummyc
   real, dimension(1,1) :: dummyr
   integer :: it,ik

   flux_by_mode_name = "es_"//flux_name//"_by_mode"
   total_flux_by_mode_name =  "total_"//flux_name//"_by_mode" 

   !flux_t_name = flux_name//"_t"


 
   if (gnostics%create) then 
    call create_variable(gnostics%sfile, gnostics%rtype, flux_by_mode_name, "XYst", &
       flux_description//" as a function of species, kx and ky", flux_units)
    call create_variable(gnostics%sfile, gnostics%rtype, total_flux_by_mode_name, "XYt", &
       flux_description//" summed over species, as a function of kx and ky" , flux_units)


   end if

   if (gnostics%create .or. .not. gnostics%wryte) return
  !if (write_flux_by_time .and. .not. variable_exists(gnostics%sfile, flux_t_name)) then 
     !call create_variable(gnostics%sfile, gnostics%rtype, flux_t_name, "rzXYt", &
       !flux_description//" as a function of time", flux_units)
   !end if
   
!>  Given that the fluxes are calculated from the fields,
!!   they will in general be distributed across processes
!!   in the same way as the fields. 
!!  In this case distributed is true
   if (.not. distributed) then

     call write_variable(gnostics%sfile, flux_by_mode_name, flux_value)
     call write_variable(gnostics%sfile, total_flux_by_mode_name, total_flux_by_mode)
     !if (write_flux_by_time) call write_variable(gnostics%sfile, flux_t_name, flux_value)

   else

     call set_count(gnostics%sfile, flux_by_mode_name, "X", 1)
     call set_count(gnostics%sfile, flux_by_mode_name, "Y", 1)
     call set_count(gnostics%sfile, total_flux_by_mode_name, "X", 1)
     call set_count(gnostics%sfile, total_flux_by_mode_name, "Y", 1)
     !if (write_flux_by_time) then
       !call set_count(gnostics%sfile, flux_t_name, "X", 1)
       !call set_count(gnostics%sfile, flux_t_name, "Y", 1)
     !end if

     ! For some reason every process has to make at least
     ! one write to a variable with an infinite dimension.
     ! Here we make some dummy writes to satisfy that
     call write_variable_with_offset(gnostics%sfile, flux_by_mode_name, dummyc)
     call write_variable_with_offset(gnostics%sfile, total_flux_by_mode_name, dummyr)

     do ik = 1,naky
       do it = 1,ntheta0
         if (field_k_local(it,ik)) then

           call set_start(gnostics%sfile, flux_by_mode_name, "X", it)
           call set_start(gnostics%sfile, flux_by_mode_name, "Y", ik)
           call write_variable_with_offset(gnostics%sfile, flux_by_mode_name, flux_value)

           call set_start(gnostics%sfile, total_flux_by_mode_name, "X", it)
           call set_start(gnostics%sfile, total_flux_by_mode_name, "Y", ik)
           call write_variable_with_offset(gnostics%sfile, total_flux_by_mode_name, total_flux_by_mode)

           !if (write_flux_by_time) then
             !call set_start(gnostics%sfile, flux_t_name, "X", it)
             !call set_start(gnostics%sfile, flux_t_name, "Y", ik)
             !call write_variable(gnostics%sfile, flux_t_name, flux_value)
           !end if
         end if
       end do
     end do
   end if
  end subroutine create_and_write_flux_by_mode

end module diagnostics_write_fluxes
