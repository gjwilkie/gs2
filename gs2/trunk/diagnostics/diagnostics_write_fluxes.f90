!> Module which contains functions for calculating
!! and writing out the fluxes of heat and momentum etc.
module diagnostics_write_fluxes

  use kt_grids, only: naky, ntheta0
  use species, only: nspec
  use theta_grid, only: ntgrid
  use simpledataio
  use simpledataio_write
  

implicit none

!> Allocate arrays
public :: init_diagnostics_write_fluxes

public :: write_fluxes
  logical :: fluxes_local = .false.

private
  real, dimension (:,:,:,:), allocatable ::  qheat, qmheat, qbheat
  real, dimension (:,:,:), allocatable ::  pflux,  vflux, vflux_par, vflux_perp
  real, dimension (:,:,:), allocatable ::  pflux_tormom
  real, dimension (:,:,:), allocatable :: vflux0, vflux1  ! low flow correction to turbulent momentum flux
  real, dimension (:,:,:), allocatable :: pmflux, vmflux
  real, dimension (:,:,:), allocatable :: pbflux, vbflux
  real, dimension (:,:,:), allocatable :: exchange
  integer, parameter :: REAL_TYPE = SDATIO_DOUBLE


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

  subroutine write_fluxes(netcdf_file, istep)
    use dist_fn, only: flux
    use dist_fn_arrays, only: g_adjust, gnew
    use species, only: nspec, spec
    use fields_arrays, only: phinew, bparnew, aparnew
    use run_parameters, only: fphi, fapar, fbpar
    type(sdatio_file), intent(in) :: netcdf_file
    integer, intent(in) :: istep
    integer :: is
    !if (istep > 0) then
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
       end do
    end if
    call write_standard_flux_properties(netcdf_file, &
      'heat_flux',  'Turbulent flux of heat', 'Q_gB = ', qheat(:,:,:,1), .true.)
    call write_standard_flux_properties(netcdf_file, &
      'heat_par',  'Turbulent flux of parallel heat', 'Q_gB = ', qheat(:,:,:,2), .true.)
    call write_standard_flux_properties(netcdf_file, &
      'heat_perp',  'Turbulent flux of perpendicular heat', 'Q_gB = ', qheat(:,:,:,3), .true.)
    call write_standard_flux_properties(netcdf_file, &
      'part_flux',  'Turbulent flux of particles', 'n_r? ', pflux, .true.)
          !call get_volume_average (qheat(:,:,is,1), heat_fluxes(is))

          !call get_volume_average (qheat(:,:,is,2), heat_par(is))

          !call get_volume_average (qheat(:,:,is,3), heat_perp(is))
          
          !call get_volume_average (pflux(:,:,is), part_fluxes(is))

!pflux_tormom(:,:,is) = pflux_tormom(:,:,is) * spec(is)%dens  
      !call get_volume_average (pflux_tormom(:,:,is), part_tormom_fluxes(is))

          !vflux(:,:,is) = vflux(:,:,is) * spec(is)%dens*sqrt(spec(is)%mass*spec(is)%temp)
          !call get_volume_average (vflux(:,:,is), mom_fluxes(is))

          !vflux_par(:,:,is) = vflux_par(:,:,is) * spec(is)%dens*sqrt(spec(is)%mass*spec(is)%temp)
          !call get_volume_average (vflux_par(:,:,is), parmom_fluxes(is))

          !vflux_perp(:,:,is) = vflux_perp(:,:,is) * spec(is)%dens*spec(is)%mass*spec(is)%stm
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
  subroutine write_standard_flux_properties(netcdf_file, flux_name, flux_description, &
    flux_units, flux_value, distributed)
    use kt_grids, only: ntheta0, naky
    use species, only: nspec
    use diagnostics_create_and_write, only: create_and_write_variable
    use volume_averages
    use fields_parallelization, only: field_k_local
    type(sdatio_file), intent(in) :: netcdf_file
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
    !!call create_and_write_flux(netcdf_file, flux_name, flux_description, flux_units, flux_value)
    call average_all(flux_value, flux_by_species, distributed) 
    call create_and_write_variable(netcdf_file, REAL_TYPE, "es_"//flux_name, "st", &
      flux_description//" averaged over kx and ky, as a function of species and time", &
      flux_units, flux_by_species)

    call create_and_write_flux_by_mode(netcdf_file, flux_name, flux_description, flux_units, &
      flux_value, total_flux_by_mode, distributed)


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


    call average_kx(total_flux_by_mode, total_flux_by_ky, .true.)
    call create_and_write_variable(netcdf_file, REAL_TYPE, "total_"//flux_name//"_by_ky", "yt", &
      flux_description//" summed over species and averaged over kx, as a function of ky and time", &
      flux_units, total_flux_by_ky)

    call average_ky(total_flux_by_mode, total_flux_by_kx, .true.)
    call create_and_write_variable(netcdf_file, REAL_TYPE, "total_"//flux_name//"_by_kx", "xt", &
      flux_description//" summed over species and averaged over ky, as a function of kx and time", &
      flux_units, total_flux_by_kx)

    call create_and_write_variable(netcdf_file, REAL_TYPE, flux_name//"_tot", "t", &
      flux_description//" summed over species and averaged over kx and ky, as a function of time", &
      flux_units, sum(total_flux_by_kx))

  end subroutine write_standard_flux_properties

  subroutine create_and_write_flux_by_mode(sfile, flux_name, flux_description, flux_units, &
      flux_value, total_flux_by_mode, write_flux_by_time)
    use fields_parallelization, only: field_k_local
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: flux_name
   character(*), intent(in) :: flux_description
   character(*), intent(in) :: flux_units
   character(len=len(flux_name)+14) :: total_flux_by_mode_name
   character(len=len(flux_name)+3+8) :: flux_by_mode_name
   real, intent(in), dimension(ntheta0, naky, nspec) :: flux_value
   real, dimension(ntheta0, naky), intent(in) :: total_flux_by_mode
   logical, intent(in) :: write_flux_by_time
   !logical, external :: k_local

   real, dimension(1, 1, nspec) :: dummyc
   real, dimension(1,1) :: dummyr
   integer :: it,ik

   flux_by_mode_name = "es_"//flux_name//"_by_mode"
   total_flux_by_mode_name =  "total_"//flux_name//"_by_mode" 

   !flux_t_name = flux_name//"_t"


	 
	 if (.not. variable_exists(sfile, flux_by_mode_name)) then 
	   call create_variable(sfile, SDATIO_DOUBLE, flux_by_mode_name, "xyst", &
       flux_description//" as a function of species, kx and ky", flux_units)
	 end if
	 if (.not. variable_exists(sfile, total_flux_by_mode_name)) then 
	   call create_variable(sfile, SDATIO_DOUBLE, total_flux_by_mode_name, "xyt", &
       flux_description//" summed over species, as a function of kx and ky" , flux_units)
	 end if
	 !if (write_flux_by_time .and. .not. variable_exists(sfile, flux_t_name)) then 
	   !call create_variable(sfile, SDATIO_DOUBLE, flux_t_name, "rzxyt", &
       !flux_description//" as a function of time", flux_units)
	 !end if
	 
   if (fluxes_local) then

     call write_variable(sfile, flux_by_mode_name, flux_value)
     call write_variable(sfile, total_flux_by_mode_name, total_flux_by_mode)
     !if (write_flux_by_time) call write_variable(sfile, flux_t_name, flux_value)

   else

     call set_count(sfile, flux_by_mode_name, "x", 1)
     call set_count(sfile, flux_by_mode_name, "y", 1)
     call set_count(sfile, total_flux_by_mode_name, "x", 1)
     call set_count(sfile, total_flux_by_mode_name, "y", 1)
     !if (write_flux_by_time) then
       !call set_count(sfile, flux_t_name, "x", 1)
       !call set_count(sfile, flux_t_name, "y", 1)
     !end if

     ! For some reason every process has to make at least
     ! one write to a variable with an infinite dimension.
     ! Here we make some dummy writes to satisfy that
     call write_variable(sfile, flux_by_mode_name, dummyc)
     call write_variable(sfile, total_flux_by_mode_name, dummyr)

     do ik = 1,naky
       do it = 1,ntheta0
         if (field_k_local(it,ik)) then

           call set_start(sfile, flux_by_mode_name, "x", it)
           call set_start(sfile, flux_by_mode_name, "y", ik)
           call write_variable(sfile, flux_by_mode_name, flux_value)

           call set_start(sfile, total_flux_by_mode_name, "x", it)
           call set_start(sfile, total_flux_by_mode_name, "y", ik)
           call write_variable(sfile, total_flux_by_mode_name, total_flux_by_mode)

           !if (write_flux_by_time) then
             !call set_start(sfile, flux_t_name, "x", it)
             !call set_start(sfile, flux_t_name, "y", ik)
             !call write_variable(sfile, flux_t_name, flux_value)
           !end if
         end if
       end do
     end do
   end if
  end subroutine create_and_write_flux_by_mode

end module diagnostics_write_fluxes
