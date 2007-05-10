module gs2_diagnostics
  implicit none

  public :: init_gs2_diagnostics
  public :: finish_gs2_diagnostics
  public :: loop_diagnostics

  ! knobs

  real :: omegatol, omegatinst
  logical :: print_line, print_old_units, print_flux_line
  logical :: write_line, write_flux_line, write_phi, write_apar, write_aperp
  logical :: write_omega, write_omavg, write_ascii
  logical :: write_qheat, write_pflux, write_vflux
  logical :: write_qmheat, write_pmflux, write_vmflux
  logical :: write_qbheat, write_pbflux, write_vbflux
  logical :: write_dmix, write_kperpnorm, write_phitot, write_epartot
  logical :: write_eigenfunc, write_final_fields, write_final_antot
  logical :: write_final_moments, write_avg_moments
  logical :: write_fcheck, write_final_epar
  logical :: write_intcheck, write_vortcheck, write_fieldcheck
  logical :: write_fieldline_avg_phi
  logical :: write_neoclassical_flux, write_nl_flux
  logical :: exit_when_converged
  logical :: dump_neoclassical_flux, dump_check1
  logical :: dump_fields_periodically
  logical :: dump_final_xfields
  logical :: use_shmem_for_xfields
  logical :: save_for_restart

  integer :: nperiod_output
  integer :: nwrite, igomega
  integer :: navg

  character(300) :: restart_file

  ! internal
  logical :: write_any, write_any_fluxes, dump_any
  logical, private :: initialized = .false.

  integer :: out_unit
  integer :: dump_neoclassical_flux_unit, dump_check1_unit

  complex, dimension (:,:,:), allocatable :: omegahist
  ! (navg,ntheta0,naky)

  real, dimension (:,:,:), allocatable :: pflux, qheat, qheat_par, qheat_perp, vflux
  real, dimension (:,:,:), allocatable :: pmflux, qmheat, qmheat_par, qmheat_perp, vmflux
  real, dimension (:,:,:), allocatable :: pbflux, qbheat, qbheat_par, qbheat_perp, vbflux
  ! (ntheta0,naky,nspec)

  integer :: ntg_out

contains

  subroutine init_gs2_diagnostics
    use theta_grid, only: init_theta_grid
    use kt_grids, only: init_kt_grids
    use run_parameters, only: init_run_parameters
    use species, only: init_species
    use dist_fn, only: init_dist_fn
    use dist_fn, only: init_intcheck, init_vortcheck, init_fieldcheck
    use gs2_io, only: init_gs2_io
    use mp, only: broadcast
    implicit none
    real :: denom
    integer :: ik, it

    if (initialized) return
    initialized = .true.

    call init_theta_grid
    call init_kt_grids
    call init_run_parameters
    call init_species
    call init_dist_fn
    call init_gs2_io 

    call real_init
    call broadcast (nwrite)
    call broadcast (write_any)
    call broadcast (write_any_fluxes)
    call broadcast (write_neoclassical_flux)
    call broadcast (write_nl_flux)
    call broadcast (write_fcheck)
    call broadcast (dump_any)
    call broadcast (dump_neoclassical_flux)
    call broadcast (dump_check1)
    call broadcast (dump_fields_periodically)
    call broadcast (dump_final_xfields)
    call broadcast (use_shmem_for_xfields)
    call broadcast (save_for_restart)
    call broadcast (restart_file)

    call broadcast (write_intcheck)
    call broadcast (write_vortcheck)
    call broadcast (write_fieldcheck)
    if (write_intcheck) call init_intcheck
    if (write_vortcheck) call init_vortcheck
    if (write_fieldcheck) call init_fieldcheck
    call broadcast (ntg_out)
    

  end subroutine init_gs2_diagnostics

  subroutine real_init
    use file_utils, only: open_output_file, get_unused_unit
    use kt_grids, only: naky, ntheta0
    use species, only: nspec
    use mp, only: proc0
    implicit none
    character(20) :: datestamp, timestamp, zone

    call read_parameters
    if (proc0) then
       if (write_ascii) call open_output_file (out_unit, ".out")
       
       if (dump_neoclassical_flux) then
          call get_unused_unit (dump_neoclassical_flux_unit)
          open (unit=dump_neoclassical_flux_unit, file="dump.neoflux", &
               status="unknown")
       end if
       
       if (dump_check1) then
          call get_unused_unit (dump_check1_unit)
          open (unit=dump_check1_unit, file="dump.check1", status="unknown")
       end if
       
       if (write_ascii) then
          write (unit=out_unit, fmt="('gs2')")
          datestamp(:) = ' '
          timestamp(:) = ' '
          zone(:) = ' '
          call date_and_time (datestamp, timestamp, zone)
          write (unit=out_unit, fmt="('Date: ',a,' Time: ',a,1x,a)") &
               trim(datestamp), trim(timestamp), trim(zone)
       end if
       
       allocate (omegahist(0:navg-1,ntheta0,naky))
       omegahist = 0.0
       
       allocate (pflux (ntheta0,naky,nspec))
       allocate (qheat (ntheta0,naky,nspec))
       allocate (qheat_par  (ntheta0,naky,nspec))
       allocate (qheat_perp (ntheta0,naky,nspec))
       allocate (vflux (ntheta0,naky,nspec))
       allocate (pmflux(ntheta0,naky,nspec))
       allocate (qmheat(ntheta0,naky,nspec))
       allocate (qmheat_par (ntheta0,naky,nspec))
       allocate (qmheat_perp(ntheta0,naky,nspec))
       allocate (vmflux(ntheta0,naky,nspec))
       allocate (pbflux(ntheta0,naky,nspec))
       allocate (qbheat(ntheta0,naky,nspec))
       allocate (qbheat_par (ntheta0,naky,nspec))
       allocate (qbheat_perp(ntheta0,naky,nspec))
       allocate (vbflux(ntheta0,naky,nspec))
    end if

  end subroutine real_init

  subroutine read_parameters
    use file_utils, only: input_unit, run_name
    use theta_grid, only: nperiod, ntheta
    use dist_fn, only: nperiod_guard
    use kt_grids, only: box
    use mp, only: proc0
    implicit none
    namelist /gs2_diagnostics_knobs/ &
         print_line, print_old_units, print_flux_line, &
         write_line, write_flux_line, write_phi, write_apar, write_aperp, &
         write_omega, write_omavg, write_ascii, &
         write_qheat, write_pflux, write_vflux, &
         write_qmheat, write_pmflux, write_vmflux, &
         write_qbheat, write_pbflux, write_vbflux, &
         write_dmix, write_kperpnorm, write_phitot, write_epartot, &
         write_eigenfunc, write_final_fields, write_final_antot, &
         write_fcheck, write_final_epar, write_final_moments, &
         write_intcheck, write_vortcheck, write_fieldcheck, &
         write_fieldline_avg_phi, write_neoclassical_flux, write_nl_flux, &
         nwrite, navg, omegatol, omegatinst, igomega, &
         exit_when_converged, write_avg_moments, &
         dump_neoclassical_flux, dump_check1, &
         dump_fields_periodically, &
         dump_final_xfields, use_shmem_for_xfields, &
         nperiod_output, &
         save_for_restart, restart_file

    if (proc0) then
       print_line = .true.
       print_old_units = .true.
       print_flux_line = .false.
       write_line = .true.
       write_flux_line = .true.
       write_phi = .true.
       write_apar = .true.
       write_aperp = .true.
       write_omega = .true.
       write_ascii = .true.
       write_omavg = .true.
       write_qheat = .true.
       write_pflux = .true.
       write_vflux = .true.
       write_qmheat = .true.
       write_pmflux = .true.
       write_vmflux = .true.
       write_qbheat = .true.
       write_pbflux = .true.
       write_vbflux = .true.
       write_dmix = .true.
       write_kperpnorm = .true.
       write_phitot = .true.
       write_epartot = .false.
       write_fieldline_avg_phi = .false.
       write_neoclassical_flux = .false.
       write_nl_flux = .false.
       write_eigenfunc = .true.
       write_final_moments = .false.
       write_avg_moments = .false.
       write_final_fields = .false.
       write_final_antot = .false.
       write_final_epar = .false.
       write_fcheck = .false.
       write_intcheck = .false.
       write_vortcheck = .false.
       write_fieldcheck = .false.
       nwrite = 100
       navg = 100
       omegatol = 1e-3
       omegatinst = 1.0
       igomega = 0
       exit_when_converged = .true.
       dump_neoclassical_flux = .false.
       dump_check1 = .false.
       dump_fields_periodically = .false.
       dump_final_xfields = .false.
       use_shmem_for_xfields = .true.
       nperiod_output = nperiod - nperiod_guard
       save_for_restart = .false.
       restart_file = trim(run_name)//".nc"
       read (unit=input_unit("gs2_diagnostics_knobs"), nml=gs2_diagnostics_knobs)

       write_avg_moments = write_avg_moments .and. box

       write_any = write_line .or. write_phi       .or. write_apar &
            .or. write_aperp  .or. write_omega     .or. write_omavg &
            .or. write_qheat  .or. write_pflux     .or. write_vflux &
            .or. write_qmheat .or. write_pmflux    .or. write_vmflux &
            .or. write_qbheat .or. write_pbflux    .or. write_vbflux &
            .or. write_dmix   .or. write_kperpnorm .or. write_phitot &
            .or. write_fieldline_avg_phi           .or. write_neoclassical_flux &
            .or. write_flux_line                   .or. write_nl_flux
       write_any_fluxes =  &
            write_qheat  .or. write_pflux     .or. write_vflux &
            .or. write_qmheat .or. write_pmflux    .or. write_vmflux &
            .or. write_qbheat .or. write_pbflux    .or. write_vbflux &
            .or. write_flux_line                   .or. print_flux_line &
            .or. write_nl_flux
       dump_any = dump_neoclassical_flux .or. dump_check1  .or. dump_fields_periodically

       nperiod_output = min(nperiod,nperiod_output)
       ntg_out = ntheta/2 + (nperiod_output-1)*ntheta
    end if
  end subroutine read_parameters

  subroutine finish_gs2_diagnostics (istep)
    use file_utils, only: open_output_file, close_output_file, get_unused_unit
    use mp, only: proc0, broadcast
    use species, only: nspec
    use run_parameters, only: tnorm
    use theta_grid, only: ntgrid, theta, delthet, jacob
    use kt_grids, only: naky, ntheta0, theta0, nx, aky_out, akx_out
    use le_grids, only: nlambda, negrid, fcheck, al
    use fields_arrays, only: phi, apar, aperp, phinew, aparnew
    use dist_fn, only: getan, get_epar, getmoms
!    use dist_fn, only: write_g
    use dist_fn, only: finish_intcheck, finish_vortcheck, finish_fieldcheck
    use dist_fn_arrays, only: g, gnew
    use gs2_layouts, only: xxf_lo
    use gs2_transforms, only: init_x_transform, transform_x
    use gs2_save, only: gs2_save_for_restart
    use constants
    use gs2_time, only: stime, simdt
    use gs2_io, only: nc_eigenfunc, nc_final_fields, nc_final_epar, nc_final_an
    use gs2_io, only: nc_final_moments, nc_finish
    implicit none
    integer, intent (in) :: istep
    integer :: ig, ik, it, il, is, unit
    complex, dimension (nlambda,ntheta0,naky,nspec) :: fcheck_f
    complex, dimension (:,:,:), allocatable :: xphi, xapar, xaperp
    complex, dimension (-ntgrid:ntgrid,ntheta0,naky) :: antot, antota, antotp, epar
    complex, dimension (-ntgrid:ntgrid,ntheta0,naky,nspec) :: ntot, density, &
         upar, tpar, tperp
    real, dimension (:), allocatable :: dl_over_b
    complex, dimension (ntheta0, naky) :: phi0
    real :: simtime, deltsave
    real, dimension (nspec) :: weights
    integer :: istatus

!    call write_g

    phi0 = 1.

    if (proc0) then
       if (write_ascii) call close_output_file (out_unit)
       if (dump_neoclassical_flux) close (dump_neoclassical_flux_unit)
       if (dump_check1) close (dump_check1_unit)

       if (write_eigenfunc) then
          if (write_ascii) call open_output_file (unit, ".eigenfunc")
          phi0 = phi(0,:,:)

          where (abs(phi0) < 10.0*epsilon(0.0)) 
             phi0 = phi(1,:,:)/(theta(1)-theta(0))
          end where

          where (abs(phi0) < 10.0*epsilon(0.0)) 
             phi0 = 1.0
          end where

          if (write_ascii) then
             do ik = 1, naky
                do it = 1, ntheta0
                   do ig = -ntg_out, ntg_out
                      write (unit, "(9(1x,e12.5))") &
                           theta(ig), theta0(it,ik), aky_out(ik), &
                             phi(ig,it,ik)/phi0(it,ik), &
                            apar(ig,it,ik)/phi0(it,ik), &
                           aperp(ig,it,ik)/phi0(it,ik)
                   end do
                   write (unit, "()")
                end do
             end do
             call close_output_file (unit)
          end if
          call nc_eigenfunc (phi0)
       end if

       if (write_final_fields) then
          if (write_ascii) then
             call open_output_file (unit, ".fields")
             do ik = 1, naky
                do it = 1, ntheta0
                   do ig = -ntg_out, ntg_out
                      write (unit, "(15(1x,e12.5))") &
!                           theta(ig), theta0(it,ik), aky_out(ik), &
                           theta(ig), aky_out(ik), akx_out(it), &
                             phi(ig,it,ik), &
                            apar(ig,it,ik), &
                           aperp(ig,it,ik), &
                           theta(ig) - theta0(it,ik), &
                           cabs(phi(ig,it,ik))
                   end do
                   write (unit, "()")
                end do
             end do
             call close_output_file (unit)
          end if
          call nc_final_fields
       end if
       
       if (write_final_epar) then
          call get_epar (phi, apar, aparnew, epar)
          if (write_ascii) then
             call open_output_file (unit, ".epar")
             do ik = 1, naky
                do it = 1, ntheta0
                   do ig = -ntg_out, ntg_out-1
                      write (unit, "(6(1x,e12.5))") &
!                           theta(ig), theta0(it,ik), aky_out(ik), &
                           theta(ig), aky_out(ik), akx_out(it), &
                           epar(ig,it,ik), &
                           theta(ig) - theta0(it,ik)
                   end do
                   write (unit, "()")
                end do
             end do
             call close_output_file (unit)
          end if
          call nc_final_epar (epar  )
       end if

    end if

    call broadcast (write_final_moments)
    if (write_final_moments) then

       call getmoms (phinew, ntot, density, upar, tpar, tperp)

       if (proc0) then
          if (write_ascii) then
             call open_output_file (unit, ".moments")
             do is  = 1, nspec
                do ik = 1, naky
                   do it = 1, ntheta0
                      do ig = -ntg_out, ntg_out
                         write (unit, "(15(1x,e12.5))") &
                              theta(ig), aky_out(ik), akx_out(it), &
                              ntot(ig,it,ik,is), &
                              density(ig,it,ik,is), &
                              upar(ig,it,ik,is), &
                              tpar(ig,it,ik,is), &
                              tperp(ig,it,ik,is), &
                              theta(ig) - theta0(it,ik), &
                              real(is)
                      end do
                      write (unit, "()")
                   end do
                end do
             end do
             call close_output_file (unit)          
          end if
          call nc_final_moments (ntot, density, upar, tpar, tperp)

          call open_output_file (unit, ".amoments")

          allocate (dl_over_b(-ntg_out:ntg_out))

          dl_over_b = delthet(-ntg_out:ntg_out)*jacob(-ntg_out:ntg_out)
          dl_over_b = dl_over_b / sum(dl_over_b)

          do is  = 1, nspec
             do it = 2, ntheta0
                write (unit, "(14(1x,e10.3))") real(is), akx_out(it), &
                     sum(phinew(-ntg_out:ntg_out,it,1)*dl_over_b), &
                     sum(ntot(-ntg_out:ntg_out,it,1,is)*dl_over_b), &
                     sum(density(-ntg_out:ntg_out,it,1,is)*dl_over_b), &
                     sum(upar(-ntg_out:ntg_out,it,1,is)*dl_over_b), &
                     sum(tpar(-ntg_out:ntg_out,it,1,is)*dl_over_b), &
                     sum(tperp(-ntg_out:ntg_out,it,1,is)*dl_over_b)
             end do
          end do
          deallocate (dl_over_b)
          call close_output_file (unit)          

       end if
    end if

    call broadcast (write_final_antot)
    if (write_final_antot) then
       call getan (antot, antota, antotp)
       if (proc0) then
          if (write_ascii) then
             call open_output_file (unit, ".antot")
             do ik = 1, naky
                do it = 1, ntheta0
                   do ig = -ntg_out, ntg_out
                      write (unit, "(10(1x,e12.5))") &
                           theta(ig), theta0(it,ik), aky_out(ik), &
                           antot(ig,it,ik), &
                           antota(ig,it,ik), &
                           antotp(ig,it,ik), &
                           theta(ig) - theta0(it,ik)
                   end do
                   write (unit, "()")
                end do
             end do
             call close_output_file (unit)
          end if
          call nc_final_an (antot, antota, antotp)
       end if
    end if

    if (dump_final_xfields) then
       call init_x_transform (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx)

       allocate (xphi(xxf_lo%nx,naky,-ntgrid:ntgrid))
       allocate (xapar(xxf_lo%nx,naky,-ntgrid:ntgrid))
       allocate (xaperp(xxf_lo%nx,naky,-ntgrid:ntgrid))
       call transform_x (phi, xphi)
       call transform_x (apar, xapar)
       call transform_x (aperp, xaperp)
       if (proc0) then
          call get_unused_unit (unit)
          open (unit=unit, file="dump.xfields", status="unknown")
          do ik = 1, naky
             do ig = -ntg_out, ntg_out
                do it = 1, xxf_lo%nx
                   write (unit, "(15(1x,e12.5))") &
                        2.0*pi/akx_out(2)*real(it-1-xxf_lo%nx/2)/real(xxf_lo%nx), &
                        theta(ig), aky_out(ik), &
                        xphi(it,ik,ig), &
                        xapar(it,ik,ig), &
                        xaperp(it,ik,ig), &
                        abs(xphi(it,ik,ig)), &
                        abs(xapar(it,ik,ig)), &
                        abs(xaperp(it,ik,ig))
                end do
                write (unit, "()")
             end do
          end do
          close (unit=unit)
       end if
       deallocate (xphi, xapar, xaperp)
    end if

    if (write_fcheck) then
       call fcheck (g, fcheck_f)
       if (proc0) then
          call open_output_file (unit, ".fcheck")
          do il = 2, nlambda
             do is = 1, nspec
                do ik = 1, naky
                   do it = 1, ntheta0
                      write (unit, "(20(1x,e12.6))") 0.5*(al(il) + al(il-1)), &
                           aky_out(ik), akx_out(it), fcheck_f(il,it,ik,is), &
                           sum((al(2:il)-al(1:il-1))*fcheck_f(2:il,it,ik,is))
                   end do
                end do
             end do
          end do
          call close_output_file (unit)
       end if
    end if

    if (write_intcheck) call finish_intcheck
    if (write_vortcheck) call finish_vortcheck
    if (write_fieldcheck) call finish_fieldcheck
    
    simtime = stime()/tnorm
    deltsave = simdt()/tnorm
    if (save_for_restart) then
       call gs2_save_for_restart (gnew, simtime, deltsave, istatus, .true.)
    end if

    call nc_finish

  end subroutine finish_gs2_diagnostics

  subroutine loop_diagnostics (istep, exit)
    use species, only: nspec, spec
    use theta_grid, only: theta, ntgrid, delthet, jacob
    use kt_grids, only: naky, ntheta0, aky_out, theta0, akx_out, aky
    use run_parameters, only: delt, woutunits, funits, fapar, fphi, faperp
    use run_parameters, only: tnorm
    use fields, only: phinorm, phinew, aparnew, aperpnew
    use fields, only: kperp, fieldlineavgphi
    use dist_fn, only: flux, intcheck, vortcheck, fieldcheck
    use dist_fn, only: neoclassical_flux, omega0, gamma0, getmoms
    use mp, only: proc0, broadcast
    use file_utils, only: get_unused_unit
    use prof, only: prof_entering, prof_leaving
    use gs2_time, only: update_time => update, stime
    use gs2_io, only: nc_qflux, nc_vflux, nc_pflux, nc_loop, nc_loop_moments
    use constants
    implicit none
    integer :: nout = 1
    integer, intent (in) :: istep
    logical, intent (out) :: exit
    complex, dimension (ntheta0, naky) :: omega, omegaavg
    real, dimension (ntheta0, naky) :: phitot, akperp
    complex, dimension (ntheta0, naky, nspec) :: pfluxneo,qfluxneo
    real :: phi2, apar2, aperp2
    real, dimension (ntheta0, naky) :: phi2_by_mode, apar2_by_mode, aperp2_by_mode
    real, dimension (ntheta0, naky, nspec) :: ntot2_by_mode
    real :: anorm, dmix, dmix4, dmixx
    real :: t, denom
    integer :: ig, ik, it, is, unit
    complex :: phiavg, sourcefac
    complex, dimension (-ntgrid:ntgrid,ntheta0,naky,nspec) :: ntot, density, &
         upar, tpar, tperp
    complex, dimension (ntheta0, nspec) :: ntot00, density00, upar00, tpar00, tperp00
    complex, dimension (ntheta0) :: phi00
    real, dimension (:), allocatable :: dl_over_b
    real, dimension (nspec) ::  heat_fluxes,  part_fluxes, mom_fluxes,  ntot2
    real, dimension (nspec) :: mheat_fluxes, mpart_fluxes, mmom_fluxes
    real, dimension (nspec) :: bheat_fluxes, bpart_fluxes, bmom_fluxes
    real, dimension (nspec) ::  heat_par,  heat_perp
    real, dimension (nspec) :: mheat_par, mheat_perp
    real, dimension (nspec) :: bheat_par, bheat_perp
    real, dimension (naky) :: fluxfac
    real :: hflux_tot, zflux_tot, vflux_tot
    character(200) :: filename

    call prof_entering ("loop_diagnostics")

    exit = .false.

    if (proc0) call get_omegaavg (istep, exit, omegaavg)
    call broadcast (exit)

    call prof_leaving ("loop_diagnostics")

    call update_time (delt)

    if (mod(istep,nwrite) /= 0 .and. .not. exit) return
    t = stime()/tnorm

    call prof_entering ("loop_diagnostics-1")

    if (proc0) then
       omega = omegahist(mod(istep,navg),:,:)
       sourcefac = exp(-zi*omega0*t+gamma0*t)
       call phinorm (phitot)
       call get_vol_average (phinew, phinew, phi2, phi2_by_mode)
       if (fapar > epsilon(0.0)) then
          call get_vol_average (aparnew, aparnew, apar2, apar2_by_mode)
          apar2 = apar2
          apar2_by_mode = apar2_by_mode
       end if
       if (faperp > epsilon(0.0)) then
          call get_vol_average (aperpnew, aperpnew, aperp2, aperp2_by_mode)
          aperp2 = aperp2
          aperp2_by_mode = aperp2_by_mode
       end if
    end if

    if (write_any_fluxes) then
       call flux (phinew, aparnew, aperpnew, &
            pflux, qheat, qheat_par, qheat_perp, vflux, &
            pmflux, qmheat, qmheat_par, qmheat_perp, vmflux, &
            pbflux, qbheat, qbheat_par, qbheat_perp, vbflux, &
            anorm)
       if (proc0) then
          if (fphi > epsilon(0.0)) then
             do is = 1, nspec
                call get_volume_average (qheat(:,:,is), heat_fluxes(is))
                heat_fluxes(is) = heat_fluxes(is)*anorm*funits &
                     *spec(is)%dens*spec(is)%temp
                
                call get_volume_average (qheat_par(:,:,is), heat_par(is))
                heat_par(is) = heat_par(is)*anorm*funits*spec(is)%dens*spec(is)%temp

                call get_volume_average (qheat_perp(:,:,is), heat_perp(is))
                heat_perp(is) = heat_perp(is)*anorm*funits*spec(is)%dens*spec(is)%temp
                
                call get_volume_average (pflux(:,:,is), part_fluxes(is))
                part_fluxes(is) = part_fluxes(is)*anorm*funits &
                     *spec(is)%dens
                
                call get_volume_average (vflux(:,:,is), mom_fluxes(is))
                mom_fluxes(is) = mom_fluxes(is)*anorm*funits &
                     *spec(is)%dens*spec(is)%mass*spec(is)%stm
             end do
          end if
          if (fapar > epsilon(0.0)) then
             do is = 1, nspec
                call get_volume_average (qmheat(:,:,is), mheat_fluxes(is))
                mheat_fluxes(is) = mheat_fluxes(is)*anorm*funits &
                     *spec(is)%dens*spec(is)%temp

                call get_volume_average (qmheat_par(:,:,is), mheat_par(is))
                mheat_par(is) = mheat_par(is)*anorm*funits*spec(is)%dens*spec(is)%temp

                call get_volume_average (qmheat_perp(:,:,is), mheat_perp(is))
                mheat_perp(is) = mheat_perp(is)*anorm*funits*spec(is)%dens*spec(is)%temp
                
                call get_volume_average (pmflux(:,:,is), mpart_fluxes(is))
                mpart_fluxes(is) = mpart_fluxes(is)*anorm*funits &
                     *spec(is)%dens

                call get_volume_average (vmflux(:,:,is), mmom_fluxes(is))
                mmom_fluxes(is) = mmom_fluxes(is)*anorm*funits &
                     *spec(is)%dens*spec(is)%mass*spec(is)%stm
             end do
          end if
          if (faperp > epsilon(0.0)) then
             do is = 1, nspec
                call get_volume_average (qbheat(:,:,is), bheat_fluxes(is))
                bheat_fluxes(is) = bheat_fluxes(is)*anorm*funits &
                     *spec(is)%dens*spec(is)%temp

                call get_volume_average (qbheat_par(:,:,is), bheat_par(is))
                bheat_par(is) = bheat_par(is)*anorm*funits*spec(is)%dens*spec(is)%temp

                call get_volume_average (qbheat_perp(:,:,is), bheat_perp(is))
                bheat_perp(is) = bheat_perp(is)*anorm*funits*spec(is)%dens*spec(is)%temp
                
                call get_volume_average (pbflux(:,:,is), bpart_fluxes(is))
                bpart_fluxes(is) = bpart_fluxes(is)*anorm*funits &
                     *spec(is)%dens

                call get_volume_average (vbflux(:,:,is), bmom_fluxes(is))
                bmom_fluxes(is) = bmom_fluxes(is)*anorm*funits &
                     *spec(is)%dens*spec(is)%mass*spec(is)%stm
             end do
          end if
       end if
    end if

    fluxfac = 0.5
    fluxfac(1) = 1.0

    if (proc0) then
       if (print_flux_line) then
          write (unit=*, fmt="('t= ',e16.10,' <phi**2>= ',e10.4, &
               & ' heat fluxes: ', 3(1x,e10.4))") &
               t, phi2, heat_fluxes(1:min(nspec,3))
          if (fapar > epsilon(0.0)) then
             write (unit=*, fmt="('t= ',e16.10,' <apar**2>= ',e10.4, &
                  & ' heat flux m: ', 3(1x,e10.4))") &
                  t, apar2, mheat_fluxes(1:min(nspec,3))
          end if
          if (faperp > epsilon(0.0)) then
             write (unit=*, fmt="('t= ',e16.10,' <aperp**2>= ',e10.4, &
                  & ' heat flux b: ', 3(1x,e10.4))") &
                  t, aperp2, bheat_fluxes(1:min(nspec,3))
          end if
       end if
       if (print_line) then
          if (print_old_units) then
             do ik = 1, naky
                do it = 1, ntheta0
                   write (unit=*, fmt="('aky=',f5.2, ' th0=',f5.2, &
                          &' om=',2f8.3,' omav=', 2f8.3,' phtot=',e10.4)") &
                        aky_out(ik), theta0(it,ik), &
                        real (omega(it,ik)), &
                        aimag(omega(it,ik)), &
                        real( omegaavg(it,ik)), &
                        aimag(omegaavg(it,ik)), &
                        phitot(it,ik)
                end do
             end do
          else
             do ik = 1, naky
                do it = 1, ntheta0
                   write (unit=*, fmt="('aky=',f5.2, ' th0=',f7.2, &
                          &' om=',2f8.3,' omav=', 2f8.3,' phtot=',e10.4)") &
                        aky_out(ik), theta0(it,ik), &
                        real( omega(it,ik)*woutunits(ik)), &
                        aimag(omega(it,ik)*woutunits(ik)), &
                        real( omegaavg(it,ik)*woutunits(ik)), &
                        aimag(omegaavg(it,ik)*woutunits(ik)), &
                        phitot(it,ik)
                end do
             end do
          end if
          write (*,*) 
       end if
    end if

    if (write_intcheck) call intcheck
    if (write_vortcheck) call vortcheck (phinew, aperpnew)
    if (write_fieldcheck) call fieldcheck (phinew, aparnew, aperpnew)

    call prof_leaving ("loop_diagnostics-1")

    if (.not. (write_any .or. dump_any)) return

    call prof_entering ("loop_diagnostics-2")

    call kperp (ntg_out, akperp)

    if (write_neoclassical_flux .or. dump_neoclassical_flux) then
       call neoclassical_flux (pfluxneo, qfluxneo)
    end if

    if (proc0 .and. write_any) then
       if (write_ascii) write (unit=out_unit, fmt=*) 'time=', t
       if (write_flux_line) then
          hflux_tot = 0.
          zflux_tot = 0.
          if (fphi > epsilon(0.0)) then
             if (write_ascii) then
                write (unit=out_unit, fmt="('t= ',e16.10,' <phi**2>= ',e10.4, &
                     & ' heat fluxes: ', 3(1x,e10.4))") &
                     t, phi2, heat_fluxes(1:min(nspec,3))
                write (unit=out_unit, fmt="('t= ',e16.10,' <phi**2>= ',e10.4, &
                     & ' part fluxes: ', 3(1x,e10.4))") &
                     t, phi2, part_fluxes(1:min(nspec,3))
             end if
             hflux_tot = sum(heat_fluxes)
             vflux_tot = sum(mom_fluxes)
             zflux_tot = sum(part_fluxes*spec%z)
          end if
          if (fapar > epsilon(0.0)) then
             if (write_ascii) then
                write (unit=out_unit, fmt="('t= ',e16.10,' <apar**2>= ',e10.4, &
                     & ' heat mluxes: ', 3(1x,e10.4))") &
                     t, apar2, mheat_fluxes(1:min(nspec,3))
                write (unit=out_unit, fmt="('t= ',e16.10,' <apar**2>= ',e10.4, &
                     & ' part mluxes: ', 3(1x,e10.4))") &
                     t, apar2, mpart_fluxes(1:min(nspec,3))
             end if
             hflux_tot = hflux_tot + sum(mheat_fluxes)
             vflux_tot = vflux_tot + sum(mmom_fluxes)
             zflux_tot = zflux_tot + sum(mpart_fluxes*spec%z)
          end if
          if (faperp > epsilon(0.0)) then
             if (write_ascii) then
                write (unit=out_unit, fmt="('t= ',e16.10,' <aperp**2>= ',e10.4, &
                     & ' heat bluxes: ', 3(1x,e10.4))") &
                     t, aperp2, bheat_fluxes(1:min(nspec,3))
                write (unit=out_unit, fmt="('t= ',e16.10,' <aperp**2>= ',e10.4, &
                     & ' part bluxes: ', 3(1x,e10.4))") &
                     t, aperp2, bpart_fluxes(1:min(nspec,3))
             end if
             hflux_tot = hflux_tot + sum(bheat_fluxes)
             vflux_tot = vflux_tot + sum(bmom_fluxes)
             zflux_tot = zflux_tot + sum(bpart_fluxes*spec%z)
          end if
          if (write_ascii) write (unit=out_unit, fmt="('t= ',e16.10,' h_tot= ',e10.4, &
               & ' z_tot= ',e10.4)") t, hflux_tot, zflux_tot
          call nc_qflux (nout, qheat, qmheat, qbheat, &
               heat_par, mheat_par, bheat_par, &
               heat_perp, mheat_perp, bheat_perp, &
               heat_fluxes, mheat_fluxes, bheat_fluxes, hflux_tot, anorm)
          call nc_vflux (nout, vflux, vmflux, vbflux, &
               mom_fluxes, mmom_fluxes, bmom_fluxes, vflux_tot)
          call nc_pflux (nout, pflux, pmflux, pbflux, &
               part_fluxes, mpart_fluxes, bpart_fluxes, zflux_tot)
          call nc_loop (nout, t, fluxfac, &
               phinew(igomega,:,:), phi2, phi2_by_mode, &
               aparnew(igomega,:,:), apar2, apar2_by_mode, &
               aperpnew(igomega,:,:), aperp2, aperp2_by_mode, &
               omega, omegaavg, woutunits, phitot)
       end if
       if (write_ascii) then
          do ik = 1, naky
             do it = 1, ntheta0
                write (out_unit,*) "aky=",aky_out(ik)," theta0=",theta0(it,ik)

                if (write_line) then
                   write (out_unit, "('t= ',e16.10,' aky= ',f5.2, ' akx= ',f5.2, &
                        &' om= ',2f8.3,' omav= ', 2f8.3,' phtot= ',e8.2)") &
                        t, aky_out(ik), akx_out(it), &
                        real( omega(it,ik)*woutunits(ik)), &
                        aimag(omega(it,ik)*woutunits(ik)), &
                        real( omegaavg(it,ik)*woutunits(ik)), &
                        aimag(omegaavg(it,ik)*woutunits(ik)), &
                        phitot(it,ik)
                end if
                
                if (write_omega) write (out_unit, *) &
                     'omega=', omega(it,ik), &
                     ' omega/(vt/a)=', real(omega(it,ik)*woutunits(ik)), &
                     aimag(omega(it,ik)*woutunits(ik))
                if (write_omavg) write (out_unit, *) &
                     'omavg=', omegaavg(it,ik), &
                     ' omavg/(vt/a)=', real(omegaavg(it,ik)*woutunits(ik)), &
                     aimag(omegaavg(it,ik)*woutunits(ik))
                
                if (write_dmix) then
                   if (abs(akperp(it,ik)) < 2.0*epsilon(0.0)) then
                      write (out_unit, *) 'dmix indeterminate'
                   else
                      dmix = 2.0*woutunits(ik)*aimag(omega(it,ik))/akperp(it,ik)**2
                      dmix4 = dmix*dmix/(woutunits(ik)*aimag(omega(it,ik)))
                      dmixx = 2.0*woutunits(ik)*aimag(omega(it,ik))/(akperp(it,ik)**2-aky(ik)**2)
                      write (out_unit, *) 'dmix=', dmix,' dmix4= ',dmix4,' dmixx= ',dmixx
                   end if
                end if
             end do
          end do
       end if
       do ik = 1, naky
          do it = 1, ntheta0
!             if (write_kperpnorm) then
!                if (aky_out(ik) /= 0.0) then
!                   write (out_unit, *) 'kperpnorm=',akperp(it,ik)/aky_out(ik)
!                else
!! huh? 
!                   write (out_unit, *) 'kperpnorm*aky=',akperp(it,ik)
!                end if
!             end if
!             if (write_phitot) write (out_unit, *) 'phitot=', phitot(it,ik)

!             if (write_fieldline_avg_phi) then
!                call fieldlineavgphi (ntg_out, it, ik, phiavg)
!                write (out_unit, *) '<phi>=',phiavg
!             end if

!             if (write_neoclassical_flux) then
!                write (out_unit, *) 'pfluxneo=', pfluxneo(it,ik,:)
!                write (out_unit, *) 'qfluxneo=', qfluxneo(it,ik,:)
!                write (out_unit, *) 'pfluxneo/sourcefac=', &
!                     pfluxneo(it,ik,:)/sourcefac
!                write (out_unit, *) 'qfluxneo/sourcefac=', &
!                     qfluxneo(it,ik,:)/sourcefac
!             end if

             if (write_nl_flux) then
                if (write_ascii) then
                   write (out_unit,"('ik,it,aky,akx,<phi**2>,t: ', &
                        & 2i5,4(1x,e12.6))") &
                        ik, it, aky_out(ik), akx_out(it), phi2_by_mode(it,ik), t
                end if
!                if (ntheta0 > 1 .and. ik == 1) then
!                   write (out_unit, "('it,akx,<phi**2>,t: ',i5,3(1x,e12.6))") &
!                        it, akx_out(it), sum(phi2_by_mode(it,:)*fluxfac(:)), t
!                end if
!                if (naky > 1 .and. it == 1) then
!                   write (out_unit, "('ik,aky,<phi**2>,t: ',i5,3(1x,e12.6))") &
!                        ik, aky_out(ik), sum(phi2_by_mode(:,ik))*fluxfac(ik), t
!                end if
             end if
          end do
       end do
    end if

    call broadcast (write_avg_moments)
    if (write_avg_moments) then

       call getmoms (phinew, ntot, density, upar, tpar, tperp)

       if (proc0) then
          allocate (dl_over_b(-ntg_out:ntg_out))

          dl_over_b = delthet(-ntg_out:ntg_out)*jacob(-ntg_out:ntg_out)
          dl_over_b = dl_over_b / sum(dl_over_b)

          do is = 1, nspec
             do it = 1, ntheta0
                ntot00(it, is)   = sum( ntot(-ntg_out:ntg_out,it,1,is)*dl_over_b)
                density00(it,is) = sum(density(-ntg_out:ntg_out,it,1,is)*dl_over_b)
                upar00(it, is)   = sum( upar(-ntg_out:ntg_out,it,1,is)*dl_over_b)
                tpar00(it, is)   = sum( tpar(-ntg_out:ntg_out,it,1,is)*dl_over_b)
                tperp00(it, is)  = sum(tperp(-ntg_out:ntg_out,it,1,is)*dl_over_b)
             end do
          end do
          
          do it = 1, ntheta0
             phi00(it)   = sum( phinew(-ntg_out:ntg_out,it,1)*dl_over_b)
          end do

          deallocate (dl_over_b)

          do is = 1, nspec
             call get_vol_average (ntot(:,:,:,is), ntot(:,:,:,is), &
                  ntot2(is), ntot2_by_mode(:,:,is))
          end do

          call nc_loop_moments (nout, ntot2, ntot2_by_mode, &
               phi00, ntot00, density00, upar00, tpar00, tperp00)

       end if
    end if

    if (proc0 .and. dump_any) then
!
! I have not checked the units in this section. BD
!       
       do ik = 1, naky
          do it = 1, ntheta0
!             it = 2
             if (dump_neoclassical_flux) then
                do is = 1, nspec
                   write (dump_neoclassical_flux_unit, "(20(1x,e12.5))") &
                        t, aky_out(ik), akx_out(it), real(is), pfluxneo(it,ik,is), &
                        qfluxneo(it,ik,is), &
                        pfluxneo(it,ik,is)/sourcefac, &
                        qfluxneo(it,ik,is)/sourcefac
                end do
             end if
             if (dump_check1) then
                denom=sum(delthet(-ntg_out:ntg_out-1)*jacob(-ntg_out:ntg_out-1))
                write (dump_check1_unit, "(20(1x,e12.5))") &
                     t, aky_out(ik), akx_out(it), &
                     sum(phinew(-ntg_out:ntg_out-1,it,ik) &
                         *delthet(-ntg_out:ntg_out-1) &
                         *jacob(-ntg_out:ntg_out-1))/denom, &
                     sum(phinew(-ntg_out:ntg_out-1,it,ik) &
                         *delthet(-ntg_out:ntg_out-1) &
                         *jacob(-ntg_out:ntg_out-1)) &
                         /denom/sourcefac
             end if
          end do
       end do
       if (dump_fields_periodically) then
          call get_unused_unit (unit)
          write (filename, "('dump.fields.t=',e12.6)") t
          open (unit=unit, file=filename, status="unknown")
          do ik = 1, naky
             do it = 1, ntheta0
                do ig = -ntg_out, ntg_out
                   write (unit=unit, fmt="(20(1x,e12.5))") &
                        theta(ig), aky_out(ik), theta0(it,ik), &
                        phinew(ig,it,ik), aparnew(ig,it,ik), &
                        aperpnew(ig,it,ik), t, akx_out(it)
                end do
             end do
          end do
          close (unit=unit)
       end if
     end if

     nout = nout + 1
    call prof_leaving ("loop_diagnostics-2")
  end subroutine loop_diagnostics

  subroutine get_omegaavg (istep, exit, omegaavg)
    use kt_grids, only: naky, ntheta0
    use fields, only: phi, apar, aperp, phinew, aparnew, aperpnew
    use run_parameters, only: delt
    use constants
    implicit none
    integer, intent (in) :: istep
    logical, intent (in out) :: exit
    complex, dimension (:,:), intent (out) :: omegaavg
    complex, dimension (navg,ntheta0,naky) :: domega
    integer :: j

    j = igomega
    where (abs(phinew(j,:,:)+aparnew(j,:,:)+aperpnew(j,:,:)) < epsilon(0.0) &
           .or. abs(phi(j,:,:)+apar(j,:,:)+aperp(j,:,:)) < epsilon(0.0))
       omegahist(mod(istep,navg),:,:) = 0.0
    elsewhere
       omegahist(mod(istep,navg),:,:) &
            = log((phinew(j,:,:) + aparnew(j,:,:) + aperpnew(j,:,:)) &
                  /(phi(j,:,:)   + apar(j,:,:)    + aperp(j,:,:)))*zi/delt
    end where

    omegaavg = sum(omegahist/real(navg),dim=1)

    if (istep > navg) then
       domega = spread(omegaavg,1,navg) - omegahist
       if (all(sqrt(sum(abs(domega)**2/real(navg),dim=1)) &
            < min(abs(omegaavg),1.0)*omegatol)) &
       then
          if (write_ascii) write (out_unit, "('*** omega converged')")
          exit = .true. .and. exit_when_converged
       end if

       if (any(abs(omegaavg)*delt > omegatinst)) then
          if (write_ascii) write (out_unit, "('*** numerical instability detected')") 
          exit = .true.
       end if
    end if
  end subroutine get_omegaavg

  subroutine get_vol_average (a, b, axb, axb_by_mode)
    use theta_grid, only: ntgrid, delthet, jacob
    use kt_grids, only: naky, ntheta0
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: a, b
    real, intent (out) :: axb
    real, dimension (:,:), intent (out) :: axb_by_mode

    integer :: ig, ik, it
    integer :: ng
    real, dimension (-ntg_out:ntg_out) :: wgt
    real :: anorm

    ng = ntg_out
    wgt = delthet(-ng:ng)*jacob(-ng:ng)
    anorm = sum(wgt)

    do ik = 1, naky
       do it = 1, ntheta0
          axb_by_mode(it,ik) &
               = sum(real(conjg(a(-ng:ng,it,ik))*b(-ng:ng,it,ik))*wgt)/anorm
       end do
    end do

    call get_volume_average (axb_by_mode, axb)
  end subroutine get_vol_average

  subroutine get_volume_average (f, favg)
    use mp, only: iproc
    use kt_grids, only: naky, ntheta0, aky
    implicit none
    real, dimension (:,:), intent (in) :: f
    real, intent (out) :: favg
!    real, dimension (naky) :: fac
    real :: fac
    integer :: ik, it

! bug fixed 3.17.00
!    fac = 1.0
!    if (naky > 1) fac(1) = 0.5
! old field array order: 
!    favg = sum(f*spread(fac,2,ntheta0))

!    fac = 0.5
!    fac(1) = 1.0

    favg = 0.
    do ik = 1, naky
       fac = 0.5
       if (aky(ik) == 0.) fac = 1.0
       do it = 1, ntheta0
          favg = favg + f(it, ik) * fac
       end do
    end do

  end subroutine get_volume_average

end module gs2_diagnostics
