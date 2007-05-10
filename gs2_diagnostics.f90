module gs2_diagnostics
  implicit none

  public :: init_gs2_diagnostics
  public :: finish_gs2_diagnostics
  public :: loop_diagnostics

  ! knobs

  logical :: print_line, print_old_units, print_flux_line
  logical :: write_line, write_flux_line, write_phi, write_apar, write_aperp
  logical :: write_omega, write_omavg
  logical :: write_qheat, write_pflux, write_vflux
  logical :: write_qmheat, write_pmflux, write_vmflux
  logical :: write_dmix, write_kperpnorm, write_phitot
  logical :: write_eigenfunc, write_final_fields, write_final_antot
  logical :: write_fcheck
  logical :: write_intcheck, write_vortcheck, write_fieldcheck
  logical :: write_fieldline_avg_phi
  logical :: write_neoclassical_flux, write_nl_flux
  integer :: nwrite, igomega
  integer :: navg
  real :: omegatol, omegatinst
  logical :: exit_when_converged

  logical :: dump_neoclassical_flux, dump_omega, dump_check1
  logical :: dump_fields_periodically

  logical :: save_for_restart
  character(300) :: restart_file

  integer :: nperiod_output

  ! internal
  logical :: write_any, write_any_fluxes, dump_any
  logical, private :: initialized = .false.

  integer :: out_unit
  integer :: dump_neoclassical_flux_unit, dump_omega_unit, dump_check1_unit

  complex, dimension (:,:,:), allocatable :: omegahist
  ! (navg,naky,ntheta0)

  real, dimension (:,:,:), allocatable :: pflux, qheat, vflux
  real, dimension (:,:,:), allocatable :: pmflux, qmheat, vmflux
  ! (naky,ntheta0,nspec)

  integer :: ntg_out

contains

  subroutine init_gs2_diagnostics
    use file_utils, only: open_output_file
    use theta_grid, only: init_theta_grid
    use theta_grid, only: theta, delthet, jacob
    use kt_grids, only: init_kt_grids, naky, ntheta0
    use fields_arrays, only: phinew
    use run_parameters, only: init_run_parameters
    use species, only: init_species
    use dist_fn, only: init_dist_fn
    use dist_fn, only: init_intcheck, init_vortcheck, init_fieldcheck
    use mp, only: proc0, broadcast
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

    if (proc0) call real_init
    call broadcast (nwrite)
    call broadcast (write_any)
    call broadcast (write_any_fluxes)
    call broadcast (write_neoclassical_flux)
    call broadcast (write_nl_flux)
    call broadcast (write_fcheck)
    call broadcast (dump_any)
    call broadcast (dump_neoclassical_flux)
    call broadcast (dump_omega)
    call broadcast (dump_check1)
    call broadcast (dump_fields_periodically)
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
    call open_output_file (out_unit, ".out")

    if (dump_neoclassical_flux) then
       call get_unused_unit (dump_neoclassical_flux_unit)
       open (unit=dump_neoclassical_flux_unit, file="dump.neoflux", &
            status="unknown")
    end if

    if (dump_omega) then
       call get_unused_unit (dump_omega_unit)
       open (unit=dump_omega_unit, file="dump.omega", status="unknown")
    end if

    if (dump_check1) then
       call get_unused_unit (dump_check1_unit)
       open (unit=dump_check1_unit, file="dump.check1", status="unknown")
    end if

    write (unit=out_unit, fmt="('gs2')")
    datestamp(:) = ' '
    timestamp(:) = ' '
    zone(:) = ' '
    call date_and_time (datestamp, timestamp, zone)
    write (unit=out_unit, fmt="('Date: ',a,' Time: ',a,x,a)") &
         trim(datestamp), trim(timestamp), trim(zone)

    allocate (omegahist(0:navg-1,naky,ntheta0))
    omegahist = 0.0

    allocate (pflux(naky,ntheta0,nspec))
    allocate (qheat(naky,ntheta0,nspec))
    allocate (vflux(naky,ntheta0,nspec))
    allocate (pmflux(naky,ntheta0,nspec))
    allocate (qmheat(naky,ntheta0,nspec))
    allocate (vmflux(naky,ntheta0,nspec))

  end subroutine real_init

  subroutine read_parameters
    use file_utils, only: input_unit, run_name
    use theta_grid, only: nperiod, ntheta
    use dist_fn, only: nperiod_guard
    implicit none
    namelist /gs2_diagnostics_knobs/ &
         print_line, print_old_units, print_flux_line, &
         write_line, write_flux_line, write_phi, write_apar, write_aperp, &
         write_omega, write_omavg, &
         write_qheat, write_pflux, write_vflux, &
         write_qmheat, write_pmflux, write_vmflux, &
         write_dmix, write_kperpnorm, write_phitot, &
         write_eigenfunc, write_final_fields, write_final_antot, &
         write_fcheck, &
         write_intcheck, write_vortcheck, write_fieldcheck, &
         write_fieldline_avg_phi, write_neoclassical_flux, write_nl_flux, &
         nwrite, navg, omegatol, omegatinst, igomega, &
         exit_when_converged, &
         dump_neoclassical_flux, dump_omega, dump_check1, &
         dump_fields_periodically, &
         nperiod_output, &
         save_for_restart, restart_file

    print_line = .true.
    print_old_units = .true.
    print_flux_line = .false.
    write_line = .true.
    write_flux_line = .true.
    write_phi = .true.
    write_apar = .true.
    write_aperp = .true.
    write_omega = .true.
    write_omavg = .true.
    write_qheat = .true.
    write_pflux = .true.
    write_vflux = .true.
    write_qmheat = .true.
    write_pmflux = .true.
    write_vmflux = .true.
    write_dmix = .true.
    write_kperpnorm = .true.
    write_phitot = .true.
    write_fieldline_avg_phi = .false.
    write_neoclassical_flux = .false.
    write_nl_flux = .false.
    write_eigenfunc = .true.
    write_final_fields = .false.
    write_final_antot = .false.
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
    dump_omega = .false.
    dump_check1 = .false.
    dump_fields_periodically = .false.
    nperiod_output = nperiod - nperiod_guard
    save_for_restart = .false.
    restart_file = trim(run_name)//".nc"
    read (unit=input_unit("gs2_diagnostics_knobs"), nml=gs2_diagnostics_knobs)

    write_any = write_line .or. write_phi       .or. write_apar &
         .or. write_aperp  .or. write_omega     .or. write_omavg &
         .or. write_qheat  .or. write_pflux     .or. write_vflux &
         .or. write_qmheat .or. write_pmflux    .or. write_vmflux &
         .or. write_dmix   .or. write_kperpnorm .or. write_phitot &
         .or. write_fieldline_avg_phi           .or. write_neoclassical_flux &
         .or. write_flux_line                   .or. write_nl_flux
    write_any_fluxes =  &
              write_qheat  .or. write_pflux     .or. write_vflux &
         .or. write_qmheat .or. write_pmflux    .or. write_vmflux &
         .or. write_flux_line                   .or. print_flux_line &
         .or. write_nl_flux
    dump_any = dump_neoclassical_flux           .or. dump_omega &
         .or. dump_check1  .or. dump_fields_periodically

    nperiod_output = min(nperiod,nperiod_output)
    ntg_out = ntheta/2 + (nperiod_output-1)*ntheta
  end subroutine read_parameters

  subroutine finish_gs2_diagnostics (istep)
    use file_utils, only: open_output_file, close_output_file, get_unused_unit
    use mp, only: proc0, broadcast
    use species, only: nspec
    use theta_grid, only: ntgrid, theta
    use kt_grids, only: naky, ntheta0, theta0, nx, aky_out, akx_out
    use le_grids, only: nlambda, negrid, fcheck, al
    use fields_arrays, only: phi, apar, aperp, phinew, aparnew, aperpnew, aminv
    use dist_fn, only: getan
    use dist_fn, only: finish_intcheck, finish_vortcheck, finish_fieldcheck
    use dist_fn_arrays, only: g, gnew
    use gs2_layouts, only: g_lo, proc_id
    use gs2_layouts, only: ik_idx, it_idx, il_idx, ie_idx, is_idx
    use gs2_save, only: gs2_save_for_restart
    use constants
    use gs2_time, only: stime
    implicit none
    integer, intent (in) :: istep
    integer :: ig, ik, it, il, is, unit
    complex, dimension (nlambda,naky,ntheta0,nspec) :: fcheck_f
    complex, dimension (:,:,:), allocatable :: xphi, xapar, xaperp
    complex, dimension (-ntgrid:ntgrid,naky,ntheta0) :: antot, antota, antotp
    complex :: phi0, phase
    real :: simtime
    integer :: istatus

    if (proc0) then
       call close_output_file (out_unit)
       if (dump_neoclassical_flux) close (dump_neoclassical_flux_unit)
       if (dump_omega) close (dump_omega_unit)
       if (dump_check1) close (dump_check1_unit)

       if (write_eigenfunc) then
          call open_output_file (unit, ".eigenfunc")
          do ik = 1, naky
             do it = 1, ntheta0
                phi0 = phi(0,ik,it)
                if (abs(phi0) < 10.0*epsilon(0.0)) then
                   phi0 = phi(1,ik,it)/(theta(1)-theta(0))
                end if
                if (abs(phi0) < 10.0*epsilon(0.0)) then
                   phi0 = 1.0
                end if
                do ig = -ntg_out, ntg_out
                   write (unit, "(9(x,e12.5))") &
                        theta(ig), theta0(ik,it), aky_out(ik), &
                        phi(ig,ik,it)/phi0, &
                        apar(ig,ik,it)/phi0, &
                        aperp(ig,ik,it)/phi0
                end do
                write (unit, "()")
             end do
          end do
          call close_output_file (unit)
       end if

       if (write_final_fields) then
          call open_output_file (unit, ".fields")
          do ik = 1, naky
!             if (phi(0,ik,1) /= 0.) then
!                phase = 1. / phi(0,ik,1)
!             else
                phase = 1.0
!             endif
             do it = 1, ntheta0
                do ig = -ntg_out, ntg_out
                   write (unit, "(10(x,e12.5))") &
                        theta(ig), theta0(ik,it), aky_out(ik), &
                        phase*phi(ig,ik,it), phase*apar(ig,ik,it), phase*aperp(ig,ik,it), &
                        theta(ig) - theta0(ik,it)
                end do
                write (unit, "()")
             end do
          end do
          call close_output_file (unit)
       end if
    end if

    call broadcast (write_final_antot)
    if (write_final_antot) then
       call getan (antot, antota, antotp)
       if (proc0) then
          call open_output_file (unit, ".antot")
          do ik = 1, naky
             do it = 1, ntheta0
                do ig = -ntg_out, ntg_out
                   write (unit, "(10(x,e12.5))") &
                        theta(ig), theta0(ik,it), aky_out(ik), &
                        antot(ig,ik,it), antota(ig,ik,it), antotp(ig,ik,it), &
                        theta(ig) - theta0(ik,it)
                end do
                write (unit, "()")
             end do
          end do
          call close_output_file (unit)
       end if
    end if

    if (write_fcheck) then
       call fcheck (g, fcheck_f)
       if (proc0) then
          call open_output_file (unit, ".fcheck")
          do il = 2, nlambda
             do is = 1, nspec
                do ik = 1, naky
                   do it = 1, ntheta0
                      write (unit, "(20(x,e12.6))") 0.5*(al(il) + al(il-1)), &
                           aky_out(ik), akx_out(it), fcheck_f(il,ik,it,is), &
                           sum((al(2:il)-al(1:il-1))*fcheck_f(2:il,ik,it,is))
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
    
    simtime = stime()
    if (save_for_restart) then
       call gs2_save_for_restart &
            (gnew, aminv, simtime, istatus)
    end if
  end subroutine finish_gs2_diagnostics

  subroutine loop_diagnostics (istep, nstep, exit)
    use species, only: nspec
    use theta_grid, only: theta, ntgrid, delthet, jacob
    use kt_grids, only: naky, ntheta0, aky_out, theta0, akx_out, aky
    use run_parameters, only: delt, woutunits
    use fields, only: phinorm, phinew, aparnew, aperpnew
    use fields, only: kperp, fieldlineavgphi
    use dist_fn, only: flux, intcheck, vortcheck, fieldcheck
    use dist_fn, only: neoclassical_flux, omega0, gamma0
    use mp, only: proc0, broadcast
    use file_utils, only: get_unused_unit
    use prof, only: prof_entering, prof_leaving
    use gs2_time, only: update_time => update, stime
    use constants
    implicit none
    integer, intent (in) :: istep, nstep
    logical, intent (out) :: exit
    complex, dimension (naky,ntheta0) :: omega, omegaavg
    real, dimension (naky,ntheta0) :: phitot, akperp
    complex, dimension (naky,ntheta0,nspec) :: pfluxneo,qfluxneo
    real :: phi2
    real, dimension (naky,ntheta0) :: phi2_by_mode
    real :: anorm
    real :: t, denom
    integer :: ig, ik, it, is, unit
    complex :: phiavg, sourcefac
    character(200) :: filename
    real, dimension (nspec) :: heat_fluxes
    real, dimension (naky) :: fluxfac

    call prof_entering ("loop_diagnostics")

    exit = .false.

    if (proc0) call get_omegaavg (istep, exit, omegaavg)
    call broadcast (exit)

    call prof_leaving ("loop_diagnostics")

    call update_time (delt)

    if (mod(istep,nwrite) /= 0 .and. .not. exit) return
    t = stime()

    call prof_entering ("loop_diagnostics-1")

    if (proc0) then
       omega = omegahist(mod(istep,navg),:,:)
       sourcefac = exp(-zi*omega0*t+gamma0*t)
       call phinorm (phitot)
       call get_average_flux (phinew, phinew, phi2, phi2_by_mode)
    end if

    if (write_any_fluxes) then
       call flux (phinew, aparnew, aperpnew, &
            pflux, qheat, vflux, pmflux, qmheat, vmflux, anorm)
       do is = 1, nspec
          call get_volume_average (qheat(:,:,is), heat_fluxes(is))
          heat_fluxes(is) = heat_fluxes(is)*anorm
       end do
    end if

    fluxfac = 1.0
    if (naky > 1 .and. aky(1) == 0.0) fluxfac(1) = 0.5

    if (proc0) then
       if (print_flux_line) then
          write (unit=*, fmt="('t=',e12.6,' <phi**2>=',e10.4, &
               & ' heat fluxes:', 3(1x,e10.4))") &
               t, phi2, heat_fluxes(1:min(nspec,3))
       end if
       if (print_line) then
          if (print_old_units) then
             do ik = 1, naky
                do it = 1, ntheta0
                   write (unit=*, fmt="('aky=',f5.2, ' th0=',f5.2, &
                          &' om=',2f8.3,' omav=', 2f8.3,' phtot=',e10.4)") &
                        aky_out(ik), theta0(ik,it), &
                        real(omega(ik,it)), &
                        aimag(omega(ik,it)), &
                        real(omegaavg(ik,it)), &
                        aimag(omegaavg(ik,it)), &
                        phitot(ik,it)
                end do
             end do
          else
             do ik = 1, naky
                do it = 1, ntheta0
                   write (unit=*, fmt="('aky=',f5.2, ' th0=',f5.2, &
                          &' om=',2f8.3,' omav=', 2f8.3,' phtot=',e10.4)") &
                        aky_out(ik), theta0(ik,it), &
                        real(omega(ik,it)*woutunits(ik)), &
                        aimag(omega(ik,it)*woutunits(ik)), &
                        real(omegaavg(ik,it)*woutunits(ik)), &
                        aimag(omegaavg(ik,it)*woutunits(ik)), &
                        phitot(ik,it)
                end do
             end do
          end if
          write (*,*) 
       end if
    end if

    if (write_intcheck) call intcheck
    if (write_vortcheck) call vortcheck (phinew, aparnew, aperpnew)
    if (write_fieldcheck) call fieldcheck (phinew, aparnew, aperpnew)

    call prof_leaving ("loop_diagnostics-1")

    if (.not. (write_any .or. dump_any)) return

    call prof_entering ("loop_diagnostics-2")

    call kperp (ntg_out, akperp)

    if (write_neoclassical_flux .or. dump_neoclassical_flux) then
       call neoclassical_flux (pfluxneo, qfluxneo)
    end if

    if (proc0 .and. write_any) then
       write (unit=out_unit, fmt=*) 'time=', t
       if (write_flux_line) then
          write (unit=out_unit, fmt="('t=',e12.6,' <phi**2>=',e10.4, &
               & ' heat fluxes:', 3(1x,e10.4))") &
               t, phi2, heat_fluxes(1:min(nspec,3))
       end if
       do ik = 1, naky
          do it = 1, ntheta0
             write (out_unit,*) "aky=",aky_out(ik)," theta0=",theta0(ik,it)

             if (write_line) then
                write (out_unit, "('aky=',f5.2, ' th0=',f4.2, &
                       &' om=',2f8.3,' omav=', 2f8.3,' phtot=',e8.2)") &
                     aky_out(ik), theta0(ik,it), &
                     real(omega(ik,it)), &
                     aimag(omega(ik,it)), &
                     real(omegaavg(ik,it)), &
                     aimag(omegaavg(ik,it)), &
                     phitot(ik,it)
             end if

             if (write_phi) write (out_unit, *) 'phi(0)=', phinew(0,ik,it)
             if (write_apar) write (out_unit, *) 'apar(0)=',aparnew(0,ik,it)
             if (write_aperp) write (out_unit, *) 'aperp(0)=',aperpnew(0,ik,it)
             if (write_omega) write (out_unit, *) &
                  'omega=', omega(ik,it), &
                  ' omega/(vt/a)=', real(omega(ik,it)*woutunits(ik)), &
                  aimag(omega(ik,it)*woutunits(ik))
             if (write_omavg) write (out_unit, *) &
                  'omavg=', omegaavg(ik,it), &
                  ' omavg/(vt/a)=', real(omegaavg(ik,it)*woutunits(ik)), &
                  aimag(omegaavg(ik,it)*woutunits(ik))

             if (write_qheat) write (out_unit, *) 'qheat=', qheat(ik,it,:)
             if (write_pflux) write (out_unit, *) 'pflux=', pflux(ik,it,:)
             if (write_vflux) write (out_unit, *) 'vflux=', vflux(ik,it,:)
             if (write_qmheat) write (out_unit, *) 'qmheat=', qmheat(ik,it,:)
             if (write_pmflux) write (out_unit, *) 'pmflux=', pmflux(ik,it,:)
             if (write_vmflux) write (out_unit, *) 'vmflux=', vmflux(ik,it,:)

             if (write_dmix) then
                if (abs(akperp(ik,it)) < 2.0*epsilon(0.0)) then
                   write (out_unit, *) 'dmix indeterminate'
                else
                   write (out_unit, *) 'dmix=', &
                        2.0*woutunits(ik)*aimag(omega(ik,it))/akperp(ik,it)**2
                end if
             end if
             if (write_kperpnorm) then
                if (aky_out(ik) /= 0.0) then
                   write (out_unit, *) 'kperpnorm=',akperp(ik,it)/aky_out(ik)
                else
! huh? 
                   write (out_unit, *) 'kperpnorm*aky=',akperp(ik,it)
                end if
             end if
             if (write_phitot) write (out_unit, *) 'phitot=', phitot(ik,it)

             if (write_fieldline_avg_phi) then
                call fieldlineavgphi (ntg_out, ik, it, phiavg)
                write (out_unit, *) '<phi>=',phiavg
             end if

             if (write_neoclassical_flux) then
                write (out_unit, *) 'pfluxneo=', pfluxneo(ik,it,:)
                write (out_unit, *) 'qfluxneo=', qfluxneo(ik,it,:)
                write (out_unit, *) 'pfluxneo/sourcefac=', &
                     pfluxneo(ik,it,:)/sourcefac
                write (out_unit, *) 'qfluxneo/sourcefac=', &
                     qfluxneo(ik,it,:)/sourcefac
             end if

             if (write_nl_flux) then
                write (out_unit,"('ik,it,aky,akx,<phi**2>,t: ', &
                     & 2i5,4(x,e12.6))") &
                     ik, it, aky_out(ik), akx_out(it), phi2_by_mode(ik,it), t
                if (ntheta0 > 1 .and. ik == 1) then
                   write (out_unit, "('it,akx,<phi**2>,t: ',i5,3(x,e12.6))") &
                        it, akx_out(it), sum(phi2_by_mode(:,it)*fluxfac(:)), t
                end if
                if (naky > 1 .and. it == 1) then
                   write (out_unit, "('ik,aky,<phi**2>,t: ',i5,3(x,e12.6))") &
                        ik, aky_out(ik), sum(phi2_by_mode(ik,:)*fluxfac(ik)), t
                end if
             end if
          end do
       end do
    end if

    if (proc0 .and. dump_any) then
       do ik = 1, naky
          do it = 1, ntheta0
             if (dump_neoclassical_flux) then
                do is = 1, nspec
                   write (dump_neoclassical_flux_unit, "(20(1x,e12.5))") &
                        t, aky_out(ik), akx_out(it), real(is), pfluxneo(ik,it,is), &
                        qfluxneo(ik,it,is), &
                        pfluxneo(ik,it,is)/sourcefac, &
                        qfluxneo(ik,it,is)/sourcefac
                end do
             end if
             if (dump_omega) then
                write (dump_omega_unit, "(20(1x,e12.5))") &
                     t, aky_out(ik), akx_out(it), omega(ik,it)*woutunits(ik), &
                     omegaavg(ik,it)*woutunits(ik)
             end if
             if (dump_check1) then
                denom=sum(delthet(-ntg_out:ntg_out-1)*jacob(-ntg_out:ntg_out-1))
                write (dump_check1_unit, "(20(1x,e12.5))") &
                     t, aky_out(ik), akx_out(it), &
                     sum(phinew(-ntg_out:ntg_out-1,ik,it) &
                         *delthet(-ntg_out:ntg_out-1) &
                         *jacob(-ntg_out:ntg_out-1))/denom, &
                     sum(phinew(-ntg_out:ntg_out-1,ik,it) &
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
          do it = 1, ntheta0
             do ik = 1, naky
                do ig = -ntg_out, ntg_out
                   write (unit=unit, fmt="(20(x,e12.5))") &
                        theta(ig), aky_out(ik), theta0(ik,it), &
                        phinew(ig,ik,it), aparnew(ig,ik,it), &
                        aperpnew(ig,ik,it), t, akx_out(it)
                end do
             end do
          end do
          close (unit=unit)
       end if
     end if
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
    complex, dimension (navg,naky,ntheta0) :: domega
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
          write (out_unit, "('*** omega converged')")
          exit = .true. .and. exit_when_converged
       end if

       if (any(abs(omegaavg)*delt > omegatinst)) then
          write (out_unit, "('*** numerical instability detected')")
          exit = .true.
       end if
    end if
  end subroutine get_omegaavg

  subroutine get_average_flux (a, b, axb, axb_by_mode)
    use theta_grid, only: ntgrid, delthet, grho, jacob
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
    wgt = delthet(-ng:ng)*grho(-ng:ng)*jacob(-ng:ng)
    anorm = sum(wgt)

    do it = 1, ntheta0
       do ik = 1, naky
          axb_by_mode(ik,it) &
               = sum(real(conjg(a(-ng:ng,ik,it))*b(-ng:ng,ik,it))*wgt)/anorm
       end do
    end do

    call get_volume_average (axb_by_mode, axb)
  end subroutine get_average_flux

  subroutine get_volume_average (f, favg)
    use kt_grids, only: naky, ntheta0
    implicit none
    real, dimension (:,:), intent (in) :: f
    real, intent (out) :: favg
    real, dimension (naky) :: fac

    fac = 1.0
    if (naky > 1) fac(1) = 0.5

    favg = sum(f*spread(fac,2,ntheta0))

  end subroutine get_volume_average

end module gs2_diagnostics
