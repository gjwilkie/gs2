module gs2_diagnostics
  implicit none

  public :: init_gs2_diagnostics
  public :: finish_gs2_diagnostics
  public :: loop_diagnostics

  ! knobs

  logical :: print_line, print_old_units
  logical :: write_line, write_phi, write_apar, write_aperp
  logical :: write_omega, write_omavg, write_eik
  logical :: write_qheat, write_pflux, write_vflux
  logical :: write_qmheat, write_pmflux, write_vmflux
  logical :: write_dmix, write_kperpnorm, write_phitot
  logical :: write_eigenfunc, write_final_fields, write_fcheck
  logical :: write_intcheck, write_vortcheck, write_fieldcheck
  logical :: write_fieldline_avg_phi
  logical :: write_neoclassical_flux
  integer :: nwrite, igomega
  integer :: navg
  real :: omegatol, omegatinst
  logical :: exit_when_converged

  logical :: dump_neoclassical_flux, dump_omega, dump_check1, dump_check2

  ! internal
  logical :: write_any, write_any_fluxes, dump_any

  integer :: out_unit
  integer :: dump_neoclassical_flux_unit, dump_omega_unit, &
	dump_check1_unit, dump_check2_unit

  complex, dimension (:,:,:), allocatable :: omegahist
  ! (navg,naky,ntheta0)

  real, dimension (:,:,:), allocatable :: pflux, qheat, vflux
  real, dimension (:,:,:), allocatable :: pmflux, qmheat, vmflux
  ! (naky,ntheta0,nspec)

contains

  subroutine init_gs2_diagnostics
    use file_utils, only: open_output_file
    use kt_grids, only: init_kt_grids
    use run_parameters, only: init_run_parameters
    use species, only: init_species
    use dist_fn, only: init_intcheck, init_vortcheck, init_fieldcheck
    use mp, only: proc0, broadcast
    implicit none
    logical, save :: initialized = .false.

    if (initialized) return
    initialized = .true.

    call init_kt_grids
    call init_run_parameters
    call init_species

    if (proc0) call real_init
    call broadcast (nwrite)
    call broadcast (write_any)
    call broadcast (write_any_fluxes)
    call broadcast (write_neoclassical_flux)
    call broadcast (write_fcheck)
    call broadcast (dump_any)
    call broadcast (dump_neoclassical_flux)
    call broadcast (dump_omega)
    call broadcast (dump_check1)
    call broadcast (dump_check2)

    call broadcast (write_intcheck)
    call broadcast (write_vortcheck)
    call broadcast (write_fieldcheck)
    if (write_intcheck) call init_intcheck
    if (write_vortcheck) call init_vortcheck
    if (write_fieldcheck) call init_fieldcheck
  end subroutine init_gs2_diagnostics

  subroutine real_init
    use file_utils, only: open_output_file, get_unused_unit
    use kt_grids, only: naky, ntheta0, aky
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

    if (dump_check2) then
       call get_unused_unit (dump_check2_unit)
       open (unit=dump_check2_unit, file="dump.check2", status="unknown")
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
    use file_utils, only: input_unit
    implicit none
    namelist /gs2_diagnostics_knobs/ &
         print_line, print_old_units, &
         write_line, write_phi, write_apar, write_aperp, &
         write_omega, write_omavg, write_eik, &
         write_qheat, write_pflux, write_vflux, &
         write_qmheat, write_pmflux, write_vmflux, &
         write_dmix, write_kperpnorm, write_phitot, &
         write_eigenfunc, write_final_fields, write_fcheck, &
         write_intcheck, write_vortcheck, write_fieldcheck, &
         write_fieldline_avg_phi, write_neoclassical_flux, &
         nwrite, navg, omegatol, omegatinst, igomega, &
         exit_when_converged, &
         dump_neoclassical_flux, dump_omega, dump_check1, dump_check2

    print_line = .true.
    print_old_units = .false.
    write_line = .true.
    write_phi = .true.
    write_apar = .true.
    write_aperp = .true.
    write_omega = .true.
    write_omavg = .true.
    write_eik = .true.
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
    write_eigenfunc = .true.
    write_final_fields = .false.
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
    dump_check2 = .false.
    read (unit=input_unit("gs2_diagnostics_knobs"), nml=gs2_diagnostics_knobs)

    write_any = write_line .or. write_phi       .or. write_apar &
         .or. write_aperp  .or. write_omega     .or. write_omavg &
         .or. write_qheat  .or. write_pflux     .or. write_vflux &
         .or. write_qmheat .or. write_pmflux    .or. write_vmflux &
         .or. write_dmix   .or. write_kperpnorm .or. write_phitot &
         .or. write_fieldline_avg_phi           .or. write_neoclassical_flux &
	 .or. write_eik
    write_any_fluxes =  &
              write_qheat  .or. write_pflux     .or. write_vflux &
         .or. write_qmheat .or. write_pmflux    .or. write_vmflux
    dump_any = dump_neoclassical_flux           .or. dump_omega &
         .or. dump_check1  .or. dump_check2
  end subroutine read_parameters

  subroutine finish_gs2_diagnostics
    use file_utils, only: open_output_file, close_output_file
    use mp, only: proc0
    use species, only: nspec
    use theta_grid, only: ntgrid, theta, bmag, gradpar
    use theta_grid, only: gbdrift, gbdrift0, cvdrift, cvdrift0
    use theta_grid, only: gds2, gds21, gds22
    use kt_grids, only: naky, ntheta0, aky, akx, theta0
    use le_grids, only: nlambda, fcheck, al
    use fields, only: phi, apar, aperp
    use dist_fn, only: finish_intcheck, finish_vortcheck, finish_fieldcheck
    use dist_fn_arrays, only: g
    implicit none
    integer :: ig, ik, it, il, is, unit
    complex, dimension (nlambda,naky,ntheta0,nspec) :: fcheck_f

    if (proc0) then
       call close_output_file (out_unit)
       if (dump_neoclassical_flux) close (dump_neoclassical_flux_unit)
       if (dump_omega) close (dump_omega_unit)
       if (dump_check1) close (dump_check1_unit)
       if (dump_check2) close (dump_check2_unit)

       if (write_eigenfunc) then
          call open_output_file (unit, ".eigenfunc")
          do ik = 1, naky
             do it = 1, ntheta0
                do ig = -ntgrid, ntgrid
                   write (unit, "(9(x,e12.5))") &
                        theta(ig), theta0(ik,it), aky(ik), &
                        phi(ig,ik,it)/phi(0,ik,it), &
                        apar(ig,ik,it)/phi(0,ik,it), &
                        aperp(ig,ik,it)/phi(0,ik,it)
                end do
                write (unit, "()")
             end do
          end do
          call close_output_file (unit)
       end if

       if (write_final_fields) then
          call open_output_file (unit, ".fields")
          do ik = 1, naky
             do it = 1, ntheta0
                do ig = -ntgrid, ntgrid
                   write (unit, "(9(x,e12.5))") &
                        theta(ig), theta0(ik,it), aky(ik), &
                        phi(ig,ik,it), apar(ig,ik,it), aperp(ig,ik,it)
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
                           aky(ik), akx(it), fcheck_f(il,ik,it,is), &
                           sum((al(2:il)-al(1:il-1))*fcheck_f(2:il,ik,it,is))
                   end do
                end do
             end do
          end do
          call close_output_file (unit)
       end if
    end if

    if (write_eik) then
       if (proc0) then
          call open_output_file (unit, ".eik")
          do ig = -ntgrid, ntgrid
             write(unit, "(11(x,e12.5))") &
	       theta(ig), gradpar(ig), bmag(ig), theta(ig), &
               gbdrift(ig), gbdrift0(ig), cvdrift(ig), cvdrift0(ig), &
               gds2(ig), gds21(ig), gds22(ig)	       
          end do
	  call close_output_file (unit)
       end if
    end if

    if (write_intcheck) call finish_intcheck
    if (write_vortcheck) call finish_vortcheck
    if (write_fieldcheck) call finish_fieldcheck
  end subroutine finish_gs2_diagnostics

  subroutine loop_diagnostics (istep, nstep, exit)
    use species, only: nspec
    use theta_grid, only: ntgrid, delthet, bmag, gradpar
    use kt_grids, only: naky, ntheta0, aky, theta0, akx
    use run_parameters, only: delt, woutunits
    use fields, only: phinorm, phinew, aparnew, aperpnew
    use fields, only: kperp, fieldlineavgphi
    use dist_fn, only: flux, intcheck, vortcheck, fieldcheck
    use dist_fn, only: neoclassical_flux, omega0, gamma0
    use mp, only: proc0, broadcast
    use prof, only: prof_entering, prof_leaving
    use constants
    implicit none
    integer, intent (in) :: istep, nstep
    logical, intent (out) :: exit
    complex, dimension (naky,ntheta0) :: omega, omegaavg
    real, dimension (naky,ntheta0) :: phitot, akperp
    complex, dimension (naky,ntheta0,nspec) :: pfluxneo,qfluxneo
    real :: t
    integer :: ik, it, is
    complex :: phiavg, sourcefac

    call prof_entering ("loop_diagnostics")

    exit = .false.

    if (proc0) call get_omegaavg (istep, exit, omegaavg)
    call broadcast (exit)

    call prof_leaving ("loop_diagnostics")

    if (mod(istep,nwrite) /= 0 .and. istep /= 1) return

    call prof_entering ("loop_diagnostics-1")

    if (proc0) then
       omega = omegahist(mod(istep,navg),:,:)
       t = delt*real(istep)
       sourcefac = exp(-zi*omega0*t+gamma0*t)
       call phinorm (phitot)
    end if

    if (proc0) then
       if (print_line) then
          if (print_old_units) then
             do ik = 1, naky
                do it = 1, ntheta0
                   write (unit=*, fmt="('aky=',f5.2, ' th0=',f4.2, &
                          &' om=',2f8.3,' omav=', 2f8.3,' phtot=',e10.4)") &
                        aky(ik), theta0(ik,it), &
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
                   write (unit=*, fmt="('aky=',f5.2, ' th0=',f4.2, &
                          &' om=',2f8.3,' omav=', 2f8.3,' phtot=',e10.4)") &
                        aky(ik), theta0(ik,it), &
                        real(omega(ik,it)*woutunits(ik)), &
                        aimag(omega(ik,it)*woutunits(ik)), &
                        real(omegaavg(ik,it)*woutunits(ik)), &
                        aimag(omegaavg(ik,it)*woutunits(ik)), &
                        phitot(ik,it)
                end do
             end do
          end if
       end if
    end if

    if (write_intcheck) call intcheck
    if (write_vortcheck) call vortcheck (phinew, aparnew, aperpnew)
    if (write_fieldcheck) call fieldcheck (phinew, aparnew, aperpnew)

    call prof_leaving ("loop_diagnostics-1")

    if (.not. (write_any .or. dump_any)) return

    call prof_entering ("loop_diagnostics-2")

    if (write_any_fluxes) then
       call flux (phinew, aparnew, aperpnew, &
            pflux, qheat, vflux, pmflux, qmheat, vmflux)
    end if

    call kperp (akperp)

    if (write_neoclassical_flux .or. dump_neoclassical_flux) then
       call neoclassical_flux (pfluxneo, qfluxneo, istep)
    end if

    if (proc0 .and. write_any) then
       write (unit=out_unit, fmt=*) 'time=', t
       do ik = 1, naky
          do it = 1, ntheta0
             write (out_unit,*) "aky=",aky(ik)," theta0=",theta0(ik,it)

             if (write_line) then
                write (out_unit, "('aky=',f5.2, ' th0=',f4.2, &
                       &' om=',2f8.3,' omav=', 2f8.3,' phtot=',e8.2)") &
                     aky(ik), theta0(ik,it), &
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
                write (out_unit, *) 'kperpnorm=',akperp(ik,it)/aky(ik)
             end if
             if (write_phitot) write (out_unit, *) 'phitot=', phitot(ik,it)

             if (write_fieldline_avg_phi) then
                call fieldlineavgphi (ik, it, phiavg)
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
          end do
       end do
    end if

    if (proc0 .and. dump_any) then
       do ik = 1, naky
          do it = 1, ntheta0
             if (dump_neoclassical_flux) then
                do is = 1, nspec
                   write (dump_neoclassical_flux_unit, "(20(1x,e12.5))") &
                        t, aky(ik), akx(it), real(is), pfluxneo(ik,it,is), &
                        qfluxneo(ik,it,is), &
                        pfluxneo(ik,it,is)/sourcefac, &
                        qfluxneo(ik,it,is)/sourcefac
                end do
             end if
             if (dump_omega) then
                write (dump_omega_unit, "(20(1x,e12.5))") &
                     t, aky(ik), akx(it), omega(ik,it)*woutunits(ik), &
                     omegaavg(ik,it)*woutunits(ik)
             end if
             if (dump_check1) then
                write (dump_check1_unit, "(20(1x,e12.5))") &
                     t, aky(ik), akx(it), &
                     sum(phinew(:ntgrid-1,ik,it)*delthet(:ntgrid-1) &
                         /gradpar(:ntgrid-1)/bmag(:ntgrid-1)) &
                     /sum(delthet(:ntgrid-1) &
                          /gradpar(:ntgrid-1)/bmag(:ntgrid-1)), &
                     sum(phinew(:ntgrid-1,ik,it)*delthet(:ntgrid-1) &
                         /gradpar(:ntgrid-1)/bmag(:ntgrid-1)) &
                     /sum(delthet(:ntgrid-1) &
                          /gradpar(:ntgrid-1)/bmag(:ntgrid-1)) &
                     /sourcefac
             end if
             if (dump_check2) then
                write (dump_check2_unit, "(20(1x,e12.5))") &
                     t, aky(ik), akx(it), &
                     sum(phinew(:ntgrid-1,ik,it)*delthet(:ntgrid-1) &
                         /gradpar(:ntgrid-1)/bmag(:ntgrid-1)) &
                     /sum(delthet(:ntgrid-1) &
                          /gradpar(:ntgrid-1)/bmag(:ntgrid-1))
             end if
          end do
       end do
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

end module gs2_diagnostics
