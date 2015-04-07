!> Contains a lot of old routines for writing stuff out at the end
!! of the simulation that I had no time to upgrade
module diagnostics_final_routines
  implicit none

  private

  public :: ntg_out, init_par_filter, do_write_gs
  public :: do_write_final_antot, do_write_final_moments, do_write_final_db
  public :: do_write_final_epar, do_write_kpar, do_write_final_fields
  public :: do_write_eigenfunc

  integer :: ntg_out

contains

  subroutine init_par_filter
    use theta_grid, only: ntgrid
    use gs2_transforms, only: init_zf
    use kt_grids, only: naky, ntheta0
    use mp, only: proc0
    implicit none
    
    if ( proc0.and.(naky*ntheta0.eq.0)) print *,"WARNING: kt_grids used in init_par_filter before initialised?"    
    call init_zf (ntgrid, ntheta0*naky)
  end subroutine init_par_filter

  subroutine par_spectrum(an, an2)
    use gs2_transforms, only: kz_spectrum
    use theta_grid, only: ntgrid
    implicit none
    complex, dimension(:,:,:), intent(in) :: an
    complex, dimension(:,:,:), intent(out) :: an2
    real :: scale
    
    call kz_spectrum (an, an2)
    scale = 1./real(4*ntgrid**2)
    an2 = an2*scale
  end subroutine par_spectrum

  subroutine do_write_eigenfunc(gnostics, phi0)
    use mp, only: proc0
    use file_utils, only: open_output_file, close_output_file
    use fields_arrays, only: phi, apar, bpar
    use kt_grids, only: ntheta0, naky, theta0, aky
    use theta_grid, only: theta
    use gs2_io, only: nc_eigenfunc
    use dist_fn, only: def_parity, even
    use run_parameters, only: fapar
    use diagnostics_config, only: diagnostics_type
    implicit none
    type(diagnostics_type), intent(in) :: gnostics
    complex, dimension (:,:), intent(inout) :: phi0
    integer :: it, ik, ig, unit
    
    if(.not.proc0) return
    
    !Should this actually use igomega instead of 0?
    !What if fphi==0? --> See where statements below
    phi0 = phi(0,:,:)

    !This looks like a hack for the case where we know we've forced phi(theta=0) to be 0
    !this could probably be better addressed by the use of igomega above
    if (def_parity .and. fapar > 0 .and. (.not. even)) phi0 = apar(0, :, :)
    
    !Address locations where phi0=0 by using next point
    where (abs(phi0) < 10.0*epsilon(0.0)) 
       phi0 = phi(1,:,:)/(theta(1)-theta(0))
    end where
    
    !Check again if any locations are 0, this could be true if fphi (or fapar)
    !is zero.
    where (abs(phi0) < 10.0*epsilon(0.0)) 
       phi0 = 1.0
    end where
    
    !Do ascii output
    if (gnostics%write_ascii) then
       call open_output_file (unit, ".eigenfunc")
       do ik = 1, naky
          do it = 1, ntheta0
             do ig = -ntg_out, ntg_out
                write (unit, "(9(1x,e12.5))") &
                     theta(ig), theta0(it,ik), aky(ik), &
                     phi(ig,it,ik)/phi0(it,ik), &
                     apar(ig,it,ik)/phi0(it,ik), &
                     bpar(ig,it,ik)/phi0(it,ik)
             end do
             write (unit, "()")
          end do
       end do
       call close_output_file (unit)
    end if
    
    !Do netcdf output
    call nc_eigenfunc (phi0)
  end subroutine do_write_eigenfunc
  
  subroutine do_write_final_fields(gnostics)
    use mp, only: proc0
    use file_utils, only: open_output_file, close_output_file
    use kt_grids, only: naky, ntheta0, aky, akx, theta0
    use theta_grid, only: theta
    use fields_arrays, only: phi, apar, bpar
    use diagnostics_config, only: diagnostics_type
    implicit none
    type(diagnostics_type), intent(in) :: gnostics
    integer :: ik, it, ig, unit

    if(.not.proc0) return
    
    !Do ascii output
    if (gnostics%write_ascii) then
       call open_output_file (unit, ".fields")
       do ik = 1, naky
          do it = 1, ntheta0
             do ig = -ntg_out, ntg_out
                write (unit, "(15(1x,e12.5))") &
                     theta(ig), aky(ik), akx(it), &
                     phi(ig,it,ik), &
                     apar(ig,it,ik), &
                     bpar(ig,it,ik), &
                     theta(ig) - theta0(it,ik), &
                     abs(phi(ig,it,ik))
             end do
             write (unit, "()")
          end do
       end do
       call close_output_file (unit)
    end if
  end subroutine do_write_final_fields

  subroutine do_write_kpar(gnostics)
    use mp, only: proc0
    use file_utils, only: open_output_file, close_output_file
    use theta_grid, only: ntgrid, gradpar, nperiod
    use kt_grids, only: naky, ntheta0, aky, akx
    use run_parameters, only: fphi,fapar,fbpar
    use fields_arrays, only: phi, apar, bpar
    use diagnostics_config, only: diagnostics_type
    implicit none
    type(diagnostics_type), intent(in) :: gnostics
    complex, dimension (:,:,:), allocatable :: phi2, apar2, bpar2
    real, dimension (2*ntgrid) :: kpar
    integer :: ig, ik, it, unit

    if(.not.proc0) return
    
    allocate (phi2(-ntgrid:ntgrid,ntheta0,naky)) 
    allocate (apar2(-ntgrid:ntgrid,ntheta0,naky))
    allocate (bpar2(-ntgrid:ntgrid,ntheta0,naky))

    if (fphi > epsilon(0.0)) then
       call par_spectrum(phi, phi2)
    else
       phi2=0.
    end if
    if (fapar > epsilon(0.0)) then
       call par_spectrum(apar, apar2)
    else
       apar2=0.
    endif
    if (fbpar > epsilon(0.0)) then
       call par_spectrum(bpar, bpar2)
    else
       bpar2=0.
    endif
    
    !Do ascii output
    if(gnostics%write_ascii)then
       call open_output_file (unit, ".kpar")
       do ig = 1, ntgrid
          kpar(ig) = (ig-1)*gradpar(ig)/real(2*nperiod-1)
          kpar(2*ntgrid-ig+1)=-(ig)*gradpar(ig)/real(2*nperiod-1)
       end do
       do ik = 1, naky
          do it = 1, ntheta0
             do ig = ntgrid+1,2*ntgrid
                write (unit, "(9(1x,e12.5))") &
                     kpar(ig), aky(ik), akx(it), &
                     phi2(ig-ntgrid-1,it,ik), &
                     apar2(ig-ntgrid-1,it,ik), &
                     bpar2(ig-ntgrid-1,it,ik)                        
             end do
             do ig = 1, ntgrid
                write (unit, "(9(1x,e12.5))") &
                     kpar(ig), aky(ik), akx(it), &
                     phi2(ig-ntgrid-1,it,ik), &
                     apar2(ig-ntgrid-1,it,ik), &
                     bpar2(ig-ntgrid-1,it,ik)
             end do
             write (unit, "()")
          end do
       end do
       call close_output_file (unit)
    endif
    
    !Currently no netcdf output for this diagnostic.
    deallocate (phi2, apar2, bpar2)
  end subroutine do_write_kpar

  subroutine do_write_final_epar(gnostics)
    use mp, only: proc0
    use file_utils, only: open_output_file, close_output_file
    use kt_grids, only: naky, ntheta0, theta0, aky, akx
    use theta_grid, only: theta, ntgrid
    use dist_fn, only: get_epar
    use fields_arrays, only: phi, apar, phinew, aparnew
    use gs2_io, only: nc_final_epar
    use diagnostics_config, only: diagnostics_type
    implicit none
    type(diagnostics_type), intent(in) :: gnostics
    complex, dimension (:,:,:), allocatable :: epar
    integer :: ik, it, ig, unit
    
    if(.not.proc0) return
    
    allocate (epar(-ntgrid:ntgrid,ntheta0,naky))
    
    !Calculate
    call get_epar (phi, apar, phinew, aparnew, epar)
    
    !Write to ascii
    if (gnostics%write_ascii) then
       call open_output_file (unit, ".epar")
       do ik = 1, naky
          do it = 1, ntheta0
             do ig = -ntg_out, ntg_out-1
                write (unit, "(6(1x,e12.5))") &
                     theta(ig), aky(ik), akx(it), &
                     epar(ig,it,ik), &
                     theta(ig) - theta0(it,ik)
             end do
             write (unit, "()")
          end do
       end do
       call close_output_file (unit)
    end if
    
    !Do netcdf output
    call nc_final_epar (epar)
    
    deallocate (epar)
  end subroutine do_write_final_epar

  subroutine do_write_final_db(gnostics)
    use file_utils, only: open_output_file, close_output_file
    use mp, only: proc0
    use theta_grid, only: ntgrid, gradpar, delthet, bmag, theta
    use kt_grids, only: ntheta0, naky, aky, akx
    use fields_arrays, only: phinew, aparnew, apar
    use gs2_time, only: code_dt
    use diagnostics_config, only: diagnostics_type
    implicit none
    type(diagnostics_type), intent(in) :: gnostics
    complex, dimension (-ntgrid:ntgrid, ntheta0, naky) :: db
    complex, dimension (ntheta0, naky) :: dbfac
    integer :: ik, it, ig, unit

    if(.not.proc0) return
    !Only write to ascii for now so if not gnostics%write_ascii return
    if(.not.gnostics%write_ascii) return
    
    !Calculate db
    do ik = 1, naky
       do it = 1, ntheta0
          dbfac(it,ik) = 1./sum(delthet/bmag/gradpar)/maxval(abs(phinew(:,it,ik)),1) &
               * abs(log(aparnew(1,it,ik)/apar(1,it,ik)))/code_dt
          ig = -ntg_out
          db(ig, it, ik) = aparnew(ig,it,ik)*delthet(ig)/bmag(ig)/gradpar(ig)*dbfac(it,ik)
          do ig = -ntg_out+1, ntg_out-1
             db(ig, it, ik) = db(ig-1, it, ik) + aparnew(ig,it,ik)*delthet(ig)/bmag(ig)/gradpar(ig)*dbfac(it,ik)
          end do
       end do
    end do
    
    if (gnostics%write_ascii) then
       call open_output_file (unit, ".db")
       do ik = 1, naky
          do it = 1, ntheta0
             do ig = -ntg_out, ntg_out-1
                write (unit, "(5(1x,e12.5))") &
                     theta(ig), aky(ik), akx(it), real(db(ig, it,ik)), aimag(db(ig, it, ik))
             end do
             write (unit, "()")
          end do
       end do
       call close_output_file (unit)
    end if
    
    !No netcdf output for this diagnostic yet
  end subroutine do_write_final_db

  subroutine do_write_final_moments(gnostics, phi0_in)
    use theta_grid, only: ntgrid, theta
    use kt_grids, only: ntheta0, naky, aky, akx, theta0
    use species, only: nspec, spec
    use dist_fn, only: getmoms
    use file_utils, only: open_output_file, close_output_file
    use mp, only: proc0
    use fields_arrays, only: phinew, bparnew
    use gs2_io, only: nc_final_moments
    use diagnostics_config, only: diagnostics_type
    use volume_averages, only: average_theta
    implicit none
    type(diagnostics_type), intent(in) :: gnostics
    complex, dimension (:,:,:,:), allocatable :: ntot, density, upar, tpar, tperp
    complex, dimension (:,:,:,:), allocatable :: qparflux, pperpj1, qpperpj1
    complex, dimension (:, :), intent(in) :: phi0_in !Phase calculated
    complex, dimension (:,:), allocatable :: phi0
    complex :: phi_tmp, ntot_tmp, density_tmp, upar_tmp, tpar_tmp, tperp_tmp
    real, dimension (:, :), allocatable :: phi02
    integer :: is, ik, it, ig, unit
    
    !Make storage
    allocate (ntot(-ntgrid:ntgrid,ntheta0,naky,nspec))
    allocate (density(-ntgrid:ntgrid,ntheta0,naky,nspec))
    allocate (upar(-ntgrid:ntgrid,ntheta0,naky,nspec))
    allocate (tpar(-ntgrid:ntgrid,ntheta0,naky,nspec))
    allocate (tperp(-ntgrid:ntgrid,ntheta0,naky,nspec))
    allocate (qparflux(-ntgrid:ntgrid,ntheta0,naky,nspec))
    allocate (pperpj1(-ntgrid:ntgrid,ntheta0,naky,nspec))
    allocate (qpperpj1(-ntgrid:ntgrid,ntheta0,naky,nspec))

    !Calculate moments
    call getmoms (phinew, bparnew, ntot, density, upar, tpar, tperp, qparflux, pperpj1, qpperpj1)
    
    if (proc0) then
       if (gnostics%write_ascii) then
          !Setup the phase factor
          allocate(phi0(ntheta0,naky),phi02(ntheta0,naky))
          if(.not.gnostics%write_eigenfunc) then
             phi0 = 1.
          else
             phi0 = phi0_in
          endif
          phi02=real(phi0*conjg(phi0))
          
          !Write out the moments normalised by phase factor
          call open_output_file (unit, ".moments")          
          do is  = 1, nspec
             do ik = 1, naky
                do it = 1, ntheta0
                   do ig = -ntg_out, ntg_out
                      write (unit, "(15(1x,e12.5))") &
                           theta(ig), aky(ik), akx(it), &
                           ntot(ig,it,ik,is)/phi0(it,ik), &
                           density(ig,it,ik,is)/phi0(it,ik), &
                           upar(ig,it,ik,is)/phi0(it,ik), &
                           tpar(ig,it,ik,is)/phi0(it,ik), &
                           tperp(ig,it,ik,is)/phi0(it,ik), &
                           theta(ig) - theta0(it,ik), &
                           real(is)
                   end do
                   write (unit, "()")
                end do
             end do
          end do
          call close_output_file (unit)          
          
          !Write out magnitude of moments
          call open_output_file (unit, ".mom2")
          do is  = 1, nspec
             do ik = 1, naky
                do it = 1, ntheta0
                   do ig = -ntg_out, ntg_out
                      write (unit, "(15(1x,e12.5))") &
                           theta(ig), aky(ik), akx(it), &
                           real(ntot(ig,it,ik,is)*conjg(ntot(ig,it,ik,is)))/phi02(it,ik), &
                           real(density(ig,it,ik,is)*conjg(density(ig,it,ik,is)))/phi02(it,ik), &
                           real(upar(ig,it,ik,is)*conjg(upar(ig,it,ik,is)))/phi02(it,ik), &
                           real(tpar(ig,it,ik,is)*conjg(tpar(ig,it,ik,is)))/phi02(it,ik), &
                           real(tperp(ig,it,ik,is)*conjg(tperp(ig,it,ik,is)))/phi02(it,ik), &
                           theta(ig) - theta0(it,ik), &
                           real(is)
                   end do
                   write (unit, "()")
                end do
             end do
          end do
          call close_output_file (unit)          
          
          !Write out the field line averaged moments
          call open_output_file (unit, ".amoments")
          write (unit,*) 'type    kx     re(phi)    im(phi)    re(ntot)   im(ntot)   ',&
               &'re(dens)   im(dens)   re(upar)   im(upar)   re(tpar)',&
               &'   im(tpar)   re(tperp)  im(tperp)'
          
          do is  = 1, nspec
             do it = 2, ntheta0/2+1
                call average_theta(phinew(-ntg_out:ntg_out,it,1),phi_tmp)
                call average_theta(ntot(-ntg_out:ntg_out,it,1,is),ntot_tmp)
                call average_theta(density(-ntg_out:ntg_out,it,1,is),density_tmp)
                call average_theta(upar(-ntg_out:ntg_out,it,1,is),upar_tmp)
                call average_theta(tpar(-ntg_out:ntg_out,it,1,is),tpar_tmp)
                call average_theta(tperp(-ntg_out:ntg_out,it,1,is),tperp_tmp)
                
                write (unit, "(i2,14(1x,e10.3))") spec(is)%type, akx(it), &
                     phi_tmp, ntot_tmp, density_tmp, upar_tmp, tpar_tmp, tperp_tmp
             end do
          end do
          call close_output_file (unit)          
          
          deallocate(phi0,phi02)
       end if
       
       !Do the netcdf output -- Note we don't have the phi0 phase factor here
       call nc_final_moments (ntot, density, upar, tpar, tperp, qparflux, pperpj1, qpperpj1)
    end if
    
    deallocate (ntot, density, upar, tpar, tperp, qparflux, pperpj1, qpperpj1)
  end subroutine do_write_final_moments
  
  subroutine do_write_final_antot(gnostics)
    use file_utils, only: open_output_file, close_output_file
    use mp, only: proc0
    use kt_grids, only: ntheta0, naky, theta0, aky
    use theta_grid, only: theta, ntgrid
    use gs2_io, only: nc_final_an
    use dist_fn, only: getan
    use diagnostics_config, only: diagnostics_type
    implicit none
    type(diagnostics_type), intent(in) :: gnostics
    complex, dimension (:,:,:), allocatable :: antot, antota, antotp
    integer :: ik, it, ig, unit
    
    allocate ( antot(-ntgrid:ntgrid,ntheta0,naky))
    allocate (antota(-ntgrid:ntgrid,ntheta0,naky))
    allocate (antotp(-ntgrid:ntgrid,ntheta0,naky))
    call getan (antot, antota, antotp)

    if (proc0) then
       !Do ascii output
       if (gnostics%write_ascii) then
          call open_output_file (unit, ".antot")
          do ik = 1, naky
             do it = 1, ntheta0
                do ig = -ntg_out, ntg_out
                   write (unit, "(10(1x,e12.5))") &
                        theta(ig), theta0(it,ik), aky(ik), &
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
       
       !Do netcdf output
       call nc_final_an (antot, antota, antotp)
    end if
    
    deallocate (antot, antota, antotp)
  end subroutine do_write_final_antot

  subroutine do_write_gs(gnostics)
    use file_utils, only: open_output_file, close_output_file
    use mp, only: proc0, sum_reduce, iproc, nproc
    use gs2_transforms, only: transform2, inverse2
    use fields_arrays, only: apar, phi
    use kt_grids, only: naky, ntheta0, akx, aky, nx, ny
    use theta_grid, only: ntgrid, delthet, jacob, gradpar, nperiod
    use constants, only: pi, zi
    use splines, only: fitp_surf1, fitp_surf2
    use diagnostics_config, only: diagnostics_type
    implicit none
    type(diagnostics_type), intent(in) :: gnostics
    real, dimension (:), allocatable :: total
    real, dimension (:,:,:), allocatable :: bxf, byf, vxf, vyf, bxfsavg, byfsavg
    real, dimension (:,:,:), allocatable :: bxfs, byfs, vxfs, vyfs, rvx, rvy, rx, ry
    real, dimension (:), allocatable :: xx4, yy4, dz
    real, dimension (:,:), allocatable :: bxs, bys, vxs, vys
    real, dimension (:,:), allocatable :: bxsavg, bysavg
    complex, dimension (:,:,:), allocatable :: bx, by, vx, vy, vx2, vy2
    real, dimension (:), allocatable :: stemp, zx1, zxm, zy1, zyn, xx, yy
    real :: zxy11, zxym1, zxy1n, zxymn, rxt, ryt, bxt, byt, L_x, L_y
    real, dimension (2*ntgrid) :: kpar
    integer :: nnx, nny, nnx4, nny4, ulim, llim, iblock, i, unit
    integer :: ik, it, ig, ierr
    
    !Note no netcdf output so if not write_ascii then exit now
    if(.not.gnostics%write_ascii) return
    
    nny = 2*ny
    nnx = 2*nx
    nnx4 = nnx+4
    nny4 = nny+4
    
    allocate (bx(-ntgrid:ntgrid,ntheta0,naky))
    allocate (by(-ntgrid:ntgrid,ntheta0,naky))
    allocate (vx(-ntgrid:ntgrid,ntheta0,naky))
    allocate (vy(-ntgrid:ntgrid,ntheta0,naky))

    do ik=1,naky
       do it=1,ntheta0
          do ig=-ntgrid, ntgrid
             bx(ig,it,ik) =  zi*aky(ik)*apar(ig,it,ik)
             by(ig,it,ik) = -zi*akx(it)*apar(ig,it,ik)
             vx(ig,it,ik) = -zi*aky(ik)*phi(ig,it,ik)
             vy(ig,it,ik) =  zi*akx(it)*phi(ig,it,ik)
          end do
       end do
    end do

    allocate (bxf(nnx,nny,-ntgrid:ntgrid))
    allocate (byf(nnx,nny,-ntgrid:ntgrid))
    allocate (vxf(nnx,nny,-ntgrid:ntgrid))
    allocate (vyf(nnx,nny,-ntgrid:ntgrid))

    call transform2 (bx, bxf, nny, nnx)
    call transform2 (by, byf, nny, nnx)
    call transform2 (vx, vxf, nny, nnx)
    call transform2 (vy, vyf, nny, nnx)

    ! fields come out as (x, y, z)

    deallocate (bx, by)

    allocate (   bxfs(nnx4, nny4, -ntgrid:ntgrid))
    allocate (   byfs(nnx4, nny4, -ntgrid:ntgrid))
    allocate (bxfsavg(nnx4, nny4, -ntgrid:ntgrid))
    allocate (byfsavg(nnx4, nny4, -ntgrid:ntgrid))
    allocate (   vxfs(nnx4, nny4, -ntgrid:ntgrid))
    allocate (   vyfs(nnx4, nny4, -ntgrid:ntgrid))

    do ig=-ntgrid,ntgrid
       do ik=1,2
          do it=3,nnx4-2
             bxfs(it,ik,ig) = bxf(it-2,nny-2+ik,ig)
             byfs(it,ik,ig) = byf(it-2,nny-2+ik,ig)
             vxfs(it,ik,ig) = vxf(it-2,nny-2+ik,ig)
             vyfs(it,ik,ig) = vyf(it-2,nny-2+ik,ig)

             bxfs(it,nny4-2+ik,ig) = bxf(it-2,ik,ig)
             byfs(it,nny4-2+ik,ig) = byf(it-2,ik,ig)
             vxfs(it,nny4-2+ik,ig) = vxf(it-2,ik,ig)
             vyfs(it,nny4-2+ik,ig) = vyf(it-2,ik,ig)
          end do
       end do
       do ik=3,nny4-2
          do it=3,nnx4-2
             bxfs(it,ik,ig) = bxf(it-2,ik-2,ig)
             byfs(it,ik,ig) = byf(it-2,ik-2,ig)
             vxfs(it,ik,ig) = vxf(it-2,ik-2,ig)
             vyfs(it,ik,ig) = vyf(it-2,ik-2,ig)
          end do
       end do
       do ik=1,nny4
          do it=1,2
             bxfs(it,ik,ig) = bxfs(nnx4-4+it,ik,ig)
             byfs(it,ik,ig) = byfs(nnx4-4+it,ik,ig)
             vxfs(it,ik,ig) = vxfs(nnx4-4+it,ik,ig)
             vyfs(it,ik,ig) = vyfs(nnx4-4+it,ik,ig)

             bxfs(nnx4-2+it,ik,ig) = bxfs(it+2,ik,ig)
             byfs(nnx4-2+it,ik,ig) = byfs(it+2,ik,ig)
             vxfs(nnx4-2+it,ik,ig) = vxfs(it+2,ik,ig)
             vyfs(nnx4-2+it,ik,ig) = vyfs(it+2,ik,ig)
          end do
       end do
    end do

    deallocate (vxf, vyf)

    allocate (xx4(nnx4), xx(nnx))
    allocate (yy4(nny4), yy(nny))

    L_x = 2.0*pi/akx(2)
    L_y = 2.0*pi/aky(2)

    do it = 1, nnx
       xx4(it+2) = real(it-1)*L_x/real(nnx)
       xx(it) = real(it-1)*L_x/real(nnx)
    end do
    do it=1,2
       xx4(it) = xx4(nnx4-4+it)-L_x
       xx4(nnx4-2+it) = xx4(it+2)+L_x
    end do

    do ik = 1, nny
       yy4(ik+2) = real(ik-1)*L_y/real(nny)
       yy(ik)    = real(ik-1)*L_y/real(nny)
    end do
    do ik=1,2
       yy4(ik) = yy4(nny4-4+ik)-L_y
       yy4(nny4-2+ik) = yy4(ik+2)+L_y
    end do

    allocate (dz(-ntgrid:ntgrid))
    dz = delthet*jacob

    allocate (   bxs(3*nnx4*nny4,-ntgrid:ntgrid)) ; bxs = 0.
    allocate (   bys(3*nnx4*nny4,-ntgrid:ntgrid)) ; bys = 0.
    allocate (   vxs(3*nnx4*nny4,-ntgrid:ntgrid)) ; vxs = 0.
    allocate (   vys(3*nnx4*nny4,-ntgrid:ntgrid)) ; vys = 0.

    allocate (bxsavg(3*nnx4*nny4,-ntgrid:ntgrid))
    allocate (bysavg(3*nnx4*nny4,-ntgrid:ntgrid))

    allocate (stemp(nnx4+2*nny4))
    allocate (zx1(nny4), zxm(nny4), zy1(nnx4), zyn(nnx4))

    do ig=-ntgrid, ntgrid
       call fitp_surf1(nnx4, nny4, xx, yy, bxfs(:,:,ig), &
            nnx4, zx1, zxm, zy1, zyn, zxy11, zxym1, zxy1n, zxymn, &
            255, bxs(:,ig), stemp, 1., ierr)

       call fitp_surf1(nnx4, nny4, xx, yy, byfs(:,:,ig), &
            nnx4, zx1, zxm, zy1, zyn, zxy11, zxym1, zxy1n, zxymn, &
            255, bys(:,ig), stemp, 1., ierr)

       call fitp_surf1(nnx4, nny4, xx, yy, vxfs(:,:,ig), &
            nnx4, zx1, zxm, zy1, zyn, zxy11, zxym1, zxy1n, zxymn, &
            255, vxs(:,ig), stemp, 1., ierr)

       call fitp_surf1(nnx4, nny4, xx, yy, vyfs(:,:,ig), &
            nnx4, zx1, zxm, zy1, zyn, zxy11, zxym1, zxy1n, zxymn, &
            255, vys(:,ig), stemp, 1., ierr)
    end do

    deallocate (zx1, zxm, zy1, zyn, stemp)

    do ig=-ntgrid, ntgrid-1
       bxsavg(:,ig) = 0.5*(bxs(:,ig)+bxs(:,ig+1))
       bysavg(:,ig) = 0.5*(bys(:,ig)+bys(:,ig+1))

       bxfsavg(:,:,ig) = 0.5*(bxfs(:,:,ig)+bxfs(:,:,ig+1))
       byfsavg(:,:,ig) = 0.5*(byfs(:,:,ig)+byfs(:,:,ig+1))
    end do

    ! now, integrate to find a field line

    allocate ( rx(nnx,nny,-ntgrid:ntgrid))
    allocate ( ry(nnx,nny,-ntgrid:ntgrid))
    allocate (rvx(nnx,nny,-ntgrid:ntgrid)) ; rvx = 0.
    allocate (rvy(nnx,nny,-ntgrid:ntgrid)) ; rvy = 0.

    do ik=1,nny
       do it=1,nnx
          rx(it,ik,-ntgrid) = xx(it)
          ry(it,ik,-ntgrid) = yy(ik)
       end do
    end do

    !Should these not come from layouts object?
    iblock = (nnx*nny-1)/nproc + 1
    llim = 1 + iblock * iproc  
    ulim = min(nnx*nny, llim+iblock-1)

    do i=llim, ulim
       it = 1 + mod(i-1, nnx)
       ik = 1 + mod((i-1)/nnx, nny)

       ig = -ntgrid

       rxt = rx(it,ik,ig)
       ryt = ry(it,ik,ig)

       rvx(it,ik,ig) = fitp_surf2(rxt, ryt, nnx4, nny4, xx4, yy4, vxfs(:,:,ig), nnx4, vxs(:,ig), 1.)
       rvy(it,ik,ig) = fitp_surf2(rxt, ryt, nnx4, nny4, xx4, yy4, vyfs(:,:,ig), nnx4, vys(:,ig), 1.)

       do ig=-ntgrid,ntgrid-1

          bxt = fitp_surf2(rxt, ryt, nnx4, nny4, xx4, yy4, bxfs(:,:,ig), nnx4, bxs(:,ig), 1.)
          byt = fitp_surf2(rxt, ryt, nnx4, nny4, xx4, yy4, byfs(:,:,ig), nnx4, bys(:,ig), 1.)

          rxt = rx(it,ik,ig) + 0.5*dz(ig)*bxt  
          ryt = ry(it,ik,ig) + 0.5*dz(ig)*byt  

          if (rxt > L_x) rxt = rxt - L_x
          if (ryt > L_y) ryt = ryt - L_y

          if (rxt < 0.) rxt = rxt + L_x
          if (ryt < 0.) ryt = ryt + L_y

          bxt = fitp_surf2(rxt, ryt, nnx4, nny4, xx4, yy4, bxfsavg(:,:,ig), nnx4, bxsavg(:,ig), 1.)
          byt = fitp_surf2(rxt, ryt, nnx4, nny4, xx4, yy4, byfsavg(:,:,ig), nnx4, bysavg(:,ig), 1.)

          rxt = rx(it,ik,ig) + dz(ig)*bxt  
          ryt = ry(it,ik,ig) + dz(ig)*byt  

          if (rxt > L_x) rxt = rxt - L_x
          if (ryt > L_y) ryt = ryt - L_y

          if (rxt < 0.) rxt = rxt + L_x
          if (ryt < 0.) ryt = ryt + L_y

          rx(it,ik,ig+1) = rxt
          ry(it,ik,ig+1) = ryt

          rvx(it,ik,ig+1) = fitp_surf2(rxt, ryt, nnx4, nny4, xx4, yy4, vxfs(:,:,ig+1), nnx4, vxs(:,ig+1), 1.)
          rvy(it,ik,ig+1) = fitp_surf2(rxt, ryt, nnx4, nny4, xx4, yy4, vyfs(:,:,ig+1), nnx4, vys(:,ig+1), 1.)
       end do
    end do

    deallocate (bxfs, byfs, bxfsavg, byfsavg, vxfs, vyfs)
    deallocate (rx, ry, bxs, bys, vxs, vys, bxsavg, bysavg)

    allocate (total(2*nnx*nny*(2*ntgrid+1)))

    i=1
    do ig=-ntgrid,ntgrid
       do ik=1,nny
          do it=1,nnx
             total(i) = rvx(it,ik,ig)
             total(i+1) = rvy(it,ik,ig)
             i = i + 2
          end do
       end do
    end do

    call sum_reduce(total, 0)

    i=1
    do ig=-ntgrid,ntgrid
       do ik=1,nny
          do it=1,nnx
             rvx(it,ik,ig) = total(i)
             rvy(it,ik,ig) = total(i+1)
             i = i + 2
          end do
       end do
    end do

    if (proc0) then
       call inverse2 (rvx, vx, nny, nnx)
       call inverse2 (rvy, vy, nny, nnx)

       allocate (vx2(-ntgrid:ntgrid,ntheta0,naky))
       allocate (vy2(-ntgrid:ntgrid,ntheta0,naky))

       call par_spectrum (vx, vx2)
       call par_spectrum (vy, vy2)

       call open_output_file (unit, ".gs")
       do ig = 1, ntgrid
          kpar(ig) = (ig-1)*gradpar(ig)/real(2*nperiod-1)
          kpar(2*ntgrid-ig+1)=-(ig)*gradpar(ig)/real(2*nperiod-1)
       end do
       do ik = 1, naky
          do it = 1, ntheta0
             do ig = ntgrid+1,2*ntgrid
                write (unit, "(9(1x,e12.5))") &
                     kpar(ig), aky(ik), akx(it), &
                     real(vx2(ig-ntgrid-1,it,ik)), &
                     real(vy2(ig-ntgrid-1,it,ik))
             end do
             do ig = 1, ntgrid
                write (unit, "(9(1x,e12.5))") &
                     kpar(ig), aky(ik), akx(it), &
                     real(vx2(ig-ntgrid-1,it,ik)), &
                     real(vy2(ig-ntgrid-1,it,ik))
             end do
             write (unit, "()")
          end do
       end do
       call close_output_file (unit)
       deallocate (vx2, vy2)
    end if

    if (allocated(bxf)) deallocate (bxf, byf, xx4, xx, yy4, yy, dz, total)
    deallocate (vx, vy, rvx, rvy)
  end subroutine do_write_gs
end module diagnostics_final_routines
