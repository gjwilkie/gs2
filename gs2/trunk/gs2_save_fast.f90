module gs2_save

  implicit none

  public :: gs2_restore, gs2_save_for_restart
  public :: init_save, init_dt, init_tstart

  interface gs2_restore
     module procedure gs2_restore_many, gs2_restore_one
  end interface
  
  private
  character(300), save :: restart_file

  double precision, allocatable, dimension(:,:,:) :: tmpr, tmpi, ftmpr, ftmpi
  integer :: ncid, thetaid, signid, gloid, kyid, kxid
  integer :: phir_id, phii_id, aparr_id, apari_id, aperpr_id, aperpi_id
  integer :: delt0id, t0id, gr_id, gi_id

contains

  subroutine gs2_save_for_restart (g, t0, delt0, istatus, fphi, fapar, faperp, exit_in)
    use theta_grid, only: ntgrid
! Must include g_layout_type here to avoid obscure bomb while compiling
! gs2_diagnostics.f90 (which uses this module) with the Compaq F90 compiler:
    use gs2_layouts, only: g_lo, g_layout_type
    use mp, only: nproc, iproc, proc0
    use fields_arrays, only: phinew, aparnew, aperpnew
    use kt_grids, only: naky, ntheta0
    use file_utils, only: error_unit
    implicit none
    character (305) :: file_proc
    character (5) :: suffix
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g
    real, intent (in) :: t0, delt0
    real, intent (in) :: fphi, fapar, faperp
    integer, intent (out) :: istatus
    logical, intent (in), optional :: exit_in
    include 'netcdf.inc'
    double precision :: tmp1
    integer :: i, n_elements
    integer :: th, h, t, u, ierr
    logical :: exit

    if (present(exit_in)) then
       exit = exit_in
    else
       exit = .false.
    end if

    n_elements = g_lo%ulim_proc-g_lo%llim_proc+1
    if (n_elements <= 0) return

    file_proc = trim(restart_file)
    
    if (nproc >= 10000) then
       ierr = error_unit()
       if (proc0) write(ierr,*) 'Too many procs for i/o to work right!'
    else
       th = iproc / 1000
       h = (iproc - th * 1000) / 100
       t = (iproc - th * 1000 - h * 100) / 10
       u = (iproc - th * 1000 - h * 100 - t * 10)
       suffix = '.'//achar(48+th)//achar(48+h)//achar(48+t)//achar(48+u)
       file_proc = trim(trim(file_proc)//suffix)
    endif
    
    istatus = nf_create (file_proc, 0, ncid)
    if (istatus /= 0) then
       ierr = error_unit()
       write(ierr,*) "nf_create error: ", nf_strerror(istatus)
       goto 1
    end if
    
    if (n_elements > 0) then
       istatus = nf_def_dim (ncid, "theta", 2*ntgrid+1, thetaid)
       if (istatus /= 0) then
          ierr = error_unit()
          write(ierr,*) "nf_def_dim theta error: ", nf_strerror(istatus)
          goto 1
       end if
       
       istatus = nf_def_dim (ncid, "sign", 2, signid)
       if (istatus /= 0) then
          ierr = error_unit()
          write(ierr,*) "nf_def_dim sign error: ", nf_strerror(istatus)
          goto 1
       end if
       
       istatus = nf_def_dim (ncid, "glo", n_elements, gloid)
       if (istatus /= 0) then
          ierr = error_unit()
          write(ierr,*) "nf_def_dim glo error: ", nf_strerror(istatus)
          goto 1
       end if
       
       istatus = nf_def_dim (ncid, "aky", naky, kyid)
       if (istatus /= 0) then
          ierr = error_unit()
          write(ierr,*) "nf_def_dim aky error: ", nf_strerror(istatus)
          goto 1
       end if
       
       istatus = nf_def_dim (ncid, "akx", ntheta0, kxid)
       if (istatus /= 0) then
          ierr = error_unit()
          write(ierr,*) "nf_def_dim akx error: ", nf_strerror(istatus)
          goto 1
       end if
    end if
    
    istatus = nf_def_var (ncid, "t0", NF_DOUBLE, 0, 0, t0id)
    if (istatus /= 0) then
       ierr = error_unit()
       write(ierr,*) "nf_def_var t0 error: ", nf_strerror(istatus)
       goto 1
    end if
    
    istatus = nf_def_var (ncid, "delt0", NF_DOUBLE, 0, 0, delt0id)
    if (istatus /= 0) then
       ierr = error_unit()
       write(ierr,*) "nf_def_var delt0 error: ", nf_strerror(istatus)
       goto 1
    end if
    
    if (n_elements > 0) then
       istatus = nf_def_var (ncid, "gr", NF_DOUBLE, 3, &
            (/ thetaid, signid, gloid /), gr_id)
       if (istatus /= 0) then
          ierr = error_unit()
          write(ierr,*) "nf_def_var g error: ", nf_strerror(istatus)
          goto 1
       end if
       
       istatus = nf_def_var (ncid, "gi", NF_DOUBLE, 3, &
            (/ thetaid, signid, gloid /), gi_id)
       if (istatus /= 0) then
          ierr = error_unit()
          write(ierr,*) "nf_def_var g error: ", nf_strerror(istatus)
          goto 1
       end if
       
       if (fphi > epsilon(0.)) then
          istatus = nf_def_var (ncid, "phi_r", NF_DOUBLE, 3, &
               (/ thetaid, kyid, kxid /), phir_id)
          if (istatus /= 0) then
             ierr = error_unit()
             write(ierr,*) "nf_def_var phi error: ", nf_strerror(istatus)
             goto 1
          end if
          
          istatus = nf_def_var (ncid, "phi_i", NF_DOUBLE, 3, &
               (/ thetaid, kyid, kxid /), phii_id)
          if (istatus /= 0) then
             ierr = error_unit()
             write(ierr,*) "nf_def_var phi error: ", nf_strerror(istatus)
             goto 1
          end if
       end if
       
       if (fapar > epsilon(0.)) then
          istatus = nf_def_var (ncid, "apar_r", NF_DOUBLE, 3, &
               (/ thetaid, kyid, kxid /), aparr_id)
          if (istatus /= 0) then
             ierr = error_unit()
             write(ierr,*) "nf_def_var apar error: ", nf_strerror(istatus)
             goto 1
          end if
          
          istatus = nf_def_var (ncid, "apar_i", NF_DOUBLE, 3, &
               (/ thetaid, kyid, kxid /), apari_id)
          if (istatus /= 0) then
             ierr = error_unit()
             write(ierr,*) "nf_def_var apar error: ", nf_strerror(istatus)
             goto 1
          end if
       end if
       
       if (faperp > epsilon(0.)) then
          istatus = nf_def_var (ncid, "aperp_r", NF_DOUBLE, 3, &
               (/ thetaid, kyid, kxid /), aperpr_id)
          if (istatus /= 0) then
             ierr = error_unit()
             write(ierr,*) "nf_def_var aperp error: ", nf_strerror(istatus)
             goto 1
          end if
          
          istatus = nf_def_var (ncid, "aperp_i", NF_DOUBLE, 3, &
               (/ thetaid, kyid, kxid /), aperpi_id)
          if (istatus /= 0) then
             ierr = error_unit()
             write(ierr,*) "nf_def_var aperp error: ", nf_strerror(istatus)
             goto 1
          end if
       end if
    end if
    
    istatus = nf_enddef (ncid)
    if (istatus /= 0) then
       ierr = error_unit()
       write(ierr,*) "nf_enddef error: ", nf_strerror(istatus)
       goto 1
    end if
 
    tmp1 = t0
    istatus = nf_put_var_double (ncid, t0id, tmp1)
    if (istatus /= 0) then
       ierr = error_unit()
       write(ierr,*) "nf_put_var_double t0 error: ", nf_strerror(istatus)
       goto 1
    end if
 
    tmp1 = delt0
    istatus = nf_put_var_double (ncid, delt0id, tmp1)
    if (istatus /= 0) then
       ierr = error_unit()
       write(ierr,*) "nf_put_var_double delt0 error: ", nf_strerror(istatus)
       goto 1
    end if
 
1   continue
    if (istatus /= 0) then
       i = nf_close (ncid)
       return
    end if

    if (n_elements > 0) then

       if (.not. allocated(tmpr)) allocate (tmpr(2*ntgrid+1,2,g_lo%llim_proc:g_lo%ulim_alloc))

       tmpr = real(g)
       istatus = nf_put_var_double (ncid, gr_id, tmpr)
       if (istatus /= 0) then
          ierr = error_unit()
          write(ierr,*) "nf_put_var_double gr error: ", nf_strerror(istatus),' ',iproc
       end if
       
       tmpr = aimag(g)
       istatus = nf_put_var_double (ncid, gi_id, tmpr)
       if (istatus /= 0) then
          ierr = error_unit()
          write(ierr,*) "nf_put_var_double gi error: ", nf_strerror(istatus),' ',iproc
       end if
       
       if (.not. allocated(ftmpr)) allocate (ftmpr(2*ntgrid+1,ntheta0,naky))
       if (.not. allocated(ftmpi)) allocate (ftmpi(2*ntgrid+1,ntheta0,naky))

       if (fphi > epsilon(0.)) then
          ftmpr = real(phinew)
          istatus = nf_put_var_double (ncid, phir_id, ftmpr)
          if (istatus /= 0) then
             ierr = error_unit()
             write(ierr,*) "nf_put_var_double phir error: ", nf_strerror(istatus),' ',iproc
          end if
          
          ftmpi = aimag(phinew)
          istatus = nf_put_var_double (ncid, phii_id, ftmpi)
          if (istatus /= 0) then
             ierr = error_unit()
             write(ierr,*) "nf_put_var_double phii error: ", nf_strerror(istatus),' ',iproc
          end if
       end if

       if (fapar > epsilon(0.)) then
          ftmpr = real(aparnew)
          istatus = nf_put_var_double (ncid, aparr_id, ftmpr)
          if (istatus /= 0) then
             ierr = error_unit()
             write(ierr,*) "nf_put_var_double aparr error: ", nf_strerror(istatus),' ',iproc
          end if
          
          ftmpi = aimag(aparnew)
          istatus = nf_put_var_double (ncid, apari_id, ftmpi)
          if (istatus /= 0) then
             ierr = error_unit()
             write(ierr,*) "nf_put_var_double apari error: ", nf_strerror(istatus),' ',iproc
          end if
       end if

       if (faperp > epsilon(0.)) then
          ftmpr = real(aperpnew)
          istatus = nf_put_var_double (ncid, aperpr_id, ftmpr)
          if (istatus /= 0) then
             ierr = error_unit()
             write(ierr,*) "nf_put_var_double aperpr error: ", nf_strerror(istatus),' ',iproc
          end if
          
          ftmpi = aimag(aperpnew)
          istatus = nf_put_var_double (ncid, aperpi_id, ftmpi)
          if (istatus /= 0) then
             ierr = error_unit()
             write(ierr,*) "nf_put_var_double aperpi error: ", nf_strerror(istatus),' ',iproc
          end if
       end if
    end if

!    if (exit) then
    i = nf_close (ncid)
!    else
!       i = nf_sync (ncid)
!    end if

  end subroutine gs2_save_for_restart

  subroutine gs2_restore_many (g, scale, istatus, fphi, fapar, faperp, many)
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo
    use mp, only: proc0, iproc, nproc
    use fields_arrays, only: phinew, aparnew, aperpnew
    use fields_arrays, only: phi, apar, aperp
    use kt_grids, only: naky, ntheta0
    use file_utils, only: error_unit
    implicit none
    character (305) :: file_proc
    character (5) :: suffix
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (out) :: g
    real, intent (in) :: scale
    integer, intent (out) :: istatus
    real, intent (in) :: fphi, fapar, faperp
    logical, intent (in) :: many
    include 'netcdf.inc'
    integer :: n_elements
    double precision :: tmp1
    integer :: iglo, i, th, h, t, u, ierr
    real :: fac

    n_elements = g_lo%ulim_proc-g_lo%llim_proc+1
    if (n_elements <= 0) return
    
    file_proc = trim(restart_file)
    
    if (nproc >= 10000) then
       if (proc0) write(*,*) 'Too many procs for i/o to work right!'
    else
       th = iproc / 1000
       h = (iproc - th * 1000) / 100
       t = (iproc - th * 1000 - h * 100) / 10
       u = (iproc - th * 1000 - h * 100 - t * 10)
       suffix = '.'//achar(48+th)//achar(48+h)//achar(48+t)//achar(48+u)
       file_proc = trim(trim(file_proc)//suffix)
    endif
    
    istatus = nf_open (file_proc, nf_write, ncid)
    if (istatus /= 0) then
       ierr = error_unit()
       write(ierr,*) "nf_open error: ", nf_strerror(istatus),' ',iproc
    end if
    
    istatus = nf_inq_dimid (ncid, "theta", thetaid)
    if (istatus /= 0) then
       ierr = error_unit()
       write(ierr,*) "nf_inq_dimid theta error: ", nf_strerror(istatus),' ',iproc
    end if
    
    istatus = nf_inq_dimid (ncid, "sign", signid)
    if (istatus /= 0) then
       ierr = error_unit()
       write(ierr,*) "nf_inq_dimid sign error: ", nf_strerror(istatus),' ',iproc
    end if
    
    istatus = nf_inq_dimid (ncid, "glo", gloid)
    if (istatus /= 0) then
       ierr = error_unit()
       write(ierr,*) "nf_inq_dimid glo error: ", nf_strerror(istatus),' ',iproc
    end if
    
    istatus = nf_inq_dimid (ncid, "aky", kyid)
    if (istatus /= 0) then
       ierr = error_unit()
       write(ierr,*) "nf_inq_dimid aky error: ", nf_strerror(istatus),' ',iproc
    end if
    
    istatus = nf_inq_dimid (ncid, "akx", kxid)
    if (istatus /= 0) then
       ierr = error_unit()
       write(ierr,*) "nf_inq_dimid akx error: ", nf_strerror(istatus),' ',iproc
    end if
    
    istatus = nf_inq_dimlen (ncid, thetaid, i)       
    if (istatus /= 0) then
       ierr = error_unit()
       write(ierr,*) "nf_inq_dimlen theta error: ", nf_strerror(istatus),' ',iproc
    end if
    if (i /= 2*ntgrid + 1) write(*,*) 'Restart error: ntgrid=? ',i,' : ',ntgrid,' : ',iproc
    
    istatus = nf_inq_dimlen (ncid, signid, i)       
    if (istatus /= 0) then
       ierr = error_unit()
       write(ierr,*) "nf_inq_dimlen sign error: ", nf_strerror(istatus),' ',iproc
    end if
    if (i /= 2) write(*,*) 'Restart error: sign=? ',i,' : ',iproc
    
    istatus = nf_inq_dimlen (ncid, gloid, i)       
    if (istatus /= 0) then
       ierr = error_unit()
       write(ierr,*) "nf_inq_dimlen glo error: ", nf_strerror(istatus),' ',iproc
    end if
    if (i /= g_lo%ulim_proc-g_lo%llim_proc+1) write(*,*) 'Restart error: glo=? ',i,' : ',iproc
    
    istatus = nf_inq_dimlen (ncid, kyid, i)       
    if (istatus /= 0) then
       ierr = error_unit()
       write(ierr,*) "nf_inq_dimlen aky error: ", nf_strerror(istatus),' ',iproc
    end if
    if (i /= naky) write(*,*) 'Restart error: naky=? ',i,' : ',naky,' : ',iproc
    
    istatus = nf_inq_dimlen (ncid, kxid, i)       
    if (istatus /= 0) then
       ierr = error_unit()
       write(ierr,*) "nf_inq_dimlen akx error: ", nf_strerror(istatus),' ',iproc
    end if
    if (i /= ntheta0) write(*,*) 'Restart error: ntheta0=? ',i,' : ',ntheta0,' : ',iproc
    
    if (fphi > epsilon(0.)) then
       istatus = nf_inq_varid (ncid, "phi_r", phir_id)
       if (istatus /= 0) then
          ierr = error_unit()
          write(ierr,*) "nf_inq_varid phir error: ", nf_strerror(istatus),' ',iproc
       end if
       
       istatus = nf_inq_varid (ncid, "phi_i", phii_id)
       if (istatus /= 0) then
          ierr = error_unit()
          write(ierr,*) "nf_inq_varid phii error: ", nf_strerror(istatus),' ',iproc
       end if
    end if
    
    if (fapar > epsilon(0.)) then
       istatus = nf_inq_varid (ncid, "apar_r", aparr_id)
       if (istatus /= 0) then
          ierr = error_unit()
          write(ierr,*) "nf_inq_varid aparr error: ", nf_strerror(istatus),' ',iproc
       end if
       
       istatus = nf_inq_varid (ncid, "apar_i", apari_id)
       if (istatus /= 0) then
          ierr = error_unit()
          write(ierr,*) "nf_inq_varid apari error: ", nf_strerror(istatus),' ',iproc
       end if
    end if
    
    if (faperp > epsilon(0.)) then
       istatus = nf_inq_varid (ncid, "aperp_r", aperpr_id)
       if (istatus /= 0) then
          ierr = error_unit()
          write(ierr,*) "nf_inq_varid aperpr error: ", nf_strerror(istatus),' ',iproc
       end if
       
       istatus = nf_inq_varid (ncid, "aperp_i", aperpi_id)
       if (istatus /= 0) then
          ierr = error_unit()
          write(ierr,*) "nf_inq_varid aperpi error: ", nf_strerror(istatus),' ',iproc
       end if
    end if
    
    istatus = nf_inq_varid (ncid, "gr", gr_id)
    if (istatus /= 0) then
       ierr = error_unit()
       write(ierr,*) "nf_inq_varid gr error: ", nf_strerror(istatus),' ',iproc
    end if
    
    istatus = nf_inq_varid (ncid, "gi", gi_id)
    if (istatus /= 0) then
       ierr = error_unit()
       write(ierr,*) "nf_inq_varid gi error: ", nf_strerror(istatus),' ',iproc
    end if
  
    if (.not. allocated(tmpr)) allocate (tmpr(2*ntgrid+1,2,g_lo%llim_proc:g_lo%ulim_alloc))
    if (.not. allocated(tmpi)) allocate (tmpi(2*ntgrid+1,2,g_lo%llim_proc:g_lo%ulim_alloc))

    tmpr = 0.; tmpi = 0.
    istatus = nf_get_var_double (ncid, gr_id, tmpr)
    if (istatus /= 0) then
       ierr = error_unit()
       write(ierr,*) "nf_get_var_double gr error: ", nf_strerror(istatus),' ',iproc
    end if

    istatus = nf_get_var_double (ncid, gi_id, tmpi)
    if (istatus /= 0) then
       ierr = error_unit()
       write(ierr,*) "nf_get_var_double gi error: ", nf_strerror(istatus),' ',iproc
    end if

    g = cmplx(tmpr, tmpi)

    if (.not. allocated(ftmpr)) allocate (ftmpr(2*ntgrid+1,ntheta0,naky))
    if (.not. allocated(ftmpi)) allocate (ftmpi(2*ntgrid+1,ntheta0,naky))

    if (fphi > epsilon(0.)) then
       istatus = nf_get_var_double (ncid, phir_id, ftmpr)
       if (istatus /= 0) then
          ierr = error_unit()
          write(ierr,*) "nf_get_var_double phir error: ", nf_strerror(istatus),' ',iproc
       end if
       
       istatus = nf_get_var_double (ncid, phii_id, ftmpi)
       if (istatus /= 0) then
          ierr = error_unit()
          write(ierr,*) "nf_get_var_double phii error: ", nf_strerror(istatus),' ',iproc
       end if
       
       phi = 0.
       phinew = cmplx(ftmpr, ftmpi)
    end if

    if (fapar > epsilon(0.)) then
       istatus = nf_get_var_double (ncid, aparr_id, ftmpr)
       if (istatus /= 0) then
          ierr = error_unit()
          write(ierr,*) "nf_get_var_double aparr error: ", nf_strerror(istatus),' ',iproc
       end if
       
       istatus = nf_get_var_double (ncid, apari_id, ftmpi)
       if (istatus /= 0) then
          ierr = error_unit()
          write(ierr,*) "nf_get_var_double apari error: ", nf_strerror(istatus),' ',iproc
       end if
       
       apar = 0.
       aparnew = cmplx(ftmpr, ftmpi)
    end if

    if (faperp > epsilon(0.)) then
       istatus = nf_get_var_double (ncid, aperpr_id, ftmpr)
       if (istatus /= 0) then
          ierr = error_unit()
          write(ierr,*) "nf_get_var_double aperpr error: ", nf_strerror(istatus),' ',iproc
       end if
       
       istatus = nf_get_var_double (ncid, aperpi_id, ftmpi)
       if (istatus /= 0) then
          ierr = error_unit()
          write(ierr,*) "nf_get_var_double aperpi error: ", nf_strerror(istatus),' ',iproc
       end if
       
       aperp = 0.
       aperpnew = cmplx(ftmpr, ftmpi)
    end if

    if (scale > 0.) then
       g = g*scale
       phinew = phinew*scale
       aparnew = aparnew*scale
       aperpnew = aperpnew*scale
    else
       fac = - scale/(maxval(abs(phinew)))
       g = g*fac
       phinew = phinew*fac
       aparnew = aparnew*fac
       aperpnew = aperpnew*fac
    end if

    istatus = nf_close (ncid)       
    if (istatus /= 0) then
    ierr = error_unit()
       write(ierr,*) "nf_close error: ", nf_strerror(istatus),' ',iproc
    end if

  end subroutine gs2_restore_many

  subroutine gs2_restore_one (g, scale, istatus, fphi, fapar, faperp) 
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo
    use mp, only: proc0, iproc, nproc
    use fields_arrays, only: phinew, aparnew, aperpnew
    use fields_arrays, only: phi, apar, aperp
    use kt_grids, only: naky, ntheta0
    use file_utils, only: error_unit
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (out) :: g
    real, intent (in) :: scale
    integer, intent (out) :: istatus
    real, intent (in) :: fphi, fapar, faperp
    include 'netcdf.inc'
    integer :: n_elements, ierr
    double precision :: tmp1
    integer :: iglo, i, ip

    istatus = nf_open (trim(restart_file), 0, ncid)
    if (istatus /= 0) then
       ierr = error_unit()
       write(ierr,*) "nf_open error: ", nf_strerror(istatus),' ',iproc
    end if

    istatus = nf_inq_dimid (ncid, "theta", thetaid)
    if (istatus /= 0) then
       ierr = error_unit()
       write(ierr,*) "nf_inq_dimid theta error: ", nf_strerror(istatus),' ',iproc
    end if

    istatus = nf_inq_dimid (ncid, "sign", signid)
    if (istatus /= 0) then
       ierr = error_unit()
       write(ierr,*) "nf_inq_dimid sign error: ", nf_strerror(istatus),' ',iproc
    end if

    istatus = nf_inq_dimid (ncid, "glo", gloid)
    if (istatus /= 0) then
       ierr = error_unit()
       write(ierr,*) "nf_inq_dimid glo error: ", nf_strerror(istatus),' ',iproc
    end if

    istatus = nf_inq_dimid (ncid, "aky", kyid)
    if (istatus /= 0) then
       ierr = error_unit()
       write(ierr,*) "nf_inq_dimid aky error: ", nf_strerror(istatus),' ',iproc
    end if

    istatus = nf_inq_dimid (ncid, "akx", kxid)
    if (istatus /= 0) then
       ierr = error_unit()
       write(ierr,*) "nf_inq_dimid akx error: ", nf_strerror(istatus),' ',iproc
    end if

    istatus = nf_inq_dimlen (ncid, thetaid, i)       
    if (istatus /= 0) then
       ierr = error_unit()
       write(ierr,*) "nf_inq_dimlen theta error: ", nf_strerror(istatus),' ',iproc
    end if
    if (i /= 2*ntgrid + 1) write(*,*) 'Restart error: ntgrid=? ',i,' : ',ntgrid,' : ',iproc

    istatus = nf_inq_dimlen (ncid, signid, i)       
    if (istatus /= 0) then
       ierr = error_unit()
       write(ierr,*) "nf_inq_dimlen sign error: ", nf_strerror(istatus),' ',iproc
    end if
    if (i /= 2) write(*,*) 'Restart error: sign=? ',i,' : ',iproc

    istatus = nf_inq_dimlen (ncid, gloid, i)       
    if (istatus /= 0) then
       ierr = error_unit()
       write(ierr,*) "nf_inq_dimlen glo error: ", nf_strerror(istatus),' ',iproc
    end if
    if (i /= g_lo%ulim_world-g_lo%llim_world+1) write(*,*) 'Restart error: glo=? ',i,' : ',iproc

    istatus = nf_inq_dimlen (ncid, kyid, i)       
    if (istatus /= 0) then
       ierr = error_unit()
       write(ierr,*) "nf_inq_dimlen aky error: ", nf_strerror(istatus),' ',iproc
    end if
    if (i /= naky) write(*,*) 'Restart error: naky=? ',i,' : ',naky,' : ',iproc

    istatus = nf_inq_dimlen (ncid, kxid, i)       
    if (istatus /= 0) then
       ierr = error_unit()
       write(ierr,*) "nf_inq_dimlen akx error: ", nf_strerror(istatus),' ',iproc
    end if
    if (i /= ntheta0) write(*,*) 'Restart error: ntheta0=? ',i,' : ',ntheta0,' : ',iproc

    if (fphi > epsilon(0.)) then
       istatus = nf_inq_varid (ncid, "phi_r", phir_id)
       if (istatus /= 0) then
          ierr = error_unit()
          write(ierr,*) "nf_inq_varid phir error: ", nf_strerror(istatus),' ',iproc
       end if
       
       istatus = nf_inq_varid (ncid, "phi_i", phii_id)
       if (istatus /= 0) then
          ierr = error_unit()
          write(ierr,*) "nf_inq_varid phii error: ", nf_strerror(istatus),' ',iproc
       end if
    end if

    if (fapar > epsilon(0.)) then
       istatus = nf_inq_varid (ncid, "apar_r", aparr_id)
       if (istatus /= 0) then
          ierr = error_unit()
          write(ierr,*) "nf_inq_varid aparr error: ", nf_strerror(istatus),' ',iproc
       end if
       
       istatus = nf_inq_varid (ncid, "apar_i", apari_id)
       if (istatus /= 0) then
          ierr = error_unit()
          write(ierr,*) "nf_inq_varid apari error: ", nf_strerror(istatus),' ',iproc
       end if
    end if

    if (faperp > epsilon(0.)) then
       istatus = nf_inq_varid (ncid, "aperp_r", aperpr_id)
       if (istatus /= 0) then
          ierr = error_unit()
          write(ierr,*) "nf_inq_varid aperpr error: ", nf_strerror(istatus),' ',iproc
       end if
       
       istatus = nf_inq_varid (ncid, "aperp_i", aperpi_id)
       if (istatus /= 0) then
          ierr = error_unit()
          write(ierr,*) "nf_inq_varid aperpi error: ", nf_strerror(istatus),' ',iproc
       end if
    end if

    istatus = nf_inq_varid (ncid, "gr", gr_id)
    if (istatus /= 0) then
       ierr = error_unit()
       write(ierr,*) "nf_inq_varid gr error: ", nf_strerror(istatus),' ',iproc
    end if

    istatus = nf_inq_varid (ncid, "gi", gi_id)
    if (istatus /= 0) then
       ierr = error_unit()
       write(ierr,*) "nf_inq_varid gi error: ", nf_strerror(istatus),' ',iproc
    end if

    n_elements = g_lo%ulim_proc - g_lo%llim_proc + 1

    if (n_elements <= 0) then
       phinew = 0.
       aparnew = 0.
       aperpnew = 0.
       return
    endif

    if (.not. allocated(tmpr)) allocate (tmpr(2*ntgrid+1, 2, n_elements))    
    if (.not. allocated(tmpi)) allocate (tmpi(2*ntgrid+1, 2, n_elements))

    tmpr = 0.;  tmpi = 0.
       
    istatus = nf_get_vara_double (ncid, gr_id, (/ 1, 1, g_lo%llim_proc /), &
         (/ 2*ntgrid + 1, 2, n_elements /), tmpr)
    if (istatus /= 0) then
       ierr = error_unit()
       write(ierr,*) "nf_get_vara_double gr error: ", nf_strerror(istatus),' ',iproc
    end if
    
    istatus = nf_get_vara_double (ncid, gi_id,  (/ 1, 1, g_lo%llim_proc /), &
         (/ 2*ntgrid + 1, 2, n_elements /), tmpi)
    if (istatus /= 0) then
       ierr = error_unit()
       write(ierr,*) "nf_get_vara_double gi error: ", nf_strerror(istatus),' ',iproc
    end if

    g = cmplx(tmpr, tmpi)*scale

    if (.not. allocated(ftmpr)) allocate (ftmpr(2*ntgrid+1,ntheta0,naky))
    if (.not. allocated(ftmpi)) allocate (ftmpi(2*ntgrid+1,ntheta0,naky))

    if (fphi > epsilon(0.)) then
       istatus = nf_get_var_double (ncid, phir_id, ftmpr)
       if (istatus /= 0) then
          ierr = error_unit()
          write(ierr,*) "nf_get_var_double phir error: ", nf_strerror(istatus),' ',iproc
       end if
       
       istatus = nf_get_var_double (ncid, phii_id, ftmpi)
       if (istatus /= 0) then
          ierr = error_unit()
          write(ierr,*) "nf_get_var_double phii error: ", nf_strerror(istatus),' ',iproc
       end if
       
       phi = 0.
       phinew = cmplx(ftmpr, ftmpi)*scale
    end if

    if (fapar > epsilon(0.)) then
       istatus = nf_get_var_double (ncid, aparr_id, ftmpr)
       if (istatus /= 0) then
          ierr = error_unit()
          write(ierr,*) "nf_get_var_double aparr error: ", nf_strerror(istatus),' ',iproc
       end if
       
       istatus = nf_get_var_double (ncid, apari_id, ftmpi)
       if (istatus /= 0) then
          ierr = error_unit()
          write(ierr,*) "nf_get_var_double apari error: ", nf_strerror(istatus),' ',iproc
       end if
       
       apar = 0.
       aparnew = cmplx(ftmpr, ftmpi)*scale
    end if

    if (faperp > epsilon(0.)) then
       istatus = nf_get_var_double (ncid, aperpr_id, ftmpr)
       if (istatus /= 0) then
          ierr = error_unit()
          write(ierr,*) "nf_get_var_double aperpr error: ", nf_strerror(istatus),' ',iproc
       end if
       
       istatus = nf_get_var_double (ncid, aperpi_id, ftmpi)
       if (istatus /= 0) then
          ierr = error_unit()
          write(ierr,*) "nf_get_var_double aperpi error: ", nf_strerror(istatus),' ',iproc
       end if
       
       aperp = 0.
       aperpnew = cmplx(ftmpr, ftmpi)*scale
    end if

    istatus = nf_close (ncid)       
    if (istatus /= 0) then
       ierr = error_unit()
       write(ierr,*) "nf_close error: ", nf_strerror(istatus),' ',iproc
    end if

  end subroutine gs2_restore_one

  subroutine init_save (file)
    character(300), intent (in) :: file
    
    restart_file = file

  end subroutine init_save

  subroutine init_dt (delt0, istatus)

    use mp, only: nproc, proc0, broadcast
    use file_utils, only: error_unit
    implicit none
    character (305) :: file_proc
    real, intent (in out) :: delt0
    integer, intent (out) :: istatus
    integer :: ierr
    include 'netcdf.inc'
    double precision :: tmp1

    if (proc0) then
       file_proc=trim(trim(restart_file)//'.0000')
       
       istatus = nf_open (file_proc, 0, ncid)
       if (istatus /= 0) then
          ierr = error_unit()
          write(ierr,*) "nf_open in init_dt error: ", nf_strerror(istatus) 
       endif

       istatus = nf_inq_varid (ncid, "delt0", delt0id)
       if (istatus /= 0) then
          ierr = error_unit()
          write(ierr,*) "nf_inq_varid for delt0 in init_dt: ", nf_strerror(istatus) 
       endif

       istatus = nf_get_var_double (ncid, delt0id, tmp1)
       
       if (istatus /= 0) then
          ierr = error_unit()
          write(ierr,*) "nf_get_var_double delt0 error: ", nf_strerror(istatus) 
          delt0 = -1.
       else
          delt0 = tmp1
       endif
       istatus = nf_close (ncid)       
    endif

    call broadcast (istatus)
    call broadcast (delt0)

  end subroutine init_dt

  subroutine init_tstart (tstart, istatus)

    use mp, only: nproc, proc0, broadcast
    use file_utils, only: error_unit
    implicit none
    character (305) :: file_proc
    real, intent (in out) :: tstart
    integer, intent (out) :: istatus
    integer :: ierr
    include 'netcdf.inc'
    double precision :: tmp1

    if (proc0) then
       file_proc=trim(trim(restart_file)//'.0000')
       
       istatus = nf_open (file_proc, 0, ncid)
       if (istatus /= 0) then
          ierr = error_unit()
          write(ierr,*) "nf_open in init_dt error: ", nf_strerror(istatus) 
       endif

       istatus = nf_inq_varid (ncid, "t0", t0id)
       if (istatus /= 0) then
          ierr = error_unit()
          write(ierr,*) "nf_inq_varid for t0 in init_tstart: ", nf_strerror(istatus) 
       endif

       istatus = nf_get_var_double (ncid, t0id, tmp1)

       if (istatus /= 0) then
          ierr = error_unit()
          write(ierr,*) "nf_get_var_double t0 error: ", nf_strerror(istatus) 
          tstart = -1.
       else
          tstart = tmp1
       endif           
       istatus = nf_close (ncid)       
    endif

    call broadcast (istatus)
    call broadcast (tstart)

  end subroutine init_tstart

end module gs2_save
