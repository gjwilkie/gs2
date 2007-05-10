module gs2_save

  implicit none

  public :: gs2_restore, gs2_save_for_restart
  public :: init_save, gs2_read_aminv

  interface gs2_restore
     module procedure gs2_restore_many, gs2_restore_one
  end interface
  
  private
  character(300), save :: restart_file

contains

  subroutine gs2_save_for_restart (g, aminv, t0, istatus)
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo, f_lo
    use mp, only: nproc, iproc, proc0
    use fields_arrays, only: phinew, aparnew, aperpnew, nidx
    use kt_grids, only: naky, ntheta0
    implicit none
    character (305) :: file_proc
    character (5) :: suffix
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,f_lo%llim_proc:), intent (in) :: aminv
    real, intent (in) :: t0
    integer, intent (out) :: istatus
    include 'netcdf.inc'
    integer :: ncid, thetaid, signid, gloid, kyid, kxid, nidxid, floid
    integer :: phir_id, phii_id, aparr_id, apari_id, aperpr_id, aperpi_id
    integer :: t0id, gr_id, gi_id
    integer :: aminvr_id, aminvi_id
    double precision :: tmp1
    double precision, dimension(2*ntgrid+1,2,g_lo%llim_proc:g_lo%ulim_alloc) :: tmp 
    double precision, dimension(2*ntgrid+1,naky,ntheta0) :: ftmpr, ftmpi
    double precision, dimension(nidx,f_lo%llim_proc:f_lo%ulim_alloc) :: atmp
    integer :: i, n_elements, n_f_elements
    integer :: th, h, t, u
    
    n_f_elements = f_lo%ulim_proc-f_lo%llim_proc+1

    n_elements = g_lo%ulim_proc-g_lo%llim_proc+1
    if (n_elements <= 0 .and. n_f_elements <= 0) return

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
    
    istatus = nf_create (file_proc, 0, ncid)
    if (istatus /= 0) then
       print *, "nf_create error: ", nf_strerror(istatus)
       goto 1
    end if

    if (n_f_elements > 0) then
       istatus = nf_def_dim (ncid, "impidx", nidx, nidxid)
       if (istatus /= 0) then
          print *, "nf_def_dim nidx error: ", nf_strerror(istatus)
          goto 1
       end if

       istatus = nf_def_dim (ncid, "flo", n_f_elements, floid)
       if (istatus /= 0) then
          print *, "nf_def_dim flo error: ", nf_strerror(istatus)
          goto 1
       end if       
    end if

    if (n_elements > 0) then
       istatus = nf_def_dim (ncid, "theta", 2*ntgrid+1, thetaid)
       if (istatus /= 0) then
          print *, "nf_def_dim theta error: ", nf_strerror(istatus)
          goto 1
       end if

       istatus = nf_def_dim (ncid, "sign", 2, signid)
       if (istatus /= 0) then
          print *, "nf_def_dim sign error: ", nf_strerror(istatus)
          goto 1
       end if
       
       istatus = nf_def_dim (ncid, "glo", n_elements, gloid)
       if (istatus /= 0) then
          print *, "nf_def_dim glo error: ", nf_strerror(istatus)
          goto 1
       end if
       
       istatus = nf_def_dim (ncid, "aky", naky, kyid)
       if (istatus /= 0) then
          print *, "nf_def_dim aky error: ", nf_strerror(istatus)
          goto 1
       end if
       
       istatus = nf_def_dim (ncid, "akx", ntheta0, kxid)
       if (istatus /= 0) then
          print *, "nf_def_dim akx error: ", nf_strerror(istatus)
          goto 1
       end if
    end if

    istatus = nf_def_var (ncid, "t0", NF_DOUBLE, 0, 0, t0id)
    if (istatus /= 0) then
       print *, "nf_def_var t0 error: ", nf_strerror(istatus)
       goto 1
    end if
    
    if (n_f_elements > 0) then
       istatus = nf_def_var (ncid, "aminvr", NF_DOUBLE, 2, &
            (/ thetaid, floid /), aminvr_id)
       if (istatus /= 0) then
          print *, "nf_def_var aminv error: ", nf_strerror(istatus)
          goto 1
       end if
    
       istatus = nf_def_var (ncid, "aminvi", NF_DOUBLE, 2, &
            (/ thetaid, floid /), aminvi_id)
       if (istatus /= 0) then
          print *, "nf_def_var aminv error: ", nf_strerror(istatus)
          goto 1
       end if
    end if

    if (n_elements > 0) then
       istatus = nf_def_var (ncid, "gr", NF_DOUBLE, 3, &
            (/ thetaid, signid, gloid /), gr_id)
       if (istatus /= 0) then
          print *, "nf_def_var g error: ", nf_strerror(istatus)
          goto 1
       end if
       
       istatus = nf_def_var (ncid, "gi", NF_DOUBLE, 3, &
            (/ thetaid, signid, gloid /), gi_id)
       if (istatus /= 0) then
          print *, "nf_def_var g error: ", nf_strerror(istatus)
          goto 1
       end if
       
       istatus = nf_def_var (ncid, "phi_r", NF_DOUBLE, 3, &
            (/ thetaid, kyid, kxid /), phir_id)
       if (istatus /= 0) then
          print *, "nf_def_var phi error: ", nf_strerror(istatus)
          goto 1
       end if
       
       istatus = nf_def_var (ncid, "phi_i", NF_DOUBLE, 3, &
            (/ thetaid, kyid, kxid /), phii_id)
       if (istatus /= 0) then
          print *, "nf_def_var phi error: ", nf_strerror(istatus)
          goto 1
       end if
       
       istatus = nf_def_var (ncid, "apar_r", NF_DOUBLE, 3, &
            (/ thetaid, kyid, kxid /), aparr_id)
       if (istatus /= 0) then
          print *, "nf_def_var apar error: ", nf_strerror(istatus)
          goto 1
       end if
       
       istatus = nf_def_var (ncid, "apar_i", NF_DOUBLE, 3, &
            (/ thetaid, kyid, kxid /), apari_id)
       if (istatus /= 0) then
          print *, "nf_def_var apar error: ", nf_strerror(istatus)
          goto 1
       end if
       
       istatus = nf_def_var (ncid, "aperp_r", NF_DOUBLE, 3, &
            (/ thetaid, kyid, kxid /), aperpr_id)
       if (istatus /= 0) then
          print *, "nf_def_var aperp error: ", nf_strerror(istatus)
          goto 1
       end if
       
       istatus = nf_def_var (ncid, "aperp_i", NF_DOUBLE, 3, &
            (/ thetaid, kyid, kxid /), aperpi_id)
       if (istatus /= 0) then
          print *, "nf_def_var aperp error: ", nf_strerror(istatus)
          goto 1
       end if
    end if

    istatus = nf_enddef (ncid)
    if (istatus /= 0) then
       print *, "nf_enddef error: ", nf_strerror(istatus)
       goto 1
    end if
    
    tmp1 = t0
    istatus = nf_put_var_double (ncid, t0id, tmp1)
    if (istatus /= 0) then
       print *, "nf_put_var_double t0 error: ", nf_strerror(istatus)
       goto 1
    end if
 
1   continue
    if (istatus /= 0) then
       i = nf_close (ncid)
       return
    end if

    if (n_f_elements > 0) then
       atmp = real(aminv)
       istatus = nf_put_var_double (ncid, aminvr_id, atmp)
       if (istatus /= 0) then
          print *, "nf_put_var_double aminvr error: ", nf_strerror(istatus),' ',iproc
       end if
       
       atmp = aimag(aminv)
       istatus = nf_put_var_double (ncid, aminvi_id, atmp)
       if (istatus /= 0) then
          print *, "nf_put_var_double aminvi error: ", nf_strerror(istatus),' ',iproc
       end if             
    end if

    if (n_elements > 0) then
       tmp = real(g)
       istatus = nf_put_var_double (ncid, gr_id, tmp)
       if (istatus /= 0) then
          print *, "nf_put_var_double gr error: ", nf_strerror(istatus),' ',iproc
       end if
       
       tmp = aimag(g)
       istatus = nf_put_var_double (ncid, gi_id, tmp)
       if (istatus /= 0) then
          print *, "nf_put_var_double gi error: ", nf_strerror(istatus),' ',iproc
       end if
       
       ftmpr = real(phinew)
       istatus = nf_put_var_double (ncid, phir_id, ftmpr)
       if (istatus /= 0) then
          print *, "nf_put_var_double phir error: ", nf_strerror(istatus),' ',iproc
       end if
       
       ftmpi = aimag(phinew)
       istatus = nf_put_var_double (ncid, phii_id, ftmpi)
       if (istatus /= 0) then
          print *, "nf_put_var_double phii error: ", nf_strerror(istatus),' ',iproc
       end if
       
       ftmpr = real(aparnew)
       istatus = nf_put_var_double (ncid, aparr_id, ftmpr)
       if (istatus /= 0) then
          print *, "nf_put_var_double aparr error: ", nf_strerror(istatus),' ',iproc
       end if
       
       ftmpi = aimag(aparnew)
       istatus = nf_put_var_double (ncid, apari_id, ftmpi)
       if (istatus /= 0) then
          print *, "nf_put_var_double apari error: ", nf_strerror(istatus),' ',iproc
       end if
       
       ftmpr = real(aperpnew)
       istatus = nf_put_var_double (ncid, aperpr_id, ftmpr)
       if (istatus /= 0) then
          print *, "nf_put_var_double aperpr error: ", nf_strerror(istatus),' ',iproc
       end if
       
       ftmpi = aimag(aperpnew)
       istatus = nf_put_var_double (ncid, aperpi_id, ftmpi)
       if (istatus /= 0) then
          print *, "nf_put_var_double aperpi error: ", nf_strerror(istatus),' ',iproc
       end if
    end if

    i = nf_close (ncid)

  end subroutine gs2_save_for_restart

  subroutine gs2_restore_many (g, t0, istatus, many)
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo
    use mp, only: proc0, iproc, nproc
    use fields_arrays, only: phinew, aparnew, aperpnew
    use fields_arrays, only: phi, apar, aperp
    use kt_grids, only: naky, ntheta0
    implicit none
    character (305) :: file_proc
    character (5) :: suffix
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (out) :: g
    real, intent (out) :: t0
    integer, intent (out) :: istatus
    logical, intent (in) :: many
    include 'netcdf.inc'
    integer :: ncid, thetaid, signid, gloid, kyid, kxid
    integer :: phir_id, phii_id, aparr_id, apari_id, aperpr_id, aperpi_id
    integer :: t0id, gr_id, gi_id
    integer :: n_elements
    double precision :: tmp1
    double precision, dimension(2*ntgrid+1,2,g_lo%llim_proc:g_lo%ulim_alloc) :: tmpr, tmpi
    double precision, dimension(2*ntgrid+1,naky,ntheta0) :: ftmpr, ftmpi
    integer :: iglo, i, th, h, t, u

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
    
    istatus = nf_open (file_proc, 0, ncid)
    if (istatus /= 0) then
       print *, "nf_open error: ", nf_strerror(istatus),' ',iproc
    end if
    
    istatus = nf_inq_dimid (ncid, "theta", thetaid)
    if (istatus /= 0) then
       print *, "nf_inq_dimid theta error: ", nf_strerror(istatus),' ',iproc
    end if

    istatus = nf_inq_dimid (ncid, "sign", signid)
    if (istatus /= 0) then
       print *, "nf_inq_dimid sign error: ", nf_strerror(istatus),' ',iproc
    end if

    istatus = nf_inq_dimid (ncid, "glo", gloid)
    if (istatus /= 0) then
       print *, "nf_inq_dimid glo error: ", nf_strerror(istatus),' ',iproc
    end if

    istatus = nf_inq_dimid (ncid, "aky", kyid)
    if (istatus /= 0) then
       print *, "nf_inq_dimid aky error: ", nf_strerror(istatus),' ',iproc
    end if

    istatus = nf_inq_dimid (ncid, "akx", kxid)
    if (istatus /= 0) then
       print *, "nf_inq_dimid akx error: ", nf_strerror(istatus),' ',iproc
    end if

    istatus = nf_inq_dimlen (ncid, thetaid, i)       
    if (istatus /= 0) then
       print *, "nf_inq_dimlen theta error: ", nf_strerror(istatus),' ',iproc
    end if
    if (i /= 2*ntgrid + 1) write(*,*) 'Restart error: ntgrid=? ',i,' : ',ntgrid,' : ',iproc

    istatus = nf_inq_dimlen (ncid, signid, i)       
    if (istatus /= 0) then
       print *, "nf_inq_dimlen sign error: ", nf_strerror(istatus),' ',iproc
    end if
    if (i /= 2) write(*,*) 'Restart error: sign=? ',i,' : ',iproc

    istatus = nf_inq_dimlen (ncid, gloid, i)       
    if (istatus /= 0) then
       print *, "nf_inq_dimlen glo error: ", nf_strerror(istatus),' ',iproc
    end if
    if (i /= g_lo%ulim_proc-g_lo%llim_proc+1) write(*,*) 'Restart error: glo=? ',i,' : ',iproc

    istatus = nf_inq_dimlen (ncid, kyid, i)       
    if (istatus /= 0) then
       print *, "nf_inq_dimlen aky error: ", nf_strerror(istatus),' ',iproc
    end if
    if (i /= naky) write(*,*) 'Restart error: naky=? ',i,' : ',naky,' : ',iproc

    istatus = nf_inq_dimlen (ncid, kxid, i)       
    if (istatus /= 0) then
       print *, "nf_inq_dimlen akx error: ", nf_strerror(istatus),' ',iproc
    end if
    if (i /= ntheta0) write(*,*) 'Restart error: ntheta0=? ',i,' : ',ntheta0,' : ',iproc

    istatus = nf_inq_varid (ncid, "t0", t0id)       
    if (istatus /= 0) then
       print *, "nf_inq_varid t0 error: ", nf_strerror(istatus),' ',iproc
    end if

    istatus = nf_inq_varid (ncid, "phi_r", phir_id)
    if (istatus /= 0) then
       print *, "nf_inq_varid phir error: ", nf_strerror(istatus),' ',iproc
    end if

    istatus = nf_inq_varid (ncid, "phi_i", phii_id)
    if (istatus /= 0) then
       print *, "nf_inq_varid phii error: ", nf_strerror(istatus),' ',iproc
    end if

    istatus = nf_inq_varid (ncid, "apar_r", aparr_id)
    if (istatus /= 0) then
       print *, "nf_inq_varid aparr error: ", nf_strerror(istatus),' ',iproc
    end if

    istatus = nf_inq_varid (ncid, "apar_i", apari_id)
    if (istatus /= 0) then
       print *, "nf_inq_varid apari error: ", nf_strerror(istatus),' ',iproc
    end if

    istatus = nf_inq_varid (ncid, "aperp_r", aperpr_id)
    if (istatus /= 0) then
       print *, "nf_inq_varid aperpr error: ", nf_strerror(istatus),' ',iproc
    end if

    istatus = nf_inq_varid (ncid, "aperp_i", aperpi_id)
    if (istatus /= 0) then
       print *, "nf_inq_varid aperpi error: ", nf_strerror(istatus),' ',iproc
    end if

    istatus = nf_inq_varid (ncid, "gr", gr_id)
    if (istatus /= 0) then
       print *, "nf_inq_varid gr error: ", nf_strerror(istatus),' ',iproc
    end if

    istatus = nf_inq_varid (ncid, "gi", gi_id)
    if (istatus /= 0) then
       print *, "nf_inq_varid gi error: ", nf_strerror(istatus),' ',iproc
    end if

    istatus = nf_get_var_double (ncid, t0id, tmp1)
    if (istatus /= 0) then
       print *, "nf_get_var_double t0 error: ", nf_strerror(istatus),' ',iproc
    end if
    t0 = tmp1

    tmpr = 0.; tmpi = 0.
    istatus = nf_get_var_double (ncid, gr_id, tmpr)
    if (istatus /= 0) then
       print *, "nf_get_var_double gr error: ", nf_strerror(istatus),' ',iproc
    end if

    istatus = nf_get_var_double (ncid, gi_id, tmpi)
    if (istatus /= 0) then
       print *, "nf_get_var_double gi error: ", nf_strerror(istatus),' ',iproc
    end if

    g = cmplx(tmpr, tmpi)

    istatus = nf_get_var_double (ncid, phir_id, ftmpr)
    if (istatus /= 0) then
       print *, "nf_get_var_double phir error: ", nf_strerror(istatus),' ',iproc
    end if

    istatus = nf_get_var_double (ncid, phii_id, ftmpi)
    if (istatus /= 0) then
       print *, "nf_get_var_double phii error: ", nf_strerror(istatus),' ',iproc
    end if

    phi = 0.
    phinew = cmplx(ftmpr, ftmpi)

    istatus = nf_get_var_double (ncid, aparr_id, ftmpr)
    if (istatus /= 0) then
       print *, "nf_get_var_double aparr error: ", nf_strerror(istatus),' ',iproc
    end if

    istatus = nf_get_var_double (ncid, apari_id, ftmpi)
    if (istatus /= 0) then
       print *, "nf_get_var_double apari error: ", nf_strerror(istatus),' ',iproc
    end if

    apar = 0.
    aparnew = cmplx(ftmpr, ftmpi)

    istatus = nf_get_var_double (ncid, aperpr_id, ftmpr)
    if (istatus /= 0) then
       print *, "nf_get_var_double aperpr error: ", nf_strerror(istatus),' ',iproc
    end if

    istatus = nf_get_var_double (ncid, aperpi_id, ftmpi)
    if (istatus /= 0) then
       print *, "nf_get_var_double aperpi error: ", nf_strerror(istatus),' ',iproc
    end if

    aperp = 0.
    aperpnew = cmplx(ftmpr, ftmpi)

    istatus = nf_close (ncid)       
    if (istatus /= 0) then
       print *, "nf_close error: ", nf_strerror(istatus),' ',iproc
    end if

  end subroutine gs2_restore_many

  subroutine gs2_restore_one (g, t0, istatus)
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo
    use mp, only: proc0, iproc, nproc
    use fields_arrays, only: phinew, aparnew, aperpnew
    use fields_arrays, only: phi, apar, aperp
    use kt_grids, only: naky, ntheta0
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (out) :: g
    real, intent (out) :: t0
    integer, intent (out) :: istatus
    include 'netcdf.inc'
    integer :: ncid, thetaid, signid, gloid, kyid, kxid
    integer :: phir_id, phii_id, aparr_id, apari_id, aperpr_id, aperpi_id
    integer :: t0id, gr_id, gi_id
    integer :: n_elements
    double precision :: tmp1
    double precision, dimension(:,:,:), allocatable :: tmpr, tmpi
    double precision, dimension(2*ntgrid+1,naky,ntheta0) :: ftmpr, ftmpi
    integer :: iglo, i, ip

    istatus = nf_open (trim(restart_file), 0, ncid)
    if (istatus /= 0) then
       print *, "nf_open error: ", nf_strerror(istatus),' ',iproc
    end if

    istatus = nf_inq_dimid (ncid, "theta", thetaid)
    if (istatus /= 0) then
       print *, "nf_inq_dimid theta error: ", nf_strerror(istatus),' ',iproc
    end if

    istatus = nf_inq_dimid (ncid, "sign", signid)
    if (istatus /= 0) then
       print *, "nf_inq_dimid sign error: ", nf_strerror(istatus),' ',iproc
    end if

    istatus = nf_inq_dimid (ncid, "glo", gloid)
    if (istatus /= 0) then
       print *, "nf_inq_dimid glo error: ", nf_strerror(istatus),' ',iproc
    end if

    istatus = nf_inq_dimid (ncid, "aky", kyid)
    if (istatus /= 0) then
       print *, "nf_inq_dimid aky error: ", nf_strerror(istatus),' ',iproc
    end if

    istatus = nf_inq_dimid (ncid, "akx", kxid)
    if (istatus /= 0) then
       print *, "nf_inq_dimid akx error: ", nf_strerror(istatus),' ',iproc
    end if

    istatus = nf_inq_dimlen (ncid, thetaid, i)       
    if (istatus /= 0) then
       print *, "nf_inq_dimlen theta error: ", nf_strerror(istatus),' ',iproc
    end if
    if (i /= 2*ntgrid + 1) write(*,*) 'Restart error: ntgrid=? ',i,' : ',ntgrid,' : ',iproc

    istatus = nf_inq_dimlen (ncid, signid, i)       
    if (istatus /= 0) then
       print *, "nf_inq_dimlen sign error: ", nf_strerror(istatus),' ',iproc
    end if
    if (i /= 2) write(*,*) 'Restart error: sign=? ',i,' : ',iproc

    istatus = nf_inq_dimlen (ncid, gloid, i)       
    if (istatus /= 0) then
       print *, "nf_inq_dimlen glo error: ", nf_strerror(istatus),' ',iproc
    end if
    if (i /= g_lo%ulim_world-g_lo%llim_world+1) write(*,*) 'Restart error: glo=? ',i,' : ',iproc

    istatus = nf_inq_dimlen (ncid, kyid, i)       
    if (istatus /= 0) then
       print *, "nf_inq_dimlen aky error: ", nf_strerror(istatus),' ',iproc
    end if
    if (i /= naky) write(*,*) 'Restart error: naky=? ',i,' : ',naky,' : ',iproc

    istatus = nf_inq_dimlen (ncid, kxid, i)       
    if (istatus /= 0) then
       print *, "nf_inq_dimlen akx error: ", nf_strerror(istatus),' ',iproc
    end if
    if (i /= ntheta0) write(*,*) 'Restart error: ntheta0=? ',i,' : ',ntheta0,' : ',iproc

    istatus = nf_inq_varid (ncid, "t0", t0id)       
    if (istatus /= 0) then
       print *, "nf_inq_varid t0 error: ", nf_strerror(istatus),' ',iproc
    end if

    istatus = nf_inq_varid (ncid, "phi_r", phir_id)
    if (istatus /= 0) then
       print *, "nf_inq_varid phir error: ", nf_strerror(istatus),' ',iproc
    end if

    istatus = nf_inq_varid (ncid, "phi_i", phii_id)
    if (istatus /= 0) then
       print *, "nf_inq_varid phii error: ", nf_strerror(istatus),' ',iproc
    end if

    istatus = nf_inq_varid (ncid, "apar_r", aparr_id)
    if (istatus /= 0) then
       print *, "nf_inq_varid aparr error: ", nf_strerror(istatus),' ',iproc
    end if

    istatus = nf_inq_varid (ncid, "apar_i", apari_id)
    if (istatus /= 0) then
       print *, "nf_inq_varid apari error: ", nf_strerror(istatus),' ',iproc
    end if

    istatus = nf_inq_varid (ncid, "aperp_r", aperpr_id)
    if (istatus /= 0) then
       print *, "nf_inq_varid aperpr error: ", nf_strerror(istatus),' ',iproc
    end if

    istatus = nf_inq_varid (ncid, "aperp_i", aperpi_id)
    if (istatus /= 0) then
       print *, "nf_inq_varid aperpi error: ", nf_strerror(istatus),' ',iproc
    end if

    istatus = nf_inq_varid (ncid, "gr", gr_id)
    if (istatus /= 0) then
       print *, "nf_inq_varid gr error: ", nf_strerror(istatus),' ',iproc
    end if

    istatus = nf_inq_varid (ncid, "gi", gi_id)
    if (istatus /= 0) then
       print *, "nf_inq_varid gi error: ", nf_strerror(istatus),' ',iproc
    end if

    istatus = nf_get_var_double (ncid, t0id, tmp1)
    if (istatus /= 0) then
       print *, "nf_get_var_double t0 error: ", nf_strerror(istatus),' ',iproc
    end if
    t0 = tmp1
           
    n_elements = g_lo%ulim_proc - g_lo%llim_proc + 1

    if (n_elements <= 0) then
       phinew = 0.
       aparnew = 0.
       aperpnew = 0.
       return
    endif

    allocate (tmpr(2*ntgrid+1, 2, n_elements))
    allocate (tmpi(2*ntgrid+1, 2, n_elements))
    tmpr = 0.;  tmpi = 0.
       
    istatus = nf_get_vara_double (ncid, gr_id, (/ 1, 1, g_lo%llim_proc /), &
         (/ 2*ntgrid + 1, 2, n_elements /), tmpr)
    if (istatus /= 0) then
       print *, "nf_get_vara_double gr error: ", nf_strerror(istatus),' ',iproc
    end if
    
    istatus = nf_get_vara_double (ncid, gi_id,  (/ 1, 1, g_lo%llim_proc /), &
         (/ 2*ntgrid + 1, 2, n_elements /), tmpi)
    if (istatus /= 0) then
       print *, "nf_get_vara_double gi error: ", nf_strerror(istatus),' ',iproc
    end if

    g = cmplx(tmpr, tmpi)

    deallocate (tmpr)
    deallocate (tmpi)

    istatus = nf_get_var_double (ncid, phir_id, ftmpr)
    if (istatus /= 0) then
       print *, "nf_get_var_double phir error: ", nf_strerror(istatus),' ',iproc
    end if

    istatus = nf_get_var_double (ncid, phii_id, ftmpi)
    if (istatus /= 0) then
       print *, "nf_get_var_double phii error: ", nf_strerror(istatus),' ',iproc
    end if

    phi = 0.
    phinew = cmplx(ftmpr, ftmpi)

    istatus = nf_get_var_double (ncid, aparr_id, ftmpr)
    if (istatus /= 0) then
       print *, "nf_get_var_double aparr error: ", nf_strerror(istatus),' ',iproc
    end if

    istatus = nf_get_var_double (ncid, apari_id, ftmpi)
    if (istatus /= 0) then
       print *, "nf_get_var_double apari error: ", nf_strerror(istatus),' ',iproc
    end if

    apar = 0.
    aparnew = cmplx(ftmpr, ftmpi)

    istatus = nf_get_var_double (ncid, aperpr_id, ftmpr)
    if (istatus /= 0) then
       print *, "nf_get_var_double aperpr error: ", nf_strerror(istatus),' ',iproc
    end if

    istatus = nf_get_var_double (ncid, aperpi_id, ftmpi)
    if (istatus /= 0) then
       print *, "nf_get_var_double aperpi error: ", nf_strerror(istatus),' ',iproc
    end if

    aperp = 0.
    aperpnew = cmplx(ftmpr, ftmpi)

    istatus = nf_close (ncid)       
    if (istatus /= 0) then
       print *, "nf_close error: ", nf_strerror(istatus),' ',iproc
    end if

  end subroutine gs2_restore_one

  subroutine gs2_read_aminv (aminv)

    use fields_arrays, only: nidx
    use mp, only: proc0, iproc, nproc
    use gs2_layouts, only: f_lo
    implicit none

    character (305) :: file_proc
    character (5) :: suffix
    complex, dimension (:,f_lo%llim_proc:) :: aminv

    include 'netcdf.inc'
    integer :: nidxid, floid
    integer :: aminvr_id, aminvi_id
    integer :: n_f_elements, istatus, ncid
    double precision, dimension(nidx,f_lo%llim_proc:f_lo%ulim_alloc) :: tmpr, tmpi
    integer :: iflo, i, th, h, t, u

    n_f_elements = f_lo%ulim_proc-f_lo%llim_proc+1
    if (n_f_elements <= 0) then
       aminv = 0.
       return
    end if
    
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
    
    istatus = nf_open (file_proc, 0, ncid)
    if (istatus /= 0) then
       print *, "nf_open error: ", nf_strerror(istatus),' ',iproc
    end if
    
    istatus = nf_inq_dimid (ncid, "impidx", nidxid)
    if (istatus /= 0) then
       print *, "nf_inq_dimid nidx error: ", nf_strerror(istatus),' ',iproc
    end if

    istatus = nf_inq_dimid (ncid, "flo", floid)
    if (istatus /= 0) then
       print *, "nf_inq_dimid flo error: ", nf_strerror(istatus),' ',iproc
    end if

    istatus = nf_inq_dimlen (ncid, nidxid, i)       
    if (istatus /= 0) then
       print *, "nf_inq_dimlen nidx error: ", nf_strerror(istatus),' ',iproc
    end if
    if (i /= nidx) write(*,*) 'aminv restart error: ',i,' ',nidx

    istatus = nf_inq_dimlen (ncid, floid, i)       
    if (istatus /= 0) then
       print *, "nf_inq_dimlen flo error: ", nf_strerror(istatus),' ',iproc
    end if
    if (i /= f_lo%ulim_proc-f_lo%llim_proc+1) write(*,*) 'Restart error: flo=? ',i,' : ',iproc

    istatus = nf_inq_varid (ncid, "aminvr", aminvr_id)
    if (istatus /= 0) then
       print *, "nf_inq_varid aminvr error: ", nf_strerror(istatus),' ',iproc
    end if

    istatus = nf_inq_varid (ncid, "aminvi", aminvi_id)
    if (istatus /= 0) then
       print *, "nf_inq_varid aminvi error: ", nf_strerror(istatus),' ',iproc
    end if

    tmpr = 0.; tmpi = 0.
    istatus = nf_get_var_double (ncid, aminvi_id, tmpi)
    if (istatus /= 0) then
       print *, "nf_get_var_double aminvi error: ", nf_strerror(istatus),' ',iproc
    end if

    istatus = nf_get_var_double (ncid, aminvr_id, tmpr)
    if (istatus /= 0) then
       print *, "nf_get_var_double aminvr error: ", nf_strerror(istatus),' ',iproc
    end if

    aminv = cmplx(tmpr, tmpi)

    istatus = nf_close (ncid)       
    if (istatus /= 0) then
       print *, "nf_close error: ", nf_strerror(istatus),' ',iproc
    end if

  end subroutine gs2_read_aminv

  subroutine init_save (file)
    character(300), intent (in) :: file
    
    restart_file = file

  end subroutine init_save

end module gs2_save
