program combine
  
  use command_line
  implicit none
  include 'netcdf.inc'
  character (500) :: run_name
  character (5) :: suff
  integer :: i, j, ierr, th, h, t, u, num_files
  integer :: istatus, ncid, gloid, g_tot, old_ncid, old_gloid
  integer :: thetaid, signid, t0id, nt, kyid, kxid, old_kxid, old_kyid
  integer :: old_thetaid, old_signid, ng, naky, ntheta0
  integer :: gr_id, gi_id, old_gr_id, old_gi_id
  integer :: phir_id, phii_id, aparr_id, apari_id, aperpr_id, aperpi_id
  integer :: ophir_id, ophii_id, oaparr_id, oapari_id, oaperpr_id, oaperpi_id
  double precision :: tmp1
  double precision, dimension(:,:,:), allocatable :: gr, gi, ftmpr, ftmpi

  if (iargc() /= 0) then
     call cl_getarg (1, run_name, 500, ierr)
     if (ierr /= 0) then
        print *, "Error getting run name."
        stop
     endif
  endif

  write(*,*) 'Run name is ',trim(run_name)

  do num_files = 0, 9999
     suff = suffix(num_files)
     open (unit=32, file=trim(run_name)//suff, &
          form='unformatted', status='old', err=100)
     close(32)
  enddo

  write(*,*) 'Too many files!'
  stop

100 write(*,*) 'Number of files = ',num_files

  g_tot = 0
  do i = 0, num_files-1
     suff = suffix(i)
     istatus = nf_open (trim(run_name)//suff, 0, ncid)
     istatus = nf_inq_dimid (ncid, "glo", gloid)
     istatus = nf_inq_dimlen (ncid, gloid, j)
     istatus = nf_close (ncid)
     g_tot = g_tot + j     
  enddo
     
  write(*,*) 'Total length of g = ',g_tot

! Get size of theta grid and time
  istatus = nf_open (trim(run_name)//'.0000', 0, ncid)
  istatus = nf_inq_dimid (ncid, "theta", thetaid)
  istatus = nf_inq_dimlen (ncid, thetaid, nt)
  
  istatus = nf_inq_dimid (ncid, "aky", kyid)
  istatus = nf_inq_dimlen (ncid, kyid, naky)

  istatus = nf_inq_dimid (ncid, "akx", kxid)
  istatus = nf_inq_dimlen (ncid, kxid, ntheta0)

  istatus = nf_inq_varid (ncid, "t0", t0id)
  istatus = nf_get_var_double (ncid, t0id, tmp1)
  istatus = nf_close (ncid)

! Start to write new file: 
  istatus = nf_create (trim(run_name), 0, ncid)
  istatus = nf_def_dim (ncid, "theta", nt, thetaid)
  istatus = nf_def_dim (ncid, "sign", 2, signid)
  istatus = nf_def_dim (ncid, "glo", g_tot, gloid)
  istatus = nf_def_dim (ncid, "aky", naky, kyid)
  istatus = nf_def_dim (ncid, "akx", ntheta0, kxid)
  
  istatus = nf_def_var (ncid, "t0", NF_DOUBLE, 0, 0, t0id)
  istatus = nf_def_var (ncid, "gr", NF_DOUBLE, 3, &
       (/ thetaid, signid, gloid /), gr_id)
  istatus = nf_def_var (ncid, "gi", NF_DOUBLE, 3, &
       (/ thetaid, signid, gloid /), gi_id)
  istatus = nf_def_var (ncid, "phi_r", NF_DOUBLE, 3, &
       (/ thetaid, kyid, kxid /), phir_id)
  istatus = nf_def_var (ncid, "phi_i", NF_DOUBLE, 3, &
       (/ thetaid, kyid, kxid /), phii_id)

  istatus = nf_def_var (ncid, "apar_r", NF_DOUBLE, 3, &
       (/ thetaid, kyid, kxid /), aparr_id)
  istatus = nf_def_var (ncid, "apar_i", NF_DOUBLE, 3, &
       (/ thetaid, kyid, kxid /), apari_id)

  istatus = nf_def_var (ncid, "aperp_r", NF_DOUBLE, 3, &
       (/ thetaid, kyid, kxid /), aperpr_id)
  istatus = nf_def_var (ncid, "aperp_i", NF_DOUBLE, 3, &
       (/ thetaid, kyid, kxid /), aperpi_id)

  istatus = nf_enddef (ncid)
  istatus = nf_put_var_double (ncid, t0id, tmp1)

! Read them in one at a time and write into the new file.

  j = 1
  do i = 0, num_files - 1
     suff = suffix(i)
     istatus = nf_open (trim(run_name)//suff, 0, old_ncid)
     istatus = nf_inq_dimid (old_ncid, "theta", old_thetaid)
     istatus = nf_inq_dimid (old_ncid, "sign", old_signid)
     istatus = nf_inq_dimid (old_ncid, "glo", old_gloid)
     istatus = nf_inq_dimid (old_ncid, "akx", old_kxid)
     istatus = nf_inq_dimid (old_ncid, "aky", old_kyid)
     istatus = nf_inq_dimlen (old_ncid, old_gloid, ng)
     
     allocate (gr(nt, 2, ng))
     allocate (gi(nt, 2, ng))

     istatus = nf_inq_varid (ncid, "gr", old_gr_id)
     istatus = nf_inq_varid (ncid, "gi", old_gi_id)

     istatus = nf_get_var_double (old_ncid, old_gr_id, gr)
     istatus = nf_get_var_double (old_ncid, old_gi_id, gi)

     if (i == 0) then  ! only do the fields once
        allocate (ftmpr(nt, naky, ntheta0))
        allocate (ftmpi(nt, naky, ntheta0))
        istatus = nf_inq_varid (ncid, "phi_r", ophir_id)
        istatus = nf_inq_varid (ncid, "phi_i", ophii_id)
        istatus = nf_inq_varid (ncid, "apar_r", oaparr_id)
        istatus = nf_inq_varid (ncid, "apar_i", oapari_id)
        istatus = nf_inq_varid (ncid, "aperp_r", oaperpr_id)
        istatus = nf_inq_varid (ncid, "aperp_i", oaperpi_id)

        istatus = nf_get_var_double (old_ncid, ophir_id, ftmpr)
        istatus = nf_get_var_double (old_ncid, ophii_id, ftmpi)

        istatus = nf_put_var_double (ncid, phir_id, ftmpr)
        istatus = nf_put_var_double (ncid, phii_id, ftmpi)

        istatus = nf_get_var_double (old_ncid, oaparr_id, ftmpr)
        istatus = nf_get_var_double (old_ncid, oapari_id, ftmpi)

        istatus = nf_put_var_double (ncid, aparr_id, ftmpr)
        istatus = nf_put_var_double (ncid, apari_id, ftmpi)

        istatus = nf_get_var_double (old_ncid, oaperpr_id, ftmpr)
        istatus = nf_get_var_double (old_ncid, oaperpi_id, ftmpi)

        istatus = nf_put_var_double (ncid, aperpr_id, ftmpr)
        istatus = nf_put_var_double (ncid, aperpi_id, ftmpi)
        deallocate (ftmpr)
        deallocate (ftmpi)
     endif
        
     istatus = nf_close (old_ncid)

     istatus = nf_put_vara_double (ncid, gr_id, (/ 1, 1, j /), &
          (/ nt, 2, ng /), gr)
     istatus = nf_put_vara_double (ncid, gi_id, (/ 1, 1, j /), &
          (/ nt, 2, ng /), gi)

     deallocate (gr)
     deallocate (gi)

     j = j + ng
  enddo

  istatus = nf_close (ncid)


contains

  function suffix (i)
    
    character (5) :: suffix
    integer :: i, th, h, t, u
    
    th = i / 1000
    h = (i - th * 1000) / 100
    t = (i - th * 1000 - h * 100) / 10
    u = (i - th * 1000 - h * 100 - t * 10)
    suffix = '.'//achar(48+th)//achar(48+h)//achar(48+t)//achar(48+u)
    
  end function suffix
    
end program combine
