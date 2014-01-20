module fields_parallelization
   use theta_grid, only: ntgrid
   use kt_grids, only: naky, ntheta0
   implicit none

contains
  !function field_k_local_allcover(ik,it)
    !use mp, only: iproc, nproc
    !integer, intent(in) :: ik, it
    !logical :: field_k_local_allcover

    !! This is temporary while the fields are being parallelised
    !!write (*,*) 'it', it, ' nproc ', nproc, 'iproc ', iproc, 'mod', mod(it,nproc)
    !field_k_local_allcover = (mod(iproc,ntheta0) == it-1)
    !!field_k_local = (mod(it,nproc) == iproc)
  !end function field_k_local_allcover

  function field_k_local(it,ik)
    use mp, only: iproc, nproc
    integer, intent(in) :: ik, it
    logical :: field_k_local

    ! This is temporary while the fields are being parallelised
    !write (*,*) 'it', it, ' nproc ', nproc, 'iproc ', iproc, 'mod', mod(it,nproc)
    !field_k_local = (mod(iproc,ntheta0) == it-1)
    field_k_local = (mod(it,nproc) == iproc)
  end function field_k_local

end module fields_parallelization
