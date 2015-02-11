module fields_parallelization
   implicit none
   private
   public :: field_k_local
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
    implicit none
    integer, intent(in) :: ik, it
    logical :: field_k_local

    ! This is temporary while the fields are being parallelised
    !write (*,*) 'it', it, ' nproc ', nproc, 'iproc ', iproc, 'mod', mod(it,nproc)
    !field_k_local = (mod(iproc,ntheta0) == it-1)
    field_k_local = (mod(it-1,nproc) == iproc)
    !field_k_local = proc0
  end function field_k_local

end module fields_parallelization
