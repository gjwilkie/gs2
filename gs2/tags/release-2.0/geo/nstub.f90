function ncopn ()
  integer ncopn

  ncopn = 1
  write(*,*) 'NETCDF not available.'
  write(*,*) 'stopping in ncopn (nstub.f90)'
  stop

end function ncopn

function ncdid ()
  integer ncdid

  ncdid = 1

end function ncdid

subroutine ncdinq

end subroutine ncdinq

function ncvid ()
  integer ncvid

  ncvid = 1

end function ncvid

subroutine ncvgt1
  
end subroutine ncvgt1

subroutine ncvgt
  
end subroutine ncvgt

subroutine ncclos

end subroutine ncclos
