module fixed

  implicit none

  integer, parameter :: nproc = 16

  integer, parameter :: ntgrid = 16
  integer, parameter :: naky = 1
  integer, parameter :: ntheta0 = 1
  integer, parameter :: nlambda = 27
  integer, parameter :: negrid = 16
  integer, parameter :: nspec = 1
  integer, parameter :: g_block = (naky*ntheta0*nlambda*negrid*nspec-1)/nproc+1
  integer, parameter :: g_llim = 1
  integer, parameter :: g_ulim = g_block

  integer, parameter :: gint_block = (naky*ntheta0*negrid*nspec-1)/nproc+1
  integer, parameter :: gint_llim = 1
  integer, parameter :: gint_ulim = gint_block

  integer, parameter :: geint_block = (naky*ntheta0*nspec-1)/nproc+1
  integer, parameter :: geint_llim = 1
  integer, parameter :: geint_ulim = geint_block

  integer, parameter :: lz_block = ((2*ntgrid+1)*naky*ntheta0*negrid*nspec-1)/nproc+1
  integer, parameter :: lz_llim = 1
  integer, parameter :: lz_ulim = lz_block

  integer, parameter :: nny = 1  ! should be consistent with naky
  integer, parameter :: nnx = 1  ! should be consistent with ntheta0
  integer, parameter :: ndky = 1 ! should be consistent with naky
  integer, parameter :: nxny = nnx*nny
  integer, parameter :: nxnky = nnx*(nny/2+1)
  integer, parameter :: nia = negrid*nlambda*nspec/nproc

  integer, parameter :: accelx_block = (nnx*nny*nlambda*negrid*nspec-1)/nproc+1
  integer, parameter :: accelx_llim = 1
  integer, parameter :: accelx_ulim = accelx_block

  integer, parameter :: accel_block = (nnx*(nny/2+1)*nlambda*negrid*nspec-1)/nproc+1
  integer, parameter :: accel_llim = 1
  integer, parameter :: accel_ulim = accel_block

end module fixed
