!> A module to calculate the properties of the
!! antenna used to drive turbulence 
module diagnostics_antenna
  use diagnostics_config, only: diagnostics_type
  use diagnostics_create_and_write
  use kt_grids, only: ntheta0, naky, kperp2
  use theta_grid, only: ntgrid
  implicit none
  real, dimension(:,:,:), allocatable ::  j_ext_hist
  contains
    subroutine init_diagnostics_antenna(gnostics)
      type(diagnostics_type), intent(in) :: gnostics
      if (gnostics%write_jext) allocate (j_ext_hist(ntheta0, naky,0:gnostics%navg-1)) 
    end subroutine init_diagnostics_antenna
    subroutine finish_diagnostics_antenna(gnostics)
      type(diagnostics_type), intent(in) :: gnostics
      if (gnostics%write_jext) deallocate(j_ext_hist)
    end subroutine finish_diagnostics_antenna
    subroutine write_jext(gnostics)
      use kt_grids, only: ntheta0, naky
      use mp, only: proc0
      implicit none
      type(diagnostics_type), intent(in) :: gnostics
      real:: t
      integer:: istep
      integer :: ik, it
      !GGH J_external
      real, dimension(:,:), allocatable ::  j_ext
       !integer :: ig,ik, it                             !Indices
    !real :: fac2                                     !Factor
    !real, dimension (:), allocatable :: wgt
    

    !Get z-centered j_ext at current time
      istep = gnostics%istep
      t = gnostics%user_time

      allocate (j_ext(ntheta0, naky)); j_ext=0.
      !if (.not.proc0) return
      call calc_jext(gnostics,j_ext)

      call create_and_write_variable(gnostics, gnostics%rtype, "antenna_j_ext", "XYt", &
        "Time averaged external current in the antenna, real(kperp^2 A_antenna), &
        &  as a function of kx and ky", "radians", j_ext)


      if (proc0 .and. gnostics%write_ascii) then 
        do ik=1,naky
           do it = 1, ntheta0
              if (j_ext(it,ik) .ne. 0.) then
                 write (unit=gnostics%ascii_files%jext, fmt="(es12.4,i4,i4,es12.4)")  &
                      t,it,ik,j_ext(it,ik)
              endif
           enddo
        enddo
      end if
      deallocate(j_ext)
    end subroutine write_jext
    !> A subroutine to calculate the time-averaged antenna current
    !! j_ext = kperp^2 A_antenna
   subroutine calc_jext (gnostics, j_ext)
      use mp, only: proc0
      use dist_fn, only: get_jext
      use antenna, only: antenna_apar
      use volume_averages, only: average_theta
      implicit none
      !Passed
      type(diagnostics_type), intent(in) :: gnostics
      integer:: istep
      real, dimension(:,:), intent(out) ::  j_ext
      complex, dimension(:,:,:), allocatable :: j_extz
      !Local 
      integer :: i

      istep = gnostics%istep

      !call get_jext(j_ext)    
      allocate (j_extz(-ntgrid:ntgrid, ntheta0, naky)) ; j_extz = 0.
      call antenna_apar (kperp2, j_extz)       
      call average_theta(real(j_extz), j_ext, gnostics%distributed)
      deallocate (j_extz)

      
      !Do averages with a history variable
      if (proc0) then
         !Save variable to history
         if (gnostics%navg > 1) then
            if (istep > 1) &
                 j_ext_hist(:,:,mod(istep,gnostics%navg))= j_ext(:,:)

            !Use average of history
            if (istep >= gnostics%navg) then
               j_ext=0.
               do i=0,gnostics%navg-1
                  j_ext(:,:) = j_ext(:,:) + j_ext_hist(:,:,i)/ real(gnostics%navg)
               end do
            end if
         end if
      end if

    end subroutine calc_jext
end module diagnostics_antenna
