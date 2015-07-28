! lowflow terms include higher-order corrections to GK equation
! such as parallel nonlinearity that require derivatives in v-space.
! most efficient way to take these derivatives is to go from g_lo to le_lo,
! i.e., bring all energies and lambdas onto each processor
# ifdef LOWFLOW

# ifndef USE_LE_LAYOUT
# define USE_LE_LAYOUT on
# endif

# endif

module le_derivatives

  implicit none

  public :: vspace_derivatives

  private

contains

  subroutine vspace_derivatives (g, gold, g1, phi, apar, bpar, phinew, aparnew, bparnew, diagnostics)

    use redistribute, only: gather, scatter
    use dist_fn_arrays, only: c_rate, g_adjust
    use gs2_layouts, only: g_lo, le_lo
    use gs2_time, only: code_dt
    use run_parameters, only: fphi, fbpar
    use theta_grid, only: ntgrid
    use le_grids, only: integrate_moment, nxi, negrid, g2le
    use collisions, only: solfp1, colls, adjust, heating, hyper_colls
# ifdef LOWFLOW
    use gs2_layouts, only: ig_idx, is_idx
    use theta_grid, only: gradpar, theta
    use le_grids, only: ixi_to_il, ixi_to_isgn, forbid, write_mpdist_le, write_mpdist, jend
    use le_grids, only: xi, speed, energy
    use species, only: spec
    use lowflow, only: dphidth
# endif

    implicit none
    
    complex, dimension (-ntgrid:,:,:,g_lo%llim_proc:), intent (in out) :: g, gold, g1
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, bpar, phinew, aparnew, bparnew
    integer, optional, intent (in) :: diagnostics

# ifdef LOWFLOW
    integer :: ile, ig, ie, ixi, isgn, il, is
    integer, dimension (2) :: i
    real :: dvp, vp0
    complex, dimension (:,:), allocatable :: gtmp
# endif
    complex, dimension (:,:,:), allocatable :: gle
    complex, dimension (:,:,:), allocatable :: gc1, gc2, gc3
    logical :: heating_flag

    heating_flag = heating .and. present(diagnostics)

# ifdef USE_LE_LAYOUT

# ifndef LOWFLOW
    if (colls) then
# endif

       if (adjust) then
          call g_adjust (g, phinew, bparnew, fphi, fbpar)
          if (heating_flag) call g_adjust (gold, phi, bpar, fphi, fbpar)
       end if
       if (heating_flag) then
          allocate (gc3(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
          gc3 = g
       end if
       
       allocate (gle(nxi+1,negrid+1,le_lo%llim_proc:le_lo%ulim_alloc)) ; gle = 0.

       ! map data from g_layout to le_layout
       call gather (g2le, g, gle)
       
# ifdef LOWFLOW

       allocate (gtmp(nxi+1,negrid+1)) ; gtmp = 0.

       ! use semi-lagrange scheme to evaluate dg/dvpa term
       do ile = le_lo%llim_proc, le_lo%ulim_proc
          gtmp = gle(:,:,ile)
          ig = ig_idx(le_lo,ile)
          is = is_idx(le_lo,ile)
          dvp = spec(is)%zstm*dphidth(ig+ntgrid+1)*gradpar(ig)*0.5*code_dt
          do ie = 1, negrid
             do ixi = 1, nxi
                il = ixi_to_il(ig,ixi)
                isgn = ixi_to_isgn(ig,ixi)
                if (.not. forbid(ig,il)) &
                     call get_gvpa (gtmp, dvp, ig, il, ixi, ie, isgn, gle(ixi,ie,ile))
             end do
          end do
       end do

       deallocate (gtmp)

       if (colls) then
# endif
       ! update distribution function to take into account collisions
       if (present(diagnostics)) then
          call solfp1 (gle, diagnostics)
       else
          call solfp1 (gle)
       end if

# ifdef LOWFLOW
       end if
# endif

       ! remap from le_layout to g_layout
       call scatter (g2le, gle, g)

       deallocate (gle)

       if (heating_flag) then
          ! form (h_i+1 + h_i)/2 * C(h_i+1) and integrate.  
          gc3 = 0.5*conjg(g+gold)*(g-gc3)/code_dt
          
          call integrate_moment (gc3, c_rate(:,:,:,:,3))
          
          deallocate (gc3)
       end if

       if (adjust) then
          call g_adjust (g, phinew, bparnew, -fphi, -fbpar)
          if (heating_flag) call g_adjust (gold, phi, bpar, -fphi, -fbpar)
       end if

# ifndef LOWFLOW
    end if
# endif

# else

    if (colls) then

       if (heating_flag) then
          allocate (gc1(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
          if (hyper_colls) then
             allocate (gc2(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
          else
             allocate (gc2(1,1,1))
          end if
          allocate (gc3(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       else
          allocate (gc1(1,1,1))
          allocate (gc2(1,1,1))
          allocate (gc3(1,1,1))
       end if
       gc1 = 0. ; gc2 = 0. ; gc3 = 0.

       if (adjust) then
          call g_adjust (g, phinew, bparnew, fphi, fbpar)
          if (heating_flag) call g_adjust (gold, phi, bpar, fphi, fbpar)
       end if

       if (heating_flag) gc3 = g

       ! update distribution function to take into account collisions
       if (present(diagnostics)) then
          call solfp1 (g, g1, gc1, gc2, diagnostics)
       else
          call solfp1 (g, g1, gc1, gc2)
       end if

       if (heating_flag) then
          call integrate_moment (gc1, c_rate(:,:,:,:,1))
          deallocate (gc1)
          
          if (hyper_colls) call integrate_moment (gc2, c_rate(:,:,:,:,2))
          deallocate (gc2)
          
          ! form (h_i+1 + h_i)/2 * C(h_i+1) and integrate.  
          gc3 = 0.5*conjg(g+gold)*(g-gc3)/code_dt
          
          call integrate_moment (gc3, c_rate(:,:,:,:,3))
          
          deallocate (gc3)
       end if

       if (adjust) then
          call g_adjust (g, phinew, bparnew, -fphi, -fbpar)
          if (heating_flag) call g_adjust (gold, phi, bpar, -fphi, -fbpar)
       end if
    end if

# endif

  contains

    subroutine get_gvpa (g_in, dv, ig0, il0, ixi0, ie0, isgn0, g_out)

      use theta_grid, only: bmag
      use le_grids, only: speed, energy, al, xi, negrid, jend, sgn

      implicit none

      complex, dimension (:,:), intent (in) :: g_in
      real, intent (in) :: dv
      integer, intent (in) :: ig0, il0, ixi0, ie0, isgn0
      complex, intent (out) :: g_out

      integer :: ie, ie_low, ie_up, il, il_low, il_up, ix_low, ix_up, isgn
      real :: vp0, vp, v, x, lam, p
      logical :: v_finished, x_finished

      if (abs(dv) < epsilon(0.)) then

         g_out = g_in(ixi0,ie0)

      else

         ! vp0 is the parallel speed at current time level
         vp0 = xi(ig0,ixi0)*speed(ie0)
         ! vp is the parallel speed at previous time level
         vp = vp0 + dv
         
         ! v is the speed corresponding to the given (vpar,vperp) point
         v = sqrt(vp**2 + bmag(ig0)*al(il0)*energy(ie0))
         ! x is magnitude of cos(pitch-angle)
         x = vp/v
         ! p is pitch-angle
         p = acos(x)
         ! lam is mu B / E
         lam = (1.0 - x**2) / bmag(ig0)
         
         v_finished = .false. ; x_finished = .false.
         
         ! find the cell in v-xi space where (vp,vperp) resides
         ! dv positive or negative influences how to search for the correct cell
         if (abs(vp) > abs(vp0)) then
            
            if (ie0 < negrid) then
               ! energy of vp will be greater than or equal to that of vp0
               do ie = ie0+1, negrid
                  ! find the two speeds on the grid that bracket (vp,vperp)
                  if (speed(ie) > v) then
                     ie_up = ie
                     ie_low = ie-1
                     v_finished = .true.
                     exit
                  end if
               end do
            end if
            
            ! if vp falls above the largest speed on the grid, set ie_up = negrid+1
            ! to let interp_g know that it must extrapolate in speed to obtain g(vp,vperp)
            if (.not. v_finished) then
               ie_up = negrid+1 ; ie_low = negrid
            end if
            
            if (il0 > 1) then
               ! pitch-angle goes to larger absolute value, corresponding to smaller il
               do il = il0-1, 1, -1
                  if (al(il) < lam) then
                     il_low = il
                     il_up = il+1
                     x_finished = .true.
                     exit
                  end if
               end do
            end if
            
            ! if vp falls below smallest lambda on grid, set il_low=0
            ! to let interp_g that it must extrapolate in lambda to obtain g(vp)
            if (.not. x_finished) then
               il_low = 0 ; il_up = 1
            end if
            
         else
            
            if (ie0 > 1) then
               ! energy of vp will be less than or equal to that of vp0
               do ie = ie0-1, 1, -1
                  ! find the two speeds on grid that bracket (vp,vperp)
                  if (speed(ie) < v) then
                     ie_low = ie
                     ie_up = ie+1
                     v_finished = .true.
                     exit
                  end if
               end do
            end if
            
            if (.not. v_finished) then
               ie_low = 0 ; ie_up = 1
            end if
            
            if (il0 == jend(ig0)) then
               write (*,*) 'Error in get_gvpa: il0=jend(ig0) should not be possible here.'
               stop
            end if
            
            ! pitch-angle goes to smaller absolute value, corresponding to larger il
            do il = il0+1, jend(ig0)
               ! find the two lambda values that bracket the lambda of (vp,vperp)
               if (al(il) > lam) then
                  il_up = il
                  il_low = il-1
                  x_finished = .true.
                  exit
               end if
            end do
            
            if (.not. x_finished) then
               write (*,*) 'Error in get_gvpa: could not bracket il'
               stop
            end if
            
         end if

         ! determine if vpar at which we want g is positive or negative
         if (abs(abs(vp)-vp) > epsilon(0.0)) then
!         if (vp < 0) then
            isgn = 2
            ix_low = 2*jend(ig0)-il_low ; ix_up = 2*jend(ig0)-il_up
         else
            isgn = 1
            ix_low = il_low ; ix_up = il_up
         end if
         
         ! bilinear interpolation using 4 grid points closest to (vp,vperp)
         call interp_g (ig0, isgn, il_low, il_up, ix_low, ix_up, ie_low, ie_up, v, x, lam, p, g_in, g_out)

      end if

    end subroutine get_gvpa

    subroutine interp_g (ig0, isgn, il_low, il_up, ix_low, ix_up, ie_low, ie_up, v0, x0, lam0, p0, g, gint)

      use theta_grid, only: ntgrid
      use le_grids, only: speed, xi, negrid, sgn, nlambda, al

      implicit none

      integer, intent (in) :: ig0, isgn, il_low, il_up, ix_low, ix_up, ie_low, ie_up
      real, intent (in) :: v0, x0, lam0, p0
      complex, dimension (:,:), intent (in) :: g
      complex, intent (out) :: gint

      real :: d1, d2, d3, d4, d5, dtot, dg, xpt
      real, dimension (:), allocatable :: var

      allocate (var(size(xi,2)))
      var = xi(ig0,:) ; xpt = x0
!      var = al ; xpt = lam0
!      var(1:nlambda) = al ; var(2*nlambda:nlambda+1:-1) = al ; xpt = lam0
!      var = acos(xi(ig0,:)) ; xpt = p0

      ! check for special cases

      ! if vp falls outside the energy grid, must extrapolate
      if (ie_up == negrid+1) then
         ! if vp also falls outside lambda grid, do 2D extrapolation
         if (il_low == 0) then
            ! use dg/dvperp = 0 as vperp -> 0
            gint = g(ix_low+sgn(isgn),negrid)*exp(-v0**2+speed(negrid)**2)
         else
            gint = (g(ix_low,negrid)*(var(ix_up)-xpt) + g(ix_up,negrid)*(xpt-var(ix_low))) / (var(ix_up)-var(ix_low)) &
                 * exp(-v0**2 + speed(negrid)**2)
         end if
      else if (il_low == 0) then
         ! handle the case where vp lies between v-grid point but outside al-grid points

         ! use dg/dvperp = 0 as vperp -> 0
         dtot = speed(ie_up) - speed(ie_low)
         d1 = speed(ie_up) - v0 ; d2 = v0 - speed(ie_low)

         gint = (g(ix_low+sgn(isgn),ie_up)*d2 + g(ix_low+sgn(isgn),ie_low)*d1) / dtot
      else if (ie_low == 0) then
         ! handle the case where vp lies between al-grid points but below smallest v-grid point

         ! use dg/dvperp = 0 as vperp -> 0
         dtot = (var(ix_up) - var(ix_low))
         d3 = var(ix_up) - xpt ; d4 = xpt - var(ix_low)

         gint = (g(ix_low,1)* d3 + g(ix_up,1)*d4) / dtot
      else
         ! linear interpolation to obtain g(vp) between four grid points
         dtot = abs(speed(ie_up) - speed(ie_low)) * abs(var(il_up) - var(il_low))
         d1 = abs(speed(ie_up) - v0) ; d2 = abs(speed(ie_low) - v0)
         d3 = abs(var(il_up) - abs(xpt)) ; d4 = abs(var(il_low) - abs(xpt))
         gint = (d4*(d2*g(ix_up,ie_up) + d1*g(ix_up,ie_low)) &
              + d3*(d2*g(ix_low,ie_up) + d1*g(ix_low,ie_low))) / dtot
      end if

      deallocate (var)

    end subroutine interp_g

  end subroutine vspace_derivatives

end module le_derivatives
