! Modifications for using FFTW version 3:
! (c) The Numerical Algorithms Group (NAG) Ltd, 2009 
!                                 on behalf of the HECToR project


# include "define.inc"

module gs2_transforms

  use redistribute, only: redist_type
  use fft_work, only: fft_type

  implicit none

  public :: init_transforms
  public :: init_x_transform, init_zf, kz_spectrum
  public :: transform_x, transform_y, transform2
  public :: inverse_x, inverse_y, inverse2
  public :: quicksort2 !<DD>
  public :: quicksort3 !<DD>
  private

  interface transform_x
     module procedure transform_x5d
!     module procedure transform_x3d
!     module procedure transform_x1d
  end interface

  interface transform_y
     module procedure transform_y5d
!     module procedure transform_y3d
!     module procedure transform_y1d
  end interface

  interface transform2
     module procedure transform2_5d_accel
     module procedure transform2_5d
     module procedure transform2_4d
     module procedure transform2_3d
     module procedure transform2_2d
  end interface

  interface inverse_x
     module procedure inverse_x5d
!     module procedure inverse_x3d
!     module procedure inverse_x1d
  end interface

  interface inverse_y
     module procedure inverse_y5d
!     module procedure inverse_y3d
!     module procedure inverse_y1d
  end interface

  interface inverse2
     module procedure inverse2_5d_accel
     module procedure inverse2_5d
     module procedure inverse2_3d
     module procedure inverse2_2d
  end interface

  ! redistribution

  type (redist_type), save :: g2x, x2y

  ! fft

  type (fft_type) :: xf_fft, xb_fft, yf_fft, yb_fft, zf_fft
  type (fft_type) :: xf3d_cr, xf3d_rc

  logical :: xfft_initted = .false.

! accel will be set to true if the v layout is used AND the number of
! PEs is such that each PE has a complete copy of the x,y space --
! in that case, no communication is needed to evaluate the nonlinear
! terms
  logical :: accel = .false.

  logical, save, dimension (:), allocatable :: aidx  ! aidx == aliased index
  integer, save, dimension (:), allocatable :: ia, iak
  complex, save, dimension (:, :), allocatable :: fft, xxf
  complex, save, dimension (:, :, :), allocatable :: ag

contains

  subroutine init_transforms &
       (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx, ny, accelerated)
    use mp, only: nproc
    use gs2_layouts, only: init_gs2_layouts, opt_redist_init
    use gs2_layouts, only: pe_layout, init_accel_transform_layouts
    use gs2_layouts, only: init_y_transform_layouts
    use gs2_layouts, only: init_x_transform_layouts
    implicit none
    integer, intent (in) :: ntgrid, naky, ntheta0, nlambda, negrid, nspec
    integer, intent (in) :: nx, ny
    logical, intent (out) :: accelerated
    
    logical, save :: initialized = .false.
    character (1) :: char

! CMR, 12/2/2010:  return correct status of "accelerated" even if already initialised
    accelerated = accel

    !Early exit if possible
    if (initialized) return
    initialized = .true.

    call init_gs2_layouts

    call pe_layout (char)

    if (char == 'v' .and. mod (negrid*nlambda*nspec, nproc) == 0) then  
       accel = .true.
       call init_accel_transform_layouts (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx, ny)
    else
       !Recommended for p+log(p)>log(N) where p is number of processors and N is total number of mesh points
       !Could automate selection, though above condition is only fairly rough
       if (opt_redist_init) then
          call init_y_redist_local (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx, ny)
       else
          call init_y_redist (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx, ny)
       endif
    end if

! need these for movies
    call init_y_transform_layouts (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx, ny)
    call init_x_transform_layouts (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx)

    call init_y_fft (ntgrid)

    accelerated = accel

  end subroutine init_transforms

  subroutine init_x_transform (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx)

  !JH As of 7th December 2011 this routine is dead code and should be considered
  !JH for removal from the source.  
  !JH This functionality is presently implemented in subroutine init_y_fft 
  !JH present later in this file.

    use mp, only: nproc
    use fft_work, only: init_ccfftw
    use gs2_layouts, only: xxf_lo, pe_layout
    integer, intent (in) :: ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx

    logical, save :: initialized = .false.
    character (1) :: char

    integer :: nb_ffts

    if (initialized) return
    initialized = .true.

    call pe_layout (char)

    if (char == 'v' .and. mod (negrid*nlambda*nspec, nproc) == 0) then
       accel = .true.
    else
       call init_x_redist (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx)

# if FFT == _FFTW_

       call init_ccfftw (xf_fft,  1, xxf_lo%nx)
       call init_ccfftw (xb_fft, -1, xxf_lo%nx)

# elif FFT == _FFTW3_

       if (.not.allocated(xxf)) then
          allocate (xxf(xxf_lo%nx,xxf_lo%llim_proc:xxf_lo%ulim_alloc))
       endif

       ! number of ffts to be calculated
       nb_ffts = xxf_lo%ulim_proc - xxf_lo%llim_proc + 1

       call init_ccfftw (xf_fft,  1, xxf_lo%nx, nb_ffts, xxf)
       call init_ccfftw (xb_fft, -1, xxf_lo%nx, nb_ffts, xxf)
# endif

       xfft_initted = .true.
    end if

  end subroutine init_x_transform

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine init_y_fft (ntgrid)

    use gs2_layouts, only: xxf_lo, yxf_lo, accel_lo, accelx_lo, dealiasing, g_lo
    use fft_work, only: init_crfftw, init_rcfftw, init_ccfftw
    implicit none
    integer, intent (in) :: ntgrid

    logical :: initialized = .false.
    integer :: idx, i

    integer :: nb_ffts

    if (initialized) return
    initialized = .true.

    if (accel) then

! prepare for dealiasing 
       allocate (ia (accel_lo%nia))
       allocate (iak(accel_lo%nia))
       allocate (              aidx (accel_lo%llim_proc:accel_lo%ulim_alloc))
       allocate (ag(-ntgrid:ntgrid,2,accel_lo%llim_proc:accel_lo%ulim_alloc))
       aidx = .true.

       do i = accel_lo%llim_proc, accel_lo%ulim_proc
          if (dealiasing(accel_lo, i)) cycle
          aidx(i) = .false.
       end do

       do idx = 1, accel_lo%nia
          ia (idx) = accelx_lo%llim_proc + (idx-1)*accelx_lo%nxny
          iak(idx) = accel_lo%llim_proc  + (idx-1)*accel_lo%nxnky
       end do

       !JH FFTW plan creation for the accelerated 2d transforms
#if FFT == _FFTW_ 

       call init_crfftw (yf_fft,  1, accel_lo%ny, accel_lo%nx)
       call init_rcfftw (yb_fft, -1, accel_lo%ny, accel_lo%nx)

#elif FFT == _FFTW3_

       call init_crfftw (yf_fft,  1, accel_lo%ny, accel_lo%nx, &
            (2*accel_lo%ntgrid+1) * 2)
       call init_rcfftw (yb_fft, -1, accel_lo%ny, accel_lo%nx, &
            (2*accel_lo%ntgrid+1) * 2)

#endif

    
    else
       ! non-accelerated
       allocate (fft(yxf_lo%ny/2+1, yxf_lo%llim_proc:yxf_lo%ulim_alloc))
       if (.not.allocated(xxf))then
          allocate (xxf(xxf_lo%nx,xxf_lo%llim_proc:xxf_lo%ulim_alloc))
       endif

       if (.not. xfft_initted) then

          !JH FFTW plan creation for transform x5d and inverse
# if FFT == _FFTW_
          call init_ccfftw (xf_fft,  1, xxf_lo%nx)
          call init_ccfftw (xb_fft, -1, xxf_lo%nx)
# elif FFT == _FFTW3_
          ! number of ffts to be calculated
          !JH 7th December 2011
	  !JH xxf_lo%ulim_alloc is used here rather than xxf_lo%lulim_proc
	  !JH because there are situations where xxf_lo%llim_proc is greater 
	  !JH than xxf_lo%ulim_proc and that would create a negative number 
	  !JH of FFTs to be calculated.  However, xxf_lo%ulim_alloc is set
	  !JH to be xxf_lo%llim_proc in this situation, and that will give 
	  !JH 1 FFT to be calculated which the code can correctly undertake.
          nb_ffts = xxf_lo%ulim_alloc - xxf_lo%llim_proc + 1
          
          call init_ccfftw (xf_fft,  1, xxf_lo%nx, nb_ffts, xxf)
          call init_ccfftw (xb_fft, -1,  xxf_lo%nx, nb_ffts, xxf)
          
# endif
          xfft_initted = .true.
       end if

       !JH FFTW plan creation for transform y5d and inverse
# if FFT == _FFTW_
       call init_crfftw (yf_fft,  1, yxf_lo%ny)
       call init_rcfftw (yb_fft, -1, yxf_lo%ny)
#elif FFT == _FFTW3_
       ! number of ffts to be calculated
       !JH 7th December 2011
       !JH yxf_lo%ulim_alloc is used here rather than yxf_lo%lulim_proc
       !JH because there are situations where yxf_lo%llim_proc is greater 
       !JH than yxf_lo%ulim_proc and that would create a negative number 
       !JH of FFTs to be calculated.  However, yxf_lo%ulim_alloc is set
       !JH to be yxf_lo%llim_proc in this situation, and that will give 
       !JH 1 FFT to be calculated which the code can correctly undertake.
       nb_ffts = yxf_lo%ulim_alloc - yxf_lo%llim_proc + 1

       call init_crfftw (yf_fft,  1, yxf_lo%ny, nb_ffts)
       call init_rcfftw (yb_fft, -1,  yxf_lo%ny, nb_ffts)

#endif

    end if

  end subroutine init_y_fft

  subroutine init_x_redist (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx)
    use gs2_layouts, only: init_x_transform_layouts
    use gs2_layouts, only: g_lo, xxf_lo, gidx2xxfidx, proc_id, idx_local
    use gs2_layouts, only: opt_local_copy, layout
    
    use mp, only: nproc
    use redistribute, only: index_list_type, init_redist, delete_list
    use redistribute, only: set_redist_character_type, set_xxf_optimised_variables
    implicit none
    type (index_list_type), dimension(0:nproc-1) :: to_list, from_list
    integer, dimension(0:nproc-1) :: nn_to, nn_from
    integer, dimension (3) :: from_low, from_high
    integer, dimension (2) :: to_high
    integer :: to_low
    integer, intent (in) :: ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx
    logical :: initialized = .false.

    integer :: iglo, isign, ig, it, ixxf
    integer :: n, ip

!<DD>This routine can probably be optimised somewhat.
!in particular the large number of calls to the index lookup routines
!and idx_local can be reduced significantly by exploiting knowledge of
!the xxf layout and similar.
!For most cases init_x_redist_local is faster
    if (initialized) return
    initialized = .true.

    call init_x_transform_layouts &
         (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx)

    ! count number of elements to be redistributed to/from each processor
    nn_to = 0
    nn_from = 0
    do iglo = g_lo%llim_world, g_lo%ulim_world
       do isign = 1, 2
          do ig = -ntgrid, ntgrid
             call gidx2xxfidx (ig, isign, iglo, g_lo, xxf_lo, it, ixxf)
             if (idx_local(g_lo,iglo)) &
                  nn_from(proc_id(xxf_lo,ixxf)) = nn_from(proc_id(xxf_lo,ixxf)) + 1
             if (idx_local(xxf_lo,ixxf)) &
                nn_to(proc_id(g_lo,iglo)) = nn_to(proc_id(g_lo,iglo)) + 1
          end do
       end do
    end do

    do ip = 0, nproc-1
       if (nn_from(ip) > 0) then
          allocate (from_list(ip)%first(nn_from(ip)))
          allocate (from_list(ip)%second(nn_from(ip)))
          allocate (from_list(ip)%third(nn_from(ip)))
       end if
       if (nn_to(ip) > 0) then
          allocate (to_list(ip)%first(nn_to(ip)))
          allocate (to_list(ip)%second(nn_to(ip)))
       end if
    end do


    ! get local indices of elements distributed to/from other processors
    nn_to = 0
    nn_from = 0
    do iglo = g_lo%llim_world, g_lo%ulim_world
       do isign = 1, 2
          do ig = -ntgrid, ntgrid
             call gidx2xxfidx (ig, isign, iglo, g_lo, xxf_lo, it, ixxf)
             if (idx_local(g_lo,iglo)) then
                ip = proc_id(xxf_lo,ixxf)
                n = nn_from(ip) + 1
                nn_from(ip) = n
                from_list(ip)%first(n) = ig
                from_list(ip)%second(n) = isign
                from_list(ip)%third(n) = iglo
             end if
             if (idx_local(xxf_lo,ixxf)) then
                ip = proc_id(g_lo,iglo)
                n = nn_to(ip) + 1
                nn_to(ip) = n
                to_list(ip)%first(n) = it
                to_list(ip)%second(n) = ixxf
             end if
          end do
       end do
    end do

    from_low (1) = -ntgrid
    from_low (2) = 1
    from_low (3) = g_lo%llim_proc

    to_low = xxf_lo%llim_proc
    
    to_high(1) = xxf_lo%nx
    to_high(2) = xxf_lo%ulim_alloc

    from_high(1) = ntgrid
    from_high(2) = 2
    from_high(3) = g_lo%ulim_alloc

    call set_redist_character_type(g2x, 'g2x')
    call set_xxf_optimised_variables(opt_local_copy, naky, ntgrid, ntheta0, &
       nlambda, nx, xxf_lo%ulim_proc, g_lo%blocksize, layout)

    call init_redist (g2x, 'c', to_low, to_high, to_list, &
         from_low, from_high, from_list)

    call delete_list (to_list)
    call delete_list (from_list)

  end subroutine init_x_redist

!<DD>
  subroutine init_x_redist_local (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx)
    use gs2_layouts, only: init_x_transform_layouts
    use gs2_layouts, only: g_lo, xxf_lo, gidx2xxfidx, proc_id, idx_local
    use gs2_layouts, only: opt_local_copy, layout
    use gs2_layouts, only: ik_idx,il_idx,ie_idx,is_idx,idx, ig_idx, isign_idx
    use mp, only: nproc
    use mp, only: iproc
    use redistribute, only: index_list_type, init_redist, delete_list
    use redistribute, only: set_redist_character_type, set_xxf_optimised_variables
    implicit none
    type (index_list_type), dimension(0:nproc-1) :: to_list, from_list, sort_list, sto_list
    integer, dimension(0:nproc-1) :: nn_to, nn_from
    integer, dimension (3) :: from_low, from_high
    integer, dimension (2) :: to_high
    integer :: to_low
    integer,dimension(:),allocatable :: nn
    integer, intent (in) :: ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx
    logical :: initialized = .false.

    integer :: iglo, isign, ig, it, ixxf, it0, count
    integer :: n, ip
    integer :: glomesg_count, xxfmesg_count
    logical, parameter :: debug=.true.

    !Early exit if possible
    if (initialized) return
    initialized = .true.

    !Setup the xxf_lo layout object
    call init_x_transform_layouts &
         (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx)

    ! count number of elements to be redistributed to/from each processor
    nn_to = 0
    nn_from = 0

    !Here we loop over the whole domain (so this doesn't get cheaper with more cores)
    !This is required to ensure that we can associate send and receive message elements
    !i.e. so we know that we put the received data in the correct place.
    !However, as we know the order the iglo,isign,ig indices should increase we could
    !loop over our ixxf local range, calculate the corresponding iglo,isign and ig indices
    !and sort the ixxf messages on this to ensure they're in the correct order.
    !Either way when just counting how much data we are going to send and receive we don't
    !care about order so we can just loop over our local range!
    !First count the sends | g_lo-->xxf_lo
    !Protect against procs with no data
    if(g_lo%ulim_proc.ge.g_lo%llim_proc)then
       do iglo = g_lo%llim_proc, g_lo%ulim_alloc
          !Convert iglo,isign=1,ig=-ntgrid into ixxf
          ixxf=idx(xxf_lo,-ntgrid,1,ik_idx(g_lo,iglo),&
               il_idx(g_lo,iglo),ie_idx(g_lo,iglo),is_idx(g_lo,iglo))

          !Now loop over other local dimensions
          do isign = 1, 2
             do ig = -ntgrid, ntgrid
                !Increase the data send count for the proc which has the ixxf
                nn_from(proc_id(xxf_lo,ixxf))=nn_from(proc_id(xxf_lo,ixxf))+1

                !Increase ixxf using knowledge of the xxf_lo layout
                ixxf=ixxf+xxf_lo%naky
             enddo
          enddo
       enddo
    endif

    !Now count the receives | xxf_lo<--g_lo
    !Protect against procs with no data
    if(xxf_lo%ulim_proc.ge.xxf_lo%llim_proc)then
       do ixxf = xxf_lo%llim_proc, xxf_lo%ulim_alloc
          !Could split it (or x) domain into two parts to account for
          !difference in it (or x) meaning and order in g_lo and xxf_lo
          !but only interested in how much data we receive and not the
          !exact indices (yet) so just do 1-->ntheta0 (this is how many 
          !non-zero x's?)
          do it=1,g_lo%ntheta0
             !Convert ixxf,it indices into iglo, ig and isign indices
             iglo=idx(g_lo,ik_idx(xxf_lo,ixxf),it,il_idx(xxf_lo,ixxf),&
                  ie_idx(xxf_lo,ixxf),is_idx(xxf_lo,ixxf))

             !Increase the data to receive count for proc with this data
             !Note, we only worry about iglo and not isign/ig because we know that
             !in g_lo each proc has all ig and isign domain.
             !The xxf_lo domain contains ig and isign so we only need to add one
             !to the count for each ixxf.
             nn_to(proc_id(g_lo,iglo))=nn_to(proc_id(g_lo,iglo))+1
          enddo
       enddo
    endif

!<DD>Debug test, are we sending and receiving the same amount to ourselves? Remove when testing completed
!if(debug.and.nn_to(iproc).ne.nn_from(iproc)) print*,"ERROR: iproc ",iproc,"nn_from",nn_from(iproc),"nn_to",nn_to(iproc)

    !Now allocate storage for data mapping indices
    do ip = 0, nproc-1
       if (nn_from(ip) > 0) then
          allocate (from_list(ip)%first(nn_from(ip)))
          allocate (from_list(ip)%second(nn_from(ip)))
          allocate (from_list(ip)%third(nn_from(ip)))
       end if
       if (nn_to(ip) > 0) then
          allocate (to_list(ip)%first(nn_to(ip)))
          allocate (to_list(ip)%second(nn_to(ip)))
          !For sorting to_list later
          allocate (sort_list(ip)%first(nn_to(ip)))
       end if
    end do

    !Reinitialise count arrays to zero
    nn_to = 0
    nn_from = 0

    !First fill in the sending indices, these define the messages data order
    !Protect against procs with no data
    if(g_lo%ulim_proc.ge.g_lo%llim_proc)then
       do iglo=g_lo%llim_proc, g_lo%ulim_alloc
          !Convert iglo,isign=1,ig=-ntgrid into ixxf
          ixxf=idx(xxf_lo,-ntgrid,1,ik_idx(g_lo,iglo),&
               il_idx(g_lo,iglo),ie_idx(g_lo,iglo),is_idx(g_lo,iglo))

          !Now loop over other local dimensions
          do isign = 1, 2
             do ig = -ntgrid, ntgrid
                !Get proc id
                ip=proc_id(xxf_lo,ixxf)

                !Increment procs message counter
                n=nn_from(ip)+1
                nn_from(ip)=n

                !Store indices
                from_list(ip)%first(n) = ig
                from_list(ip)%second(n) = isign
                from_list(ip)%third(n) = iglo

                !We could send this information, transformed to the xxf layout to the proc.

                !Increment counter
                ixxf=ixxf+xxf_lo%naky
             enddo
          enddo
       enddo
    endif

    !Now lets fill in the receiving indices, these must match the messages data order
    !Protect against procs with no data
    if(xxf_lo%ulim_proc.ge.xxf_lo%llim_proc)then
       do ixxf = xxf_lo%llim_proc, xxf_lo%ulim_alloc
          !Get ig and isign indices
          ig=ig_idx(xxf_lo,ixxf)
          isign=isign_idx(xxf_lo,ixxf)

          !Loop over receiving "it" indices 
          do it=1,g_lo%ntheta0
             !Convert from g_lo%it to xxf_lo%it
             if(it>(xxf_lo%ntheta0+1)/2) then
                it0=it+xxf_lo%nx-xxf_lo%ntheta0
             else
                it0=it
             endif

             !Convert ixxf,it indices into iglo indices
             iglo=idx(g_lo,ik_idx(xxf_lo,ixxf),it,il_idx(xxf_lo,ixxf),&
                  ie_idx(xxf_lo,ixxf),is_idx(xxf_lo,ixxf))

             !Get proc id which has this data
             ip=proc_id(g_lo,iglo)

             !Determine message position index
             n=nn_to(ip)+1
             nn_to(ip)=n

             !Store receive indices
             to_list(ip)%first(n)=it0
             to_list(ip)%second(n)=ixxf

             !Store index for sorting
             sort_list(ip)%first(n) = ig+ntgrid-1+xxf_lo%ntgridtotal*(isign-1+2*(iglo-g_lo%llim_world)) 
          enddo
       enddo
    endif

    !Now we need to sort the to_list message indices based on sort_list | This seems potentially slow + inefficient
    do ip=0,nproc-1
       !Only need to worry about procs which we are receiving from
       if(nn_to(ip)>0) then
          !Sort using quicksort
          CALL QUICKSORT2(nn_to(ip),sort_list(ip)%first,to_list(ip)%first,to_list(ip)%second)
       endif
    enddo

    !Setup array range values
    from_low (1) = -ntgrid
    from_low (2) = 1
    from_low (3) = g_lo%llim_proc

    to_low = xxf_lo%llim_proc
    
    to_high(1) = xxf_lo%nx
    to_high(2) = xxf_lo%ulim_alloc

    from_high(1) = ntgrid
    from_high(2) = 2
    from_high(3) = g_lo%ulim_alloc

    call set_redist_character_type(g2x, 'g2x')
    call set_xxf_optimised_variables(opt_local_copy, naky, ntgrid, ntheta0, &
       nlambda, nx, xxf_lo%ulim_proc, g_lo%blocksize, layout)

    !Create g2x redistribute object
    call init_redist (g2x, 'c', to_low, to_high, to_list, &
         from_low, from_high, from_list)

    !Deallocate list objects
    call delete_list (to_list)
    call delete_list (from_list)
    call delete_list(sort_list)

  end subroutine init_x_redist_local
!</DD>

  subroutine init_y_redist (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx, ny)
    use gs2_layouts, only: init_y_transform_layouts
    use gs2_layouts, only: xxf_lo, yxf_lo, xxfidx2yxfidx, proc_id, idx_local
    use mp, only: nproc
    use redistribute, only: index_list_type, init_redist, delete_list
    use redistribute, only: set_yxf_optimised_variables, set_redist_character_type
    implicit none
    type (index_list_type), dimension(0:nproc-1) :: to_list, from_list
    integer, dimension(0:nproc-1) :: nn_to, nn_from
    integer, dimension (2) :: from_low, from_high, to_high
    integer :: to_low
    integer, intent (in) :: ntgrid, naky, ntheta0, nlambda, negrid, nspec
    integer, intent (in) :: nx, ny

    integer :: it, ixxf, ik, iyxf
    integer :: n, ip
    logical :: initialized = .false.
!<DD>This routine can probably be optimised somewhat.
!in particular the large number of calls to the index lookup routines
!and idx_local can be reduced significantly by exploiting knowledge of
!the xxf layout and similar.
!For most cases init_y_redist_local is faster

    if (initialized) return
    initialized = .true.

    call init_x_redist (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx)
    call init_y_transform_layouts &
         (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx, ny)

    ! count number of elements to be redistributed to/from each processor
    nn_to = 0
    nn_from = 0
    do ixxf = xxf_lo%llim_world, xxf_lo%ulim_world
       do it = 1, yxf_lo%nx
          call xxfidx2yxfidx (it, ixxf, xxf_lo, yxf_lo, ik, iyxf)
          if (idx_local(xxf_lo,ixxf)) &
             nn_from(proc_id(yxf_lo,iyxf)) = nn_from(proc_id(yxf_lo,iyxf)) + 1
          if (idx_local(yxf_lo,iyxf)) &
             nn_to(proc_id(xxf_lo,ixxf)) = nn_to(proc_id(xxf_lo,ixxf)) + 1
       end do
    end do

    do ip = 0, nproc-1
       if (nn_from(ip) > 0) then
          allocate (from_list(ip)%first(nn_from(ip)))
          allocate (from_list(ip)%second(nn_from(ip)))
       end if
       if (nn_to(ip) > 0) then
          allocate (to_list(ip)%first(nn_to(ip)))
          allocate (to_list(ip)%second(nn_to(ip)))
       end if
    end do

    ! get local indices of elements distributed to/from other processors
    nn_to = 0
    nn_from = 0
!
!CMR: loop over all xxf indices, find corresponding yxf indices
!
    do ixxf = xxf_lo%llim_world, xxf_lo%ulim_world
       do it = 1, yxf_lo%nx
!
!CMR obtain corresponding yxf indices
!
          call xxfidx2yxfidx (it, ixxf, xxf_lo, yxf_lo, ik, iyxf)
          if (idx_local(xxf_lo,ixxf)) then
!CMR: if xxf index local, set:
!         ip = corresponding yxf processor
!        from_list%first,second arrays = it,ixxf  (ie xxf indices)
!     later will send from_list to proc ip
             ip = proc_id(yxf_lo,iyxf)
             n = nn_from(ip) + 1
             nn_from(ip) = n
             from_list(ip)%first(n) = it
             from_list(ip)%second(n) = ixxf
          end if
          if (idx_local(yxf_lo,iyxf)) then
!CMR: if yxf index local, set ip to corresponding xxf processor
!     set to_list%first,second arrays = ik,iyxf  (ie yxf indices)
!     will receive to_list from ip
             ip = proc_id(xxf_lo,ixxf)
             n = nn_to(ip) + 1
             nn_to(ip) = n
             to_list(ip)%first(n) = ik
             to_list(ip)%second(n) = iyxf
          end if
       end do
    end do

    from_low(1) = 1
    from_low(2) = xxf_lo%llim_proc

    to_low = yxf_lo%llim_proc

    to_high(1) = yxf_lo%ny/2+1
    to_high(2) = yxf_lo%ulim_alloc

    from_high(1) = xxf_lo%nx
    from_high(2) = xxf_lo%ulim_alloc

    call set_redist_character_type(x2y, 'x2y')
    call set_yxf_optimised_variables(yxf_lo%ulim_proc)

    call init_redist (x2y, 'c', to_low, to_high, to_list, &
         from_low, from_high, from_list)

    call delete_list (to_list)
    call delete_list (from_list)
    
  end subroutine init_y_redist

!<DD>
  subroutine init_y_redist_local (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx, ny)
    use gs2_layouts, only: init_y_transform_layouts
    use gs2_layouts, only: xxf_lo, yxf_lo, xxfidx2yxfidx, proc_id, idx_local
    use gs2_layouts, only: ik_idx,it_idx,il_idx,ie_idx,is_idx,idx,ig_idx,isign_idx
    use mp, only: nproc,iproc
    use redistribute, only: index_list_type, init_redist, delete_list
    use redistribute, only: set_yxf_optimised_variables, set_redist_character_type
    implicit none
    type (index_list_type), dimension(0:nproc-1) :: to_list, from_list,sort_list,sto_list
    integer, dimension(0:nproc-1) :: nn_to, nn_from
    integer, dimension (2) :: from_low, from_high, to_high
    integer :: to_low
    integer, intent (in) :: ntgrid, naky, ntheta0, nlambda, negrid, nspec
    integer, intent (in) :: nx, ny

    integer :: it, ixxf, ik, iyxf
    integer :: ixxf_start, iyxf_start, count
    integer :: n, ip
    logical :: initialized = .false.
    logical :: debug=.true.

    !Early exit if possible
    if (initialized) return
    initialized = .true.

    !Setup g_lo-->xxf_lo redist object first
    call init_x_redist_local (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx)

    !Setup the yxf layout object
    call init_y_transform_layouts &
         (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx, ny)

    !Initialise counts to zero
    nn_to = 0
    nn_from = 0

    !First count data to send | xxf_lo-->yxf_lo
    !Protect against procs with no data
    if(xxf_lo%ulim_proc.ge.xxf_lo%llim_proc)then
       do ixxf=xxf_lo%llim_proc,xxf_lo%ulim_alloc
          !Get iyxf index for "it"=1 
          iyxf_start=idx(yxf_lo,ig_idx(xxf_lo,ixxf),&
               isign_idx(xxf_lo,ixxf),1,il_idx(xxf_lo,ixxf),&
               ie_idx(xxf_lo,ixxf),is_idx(xxf_lo,ixxf))

          !Loop over "it" range, note that we actually only want to know
          !iyxf in this range and we know that it-->it+1 => iyxf-->iyxf+1
          !so replace loop with one over iyxf
          do iyxf=iyxf_start,iyxf_start+yxf_lo%nx-1
             !Increase the appropriate procs send count
             nn_from(proc_id(yxf_lo,iyxf))=nn_from(proc_id(yxf_lo,iyxf))+1
          enddo
       enddo
    endif

    !Now count data to receive | yxf_lo<--xxf_lo
    !Protect against procs with no data
    if(yxf_lo%ulim_proc.ge.yxf_lo%llim_proc)then
       do iyxf=yxf_lo%llim_proc,yxf_lo%ulim_alloc
          !Get ixxf index for "ik"=1 
          ixxf_start=idx(xxf_lo,ig_idx(yxf_lo,iyxf),&
               isign_idx(yxf_lo,iyxf),1,il_idx(yxf_lo,iyxf),&
               ie_idx(yxf_lo,iyxf),is_idx(yxf_lo,iyxf))

          !Loop over "ik" range, note that we actually only want to know
          !ixxf in this range and we know that ik-->ik+1 => ixxf-->ixxf+1
          !so replace loop with one over ixxf
          do ixxf=ixxf_start,ixxf_start+xxf_lo%naky-1
             !Increase the appropriate procs recv count
             nn_to(proc_id(xxf_lo,ixxf))=nn_to(proc_id(xxf_lo,ixxf))+1
          enddo
       enddo
    endif

!<DD>Debug test, are we sending and receiving the same amount to ourselves? Remove when testing completed
!if(debug.and.nn_to(iproc).ne.nn_from(iproc)) print*,"ERROR: iproc ",iproc,"nn_from",nn_from(iproc),"nn_to",nn_to(iproc)

    !Now allocate storage for data mapping structures
    do ip = 0, nproc-1
       if (nn_from(ip) > 0) then
          allocate (from_list(ip)%first(nn_from(ip)))
          allocate (from_list(ip)%second(nn_from(ip)))
       end if
       if (nn_to(ip) > 0) then
          allocate (to_list(ip)%first(nn_to(ip)))
          allocate (to_list(ip)%second(nn_to(ip)))
          !For sorting to_list later
          allocate (sort_list(ip)%first(nn_to(ip)))
       end if
    end do

    !Reinitialise count arrays to zero
    nn_to = 0
    nn_from = 0

    !First fill in the sending indices, these define the messages data order
    !Protect against procs with no data
    if(xxf_lo%ulim_proc.ge.xxf_lo%llim_proc)then
       do ixxf=xxf_lo%llim_proc,xxf_lo%ulim_alloc
          !Get iyxf for "it"=1
          iyxf=idx(yxf_lo,ig_idx(xxf_lo,ixxf),&
               isign_idx(xxf_lo,ixxf),1,il_idx(xxf_lo,ixxf),&
               ie_idx(xxf_lo,ixxf),is_idx(xxf_lo,ixxf))

          !Now loop over other local dimension. Note we need "it" here
          !so don't replace this with loop over iyxf
          do it=1,yxf_lo%nx
             !Get the processor id
             ip=proc_id(yxf_lo,iyxf)

             !Increment the procs message counter
             n=nn_from(ip)+1
             nn_from(ip)=n

             !Store indices
             from_list(ip)%first(n)=it
             from_list(ip)%second(n)=ixxf

             !Increment iyxf
             iyxf=iyxf+1
          enddo
       enddo
    endif

    !Now fill in the receiving indices, these must match the message data order, achieved by later sorting
    !Protect against procs with no data
    if(yxf_lo%ulim_proc.ge.yxf_lo%llim_proc)then
       do iyxf=yxf_lo%llim_proc,yxf_lo%ulim_alloc
          !Get ixxf for "ik"=1
          ixxf=idx(xxf_lo,ig_idx(yxf_lo,iyxf),&
               isign_idx(yxf_lo,iyxf),1,il_idx(yxf_lo,iyxf),&
               ie_idx(yxf_lo,iyxf),is_idx(yxf_lo,iyxf))

          !Now loop over other local dimension. Note we need "ik" here
          !so don't replace this with loop over ixxf
          do ik=1,xxf_lo%naky
             !Get the processor id
             ip=proc_id(xxf_lo,ixxf)

             !Increment the procs message counter
             n=nn_to(ip)+1
             nn_to(ip)=n

             !Store indices
             to_list(ip)%first(n)=ik
             to_list(ip)%second(n)=iyxf

             !Store index for sorting
             sort_list(ip)%first(n)=it_idx(yxf_lo,iyxf)+ixxf*yxf_lo%nx

             !Increment ixxf
             ixxf=ixxf+1
          enddo
       enddo
    endif

    !Now we need to sort the to_list message indices based on sort_list
    !This could be slow and inefficient
    do ip=0,nproc-1
       !Only need to worry about procs which we are receiving from
       if(nn_to(ip)>0) then
          !Use quicksort based on compound index
          CALL QUICKSORT2(nn_to(ip),sort_list(ip)%first,to_list(ip)%first,to_list(ip)%second)
       endif
    enddo

    !Setup array bound values
    from_low(1) = 1
    from_low(2) = xxf_lo%llim_proc

    to_low = yxf_lo%llim_proc

    to_high(1) = yxf_lo%ny/2+1
    to_high(2) = yxf_lo%ulim_alloc

    from_high(1) = xxf_lo%nx
    from_high(2) = xxf_lo%ulim_alloc

    call set_redist_character_type(x2y, 'x2y')
    call set_yxf_optimised_variables(yxf_lo%ulim_proc)

    !Create x2y redist object
    call init_redist (x2y, 'c', to_low, to_high, to_list, &
         from_low, from_high, from_list)
    
    !Deallocate list objects
    call delete_list (to_list)
    call delete_list (from_list)
    call delete_list (sort_list)

  end subroutine init_y_redist_local
!</DD>

!<DD> Numerical recipes quicksort algorithm 
  !e.g. see Section 8.2 of Numerical recipes for fortran
  !We should perhaps move this routine to some other module as
  !it will be useful for other cases (e.g. collision layout objects)
  !Goto page 326 of numerical recipes in f90 available at http://apps.nrbook.com/fortran/index.html

  SUBROUTINE QUICKSORT2 (n,key,arr1,arr2)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: n !How long is the array?
    !Arrays to be sorted into ascending order based on key arr order
    INTEGER, INTENT(inout), DIMENSION(1:n) :: key, arr1, arr2 
    INTEGER, PARAMETER :: m=7 !How small does list segment have to be to use insertion sort
    INTEGER :: nstack !What is the best way to set this? For now just relate it to message size.
    INTEGER, DIMENSION(:),ALLOCATABLE::istack
    INTEGER :: i,ir,j,jstack,k,l
    INTEGER :: ak,a1,a2, temp

    !Initialise vars
    nstack=MAX(n/2,50)
    ALLOCATE(istack(nstack))
    jstack=0
    l=1
    ir=n

!Insertion sort when list segment becomes small enough
1   if(ir-l.lt.m) then
       do j=l+1,ir
          !Get current values
          ak=key(j)
          a1=arr1(j)
          a2=arr2(j)
          !Loop over left values
          do i=j-1,l,-1
             if(key(i).le.ak) goto 2
             key(i+1)=key(i)
             arr1(i+1)=arr1(i)
             arr2(i+1)=arr2(i)
          enddo
          i=l-1
2         key(i+1)=ak
          arr1(i+1)=a1
          arr2(i+1)=a2
       enddo
       !Exit if nothing left on stack to sort
       if(jstack.eq.0) then
          deallocate(istack)
          return
       endif

       !Move on to next section of stack, ir and l set range to look at
       ir=istack(jstack)
       l=istack(jstack-1)
       jstack=jstack-2
!Partioning and quicksort part
    else
       !Choose partion value as median of left, centre and right values
       !Which index are we looking at --> Middle of sub-array
       k=(l+ir)/2

       !Swap k and l+1 values | Key
       temp=key(k)
       key(k)=key(l+1)
       key(l+1)=temp
       !Swap k and l+1 values | Arr1
       temp=arr1(k)
       arr1(k)=arr1(l+1)
       arr1(l+1)=temp
       !Swap k and l+1 values | Arr2
       temp=arr2(k)
       arr2(k)=arr2(l+1)
       arr2(l+1)=temp

       !Put in order so key(l)<=key(l+1)<=key(ir)
       if(key(l).gt.key(ir)) then
          !Key
          temp=key(l)
          key(l)=key(ir)
          key(ir)=temp
          !Arr1
          temp=arr1(l)
          arr1(l)=arr1(ir)
          arr1(ir)=temp
          !Arr2
          temp=arr2(l)
          arr2(l)=arr2(ir)
          arr2(ir)=temp
       endif
       if(key(l+1).gt.key(ir)) then
          !Key
          temp=key(l+1)
          key(l+1)=key(ir)
          key(ir)=temp
          !Arr1
          temp=arr1(l+1)
          arr1(l+1)=arr1(ir)
          arr1(ir)=temp
          !Arr2
          temp=arr2(l+1)
          arr2(l+1)=arr2(ir)
          arr2(ir)=temp
       endif
       if(key(l).gt.key(l+1)) then
          !Key
          temp=key(l)
          key(l)=key(l+1)
          key(l+1)=temp
          !Arr1
          temp=arr1(l)
          arr1(l)=arr1(l+1)
          arr1(l+1)=temp
          !Arr2
          temp=arr2(l)
          arr2(l)=arr2(l+1)
          arr2(l+1)=temp
       endif

       !Now get ready for partitioning
       i=l+1
       j=ir

       !Pick partion value
       ak=key(l+1)
       a1=arr1(l+1)
       a2=arr2(l+1)

3      continue !Effectively a while loop until we find first element bigger than ak
       i=i+1 !Scan until we find key(i)>ak
       if(key(i).lt.ak) GOTO 3

4      continue !Effectively a while loop until we find first element smaller than ak
       j=j-1
       if(key(j).gt.ak) GOTO 4

       if(j.lt.i) GOTO 5 !If j<i then partioning has completed so move onto sorting

       !Now swap elements | Key
       temp=key(i)
       key(i)=key(j)
       key(j)=temp
       !Now swap elements | Arr1
       temp=arr1(i)
       arr1(i)=arr1(j)
       arr1(j)=temp
       !Now swap elements | Arr2
       temp=arr2(i)
       arr2(i)=arr2(j)
       arr2(j)=temp
       !Now move onto next elements
       GOTO 3
       !Put partiioning element in correct place
5      continue
       key(l+1)=key(j)
       key(j)=ak
       arr1(l+1)=arr1(j)
       arr1(j)=a1
       arr2(l+1)=arr2(j)
       arr2(j)=a2

       !Increase stack counter (move onto next section of array?)
       jstack=jstack+2
       if (jstack.gt.nstack) print*,"NSTACK too small in quicksort" !Would like to remove this
       if (ir-i+1.ge.j-l) then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
       else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
       endif
    endif

    !Repeat process for new sublist, i.e. like recursive function
    GOTO 1

  END SUBROUTINE QUICKSORT2

  SUBROUTINE QUICKSORT3 (n,key,arr1,arr2,arr3)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: n !How long is the array?
    !Arrays to be sorted into ascending order based on key arr order
    INTEGER, INTENT(inout), DIMENSION(1:n) :: key, arr1, arr2, arr3
    INTEGER, PARAMETER :: m=7 !How small does list segment have to be to use insertion sort
    INTEGER :: nstack !What is the best way to set this? For now just relate it to message size.
    INTEGER, DIMENSION(:),ALLOCATABLE::istack
    INTEGER :: i,ir,j,jstack,k,l
    INTEGER :: ak,a1,a2,a3, temp

    !Initialise vars
    nstack=MAX(n/2,50)
    ALLOCATE(istack(nstack))
    jstack=0
    l=1
    ir=n

!Insertion sort when list segment becomes small enough
1   if(ir-l.lt.m) then
       do j=l+1,ir
          !Get current values
          ak=key(j)
          a1=arr1(j)
          a2=arr2(j)
          a3=arr3(j)
          !Loop over left values
          do i=j-1,l,-1
             if(key(i).le.ak) goto 2
             key(i+1)=key(i)
             arr1(i+1)=arr1(i)
             arr2(i+1)=arr2(i)
             arr3(i+1)=arr3(i)
          enddo
          i=l-1
2         key(i+1)=ak
          arr1(i+1)=a1
          arr2(i+1)=a2
          arr3(i+1)=a3
       enddo
       !Exit if nothing left on stack to sort
       if(jstack.eq.0) then
          deallocate(istack)
          return
       endif

       !Move on to next section of stack, ir and l set range to look at
       ir=istack(jstack)
       l=istack(jstack-1)
       jstack=jstack-2
!Partioning and quicksort part
    else
       !Choose partion value as median of left, centre and right values
       !Which index are we looking at --> Middle of sub-array
       k=(l+ir)/2

       !Swap k and l+1 values | Key
       temp=key(k)
       key(k)=key(l+1)
       key(l+1)=temp
       !Swap k and l+1 values | Arr1
       temp=arr1(k)
       arr1(k)=arr1(l+1)
       arr1(l+1)=temp
       !Swap k and l+1 values | Arr2
       temp=arr2(k)
       arr2(k)=arr2(l+1)
       arr2(l+1)=temp
       !Swap k and l+1 values | Arr3
       temp=arr3(k)
       arr3(k)=arr3(l+1)
       arr3(l+1)=temp

       !Put in order so key(l)<=key(l+1)<=key(ir)
       if(key(l).gt.key(ir)) then
          !Key
          temp=key(l)
          key(l)=key(ir)
          key(ir)=temp
          !Arr1
          temp=arr1(l)
          arr1(l)=arr1(ir)
          arr1(ir)=temp
          !Arr2
          temp=arr2(l)
          arr2(l)=arr2(ir)
          arr2(ir)=temp
          !Arr3
          temp=arr3(l)
          arr3(l)=arr3(ir)
          arr3(ir)=temp
       endif
       if(key(l+1).gt.key(ir)) then
          !Key
          temp=key(l+1)
          key(l+1)=key(ir)
          key(ir)=temp
          !Arr1
          temp=arr1(l+1)
          arr1(l+1)=arr1(ir)
          arr1(ir)=temp
          !Arr2
          temp=arr2(l+1)
          arr2(l+1)=arr2(ir)
          arr2(ir)=temp
          !Arr3
          temp=arr3(l+1)
          arr3(l+1)=arr3(ir)
          arr3(ir)=temp
       endif
       if(key(l).gt.key(l+1)) then
          !Key
          temp=key(l)
          key(l)=key(l+1)
          key(l+1)=temp
          !Arr1
          temp=arr1(l)
          arr1(l)=arr1(l+1)
          arr1(l+1)=temp
          !Arr2
          temp=arr2(l)
          arr2(l)=arr2(l+1)
          arr2(l+1)=temp
          !Arr3
          temp=arr3(l)
          arr3(l)=arr3(l+1)
          arr3(l+1)=temp
       endif

       !Now get ready for partitioning
       i=l+1
       j=ir

       !Pick partion value
       ak=key(l+1)
       a1=arr1(l+1)
       a2=arr2(l+1)
       a3=arr3(l+1)

3      continue !Effectively a while loop until we find first element bigger than ak
       i=i+1 !Scan until we find key(i)>ak
       if(key(i).lt.ak) GOTO 3

4      continue !Effectively a while loop until we find first element smaller than ak
       j=j-1
       if(key(j).gt.ak) GOTO 4

       if(j.lt.i) GOTO 5 !If j<i then partioning has completed so move onto sorting

       !Now swap elements | Key
       temp=key(i)
       key(i)=key(j)
       key(j)=temp
       !Now swap elements | Arr1
       temp=arr1(i)
       arr1(i)=arr1(j)
       arr1(j)=temp
       !Now swap elements | Arr2
       temp=arr2(i)
       arr2(i)=arr2(j)
       arr2(j)=temp
       !Now swap elements | Arr3
       temp=arr3(i)
       arr3(i)=arr3(j)
       arr3(j)=temp

       !Now move onto next elements
       GOTO 3
       !Put partiioning element in correct place
5      continue
       key(l+1)=key(j)
       key(j)=ak
       arr1(l+1)=arr1(j)
       arr1(j)=a1
       arr2(l+1)=arr2(j)
       arr2(j)=a2
       arr3(l+1)=arr3(j)
       arr3(j)=a3

       !Increase stack counter (move onto next section of array?)
       jstack=jstack+2
       if (jstack.gt.nstack) print*,"NSTACK too small in quicksort" !Would like to remove this
       if (ir-i+1.ge.j-l) then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
       else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
       endif
    endif

    !Repeat process for new sublist, i.e. like recursive function
    GOTO 1

  END SUBROUTINE QUICKSORT3

!</DD>

  subroutine transform_x5d (g, xxf)
    use gs2_layouts, only: xxf_lo, g_lo
    use prof, only: prof_entering, prof_leaving
    use redistribute, only: gather
    implicit none
    complex, dimension (-xxf_lo%ntgrid:,:,g_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,xxf_lo%llim_proc:), intent (out) :: xxf
    complex, dimension(:), allocatable :: aux
    integer :: i

    call prof_entering ("transform_x5d", "gs2_transforms")

    ! intent statement in gather actually makes this next line non-standard: 
    xxf = 0.
!CMR, 7/3/2011: gather pulls appropriate pieces of g onto this processor for
!    local Fourier transform in x, and may also pad with zeros for dealiasing
!
    call gather (g2x, g, xxf)

    ! do ffts
# if FFT == _FFTW_
    !JH 7th December 2011
    !JH xxf_lo%ulim_alloc is used here rather than xxf_lo%lulim_proc
    !JH because there are situations where xxf_lo%llim_proc is greater 
    !JH than xxf_lo%ulim_proc and that would create a negative number 
    !JH of FFTs to be calculated.  However, xxf_lo%ulim_alloc is set
    !JH to be xxf_lo%llim_proc in this situation, and that will give 
    !JH 1 FFT to be calculated which the code can correctly undertake.
    i = xxf_lo%ulim_alloc - xxf_lo%llim_proc + 1

    allocate (aux(xf_fft%n))
    call fftw_f77 (xf_fft%plan, i, xxf, 1, xxf_lo%nx, aux, 0, 0)
    deallocate (aux)
# elif FFT == _FFTW3_
    call dfftw_execute(xf_fft%plan)
# endif

    call prof_leaving ("transform_x5d", "gs2_transforms")
  end subroutine transform_x5d

  subroutine inverse_x5d (xxf, g)
    use gs2_layouts, only: xxf_lo, g_lo
    use prof, only: prof_entering, prof_leaving
    use redistribute, only: scatter
    implicit none
    complex, dimension (:,xxf_lo%llim_proc:), intent (in out) :: xxf
    complex, dimension (-xxf_lo%ntgrid:,:,g_lo%llim_proc:), intent (out) :: g
    complex, dimension(:), allocatable :: aux
    integer :: i

    call prof_entering ("inverse_x5d", "gs2_transforms")

    ! do ffts
# if FFT == _FFTW_
    !JH 7th December 2011
    !JH xxf_lo%ulim_alloc is used here rather than xxf_lo%lulim_proc
    !JH because there are situations where xxf_lo%llim_proc is greater 
    !JH than xxf_lo%ulim_proc and that would create a negative number 
    !JH of FFTs to be calculated.  However, xxf_lo%ulim_alloc is set
    !JH to be xxf_lo%llim_proc in this situation, and that will give 
    !JH 1 FFT to be calculated which the code can correctly undertake.
    i = xxf_lo%ulim_alloc - xxf_lo%llim_proc + 1

    allocate (aux(xb_fft%n))
    call fftw_f77 (xb_fft%plan, i, xxf, 1, xxf_lo%nx, aux, 0, 0)
    deallocate (aux)
# elif FFT == _FFTW3_
    call dfftw_execute(xb_fft%plan)
# endif

    call scatter (g2x, xxf, g)

    call prof_leaving ("inverse_x5d", "gs2_transforms")
  end subroutine inverse_x5d

  subroutine transform_y5d (xxf, yxf)
    use gs2_layouts, only: xxf_lo, yxf_lo
    use prof, only: prof_entering, prof_leaving
    use redistribute, only: gather
    implicit none
    complex, dimension (:,xxf_lo%llim_proc:), intent (in) :: xxf
# ifdef FFT
    real, dimension (:,yxf_lo%llim_proc:), intent (out) :: yxf
# else
    real, dimension (:,yxf_lo%llim_proc:) :: yxf
# endif
    integer :: i

    call prof_entering ("transform_y5d", "gs2_transforms")

    fft = 0.
    call gather (x2y, xxf, fft)

    ! do ffts
# if FFT == _FFTW_
    !JH 7th December 2011
    !JH yxf_lo%ulim_alloc is used here rather than yxf_lo%lulim_proc
    !JH because there are situations where yxf_lo%llim_proc is greater 
    !JH than yxf_lo%ulim_proc and that would create a negative number 
    !JH of FFTs to be calculated.  However, yxf_lo%ulim_alloc is set
    !JH to be yxf_lo%llim_proc in this situation, and that will give 
    !JH 1 FFT to be calculated which the code can correctly undertake.
    i = yxf_lo%ulim_alloc - yxf_lo%llim_proc + 1
    call rfftwnd_f77_complex_to_real (yf_fft%plan, i, fft, 1, yxf_lo%ny/2+1, yxf, 1, yxf_lo%ny)
# elif FFT == _FFTW3_

    call dfftw_execute_dft_c2r (yf_fft%plan, fft, yxf)
# endif

    call prof_leaving ("transform_y5d", "gs2_transforms")
  end subroutine transform_y5d

  subroutine inverse_y5d (yxf, xxf)
    use gs2_layouts, only: xxf_lo, yxf_lo
    use prof, only: prof_entering, prof_leaving
    use redistribute, only: scatter
    implicit none
    real, dimension (:,yxf_lo%llim_proc:), intent (in out) :: yxf
    complex, dimension (:,xxf_lo%llim_proc:), intent (out) :: xxf
    integer :: i 

    call prof_entering ("inverse_y5d", "gs2_transforms")

    ! do ffts
# if FFT == _FFTW_
    !JH 7th December 2011
    !JH yxf_lo%ulim_alloc is used here rather than yxf_lo%lulim_proc
    !JH because there are situations where yxf_lo%llim_proc is greater 
    !JH than yxf_lo%ulim_proc and that would create a negative number 
    !JH of FFTs to be calculated.  However, yxf_lo%ulim_alloc is set
    !JH to be yxf_lo%llim_proc in this situation, and that will give 
    !JH 1 FFT to be calculated which the code can correctly undertake.
    i = yxf_lo%ulim_alloc - yxf_lo%llim_proc + 1
    call rfftwnd_f77_real_to_complex (yb_fft%plan, i, yxf, 1, yxf_lo%ny, fft, 1, yxf_lo%ny/2+1)
# elif FFT == _FFTW3_
    call dfftw_execute_dft_r2c (yb_fft%plan, yxf, fft)
# endif

    call scatter (x2y, fft, xxf)

    call prof_leaving ("inverse_y5d", "gs2_transforms")
  end subroutine inverse_y5d

  subroutine transform2_5d (g, yxf)
    use gs2_layouts, only: g_lo, yxf_lo, ik_idx
    implicit none
    complex, dimension (:,:,g_lo%llim_proc:), intent (in out) :: g
    real, dimension (:,yxf_lo%llim_proc:), intent (out) :: yxf
    integer :: iglo

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       if (ik_idx(g_lo, iglo) == 1) cycle
       g(:,:,iglo) = g(:,:,iglo) / 2.0
    end do

    call transform_x (g, xxf)
    call transform_y (xxf, yxf)
  end subroutine transform2_5d

  subroutine inverse2_5d (yxf, g)
    use gs2_layouts, only: g_lo, yxf_lo, ik_idx
    implicit none
    real, dimension (:,yxf_lo%llim_proc:), intent (in out) :: yxf
    complex, dimension (:,:,g_lo%llim_proc:), intent (out) :: g
    integer :: iglo

    call inverse_y (yxf, xxf)
    call inverse_x (xxf, g)

    g = g * xb_fft%scale * yb_fft%scale

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       if (ik_idx(g_lo, iglo) == 1) cycle
       g(:,:,iglo) = g(:,:,iglo) * 2.0
    end do

  end subroutine inverse2_5d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine transform2_5d_accel (g, axf, i)
    use gs2_layouts, only: g_lo, accel_lo, accelx_lo, ik_idx
    implicit none
    complex, dimension (:,:,g_lo%llim_proc:), intent (in out) :: g

# ifdef FFT

    real, dimension (:,:,accelx_lo%llim_proc:), intent (out) :: axf

# else

    real, dimension (:,:,accelx_lo%llim_proc:) :: axf

# endif


    integer :: iglo, k, i, idx
    integer :: itgrid, iduo
    integer :: ntgrid

    ntgrid = accel_lo%ntgrid

    ! scale the g and copy into the anti-aliased array ag
    ! zero out empty ag
    ! touch each g and ag only once
    idx = g_lo%llim_proc
    do k = accel_lo%llim_proc, accel_lo%ulim_proc
! CMR: aidx is true for modes killed by the dealiasing mask
! so following line removes ks in dealiasing region
       if (aidx(k)) then
          ag(:,:,k) = 0.0
       else
          ! scaling only for k_y not the zero mode
          if (ik_idx(g_lo, idx) .ne. 1) then
             do iduo = 1, 2
                do itgrid = 1, 2*ntgrid +1
                   g(itgrid, iduo, idx) &
                        = 0.5 * g(itgrid, iduo, idx)
                   ag(itgrid - (ntgrid+1), iduo, k) &
                        = g(itgrid, iduo, idx)
                enddo
             enddo
             ! in case of k_y being zero: just copy
          else
             do iduo = 1, 2
                do itgrid = 1, 2*ntgrid+1
                   ag(itgrid-(ntgrid+1), iduo, k) &
                        = g(itgrid, iduo, idx)
                enddo
             enddo
          endif
          idx = idx + 1
       endif
    enddo

    ! we might not have scaled all g
    Do iglo = idx, g_lo%ulim_proc
       if (ik_idx(g_lo, iglo) .ne. 1) then
          g(:,:, iglo) = 0.5 * g(:,:,iglo)
       endif
    enddo
       
       
    ! transform
    i = (2*accel_lo%ntgrid+1)*2
    idx = 1
    do k = accel_lo%llim_proc, accel_lo%ulim_proc, accel_lo%nxnky
# if FFT == _FFTW_
       call rfftwnd_f77_complex_to_real (yf_fft%plan, i, ag(:,:,k:), i, 1, &
            axf(:,:,ia(idx):), i, 1)
# elif FFT == _FFTW3_
       ! remember FFTW3 for c2r destroys the contents of ag
       call dfftw_execute_dft_c2r (yf_fft%plan, ag(:, :, k:), &
            axf(:, :, ia(idx):))
# endif
       idx = idx + 1
    end do

  end subroutine transform2_5d_accel

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inverse2_5d_accel (axf, g, i)
    use gs2_layouts, only: g_lo, accel_lo, accelx_lo, ik_idx
    implicit none
    real, dimension (:,:,accelx_lo%llim_proc:), intent (in out) :: axf
    complex, dimension (:,:,g_lo%llim_proc:), intent (out) :: g
    integer :: iglo, i, idx, k
    integer :: itgrid, iduo
    integer :: ntgrid

    ntgrid = accel_lo%ntgrid

! transform
    i = (2*accel_lo%ntgrid+1)*2
    idx = 1
    do k = accelx_lo%llim_proc, accelx_lo%ulim_proc, accelx_lo%nxny
# if FFT == _FFTW_
       call rfftwnd_f77_real_to_complex (yb_fft%plan, i, axf(:,:,k:), i, 1, &
            ag(:,:,iak(idx):), i, 1)
# elif FFT == _FFTW3_
       call dfftw_execute_dft_r2c(yb_fft%plan, axf(:, :, k:), &
            ag(:, :, iak(idx):))
# endif
       idx = idx + 1
    end do

    idx = g_lo%llim_proc
    do k = accel_lo%llim_proc, accel_lo%ulim_proc
       ! ignore the large k (anti-alias)
       if ( .not.aidx(k)) then
          ! different scale factors depending on ky == 0
          if (ik_idx(g_lo, idx) .ne. 1) then
             do iduo = 1, 2 
                do itgrid = 1, 2*ntgrid+1
                   g(itgrid, iduo, idx) &
                        = 2.0 * yb_fft%scale * ag(itgrid-(ntgrid+1), iduo, k)
                enddo
             enddo
          else
             do iduo = 1, 2 
                do itgrid = 1, 2*ntgrid+1
                   g(itgrid, iduo, idx) &
                        = yb_fft%scale * ag(itgrid-(ntgrid+1), iduo, k)
                enddo
             enddo
          endif
          idx = idx + 1
       endif
    enddo
    
    ! we might not have scaled all g
    Do iglo = idx, g_lo%ulim_proc
       if (ik_idx(g_lo, iglo) .ne. 1) then
          g(:,:, iglo) = 2.0 * g(:,:,iglo)
       endif
    enddo

  end subroutine inverse2_5d_accel

  subroutine init_3d (nny_in, nnx_in, how_many_in)

    use fft_work, only: init_crfftw, init_rcfftw, delete_fft
    logical :: initialized = .false.
    integer :: nny_in, nnx_in, how_many_in
    integer, save :: nnx, nny, how_many

    if (initialized) then
       if (nnx /= nnx_in .or. nny /= nny_in) then
          call delete_fft(xf3d_cr)
          call delete_fft(xf3d_rc)
# if FFT == _FFTW3_
       elseif ( how_many /= how_many_in) then
          call delete_fft(xf3d_cr)
          call delete_fft(xf3d_rc)
# endif
       else
          return
       end if
    end if
    initialized = .true.
    nny = nny_in
    nnx = nnx_in
    how_many = how_many_in

#if FFT == _FFTW_    

    call init_crfftw (xf3d_cr,  1, nny, nnx)
    call init_rcfftw (xf3d_rc, -1, nny, nnx)

#elif FFT == _FFTW3_

    call init_crfftw (xf3d_cr,  1, nny, nnx, how_many)
    call init_rcfftw (xf3d_rc, -1, nny, nnx, how_many)

#endif

  end subroutine init_3d

  subroutine transform2_3d (phi, phixf, nny, nnx)
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0, nx, aky
    implicit none
    integer :: nnx, nny
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi
    real, dimension (:,:,-ntgrid:), intent (out) :: phixf  
    real, dimension (:,:,:), allocatable :: phix
    complex, dimension (:,:,:), allocatable :: aphi
    real :: fac
    integer :: ig, ik, it, i

! scale, dealias and transpose

    call init_3d (nny, nnx, 2*ntgrid+1)

    allocate (phix (-ntgrid:ntgrid, nny, nnx))
    allocate (aphi (-ntgrid:ntgrid, nny/2+1, nnx))
    aphi = 0.

    do ik=1,naky
       fac = 0.5
       if (aky(ik) < epsilon(0.)) fac = 1.0
       do it=1,(ntheta0+1)/2
          do ig=-ntgrid, ntgrid
             aphi(ig,ik,it) = phi(ig,it,ik)*fac
          end do
       end do
       do it=(ntheta0+1)/2+1,ntheta0
          do ig=-ntgrid, ntgrid
!CMR, 30/3/2010: bug fix to replace nx by nnx on next line
             aphi(ig,ik,it-ntheta0+nnx) = phi(ig,it,ik)*fac
          end do
       end do
    end do

! transform
    i = 2*ntgrid+1
# if FFT == _FFTW_
    call rfftwnd_f77_complex_to_real (xf3d_cr%plan, i, aphi, i, 1, phix, i, 1)
# elif FFT == _FFTW3_
    call dfftw_execute_dft_c2r (xf3d_cr%plan, aphi, phix)
# endif

    do it=1,nnx
       do ik=1,nny
          do ig=-ntgrid, ntgrid
             phixf (it,ik,ig) = phix (ig,ik,it)
          end do
       end do
    end do

    deallocate (aphi, phix)

  end subroutine transform2_3d
  
  subroutine inverse2_3d (phixf, phi, nny, nnx)
!CMR, 30/4/2010:  
! Fixed up previously buggy handling of dealiasing.
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0, aky
    implicit none
    real, dimension (:,:,-ntgrid:):: phixf
    complex, dimension (-ntgrid:,:,:) :: phi
    integer :: nnx, nny
    complex, dimension (:,:,:), allocatable :: aphi
    real, dimension (:,:,:), allocatable :: phix
    real :: fac
    integer :: i, ik, it, ig

    allocate (aphi (-ntgrid:ntgrid, nny/2+1, nnx))
    allocate (phix (-ntgrid:ntgrid, nny, nnx))

    do it=1,nnx
       do ik=1,nny
          do ig=-ntgrid, ntgrid
             phix (ig,ik,it) = phixf (it,ik,ig)
          end do
       end do
    end do

! transform
    i = 2*ntgrid+1
# if FFT == _FFTW_
    call rfftwnd_f77_real_to_complex (xf3d_rc%plan, i, phix, i, 1, aphi, i, 1)
# elif FFT == _FFTW3_
    call dfftw_execute_dft_r2c (xf3d_rc%plan, phix, aphi)
# endif

! dealias and scale
    do it=1,(ntheta0+1)/2
       do ik=1,naky
          fac = 2.0
          if (aky(ik) < epsilon(0.0)) fac = 1.0
          do ig=-ntgrid, ntgrid
             phi (ig,it,ik) = aphi (ig,ik,it)*fac*xf3d_rc%scale
          end do
       end do
    end do

!CMR, 30/4/2010:  fixed up previously buggy handling of dealiasing
    do it=(ntheta0+1)/2+1,ntheta0
       do ik=1,naky
          fac = 2.0
          if (aky(ik) < epsilon(0.0)) fac = 1.0
          do ig=-ntgrid, ntgrid
             phi (ig,it,ik) = aphi (ig,ik,it-ntheta0+nnx)*fac*xf3d_rc%scale
          end do
       end do
    end do

    deallocate (aphi, phix)

  end subroutine inverse2_3d

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine transform2_2d (phi, phixf, nny, nnx)
    use fft_work, only: FFTW_BACKWARD, delete_fft, init_crfftw
    use kt_grids, only: naky, nakx => ntheta0, nx, aky, akx
    implicit none
    integer :: nnx, nny
    complex, intent (in) :: phi(:,:)
    real, intent (out) :: phixf  (:,:)
    real, allocatable :: phix(:,:)
    complex, allocatable :: aphi(:,:)
    real :: fac
    integer :: ik, it
    type (fft_type) :: xf2d

#if FFT == _FFTW_

    call init_crfftw (xf2d, FFTW_BACKWARD, nny, nnx)

#elif FFT == _FFTW3_

    call init_crfftw (xf2d, FFTW_BACKWARD, nny, nnx, 1)

#endif

    allocate (phix (nny, nnx))
    allocate (aphi (nny/2+1, nnx))
    phix(:,:)=0.; aphi(:,:)=cmplx(0.,0.)

! scale, dealias and transpose
    do ik=1,naky
       fac = 0.5
       if (aky(ik) < epsilon(0.)) fac = 1.0
       do it=1,(nakx+1)/2
          aphi(ik,it) = phi(it,ik)*fac
       end do
       do it=(nakx+1)/2+1,nakx
          aphi(ik,it-nakx+nx) = phi(it,ik)*fac
       end do
    end do

! transform
# if FFT == _FFTW_
    call rfftwnd_f77_complex_to_real (xf2d%plan, 1, aphi, 1, 1, phix, 1, 1)
# elif FFT == _FFTW3_
    call dfftw_execute_dft_c2r (xf2d%plan, aphi, phix)
# endif

    phixf(:,:)=transpose(phix(:,:))

    deallocate (aphi, phix)
!RN> this statement causes error for lahey with DEBUG. I don't know why
!    call delete_fft(xf2d)
  end subroutine transform2_2d


  subroutine inverse2_2d (phixf, phi, nny, nnx)
    use fft_work, only: FFTW_FORWARD, delete_fft, init_rcfftw
    use kt_grids, only: naky, nakx => ntheta0, aky
    implicit none
    real, intent(in) :: phixf(:,:)
    complex, intent(out) :: phi(:,:)
    integer :: nnx, nny
    complex, allocatable :: aphi(:,:)
    real, allocatable :: phix(:,:)
    real :: fac
    integer :: ik, it
    type (fft_type) :: xf2d

#if FFT == _FFTW_

    call init_rcfftw (xf2d, FFTW_FORWARD, nny, nnx)

#elif FFT == _FFTW3_

    call init_rcfftw (xf2d, FFTW_FORWARD, nny, nnx, 1)

#endif    

    allocate (aphi (nny/2+1, nnx))
    allocate (phix (nny, nnx))
    phix(:,:)=cmplx(0.,0.); aphi(:,:)=cmplx(0.,0.)

    phix(:,:)=transpose(phixf(:,:))

! transform
# if FFT == _FFTW_
    call rfftwnd_f77_real_to_complex (xf2d%plan, 1, phix, 1, 1, aphi, 1, 1)
# elif FFT == _FFTW3_
    call dfftw_execute_dft_r2c (xf2d%plan, phix, aphi)
# endif

! scale, dealias and transpose
    do it=1,nakx
       do ik=1,naky
          fac = 2.0
          if (aky(ik) < epsilon(0.0)) fac = 1.0
          phi(it,ik) = aphi(ik,it)*fac*xf2d%scale
       end do
    end do

    deallocate (aphi, phix)
!RN> this statement causes error for lahey with DEBUG. I don't know why
!    call delete_fft(xf2d)
  end subroutine inverse2_2d

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine transform2_4d (den, phixf, nny, nnx)
!CMR, 30/3/2010: den input has 4th index being species, 
!     but transform only operates for species=1
!     anyone who uses this routine should be aware of/fix this!
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0, nx, aky
    implicit none
    integer :: nnx, nny
    complex, dimension (-ntgrid:,:,:,:), intent (in) :: den
    real, dimension (:,:,-ntgrid:), intent (out) :: phixf  
    real, dimension (:,:,:), allocatable :: phix
    complex, dimension (:,:,:), allocatable :: aphi
    real :: fac
    integer :: ig, ik, it, i

! scale, dealias and transpose

    call init_3d (nny, nnx, 2*ntgrid+1)

    allocate (phix (-ntgrid:ntgrid, nny, nnx))
    allocate (aphi (-ntgrid:ntgrid, nny/2+1, nnx))
    aphi = 0.

    do ik=1,naky
       fac = 0.5
       if (aky(ik) < epsilon(0.)) fac = 1.0
       do it=1,(ntheta0+1)/2
          do ig=-ntgrid, ntgrid
             aphi(ig,ik,it) = den(ig,it,ik,1)*fac
          end do
       end do
       do it=(ntheta0+1)/2+1,ntheta0
          do ig=-ntgrid, ntgrid
!CMR, 30/3/2010: bug fix to replace nx by nnx on next line
             aphi(ig,ik,it-ntheta0+nnx) = den(ig,it,ik,1)*fac
          end do
       end do
    end do

! transform
    i = 2*ntgrid+1
# if FFT == _FFTW_
    call rfftwnd_f77_complex_to_real (xf3d_cr%plan, i, aphi, i, 1, phix, i, 1)
# elif FFT == _FFTW3_
    call dfftw_execute_dft_c2r (xf3d_cr%plan, aphi, phix)
# endif

    do it=1,nnx
       do ik=1,nny
          do ig=-ntgrid, ntgrid
             phixf (it,ik,ig) = phix (ig,ik,it)
          end do
       end do
    end do

    deallocate (aphi, phix)

  end subroutine transform2_4d
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_zf (ntgrid, nperiod, howmany)

    use fft_work, only: init_z
    implicit none
    integer, intent (in) :: ntgrid, nperiod, howmany
    logical :: done = .false.

    if (done) return
    done = .true.

    call init_z (zf_fft, 1, 2*ntgrid, howmany)
    
  end subroutine init_zf

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
  subroutine kz_spectrum (an, an2, ntgrid, ntheta0, naky)

    complex, dimension (:,:,:), intent(in)  :: an
    complex, dimension (:,:,:), intent(out) :: an2
    integer, intent (in) :: ntheta0, naky, ntgrid

# if FFT == _FFTW_    
    call fftw_f77 (zf_fft%plan, ntheta0*naky, an, 1, zf_fft%n+1, an2, 1, zf_fft%n+1)
# elif FFT == _FFTW3_
    call dfftw_execute_dft(zf_fft%plan, an, an2)
# endif
    an2 = conjg(an2)*an2

  end subroutine kz_spectrum

end module gs2_transforms