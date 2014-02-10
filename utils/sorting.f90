!>This module provides various routines for sorting arrays
!This isn't something we commonly do but these routines are
!useful in a few places.
!
!As with most algorithms different options can be made depending
!on if we want to optimise memory usage, speed etc.
module sorting
  implicit none

  public :: quicksort
  private 

!Interfaces
  interface quicksort
     ! module procedure i_1_quicksort !Not yet implemented
     module procedure i_2_quicksort
     module procedure i_3_quicksort
  end interface

contains
!////////////////////
!// QUICKSORT
!////////////////////

  !<DD> Numerical recipes quicksort algorithm 
  !e.g. see Section 8.2 of Numerical recipes for fortran
  !Goto page 326 of numerical recipes in f90 available at http://apps.nrbook.com/fortran/index.html

  !Apologies for rather nasty style of these routines, feel free
  !to update to a more modern style but please check correctness and
  !performance!

  !NOTE: Currently we provide several routines for
  !sorting different numbers of arrays based on the same key.
  !It is likely to be faster to sort an "index" array (e.g.
  !an array of integers initially (/1..n/)) and use this to
  !address as many arrays (of any type) that we want to sort.
  !This will however have some additional memory overhead so
  !for now do things this way.

  !Sort two integer arrays based on integer KEY
  SUBROUTINE i_2_quicksort (n,key,arr1,arr2)
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
       !<DD>Should probably uncomment below but change it to mp_abort
!       if (jstack.gt.nstack) print*,"NSTACK too small in quicksort" !Would like to remove this
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

  END SUBROUTINE i_2_quicksort
  
  !Sort three integer arrays based on integer KEY  
  SUBROUTINE i_3_quicksort (n,key,arr1,arr2,arr3)
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
       !<DD>Should probably uncomment below but change it to mp_abort
!       if (jstack.gt.nstack) print*,"NSTACK too small in quicksort" !Would like to remove this
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

  END SUBROUTINE i_3_quicksort
  !</DD>
end module sorting
