      integer function lnblank(string)
c
c Finds the "Last non-blank character" before the first blank in a string
c
      character string*(*)

      lnblank=index(string,' ')-1
      if(lnblank .eq. -1) lnblank=len(string)

      return 
      end
