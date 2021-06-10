      subroutine getindx(gn,gmin,gsiz,gloc,gindex,ginflag)
c-----------------------------------------------------------------------
c
c     Gets the coordinate index location of a point within a grid
c     ***********************************************************
c
c
c n       number of "nodes" or "cells" in this coordinate direction
c min     origin at the center of the first cell
c siz     size of the cells
c loc     location of the point being considered
c index   output index within [1,n]
c inflag  true if the location is actually in the grid (false otherwise
c         e.g., if the location is outside then index will be set to
c         nearest boundary
c
c
c
c-----------------------------------------------------------------------
      integer   gn,gindex
      real*8    gmin,gsiz,gloc
      logical   ginflag
c
c Compute the index of "loc":
c
c      write(*,*)'daty1:',gn,gmin,gsiz,gloc
      gindex = int( (gloc-gmin)/gsiz + 1.5 )
c      write(*,*)'daty2:',gindex
c
c Check to see if in or out:
c
      if(gindex.lt.1) then
            gindex  = 1
            ginflag = .false.
      else if(gindex.gt.gn) then
            gindex  = gn
            ginflag = .false.
      else
            ginflag = .true.
      end if
c
c Return to calling program:
c
      return
      end
