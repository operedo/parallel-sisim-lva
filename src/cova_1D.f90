      subroutine cova3_1D_sp(dist,ivarg,nst,MAXNST,MAXCUT,c0,it,cc,aa,cmax,cova)
!-----------------------------------------------------------------------
!
!                    Covariance 
!                    ***********
! determine the cov given a distance
!
! INPUT VARIABLES:
!
!   dist           1D distance between points
!   nst(ivarg)       number of nested structures (maximum of 4)
!   ivarg            variogram number (set to 1 unless doing cokriging
!                       or indicator kriging)
!   MAXNST           size of variogram parameter arrays
!   c0(ivarg)        isotropic nugget constant
!   it(i)            type of each nested structure:
!                      1. spherical model of range a;
!                      2. exponential model of parameter a;
!                           i.e. practical range is 3a
!                      3. gaussian model of parameter a;
!                           i.e. practical range is a*sqrt(3)
!                      4. power model of power a (a must be gt. 0  and
!                           lt. 2).  if linear model, a=1,c=slope.
!                      5. hole effect model
!   cc(i)            multiplicative factor of each nested structure.
!                      (sill-c0) for spherical, exponential,and gaussian
!                      slope for linear model.
!   aa(i)            parameter "a" of each nested structure.
!   irot             index of the rotation matrix for the first nested 
!                    structure (the second nested structure will use
!                    irot+1, the third irot+2, and so on)
!   MAXROT           size of rotation matrix arrays
!   rotmat           rotation matrices
!
!
! OUTPUT VARIABLES:
!
!   cmax             maximum covariance
!   cova             covariance between (x1,y1,z1) and (x2,y2,z2)
!
!
!
! EXTERNAL REFERENCES: sqdist    computes anisotropic squared distance
!                      rotmat    computes rotation matrix for distance
!-----------------------------------------------------------------------

!      USE DFLIB 
      implicit none

      real,intent(in) :: dist
      real,intent(inout) :: cmax
      integer,intent(in) :: ivarg,MAXNST,MAXCUT
      integer,intent(in) :: nst(MAXCUT),it(MAXCUT*MAXNST)
      real,intent(in) :: c0(MAXCUT),cc(MAXCUT*MAXNST),aa(MAXCUT*MAXNST)
      real,intent(inout) :: cova

      real :: PI,PMX,EPSLON,h,hr
      integer :: istart,is,ist
      !integer   nst(*),it(*)
      !real      c0(*),cc(*),aa(*)
      real ::   sqdist
      real :: hsqd
      
      PI=3.14159
      PMX=999.0
      EPSLON=1.e-5
      !write(*,*)'inside-cova3',dist

      
      
!
! Calculate the maximum covariance value (used for zero distances and
! for power model covariance):
!
      istart = 1 + (ivarg-1)*MAXNST
      cmax   = c0(ivarg)
      do is=1,nst(ivarg)
            ist = istart + is - 1
            if(it(ist).eq.4) then
                  cmax = cmax + 999.0
            else
                  cmax = cmax + cc(ist)
            endif
      end do


!
! Check for "zero" distance, return with cmax if so:

      if(real(dist).lt.EPSLON) then !need cmax
            cova = cmax
            return
      endif
      cova = 0.0
      h = dist  
!
! Loop over all the structures:
!
      
      do is=1,nst(ivarg)
            ist = istart + is - 1
 
!
! Spherical Variogram Model?
!
            if(it(ist).eq.1) then
                  hr = h/aa(ist)
                  if(hr.lt.1.) cova=cova+cc(ist)*(1.-hr*(1.5-.5*hr*hr))
!
! Exponential Variogram Model?
!
            else if(it(ist).eq.2) then
                  cova = cova + cc(ist)*exp(-3.0*h/aa(ist))
!                  write(*,*)'cov-exp:',cova,h,cc(ist),aa(ist)
!
! Gaussian Variogram Model?
!
            else if(it(ist).eq.3) then
                  cova = cova + cc(ist)*exp(-3.*(h/aa(ist))*(h/aa(ist)))
!
! Power Variogram Model?
!
            else if(it(ist).eq.4) then
          
              
              cova = cova + cmax - cc(ist)*(h**aa(ist))
!
! Hole Effect Model?
!
            else if(it(ist).eq.5) then
!                 d = 10.0 * aa(ist)
!                 cova = cova + cc(ist)*exp(-3.0*h/d)*cos(h/aa(ist)*PI)
                  cova = cova + cc(ist)*cos(h/aa(ist)*PI)
            else if(it(ist).eq.6) then
                  cova = cova + cc(ist)-h/aa(ist)*cc(ist)
            endif
!         write(*,*)'ent-cova3-0',it(ist),dist,cova
      end do
!
! Finished:
!
      return
      end


      subroutine cova3_1D(dist,ivarg,nst,MAXNST,c0,it,cc,aa,cmax,cova)
!-----------------------------------------------------------------------
!
!                    Covariance 
!                    ***********
! determine the cov given a distance
!
! INPUT VARIABLES:
!
!   dist           1D distance between points
!   nst(ivarg)       number of nested structures (maximum of 4)
!   ivarg            variogram number (set to 1 unless doing cokriging
!                       or indicator kriging)
!   MAXNST           size of variogram parameter arrays
!   c0(ivarg)        isotropic nugget constant
!   it(i)            type of each nested structure:
!                      1. spherical model of range a;
!                      2. exponential model of parameter a;
!                           i.e. practical range is 3a
!                      3. gaussian model of parameter a;
!                           i.e. practical range is a*sqrt(3)
!                      4. power model of power a (a must be gt. 0  and
!                           lt. 2).  if linear model, a=1,c=slope.
!                      5. hole effect model
!   cc(i)            multiplicative factor of each nested structure.
!                      (sill-c0) for spherical, exponential,and gaussian
!                      slope for linear model.
!   aa(i)            parameter "a" of each nested structure.
!   irot             index of the rotation matrix for the first nested 
!                    structure (the second nested structure will use
!                    irot+1, the third irot+2, and so on)
!   MAXROT           size of rotation matrix arrays
!   rotmat           rotation matrices
!
!
! OUTPUT VARIABLES:
!
!   cmax             maximum covariance
!   cova             covariance between (x1,y1,z1) and (x2,y2,z2)
!
!
!
! EXTERNAL REFERENCES: sqdist    computes anisotropic squared distance
!                      rotmat    computes rotation matrix for distance
!-----------------------------------------------------------------------

!      USE DFLIB 
      implicit none

      real*8 PI,PMX,EPSLON,dist,cova,cmax,h,hr
      integer istart,ivarg,maxnst,is,ist
      parameter(PI=3.14159265,PMX=999.,EPSLON=1.e-5)
      integer   nst(*),it(*)
      real*8      c0(*),cc(*),aa(*)
      real*8    sqdist
      real*8 :: hsqd
      

      
      
!
! Calculate the maximum covariance value (used for zero distances and
! for power model covariance):
!
      istart = 1 + (ivarg-1)*MAXNST
!
! Check for "zero" distance, return with cmax if so:

      if(real(dist).lt.EPSLON) then !need cmax
            cova = cmax
            return
      endif
      cova = 0.0
      h = dist  
!
! Loop over all the structures:
!
      
      do is=1,nst(ivarg)
            ist = istart + is - 1
 
!
! Spherical Variogram Model?
!
            if(it(ist).eq.1) then
                  hr = h/aa(ist)
                  if(hr.lt.1.) cova=cova+cc(ist)*(1.-hr*(1.5-.5*hr*hr))
!
! Exponential Variogram Model?
!
            else if(it(ist).eq.2) then
                  cova = cova + cc(ist)*exp(-3.0*h/aa(ist))
!
! Gaussian Variogram Model?
!
            else if(it(ist).eq.3) then
                  cova = cova + cc(ist)*exp(-3.*(h/aa(ist))*(h/aa(ist)))
!
! Power Variogram Model?
!
            else if(it(ist).eq.4) then
          
              
              cova = cova + cmax - cc(ist)*(h**aa(ist))
!
! Hole Effect Model?
!
            else if(it(ist).eq.5) then
!                 d = 10.0 * aa(ist)
!                 cova = cova + cc(ist)*exp(-3.0*h/d)*cos(h/aa(ist)*PI)
                  cova = cova + cc(ist)*cos(h/aa(ist)*PI)
            else if(it(ist).eq.6) then
                  cova = cova + cc(ist)-h/aa(ist)*cc(ist)
            endif
      end do
!
! Finished:
!
      return
      end
