!
! Module to declare dynamic arrays in multiple subroutines:
!
      module      geostat        
! Editado por RGE, (global tiene landpts)
      use  global
      use  graph_vars
      use  graph
      
      real,allocatable   :: x(:),y(:),z(:),vr(:,:),close(:),actloc(:),& 
                            tmpdat(:),order(:),gcut(:),sim(:),tmp(:),& 
                            gcdf(:),covtab(:,:,:,:),cnodex(:),cnodey(:),&
                            cnodez(:),cnodev(:),cnodet(:),vra(:),&
                            thres(:),cdf(:),ccdf(:),ccdfo(:),beez(:),&
                            c0(:),cmax(:),cc(:),aa(:),ang1(:),ang2(:),&
                            ang3(:),anis1(:),anis2(:),aviol(:),xviol(:),&
                            sec1(:),sec2(:),sec3(:)
      real*8,allocatable  :: r(:),rr(:),s(:),a(:),distance(:)
      integer,allocatable :: ixnode(:),iynode(:),iznode(:),nisb(:),&
                             icnode(:),ixsbtosr(:),iysbtosr(:),&
                             izsbtosr(:),it(:),nst(:),nviol(:),aclose(:)
      real*8,allocatable  :: rotmat(:,:,:)
      logical,allocatable :: atnode(:),softdat(:)



      integer,allocatable :: ncnodeIndex(:),icnodeIndex(:,:),neighbours(:)
      integer,allocatable :: cnodeid(:)
      integer,allocatable :: cnodeidIndex(:,:)
      real,allocatable :: cdfvalIndex(:)
      integer,allocatable :: resultsIdxIndex(:,:)      

      end module
!
!
!
      program main
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                                                                      %
! Copyright (C) 1996, The Board of Trustees of the Leland Stanford     %
! Junior University.  All rights reserved.                             %
!                                                                      %
! The programs in GSLIB are distributed in the hope that they will be  %
! useful, but WITHOUT ANY WARRANTY.  No author or distributor accepts  %
! responsibility to anyone for the consequences of using them or for   %
! whether they serve any particular purpose or work at all, unless he  %
! says so in writing.  Everyone is granted permission to copy, modify  %
! and redistribute the programs in GSLIB, but only under the condition %
! that this notice and the above copyright notice remain intact.       %
!                                                                      %
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
!           Conditional Simulation of a 3-D Rectangular Grid
!           ************************************************
!
! The output file will be a GEOEAS file containing the simulated values
! The file is ordered by x,y,z, and then simulation (i.e., x cycles
! fastest, then y, then z, then simulation number).
!
!
!
!-----------------------------------------------------------------------
      use geostat
      use global
      !use omp_lib
      !include  'sisim.inc'
      implicit none

      integer numThreads

      VERSION = 1.000

      lin  = 1
      lout = 2
      ldbg = 3


      numThreads=1
!!$omp parallel default(shared)
!      numThreads = OMP_get_num_threads()
!      print *,numThreads
!!$omp end parallel





!      include  'sisim_new.inc'
!
! Read the Parameter File and the Data:
!
      call readparm
      !call readparm(MAXNOD,MAXCXY,MAXCTX,MAXCTY,MAXCTZ,MAXKR1,MAXSBX, &
      !              MAXSBY,MAXSBZ,MAXCUT,MXCUT,MAXROT)
!
! Call sisim for the simulation:
!
      call sisim
      !call sisim   (MAXNOD,MAXCXY,MAXCTX,MAXCTY,MAXCTZ,MAXKR1,MAXSBX, &
      !              MAXSBY,MAXSBZ,MAXCUT,MXCUT,MAXROT)
!
! Finished:
!
      close(lout)
      close(ldbg)
      write(*,9998) VERSION
 9998 format(/' SISIM Version: ',f5.3, ' Finished'/)
      stop
      end
 
 

      subroutine readparm
      !subroutine readparm(MAXNOD,MAXCXY,MAXCTX,MAXCTY,MAXCTZ,MAXKR1, &
      !                    MAXSBX,MAXSBY,MAXSBZ,MAXCUT,MXCUT,MAXROT)
!-----------------------------------------------------------------------
!
!                  Initialization and Read Parameters
!                  **********************************
!
! The input parameters and data are read in from their files. Some quick
! error checking is performed and the statistics of all the variables
! being considered are written to standard output.
!
! NOTE: 1. The variables and the data are allocated in common blocks
!          (sisim.inc)
!
!
!-----------------------------------------------------------------------
      !use       msflib
      use       geostat
      use global
      use       grid_info
      use        random2
      implicit none
      !include  'sisim.inc'
!      include  'sisim_new.inc'
      integer,parameter :: MV=20
      real      var(MV),test1
      real*8    p,acorni
      integer,allocatable :: ivrs(:)
!      integer   test
      character datafl*512,tabfl*512,softfl*512,outfl*512,dbgfl*512, &
                str*512,title*80,grid_fl*512
      logical   testfl
      real(kind=8), dimension ( 3 ) :: angles
      integer i,j,k,jj,index,nvari,ic
      real*8 xmax1,ymax1,zmax1,xmax2,ymax2,zmax2
      real*8 xx,yy,zz
      integer gin,ierr,icut,istart1,istarti,ist,index1,indexi
      real cdfval,clos
      real*8 zval
      real*4 t1,t2

      t1=secnds(0.0)
!
! Fortran unit numbers needed:
!
      lin  = 1
      lout = 2
      ldbg = 3
!
! Note VERSION number:
!
      write(*,9999) VERSION
 9999 format(/' SISIM_LVA Version: ',f5.3/)
!
! Get the name of the parameter file - try the default name if no input:
!
      do i=1,512
            str(i:i) = ' '
      end do
      call getarg(1,str)
      if(str(1:1).eq.' ')then
            write(*,*) 'Which parameter file do you want to use?'
            read (*,'(a)') str
      end if
      if(str(1:1).eq.' ') str(1:20) = 'sisim_lva.par           '
      inquire(file=str,exist=testfl)
      if(.not.testfl) then
            write(*,*) 'ERROR - the parameter file does not exist,'
            write(*,*) '        check for the file and try again  '
            write(*,*)
            if(str(1:20).eq.'sisim_lva.par           ') then
                  write(*,*) '        creating a blank parameter file'
                  call makepar
                  write(*,*)
            end if
            stop
      endif
      open(lin,file=str,status='OLD')
!
! Find Start of Parameters:
!
 1    read(lin,'(a4)',end=98) str(1:4)
      if(str(1:4).ne.'STAR') go to 1
!
! Read Input Parameters:
!

      read(lin,*,err=98) ivtype
      write(*,*) ' variable type (1=continuous, 0=categorical)= ',ivtype

      read(lin,*,err=98) ncut
      write(*,*) ' number of thresholds / categories = ',ncut
!
! Find the needed parameters:
!
      MAXCUT = ncut
      MAXROT = MAXCUT * MAXNST + 1
      MXCUT = MAXCUT + 1
!
! Allocate the needed memory:
!34
      allocate(thres(MAXCUT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 34: Allocation failed due to ', &
                        'insuffiecent memory.'
                  stop
            end if
!35
      allocate(cdf(MAXCUT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 35: Allocation failed due to ', &
                        'insuffiecent memory.'
                  stop
            end if
!36
      allocate(ccdf(MAXCUT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 36: Allocation failed due to ', &
                        'insuffiecent memory.'
                  stop
            end if
!37
      allocate(ccdfo(MAXCUT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 37: Allocation failed due to ', &
                        'insuffiecent memory.'
                  stop
            end if
!38
      allocate(beez(MAXCUT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 38: Allocation failed due to ', &
                        'insuffiecent memory.'
                  stop
            end if
!39
      allocate(c0(MAXCUT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 39: Allocation failed due to ', &
                        'insuffiecent memory.'
                  stop
            end if
!40
      allocate(cmax(MAXCUT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 40: Allocation failed due to ', &
                        'insuffiecent memory.'
                  stop
            end if
!41
      allocate(cc(MAXCUT*MAXNST),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 41: Allocation failed due to ', &
                        'insuffiecent memory.'
                  stop
            end if
!42
      allocate(aa(MAXCUT*MAXNST),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 42: Allocation failed due to ', &
                        'insuffiecent memory.'
                  stop
            end if
!43
      allocate(ang1(MAXCUT*MAXNST),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 43: Allocation failed due to ', &
                        'insuffiecent memory.'
                  stop
            end if
!44
      allocate(ang2(MAXCUT*MAXNST),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 44: Allocation failed due to ', &
                        'insuffiecent memory.'
                  stop
            end if
!45
      allocate(ang3(MAXCUT*MAXNST),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 45: Allocation failed due to ', &
                        'insuffiecent memory.'
                  stop
            end if
!46
      allocate(anis1(MAXCUT*MAXNST),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 46: Allocation failed due to ', &
                        'insuffiecent memory.'
                  stop
            end if
!47
      allocate(anis2(MAXCUT*MAXNST),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 47: Allocation failed due to ', &
                        'insuffiecent memory.'
                  stop
            end if
!48
      allocate(aviol(MAXCUT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 48: Allocation failed due to ', &
                        'insuffiecent memory.'
                  stop
            end if
!49
      allocate(xviol(MAXCUT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 49: Allocation failed due to ', &
                        'insuffiecent memory.'
                  stop
            end if
!50
      allocate(it(MAXCUT*MAXNST),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 50: Allocation failed due to ', &
                        'insuffiecent memory.'
                  stop
            end if
!51
      allocate(nst(MAXCUT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 51: Allocation failed due to ', &
                        'insuffiecent memory.'
                  stop
            end if
!52
      allocate(nviol(MAXCUT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 52: Allocation failed due to ', &
                        'insuffiecent memory.'
                  stop
            end if
!53
      allocate(rotmat(MAXROT,3,3),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 53: Allocation failed due to ', &
                        'insuffiecent memory.'
                  stop
            end if
!54
      allocate(ivrs(MAXCUT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 54: Allocation failed due to ', &
                        'insuffiecent memory.'
                  stop
            end if
!
      read(lin,*,err=98) (thres(i),i=1,ncut)
      write(*,*) ' thresholds / categories = ',(thres(i),i=1,ncut)

      read(lin,*,err=98) (cdf(i),i=1,ncut)
      write(*,*) ' global cdf / pdf        = ',(cdf(i),i=1,ncut)

      read(lin,'(a512)',err=98) datafl
      call chknam(datafl,512)
      write(*,*) ' data file = ',datafl(1:40)

      read(lin,*,err=98) ixl,iyl,izl,ivrl
      write(*,*) ' input columns = ',ixl,iyl,izl,ivrl

      read(lin,'(a512)',err=98) softfl
      call chknam(softfl,512)
      write(*,*) ' soft data file = ',softfl(1:40)
      inquire(file=softfl,exist=testfl)

      if(testfl) then
            read(lin,*,err=98) ixs,iys,izs,(ivrs(i),i=1,ncut)
            write(*,*) ' columns = ',ixs,iys,izs,(ivrs(i),i=1,ncut)
            read(lin,*,err=98) imbsim
            write(*,*) ' Markov-Bayes simulation = ',imbsim
            if(imbsim.eq.1) then
                  read(lin,*,err=98) (beez(i),i=1,ncut)
            else
                  read(lin,*,err=98)
            end if
      else
            read(lin,*,err=98)
            read(lin,*,err=98)
            read(lin,*,err=98)
      end if

      read(lin,*,err=98) tmin,tmax
      write(*,*) ' trimming limits      ',tmin,tmax

      read(lin,*,err=98) zmin,zmax
      write(*,*) ' data limits (tails)  ',zmin,zmax

      read(lin,*,err=98) ltail,ltpar
      write(*,*) ' lower tail = ',ltail,ltpar

      read(lin,*,err=98) middle,mpar
      write(*,*) ' middle = ',middle,mpar

      read(lin,*,err=98) utail,utpar
      write(*,*) ' upper tail = ',utail,utpar

      read(lin,'(a512)',err=98) tabfl
      call chknam(tabfl,512)
      write(*,*) ' file for tab. quant. ',tabfl(1:40)

      read(lin,*,err=98) itabvr,itabwt
      write(*,*) ' columns for vr wt = ',itabvr,itabwt

      read(lin,*,err=98) idbg
      write(*,*) ' debugging level = ',idbg

      read(lin,'(a512)',err=98) dbgfl
      call chknam(dbgfl,512)
      write(*,*) ' debugging file = ',dbgfl(1:40)

      read(lin,'(a512)',err=98) outfl
      call chknam(outfl,512)
      write(*,*) ' output file = ',outfl(1:40)

      read(lin,*,err=98) nsim
      write(*,*) ' number of simulations = ',nsim

      read(lin,*,err=98) nx,xmn,xsiz
      write(*,*) ' X grid specification = ',nx,xmn,xsiz

      read(lin,*,err=98) ny,ymn,ysiz
      write(*,*) ' Y grid specification = ',ny,ymn,ysiz

      read(lin,*,err=98) nz,zmn,zsiz
      write(*,*) ' Z grid specification = ',nz,zmn,zsiz
      nxy  = nx*ny
      nxyz = nx*ny*nz

      read(lin,*,err=98) ixv(1)
      write(*,*) ' random number seed = ',ixv(1)
      do i=1,1000
             p = acorni(idum)
      end do

      read(lin,*,err=98) ndmax
      write(*,*) ' ndmax = ',ndmax

      read(lin,*,err=98) nodmax
      write(*,*) ' max prev sim nodes = ',nodmax

      read(lin,*,err=98) maxsec
      write(*,*) ' max soft indicator data = ',maxsec

      read(lin,*,err=98) sstrat
      write(*,*) ' search strategy = ',sstrat

      read(lin,*,err=98) mults,nmult
      write(*,*) ' multiple grid search flag = ',mults,nmult

      read(lin,*,err=98) noct
      write(*,*) ' max per octant = ',noct

      read(lin,*,err=98) radius,radius1,radius2
      write(*,*) ' search radii = ',radius,radius1,radius2
      if(radius.lt.EPSLON) stop 'radius must be greater than zero'
      radsqd = radius  * radius
      sanis1 = radius1 / radius
      sanis2 = radius2 / radius

      read(lin,*,err=98) sang1,sang2,sang3
      write(*,*) ' search anisotropy angles = ',sang1,sang2,sang3

      read(lin,*,err=98) mxctx,mxcty,mxctz
      write(*,*) ' size of covariance lookup = ',mxctx,mxcty,mxctz
      
      read(lin,*,err=98) mik,cutmik
      write(*,*) ' median IK switch = ',mik,cutmik

      read(lin,*,err=98) ktype
      write(*,*) ' kriging type switch = ',ktype

!
! Output now goes to debugging file:
!
      open(ldbg,file=dbgfl,status='UNKNOWN')
      do i=1,ncut
            read(lin,*,err=98) nst(i),c0(i)
            if(ivtype.eq.0) &
            write(ldbg,100)  i,thres(i),cdf(i),nst(i),c0(i)
            if(ivtype.eq.1) &
            write(ldbg,101)  i,thres(i),cdf(i),nst(i),c0(i)
            write(*,*)nst(i),MAXNST 
            if(nst(i).gt.MAXNST) stop 'nst is too big'
            istart = 1 + (i-1)*MAXNST
            do j=1,nst(i)
                  index = istart + j - 1
                  read(lin,*,err=98) it(index),cc(index),ang1(index),&
                                     ang2(index),ang3(index)
                  if(it(index).eq.3) STOP 'Gaussian Model Not Allowed!'
                  read(lin,*,err=98) aa(index),aa1,aa2
                  !write(ldbg,102)  j,it(index),aa(index),cc(index)
                  anis1(index) = aa1 / max(EPSLON,aa(index))
                  anis2(index) = aa2 / max(EPSLON,aa(index))
                  !write(ldbg,103) ang1(index),ang2(index),ang3(index), &
                  !                anis1(index),anis2(index)
            end do
      end do
!     close(lin)
 
     
!   Inicio de Lectura de parámetros de LVA.
!
      read(lin,'(a512)',err=98) grid_fl
      call chknam(grid_fl,512)
      write(*,*) ' grid file to use:',grid_fl(1:40)
            
      read(lin,*,err=98) ia1g,ia2g,ia3g,ir1,ir2 !read in col for angles (1-3) and ratios (1-2)
      write(*,'(a20,5(I4))') ' columns for LVA grid = ',ia1g,ia2g,ia3g, &
           ir1,ir2

      read(lin,*,err=98) n(1),o(1),sss(1)
      write(*,'(a20,I4,f12.4,xx,f12.4)') ' ngx, xgmn, xgsiz = ',n(1), &
           o(1),sss(1)

      read(lin,*,err=98) n(2),o(2),sss(2)
      write(*,'(a20,I4,f12.4,xx,f12.4)') ' ngy, ygmn, ygsiz = ',n(2), &
           o(2),sss(2)

      read(lin,*,err=98) n(3),o(3),sss(3)
      write(*,'(a20,I4,f12.4,xx,f12.4)') ' ngz, zgmn, zgsiz = ',n(3), &
           o(3),sss(3)

      refine=1
      
      read(lin,*,err=98) graph_offset
      if(graph_offset/=int(graph_offset) .or. graph_offset<1) then
          write(*,*) 'ERROR WITH THE graph_offset VALUE', &
                'SETTING graph_offset TO 1'
          graph_offset=1
      end if
      write(*,*) ' offset parameter',graph_offset
      
      read(lin,*,err=98) MDS_opt
      
      
      cal_stress=0
      if(MDS_opt <0) then
      
      write(*,*) '****will calculate stress******'
      cal_stress = 1
      MDS_opt=-MDS_opt
      end if
      
      
      write(*,*) ' Use MDS:',MDS_opt
      
      read(lin,*,err=98) xland,yland,zland
      write(*,*) ' spacing of landmark points:',xland,yland,zland
      
      read(lin,*,err=98) dim
      write(*,*) ' max number of dimens to use:',dim
      
      !quick check on grid coverage:
      xmax1=xmn+nx*xsiz;ymax1=ymn+ny*ysiz;zmax1=zmn+nz*zsiz
      xmax2=o(1)+n(1)*sss(1)
      ymax2=o(2)+n(2)*sss(2)
      zmax2=o(3)+n(3)*sss(3)
!  if(xmax1>xmax2 .or. ymax1>ymax2 .or. zmax1>zmax2 .or. o(1)>xmn .or. o(2)>ymn .or. o(3)>zmn) &
!  write(*,*)'**WARNING** LVA GRID MAY NOT COVER KRIGING AREA, MAY CAUSE ERROR'
   

      !read in the grid and store the necessary angles etc.
      inquire(file=grid_fl,exist=testfl)
      if(.not.testfl) then
            write(*,*) 'ERROR - the grid file does not exist,'
            write(*,*) '        check for the file and try again  '
            write(*,*)
            stop
      endif
      open(gin,file=grid_fl,status='OLD')

      read(gin,*)
      read(gin,*,err=99)       nvari
      do i=1,nvari
            read(gin,*)
      end do
      
       allocate(grid(nx,ny,nz,5))


      if(nx /=n(1) .or. ny /= n(2) .or. nz /= n(3).or. &
        xmn.ne.o(1).or.ymn.ne.o(2).or.zmn.ne.o(3).or. &
        xsiz /= sss(1) .or. ysiz /= sss(2) .or. zsiz /= sss(3) )then
       
        write(*,*) 
        write(*,*) 
        write(*,*) '*******WARNING*******'
        write(*,*) 'simulation grid does not match LVA grid'
        write(*,*) 'will regrid the LVA grid to match  simulation grid'
        write(*,*) 
        write(*,*)
      
      
      
        temp_o = o
        temp_sss = sss
        temp_n = n
        allocate(temp_grid(n(1),n(2),n(3),5))
        n(1)=nx ; n(2)=ny ; n(3)=nz
        o(1)=xmn ; o(2)=ymn ; o(3)=zmn
        sss(1) = xsiz ;sss(2) = ysiz ;sss(3) = zsiz
      
      
     
         do k=1,temp_n(3)
          do j=1,temp_n(2)
           do i=1,temp_n(1)
            read(gin,*,err=98)  (var(jj),jj=1,nvari)
            temp_grid(i,j,k,1)=var(ia1g)
            temp_grid(i,j,k,2)=var(ia2g)
            temp_grid(i,j,k,3)=var(ia3g)
            temp_grid(i,j,k,4)=var(ir1)
            temp_grid(i,j,k,5)=var(ir2)
           end do
          end do
         end do
         
         do k=1,n(3)
          do j=1,n(2)
           do i=1,n(1)
           
           !where are we in the model (xx,yy,zz)
           
           xx = xmn + real(i-1)*xsiz
           yy = ymn + real(j-1)*ysiz
            zz = zmn + real(k-1)*zsiz
            !where is this in the LVA grid
            call getindx(temp_n(1),temp_o(1),temp_sss(1),xx,locx,inflag)
            call getindx(temp_n(2),temp_o(2),temp_sss(2),yy,locy,inflag)
            call getindx(temp_n(3),temp_o(3),temp_sss(3),zz,locz,inflag)
            
           !fill in the grid
            grid(i,j,k,:) = temp_grid(locx,locy,locz,:)
            
           end do
          end do
         end do
        deallocate(temp_grid)
     
      else !read in per normal
     
     
      do k=1,n(3)
      do j=1,n(2)
       do i=1,n(1)
        read(gin,*,err=98)  (var(jj),jj=1,nvari)
        grid(i,j,k,1)=var(ia1g)
        grid(i,j,k,2)=var(ia2g)
        grid(i,j,k,3)=var(ia3g)
        grid(i,j,k,4)=var(ir1)
        grid(i,j,k,5)=var(ir2)
       end do
      end do
      end do
     
      end if !loop over regrid
     
      do k=1,n(3)
      do j=1,n(2)
       do i=1,n(1)
        !if necessary, fix the three angles read in
        if(k==1 .and. j==1 .and. i==1  .and. grid(i,j,k,4)==0) then
        write(*,*) 'WARNING***** you cannot have an anisotropy ratio=0', &
       '  it will be reset to 1:1****'
        end if
        
        if(k==1 .and. j==1 .and. i==1  .and. grid(i,j,k,5)==0) then
        write(*,*) 'WARNING***** you cannot have an anisotropy ratio=0', &
       'it will be reset to 1:1****'
        end if
        
        if(grid(i,j,k,4)==0) grid(i,j,k,4)=1
        if(grid(i,j,k,5)==0) grid(i,j,k,5)=1
        
        angles(1) = grid(i,j,k,1)  
        angles(2) = grid(i,j,k,2)  
        angles(3) = grid(i,j,k,3) 
        call fix_angles  ( angles )
        grid(i,j,k,1) = angles(1) 
        grid(i,j,k,2) = angles(2)  
        grid(i,j,k,3) = angles(3)  
       end do
      end do
      
      end do
      close (gin)
      assign_nodes=1
     
      read(lin,*,err=98) radius_lva
      write(*,'(a,f18.4)') ' search radii = ',radius_lva
      
      read(lin,*,err=98) d_tree
      write(*,'(a,I4)') ' search dim to use = ',d_tree
      
      if(radius.lt.EPSLON) stop 'radius must be greater than zero'
      radsqd_lva = radius_lva  * radius_lva
      sanis1_lva = radius_lva / radius_lva
      sanis2_lva = radius_lva / radius_lva
      sang1_lva=0
      sang2_lva=0
      sang3_lva=0
!      
!      
!      


      
!
!     Fin Lectura de parámetros de LVA.  
 
 
!
! Find the needed parameters:
!
      MAXX = nx
      MAXY = ny
      MAXZ = nz
      MXYZ = MAXX * MAXY * MAXZ
      MAXCTX = mxctx
      MAXCTY = mxcty
      MAXCTZ = mxctz
      MAXCXY = MAXCTX * MAXCTY
      MAXXYZ = MAXCTX * MAXCTY * MAXCTZ
      MAXSAM = ndmax 
      MAXNOD = nodmax
      MAXKR1 = 2 * MAXNOD + 2 * MAXSAM + 1
      MAXKR2 = MAXKR1 * MAXKR1
      MAXSBX = 1
      if(nx.gt.1)then
            MAXSBX = int(nx/2)
            if(MAXSBX.gt.50)MAXSBX=50
      end if
!
      MAXSBY = 1
      if(ny.gt.1)then
            MAXSBY = int(ny/2)
            if(MAXSBY.gt.50)MAXSBY=50
      end if
!
      MAXSBZ = 1
      if(nz.gt.1)then
            MAXSBZ = int(nz/2)
            if(MAXSBZ.gt.50)MAXSBZ=50
      end if
!
      MAXSB = MAXSBX * MAXSBY * MAXSBZ
      MXSXY = 4 * MAXSBX * MAXSBY
      MXSX = 2 * MAXSBX
      av = 0.0
      ss = 0.0
!
! Find the paramater MAXDAT:
!
      MAXDAT = 1
      inquire(file=datafl,exist=testfl)
      if(testfl) then
            open(lin,file=datafl,status='OLD')
            read(lin,*,err=99)
            read(lin,*,err=99)       nvari
            do i=1,nvari
            read(lin,'()',err=99)
            end do
            MAXDAT = 0
 55         read(lin,*,end=66,err=99) (var(j),j=1,nvari)
            if(var(ivrl).lt.tmin.or.var(ivrl).ge.tmax)go to 55
            MAXDAT = MAXDAT + 1
            go to 55
 66         continue
            rewind(lin)
            close(lin)
      end if
!
! Find the paramater MAXTAB:
!
      MAXTAB = ncut
      inquire(file=tabfl,exist=testfl)
      if(testfl) then
            open(lin,file=tabfl,status='OLD')
            read(lin,*,err=99)
            read(lin,*,err=99)       nvari
            do i=1,nvari
                  read(lin,'()',err=99)
            end do
            MAXTAB = 0
 77         read(lin,*,end=88,err=99) (var(j),j=1,nvari)
            if(var(ivrl).lt.tmin.or.var(ivrl).ge.tmax)go to 77
            MAXTAB = MAXTAB + 1
            go to 77
 88         continue
            rewind(lin)
            close(lin)
      end if
!
! Allocate the needed memory:
write(*,*)'MAXDAT',MAXDAT
      allocate(cdfvalIndex(nxyz),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 1: Allocation failed due to ', &
                        'insufficient memory.'
                  stop
            end if
!1
      allocate(x(MAXDAT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 1: Allocation failed due to ', &
                        'insufficient memory.'
                  stop
            end if
!2
      allocate(y(MAXDAT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 2: Allocation failed due to ', &
                        'insufficient memory.'
                  stop
            end if
!3
      allocate(z(MAXDAT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 3: Allocation failed due to ', &
                        'insufficient memory.'
                  stop
            end if
!4
      allocate(vr(MAXDAT,MXCUT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 4: Allocation failed due to ', &
                        'insufficient memory.'
                  stop
            end if
!5
      allocate(aclose(MAXKR1),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 5: Allocation failed due to ', &
                        'insufficient memory.'
                  stop
            end if
!
      allocate(close(MAXDAT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 5: Allocation failed due to ', &
                        'insufficient memory.'
                  stop
            end if
      allocate(neighbours(MAXKR1),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 5: Allocation failed due to ', &
                        'insufficient memory.'
                  stop
            end if


      allocate(distance(MAXKR1),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 5: Allocation failed due to ', &
                        'insufficient memory.'
                  stop
            end if
!6
      allocate(actloc(MAXDAT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 6: Allocation failed due to ', &
                        'insufficient memory.'
                  stop
            end if
!7
      allocate(tmpdat(MAXDAT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 7: Allocation failed due to ', &
                        'insufficient memory.'
                  stop
            end if
!8
      allocate(sim(MXYZ),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 8: Allocation failed due to ', &
                        'insufficient memory.'
                  stop
            end if
!9
      MAXORD = MXYZ
      if(MXYZ.lt.MAXCXY) MAXORD=MAXCXY
      allocate(order(MAXORD),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 9: Allocation failed due to ', &
                        'insufficient memory.'
                  stop
            end if
!10
      allocate(tmp(MAXORD),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 10: Allocation failed due to ', &
                        'insufficient memory.'
                  stop
            end if
!11
      allocate(gcut(MAXTAB),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 11: Allocation failed due to ', &
                        'insufficient memory.'
                  stop
            end if
!12
      allocate(gcdf(MAXTAB),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 12: Allocation failed due to ', &
                        'insufficient memory.'
                  stop
            end if
!13
      allocate(covtab(MAXCTX,MAXCTY,MAXCTZ,MAXCUT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 13: Allocation failed due to ', &
                        'insuffiecent memory.'
                  stop
            end if
!13.5
      allocate(cnodeid(MAXNOD),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 14: Allocation failed due to ', &
                       'insuffiecent memory.'
                  stop
            end if
!13.75            
      allocate(cnodeidIndex(MAXNOD,MXYZ),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 14: Allocation failed due to ', &
                       'insuffiecent memory.'
                  stop
            end if

!14
      allocate(cnodex(MAXNOD),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 14: Allocation failed due to ', &
                        'insuffiecent memory.'
                  stop
            end if
!15
      allocate(cnodey(MAXNOD),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 15: Allocation failed due to ', &
                        'insuffiecent memory.'
                  stop
            end if
!16
      allocate(cnodez(MAXNOD),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 16: Allocation failed due to ', &
                        'insuffiecent memory.'
                  stop
            end if
!17
      allocate(cnodev(MAXNOD),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 17: Allocation failed due to ', &
                        'insuffiecent memory.'
                  stop
            end if
!18
      allocate(cnodet(MAXNOD),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 18: Allocation failed due to ', &
                        'insuffiecent memory.'
                  stop
            end if
!19
      allocate(vra(MAXKR1),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 19: Allocation failed due to ', &
                        'insuffiecent memory.'
                  stop
            end if
!20
      allocate(r(MAXKR1),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 20: Allocation failed due to ', &
                        'insuffiecent memory.'
                  stop
            end if
!21
      allocate(rr(MAXKR1),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 21: Allocation failed due to ', &
                        'insuffiecent memory.'
                  stop
            end if
!22
      allocate(s(MAXKR1),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 22: Allocation failed due to ', &
                        'insuffiecent memory.'
                  stop
            end if
!23
      allocate(a(MAXKR2),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 23: Allocation failed due to ', &
                        'insuffiecent memory.'
                  stop
            end if
!24
      allocate(ixnode(MAXXYZ),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 24: Allocation failed due to ', &
                        'insuffiecent memory.'
                  stop
            end if
!25
      allocate(iynode(MAXXYZ),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 25: Allocation failed due to ', &
                        'insuffiecent memory.'
                  stop
            end if
!26
      allocate(iznode(MAXXYZ),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 26: Allocation failed due to ', &
                        'insuffiecent memory.'
                  stop
            end if
!27
      allocate(nisb(MAXSB),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 27: Allocation failed due to ', &
                        'insuffiecent memory.'
                  stop
            end if
!28
      allocate(icnode(MAXNOD),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 28: Allocation failed due to ', &
                        'insuffiecent memory.'
                  stop
            end if
!28.5
      allocate(icnodeIndex(MAXNOD,MXYZ),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 28: Allocation failed due to ', &
                       'insuffiecent memory.'
                  stop
            end if
      allocate(ncnodeIndex(MXYZ),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 28: Allocation failed due to ', &
                       'insuffiecent memory.'
                  stop
            end if




!29
      allocate(ixsbtosr(8*MAXSB),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 29: Allocation failed due to ', &
                        'insuffiecent memory.'
                  stop
            end if
!30
      allocate(iysbtosr(8*MAXSB),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 30: Allocation failed due to ', &
                        'insuffiecent memory.'
                  stop
            end if
!31
      allocate(izsbtosr(8*MAXSB),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 31: Allocation failed due to ', &
                        'insuffiecent memory.'
                  stop
            end if
!32
      allocate(atnode(MAXDAT),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 32: Allocation failed due to ', &
                        'insuffiecent memory.'
                  stop
            end if
!33
      allocate(softdat(MAXKR1),stat=test)
            if(test.ne.0)then
                  write(*,*)'ERROR 33: Allocation failed due to ', &
                        'insuffiecent memory.'
                  stop
            end if
!
 100  format(/,' Category  number ',i2,' = ',f12.3,/, &
               '           global prob value = ',f8.4,/, &
               '           number of structures = ',i3,/, &
               '           nugget effect        = ',f8.4)
 101  format(/,' Threshold number ',i2,' = ',f12.3,/, &
               '           global prob value = ',f8.4,/, &
               '           number of structures = ',i3,/, &
               '           nugget effect        = ',f8.4)
 102  format(  '           type of structure ',i3,' = ',i3,/, &
               '           aa parameter         = ',f12.4,/, &
               '           cc parameter         = ',f12.4)
 103  format(  '           ang1, ang2, ang3     = ',3f6.2,/, &
               '           anis1, anis2         = ',2f12.4)
!
! Perform some quick error checking:
!
      if(nx.gt.MAXX) stop 'nx is too big - modify .inc file'
      if(ny.gt.MAXY) stop 'ny is too big - modify .inc file'
      if(nz.gt.MAXZ) stop 'nz is too big - modify .inc file'
!
! Check to make sure the data file exists, then either read in the
! data or write a warning:
!
      title = 'SISIM SIMULATIONS:                      '//&
              '                                        '
      nd = 0
      inquire(file=datafl,exist=testfl)
      if(.not.testfl) then
            write(*,113) datafl
 113        format('WARNING data file ',a40,' does not exist!',/, &
                   ' Hope your intention was to create an', &
                   ' unconditional simulation.')
      else
!
! The data file exists so open the file and read in the header
! information.
!
            write(*,*) 'Reading input data'
            av = 0.0
            ss = 0.0
            open(lin,file=datafl,status='OLD')
            read(lin,'(a60)',err=99) title(21:80)
            read(lin,*,err=99)       nvari
            do i=1,nvari
                  read(lin,'()',err=99)
            end do
!
! Read all the data until the end of the file:
!
 5          read(lin,*,end=6,err=99) (var(j),j=1,nvari)
            vrt = var(ivrl)
            if(vrt.lt.tmin.or.vrt.ge.tmax) go to 5
            nd  = nd + 1
            x(nd) = xmn
            y(nd) = ymn
            z(nd) = zmn
            if(ixl.gt.0) x(nd) = var(ixl)
            if(iyl.gt.0) y(nd) = var(iyl)
            if(izl.gt.0) z(nd) = var(izl)
            av = av + vrt
            ss = ss + vrt*vrt
!
! The indicator data are constructed knowing the thresholds and the
! data value.
!
            if(ivtype.eq.0) then
                  do ic=1,ncut
                        vr(nd,ic) = 0.0
                        if(int(vrt+0.5).eq.int(thres(ic)+0.5))& 
                        vr(nd,ic) = 1.0
                  end do
            else
                  do ic=1,ncut
                        vr(nd,ic) = 1.0
                        if(vrt.gt.thres(ic)) vr(nd,ic) = 0.0
                  end do
            end if
            vr(nd,MXCUT) = vrt
            go to 5
6           close(lin)
           
!
!Edicion RGE
      allocate(data_ind(MAXDAT))
      !write(*,*)'x = ',x
      !write(*,*)'y = ',y
      !write(*,*)'nd = ',nd
      allocate(x_loc(MAXDAT))
      allocate(y_loc(MAXDAT))
      allocate(z_loc(MAXDAT))

      do i=1,nd
      

         xind = int( (x(i)-xmn)/xsiz + 1.5 )
         yind = int( (y(i)-ymn)/ysiz + 1.5 )
         zind = int( (z(i)-zmn)/zsiz + 1.5 )
         
         !write(*,*)'xind = ',xind
         !write(*,*)'yind = ',yind

    
         data_ind(i)=get_node(xind,yind,zind,nx,ny,nz)
         
      end do
        
        !write(*,*)'data_ind = ',data_ind

 
!fIN EDICION RGE
!
           
 
!     
! Compute the averages and variances as an error check for the user:
!
            xd = max(real(nd),1.0)
            av = av / xd
            ss =(ss / xd ) - av * av
            write(*,120)    ivrl,nd,av,ss
            write(ldbg,120) ivrl,nd,av,ss
 120        format(/,'  Data for SISIM: Variable number ',i2, &
                   /,'  Number of acceptable data  = ',i8, &
                   /,'  Equal Weighted Average     = ',f12.4, &
                   /,'  Equal Weighted Variance    = ',f12.4,/)
!
! Check to make sure that the grid is compatible with the data:
!
            if(ixl.le.0.and.nx.gt.1) then
               write(*,*) 'ERROR there is no X coordinate in data file'
               write(*,*) '      nx must be set to 1'
               stop
            end if
            if(iyl.le.0.and.ny.gt.1) then
               write(*,*) 'ERROR there is no Y coordinate in data file'
               write(*,*) '      ny must be set to 1'
               stop
            end if
            if(izl.le.0.and.nz.gt.1) then
               write(*,*) 'ERROR there is no Z coordinate in data file'
               write(*,*) '      nz must be set to 1'
               stop
            end if
      endif
!
! Now, if required, read in the tabulated values for details of the dist
!
      if(ltail.eq.3.or.middle.eq.3.or.utail.eq.3) then
            ng = 0
            inquire(file=tabfl,exist=testfl)
            if(.not.testfl) stop 'ERROR tabfl does not exist'
            open(lin,file=tabfl,status='OLD')
            read(lin,*,err=97)
            read(lin,*,err=97) nvari
            do i=1,nvari
                  read(lin,*,err=97)
            end do
            tcdf = 0.0
            ng   = 0
 21         read(lin,*,end=22,err=97) (var(j),j=1,nvari)
            if(var(itabvr).lt.tmin.or.var(itabvr).ge.tmax) go to 21
            ng = ng + 1
            gcut(ng) = var(itabvr)
            gcdf(ng) = 1.0
            if(itabwt.gt.0) gcdf(ng) = var(itabwt)
            tcdf = tcdf + gcdf(ng)
            go to 21
 22         close(lin)
!
! Sort in ascending order and keep track of where the tabulated values
! switch classes:
!
            if(tcdf.le.0.0) then
                  write(*,*) 'ERROR: either the weights are zero or'
                  write(*,*) '       there are no tabulated data.'
                  stop
            endif
            !call sortem(1,ng,gcut,1,gcdf,c,d,e,f,g,h)
            call sortem(1,ng,gcut,1,gcdf)
!
! Set up gcdf for tabulated quantiles:
!
            oldcp = 0.0
            cp    = 0.0
            tcdf  = 1.0 / tcdf
            do i=1,ng
                  cp      = cp + gcdf(i) * tcdf
                  gcdf(i) =(cp + oldcp) * 0.5
                  oldcp   = cp
            end do
      end if
!
! Direct input of indicator data:
!
      nhd = nd
      inquire(file=softfl,exist=testfl)
      if(testfl) then
            write(*,*)
            write(*,*) 'Reading soft indicator data'
            open(lin,file=softfl,status='OLD')
            read(lin,*,err=97)
            read(lin,*,err=97) nvari
            if(ivrs(ncut).gt.nvari) then
                  write(*,*) ' ERROR: too few variables in ',softfl
                  write(*,*) '        inconsistent with parameters'
                  stop
            end if
            do i=1,nvari
                  read(lin,*,err=97)
            end do
 12         read(lin,*,end=13,err=96) (var(j),j=1,nvari)
!
! Don't keep soft data co-located with hard data:
!
            xx = xmn
            yy = ymn
            zz = zmn
            if(ixs.gt.0) xx = var(ixs)
            if(iys.gt.0) yy = var(iys)
            if(izs.gt.0) zz = var(izs)
            do i=1,nhd
                  test1 = abs(xx-x(i)) + abs(yy-y(i)) + abs(zz-z(i))
                  if(test1.le.EPSLON) go to 12
            end do
!
! Accept this data:
!
            nd = nd + 1
            x(nd) = xx
            y(nd) = yy
            z(nd) = zz
            do j=1,ncut
                  i = ivrs(j)
                  vr(nd,j) = var(i)
                  ccdf(j)  = var(i)
            end do
!
! Draw a value for this soft distribution (in case the distribution is
! co-located with a grid node and Markov-Bayes is not used):
!
            call init_genrand(ixv(1))  !intial the random generator

            !cdfval = real(acorni(idum))
            cdfval = real(grnd())
            call ordrel(ivtype,ncut,ccdf,ccdfo,nviol,aviol,xviol)
            zval = UNEST
            call beyond(ivtype,ncut,thres,ccdfo,ng,gcut,gcdf,zmin,zmax,&
                        ltail,ltpar,middle,mpar,utail,utpar,zval,&
                        cdfval,ierr)
            vr(nd,MXCUT) = zval
write(*,*)'vr-size:',size(vr)
!
! If performing median IK then check for missing values:
!
            if(mik.eq.1) then
                  do ic=1,ncut
                        if(vr(nd,ic).lt.0.0) then
                              write(*,150) softfl
                              stop
                        endif
                  end do
 150              format(' Since the median IK approach is being',&
                         ' considered no missing values are',&
                         ' allowed',/,' Check file ',a40)
            endif
            go to 12
 13         close(lin)
      endif
!
! Load the right variogram as the first one if performing median IK:
!
      if(mik.eq.1) then
            icut = 1
            clos = abs(cutmik-thres(1))
            do ic=2,ncut
                  test1 = abs(cutmik-thres(ic))
                  if(test1.lt.clos) then
                        icut = ic
                        clos = test1
                  end if
            end do
            c0(1)   = c0(icut)
            nst(1)  = nst(icut)
            istart1 = 1
            istarti = 1 + (icut-1)*MAXNST
            do ist=1,nst(1)
                  index1        = istart1 + ist - 1
                  indexi        = istarti + ist - 1
                  it(index1)    = it(indexi)
                  aa(index1)    = aa(indexi)
                  cc(index1)    = cc(indexi)
                  ang1(index1)  = ang1(indexi)
                  ang2(index1)  = ang2(indexi)
                  ang3(index1)  = ang3(indexi)
                  anis1(index1) = anis1(indexi)
                  anis2(index1) = anis2(indexi)
            end do
      end if
      elapsed=secnds(t1)
      write(*,'(a,f10.4,a)')'[FOR](sgs) ',elapsed, &
      's - time to read pars'
!
! Open the output file and return:
!
      open(lout,file=outfl,status='UNKNOWN')
      write(lout,104) title
 104  format(a80)
      write(lout,105) 1,nx,ny,nz
 105  format(4(1x,i4))
      write(lout,106)
 106  format('Simulated Value')
      return
!
! Error in an Input File Somewhere:
!
 96   stop 'ERROR in soft data file!'
 97   stop 'ERROR in table look up file!'
 98   stop 'ERROR in parameter file!'
 99   stop 'ERROR in data file!'
      end



      subroutine sisim
      !subroutine sisim(MAXNOD,MAXCXY,MAXCTX,MAXCTY,MAXCTZ,MAXKR1,&
      !                  MAXSBX,MAXSBY,MAXSBZ,MAXCUT,MXCUT,MAXROT)
!-----------------------------------------------------------------------
!
!           Conditional Simulation of a 3-D Rectangular Grid
!           ************************************************
!
! This subroutine generates 3-D conditional simulations of a continuous
! variable with sequential indicator simulation.
!
!
!
! PROGRAM NOTES:
!
!  1. The three dimensional anisotropy parameters of the search ellipse
!     and the variogram ranges are described in section 2.3 of the
!     manual.   The variogram parameters are described in the same place
!
!  2. The conditioning data and previously simulated grid nodes can be
!     searched separately.  There can be a different maximum number of 
!     each and a minimum number of conditioning data can be specified 
!     to restrict simulation beyond the limits of the data.  The 
!     closeness of previously simulated grid nodes is measured according
!     to the variogram structural distance.
!
!  
!
!
!
! Based on the 1990 version of IK3D and the SIS program
!
!-----------------------------------------------------------------------
      use geostat
      !  Editado por RGE
      use       global
      !use       dfport
      use       srch
      use       graph_vars
      use       kdtree2_module
      use        random2
! -----------------------     
      implicit none
      !include  'sisim.inc'
!      include  'sisim_new.inc'


      real      ntviol,atviol,big(20),line(20),cix(20),ciy(20),ciz(20),&
                cdif(20),cccdf(20),cccdfo(20),cloc(20) 
      real*8    acorni
      integer   query
      



      integer :: NODES_LENGTH 
      integer :: GRID_OUT_LENGTH 
      integer :: NODES2CAL_LENGTH 
      integer, allocatable :: cur_edge_node_array1(:)
      integer, allocatable :: cur_edge_node_array2(:)
      real*8, allocatable :: edge_dist_array(:)
      integer, allocatable :: nodes2cal_array(:)

!levels code
      integer :: numberOfLevels,maxLevel,levelThreshold,countLev,lastCount,lev
      !integer, allocatable :: level(:), ncloseIndex(:), resultsIdxIndex(:,:),indexSort(:),levelCount(:),levelStart(:)
      integer, allocatable :: level(:), ncloseIndex(:), indexSort(:),levelCount(:),levelStart(:)
      integer, allocatable, target :: lock(:)
      integer :: plock
      real*8, allocatable :: resultsDisIndex(:,:),xppIndex(:,:)
      integer :: levIni, levFin, levIniLocal, levFinLocal, ilock
      integer :: threadId,numThreads,blocknumber
      real :: invNumThreads
      INTEGER, EXTERNAL :: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS
!levels code end

!kdtree code begin
      integer, allocatable :: useablemap(:)
      real(kdkind), dimension(:,:), allocatable :: my_array
      real(kdkind), allocatable :: query_vec(:)
      type(kdtree2), pointer :: tree  ! this is how you declare a tree in your main program
      type(kdtree2_result),allocatable :: results(:), resultsb(:)
      real(kdkind) :: rv 


      integer nnx, nny,nnz,nnclose,nMAXORD

      integer MAXNSTvar,NCLOSEvar

      real test1,test2
      integer no_cpp,cnt,i,j,k,ind,indi,indj,indk,mdt,nxdis,nydis,nzdis
      integer ndb,nk,koption,nloop,irepo
      real xk,vk,xkmae,xkmse,ddh,TINY
      integer icut,ic,is,isrot,nclose,isim,imult,jz,jy,jx,ix,iy,iz,index,id,id2,testind
      real*8 rand_num
      real xx,yy,zz
      integer in,ncnode,nclosetmp,ierr,nctx,ncty,nctz,nxysim
      real cdfval
      real*8 zval
      real*4 t1,t2
      integer :: nloopchunksize,nloopmin,nloopmax
      integer :: blocksize,blockid,blocknum,nlast  

      koption=0

      MAXNSTvar = MAXNST
      allocate(results(ndmax))
!kdtree code end
      
      t1=secnds(0.0)


!      
! Edicion Rodrigo Gutierrez
!
!
! Se realizan configuraciones necesarias para el L-ISOMAP
! 
      write(*,*)'nodmax inicio sisim  = ',nodmax
      xyzland=xland*yland*zland 
      
      write(*,'(a,I10,a,I10,a)'),'valor xyzland =',xyzland
      
      allocate(landpts(xyzland))
         if(test.ne.0)then
                 write(*,*)'ERROR: Allocation failed due to'
                 stop
         end if
         
      allocate(evalues(xyzland))
         if(test.ne.0)then
                write(*,*)'ERROR: Allocation failed due to'
                stop
         end if
         
      allocate(vectors(xyzland,xyzland))
          if(test.ne.0)then
                write(*,*)'ERROR: Allocation failed due to'
                stop
         end if
      
      !write(*,*)'nx = ',nx,' xland = ',xland
      !write(*,*)'int(nx/xland) = ',int(nx/xland)
      xspace=int(real(nx/xland))
      yspace=int(real(ny/yland))
      zspace=int(real(nz/zland))
      
      write(*,*),'xpace, yspace,zspace = ',xspace,yspace,zspace
      no_cpp=0
      cnt=0
      
      do k=1,zland
         indk=zspace*k - 0.5*zspace + 1
         do j=1,yland
            indj=yspace*j - 0.5*yspace + 1
            do i=1,xland
               cnt=cnt+1 
               indi=xspace*i - 0.5*xspace + 1
               ind=get_node(indi,indj,indk,nx,ny,nz)
               landpts(cnt)=ind
            end do
         end do
      end do   
      
      !call set_graph(nx,ny,nz,xmn,ymn,zmn,xsiz,ysiz,zsiz,graph_offset) 
      !set up the graph for the input grid, calculates edge lengths
      write(*,*)nz,ny,nx
      call set_graph_onlymem(MDS_opt,NODES_LENGTH,GRID_OUT_LENGTH,cur_edge_node_array1,cur_edge_node_array2,edge_dist_array) 
      
    
      t2=secnds(t1)
      elapsed=t2
      write(*,'(a,f10.4,a)') '[FOR](sgs) ',elapsed, &
      's - time for calculate edge lengths'
      write(time_out2,'(f10.4,a)')elapsed, 's - time to cal edge legths'

      nsec = 2
      mdt = 1
      if(ktype.eq.3) mdt = mdt + 1
      if(ktype.eq.0) mdt = 0
      if(ktype.eq.2) mdt = 0

! In all cases the offsets are relative to the lower left corner.
! This is done for rescaling the drift terms in the kriging matrix.
!
      if(nxdis.lt.1) nxdis = 1
      if(nydis.lt.1) nydis = 1
      if(nzdis.lt.1) nzdis = 1
      ndb=1
      
      
      !xdis = xsiz  / max(real(nxdis),1.0)
      !ydis = ysiz  / max(real(nydis),1.0)
      !zdis = zsiz  / max(real(nzdis),1.0)
      !i    = 0
      !xloc = -0.5*(xsiz+xdis)
     
      ! do ix =1,nxdis
       !     xloc = xloc + xdis
       !     yloc = -0.5*(ysiz+ydis)
       !!     write(*,*) 'Aca vamos perdidos...2.5'
        !    do iy=1,nydis
        !!    write(*,*) 'Aca vamos perdidos...2.55'
         !         yloc = yloc + ydis
        !         zloc = -0.5*(zsiz+zdis)
         !         write(*,*) 'Aca vamos perdidos...2.6'
         !         do iz=1,nzdis
         !         write(*,*) 'Aca vamos perdidos...2.64'
         !               zloc = zloc + zdis
         !              write(*,*) 'Aca vamos perdidos...2.65'
         !               i = i+1
         !               write(*,*) 'Aca vamos perdidos...2.68'
         !               xdb(i) = xloc + 0.5*xsiz
         !               write(*,*) 'Aca vamos perdidos...2.7'
         !               ydb(i) = yloc + 0.5*ysiz
         !               zdb(i) = zloc + 0.5*zsiz
         !               write(*,*) 'Aca vamos perdidos...3'
         !         end do
         !   end do
      ! end do
! Initialize accumulators:
!
      nk    = 0
      xk    = 0.0
      vk    = 0.0
      xkmae = 0.0
      xkmse = 0.0
! Report on progress from time to time:
      !write(*,*)'nodmax inicio sisim 3  = ',nodmax
      if(koption.eq.0) then
            nxy   = nx*ny
            nxyz  = nx*ny*nz
            nloop = nxyz
            irepo = max(1,min((nxyz/10),10000))
      else
            nloop = nx*ny*nz
            irepo = max(1,min((nd/10),10000))
      end if
      ddh = 0.0      
      
      allocate(ex_results_d(nodmax))
      
      allocate(ex_sim_array(nloop+nd))
      ex_sim_array=-999
      
      ex_sim_array(1:nd) = data_ind(1:nd)
      
            
      ex_results_d = 0
      
     
      ! allocate(ex_results_d(nclose))
      ! allocate(ex_results(nclose))
      
      !  allocate(ex_sim_array(nloop+nd))
!       allocate(ex_sim_array_nodos(nloop+nd))
      
      ! n_searched=0    
      ! ex_sim_array=-999
      ! ex_sim_array_nodos=-999
      ! write(*,*) 'Aca vamos ...8'
       !write(*,*) 'data_ind',data_ind(4)
       
       !write(*,*) 'data_ind = ',data_ind
      ! 
      write(*,'(a,f10.4,a)') '[FOR](sgs) ',elapsed, &
      's - setup for ISOMAP'
      
      
      !write(*,*) 'Llegamos al calculo del MDS isomap'
          if(MDS_opt==2 .or. MDS_opt==3) then
               write(*,*)'G:',ndmax,NODES_LENGTH,GRID_OUT_LENGTH
                !call get_landmark_pts_onlymem(ndmax,nd,nx,ny,nz,NODES_LENGTH,GRID_OUT_LENGTH,cur_edge_node_array1,cur_edge_node_array2,edge_dist_array) 
                t1=secnds(0.0)
                call get_landmark_pts_onlymem(NODES_LENGTH,GRID_OUT_LENGTH,cur_edge_node_array1,cur_edge_node_array2,edge_dist_array) 
                t2=secnds(t1)
                elapsed=t2
      write(*,'(a,f10.4,a)') '[FOR](sgs) ',elapsed, &
      's - time to finish get_landamark_pts_onlymem (dijkstra)'
                !use c++ program to get distances to landmark poitns
                !call MDS_ISOMAP(ndmax,nd,nx,ny,nz)
                t1=secnds(0.0)
                call MDS_ISOMAP
                t2=secnds(t1)
                elapsed=t2
      write(*,'(a,f10.4,a)') '[FOR](sgs) ',elapsed, &
      's - time to finish MDS_ISOMAP'
!write(*,*)'TEST-1',coord_ISOMAP_trans(1,1:4)
                 !do the multidimensinoal scaling
!now have coord of all grid points in coord_ISOMAP(NODES,dim)
          end if
      !write(*,*) 'Terminamos al calculo del MDS isomap'    
      !write(*,*) 'coord_ISOMAP_trans(1,1572864,1)=',coord_ISOMAP_trans(1,1572864) 
      ! write(*,*) ' coord_ISOMAP() ',coord_ISOMAP
      t1=secnds(0.0)
      !elapsed=secnds(total)
      !write(*,'(a,f10.4,a)') '[FOR](sgs) ',elapsed, &
      !'s - time to finish multidimensional scaling'
      !write(*,'(f10.4,a)')      elapsed, &
      !      's - time to finish multidimensional scaling'
      !write(time_out2,'(f10.4,a)')      elapsed, &
      !      's - time to finish multidimensional scaling'

! fin edicion Rodrigo Gutierrez
!      
!
      
! kdtree code begin
!set up the kdtree for searching:
      if(d_tree>dim .or. d_tree<1) then
          d_tree=dim
          write(*,'(a,I5)') &
      '[FOR](sgs) actual number of dim to use in the search: ',d_tree
      end if
  
!set up the logical array for searching the tree:
      allocate(is_usable(nx*ny*nz))
      allocate(useablemap(nx*ny*nz))

      allocate(coord_ISOMAP_trans_tmp(d_tree,NODES))
      allocate(coord_ISOMAP_trans_rearranged(d_tree,NODES))

      !allocate(input_usable(nx*ny*nz))
      !allocate(input_data(d_tree,NODES) )

!now set all the data locations to usable:
      is_usable=.false.
      is_usable(data_ind(1:nd)) = .true.

      !tree => kdtree2_create(nx,ny,nz,NODES,sort=.true.,rearrange=.true.)  ! this is how you create a tree. 
      !tree => kdtree2_create(NODES,sort=.true.,rearrange=.true.)  ! this is how you create a tree. 
 
    
!Build a map array so that the is_usable array can be updated
!in the case where the kdtree rearranges the data
      !if(tree%rearrange)then
      !   do i=1, nx*ny*nz
      !      useablemap(tree%ind(i)) = i
      !   enddo
      !endif

      write(*,*) '[FOR](sgs) Preparing the simulation '

      use_kd_tree=.true.  ! if =.false. will start with an exhaustive search, and used kd when it is more efficient

! kdtree code end
!stop
     
!
! Set up the rotation/anisotropy matrices that are needed for the
! variogram and search:
!
      write(*,*) 'Setting up rotation matrices for variogram and search'
      do ic=1,ncut
      do is=1,nst(ic)
            ind = is + (ic-1)*MAXNST
            call setrot(ang1(ind),ang2(ind),ang3(ind),anis1(ind),&
                        anis2(ind),ind,MAXROT,rotmat)
      end do
      end do
      isrot = MAXNST*MAXCUT + 1
      call setrot(sang1,sang2,sang3,sanis1,sanis2,isrot,MAXROT,rotmat)
!
! Set up for super block searching:
!
      !if(sstrat.eq.0) then
      !      write(*,*) 'Setting up super block search strategy'
      !      do i=1,nd
      !            actloc(i) = real(i)
      !      end do
      !      nsec = 0
      !      call setsupr(nx,xmn,xsiz,ny,ymn,ysiz,nz,zmn,zsiz,nd,x,y,z,&
      !             actloc,tmp,nsec,sec1,sec2,sec3,MAXSBX,MAXSBY,MAXSBZ,&
      !             nisb,nxsup,xmnsup,xsizsup,nysup,ymnsup,ysizsup,nzsup,&
      !             zmnsup,zsizsup)
      !      call picksup(nxsup,xsizsup,nysup,ysizsup,nzsup,zsizsup,&
      !             isrot,MAXROT,rotmat,radsqd,nsbtosr,ixsbtosr,&
      !             iysbtosr,izsbtosr)
      !end if
!
! Set up the covariance table and the spiral search:
!
      !call ctable
      !call ctable(MAXNOD,MAXCXY,MAXCTX,MAXCTY,MAXCTZ,MAXROT)
!write(*,*)'covtab-1=',covtab(49,26,52,1)
!write(*,*)'covtab-2=',covtab(49,26,52,2)
!write(*,*)'covtab-3=',covtab(49,26,52,3)
!write(*,*)'covtab-4=',covtab(49,26,52,4)

!levels code

      t2=secnds(t1)
      elapsed=t2
      write(*,'(a,f10.4,a)')  '[FOR](sgs) ',elapsed, &
      's - time to setting arrays and variables for simulation'

      t1=secnds(0.0)

      allocate(level(nloop))
      allocate(ncloseIndex(nloop))
      nclose = min(nd,ndmax)
      allocate(resultsDisIndex(nclose,nloop))
      allocate(resultsIdxIndex(nclose,nloop))
      !allocate(xppIndex(nloop,nreal))
      allocate(indexSort(nloop))
      allocate(lock(nloop))
  
      do i=1,nloop
          level(i) = 0
      end do
      !level(data_ind(1:nd))=0
      numberOfLevels=0
      lock=0 
  
      lock(data_ind(1:nd)) = 1



!levels code end





!
! MAIN LOOP OVER ALL THE SIMULAUTIONS:
!
      write(*,*) '[FOR](sgs) Working on the simulation '
      !write(*,*)'llegamos a do de simulaciones'
      do isim=1,nsim
             !write(*,*)'Starting simulation ',isim
      
      !Arreglo post memoria simulaciones
             ex_sim_array = -999
             ex_sim_array(1:nd) = data_ind(1:nd)
             
      ! Fin arreglo
!
! Work out a random path for this realization:
!
            !write(*,*)'Start random path'
!            do ind=1,nxyz
!                  sim(ind)   = real(acorni(idum))
!                  order(ind) = ind
!write(*,*)'acorni:',ind,sim(ind)
!            end do

            call init_genrand(ixv(1))  !intial the random generator
            do i=1,10000
                rand_num = grnd()
            end do
            
            !use sim array as a temp array for the order
            do i=1,nxyz
                sim(i) = real(grnd())
                order(i) = i
!write(*,*)'acorni:',i,sim(i)
            end do

 



            !write(*,*)'Finished random path'
!
! The multiple grid search works with multiples of 4 (yes, that is
! somewhat arbitrary):
!
            if(mults.eq.1) then
                  do imult=1,nmult
                        nnz = max(1,nz/(imult*4))
                        nny = max(1,ny/(imult*4))
                        nnx = max(1,nx/(imult*4))
                        jz  = 1
                        jy  = 1
                        jx  = 1
                        do iz=1,nnz
                           if(nnz.gt.1) jz = iz*imult*4
                           do iy=1,nny
                              if(nny.gt.1) jy = iy*imult*4
                              do ix=1,nnx
                                 if(nnx.gt.1) jx = ix*imult*4
                                 index = jx + (jy-1)*nx + (jz-1)*nxy
                                 sim(index) = sim(index) - imult
                              end do
                           end do
                        end do
                  end do
            end if
            !write(*,*)'Start sortem'
            !call sortem(1,nxyz,sim,1,order,c,d,e,f,g,h)
            !call sortem2dp(1,nxyz,sim,1,order)
            call sortem2(1,nxyz,sim,1,order)
            !write(*,*)'Finished sortem'
!
! Initialize the simulation:
!
            do i=1,nxyz
                  sim(i) = UNEST
                  tmp(i) = 0.0
            end do
            write(*,*)
            write(*,*) ' Working on realization number: ',isim
            
!
! Assign the data to the closest grid node:
!
!write(*,*)'check:',nd
!            TINY = 0.0001
            TINY = 5.1
            do id=1,nd
!write(*,*)'datx0:',id,x(id),y(id),z(id)
!write(*,*)'datx1:',id,nx,ny,nz
!write(*,*)'datx2:',id,xmn,ymn,zmn
!write(*,*)'datx3:',id,xsiz,ysiz,zsiz
                  call getindx(nx,xmn,xsiz,dble(x(id)),ix,testind)
                  call getindx(ny,ymn,ysiz,dble(y(id)),iy,testind)
                  call getindx(nz,zmn,zsiz,dble(z(id)),iz,testind)
                  ind  = ix + (iy-1)*nx + (iz-1)*nxy
!write(*,*)'datt:',id,ind,ix,iy,iz,nx,nxy
                  xx   = xmn + real(ix-1)*xsiz
                  yy   = ymn + real(iy-1)*ysiz
                  zz   = zmn + real(iz-1)*zsiz
                  test1 = abs(xx-x(id)) + abs(yy-y(id)) + abs(zz-z(id))
!
! Assign this data to the node (unless there is a closer data):
!
                  atnode(id) = .false.
                  if(sstrat.eq.1)                  atnode(id) = .true.
                  if(sstrat.eq.0.and.test1.le.TINY) atnode(id) = .true.
!write(*,*)'dat0:',id,sstrat,test1,TINY,atnode(id)
                  if(atnode(id)) then
                        if(sim(ind).ge.0.0) then
                              id2 = int(sim(ind)+0.5)
                              test2 = abs(xx-x(id2)) + abs(yy-y(id2)) + abs(zz-z(id2))
!write(*,*)'dat1:',id,id2,test1,test2
                              if(test1.le.test2) then
                                  sim(ind) = real(id)
!write(*,*)'dat2:',id,ind,sim(ind)
                              end if
                              if(idbg.ge.2) write(ldbg,102) id,id2
                        else
                              sim(ind) = real(id)
!write(*,*)'dat3:',id,ind,sim(ind)
                        end if
                  end if
            end do
 102        format(' WARNING data values ',2i5,' are both assigned to',&
                 /,'         the same node - taking the closest')
!
! Now, enter the hard data values into the "sim" array and keep the
! data number in the "tmp" array (to be reset when a hard value
! is assigned to that node):
            !write(*,*) ' Enter hard data to sim array: ',isim
            do i=1,nxyz
                  id = int(sim(i)+0.5)
                  if(id.gt.0) then
                        if(id.le.nhd) then
                                sim(i) = vr(id,MXCUT)
!write(*,*)'data:',i,sim(i)
                        else
                              tmp(i) = sim(i)
                              sim(i) = UNEST
                        end if
                  end if
            end do
!
! Accumulate the number and magnitude of order relations violations:
!
            !write(*,*) ' Accumulate the number and magnitude of order: ',isim
            nclose = 0
            !irepo  = max(1,min((nxyz/10),10000))
            irepo = 1000
            ntviol = 0.0
            atviol = 0.0
            do icut=1,ncut
                  nviol(icut) =  0
                  aviol(icut) =  0.0
                  xviol(icut) = -1.0
            end do
!
! MAIN LOOP OVER ALL THE NODES:
!
query=-100 
            j = 0
            !write(*,*) ' Main loop over all nodes: ',isim



!!$omp parallel default(shared) private(threadId,numThreads,i,in,id,index,iz,iy,ix,xx,yy,zz,icnode,cnodeid,ncnode,nclose,nclosetmp,results,cdfval,tree,input_data,input_usable,nloopchunksize,nloopmin,nloopmax,maxLevel,numberOfLevels)
!$omp parallel default(shared) private(threadId,numThreads,i,in,id,index,iz,iy,ix,xx,yy,zz,icnode,cnodeid,ncnode,nclose,nclosetmp,results,cdfval,tree,nloopchunksize,nloopmin,nloopmax,maxLevel,numberOfLevels,blocksize,blockid,blocknum,nlast)


      threadId=omp_get_thread_num()+1
      numThreads = omp_get_num_threads()
      if(threadId.eq.1) write(*,*) &
      '[FOR](sgs) Starting neighbours calculation numThreads=',&
      numThreads
  
      tree => kdtree2_create_v2(nx,ny,nz,NODES,sort=.true.,rearrange=.true.)  ! this is how you create a tree. 
      !tree => kdtree2_create_v2(nx,ny,nz,NODES,sort=.true.,rearrange=.false.)  ! this is how you create a tree. 
      !write(*,*)threadId,'tree%dimen',tree%dimen
      if(threadId.eq.1)write(*,*) '[FOR](sgs) Trees computed '
      if(tree%rearrange)then
         do i=1, nx*ny*nz
             useablemap(tree%ind(i)) = i
         enddo
      endif
   
      blocksize=250
   
      blocknum=int(ceiling(real(nloop)/real(blocksize)))
      write(*,*)'blocknum=',blocknum,'nloop=',nloop,&
      'blocksize=',blocksize
      nlast=1
  
      do blockid=1,blocknum
  
      if(mod(blockid,numThreads).eq.(threadId-1))then
   
      !nloopchunksize=nxyz/numThreads
      !nloopmin=(threadId-1)*nloopchunksize + 1
      !nloopmax=min(threadId*nloopchunksize,nxyz)
      nloopmin=(blockid-1)*blocksize + 1
      nloopmax=min(blockid*blocksize,nloop)
  
  
      !write(*,*)'[FOR](sgs) threadId=',threadId,&
      !' processing from ',nloopmin,' to ',nloopmax


      !if(threadId>1)then
          !do i=1,threadId-1
              !write(*,*)'[FOR](sgs) threadId=',threadId,&
      !' marking previous nodes 1 to ',nloopmin-1
              do in=nlast,nloopmin-1
                  index = order(in)
                  tree%REARRANGED_IS_USABLE(useablemap(index)) = .true.  !need to tell the tree that this point is now simulated:
                  !do ireal=1,nreal
                  !   pp = grnd()
                  !   call gauinv(pp,xpp,ierr)
                  !   xppIndex(index,ireal) = xpp
                  !   !if(sim(index,1).eq.-999.0) sim(index,ireal) = -998.0
                  !end do
              end do
              nlast=nloopmax+1
              !write(*,*)'[FOR](sgs) threadId=',threadId,' marking done'
          !end do
      !end if

      !do index2=1,nloop
      !if(threadId.eq.1)write(*,*) '[FOR](sgs) Running tree search '
      do in=nloopmin,nloopmax

          if((int(in/irepo)*irepo).eq.in) write(*,104) in
 104      format('[FOR](sgs)   currently on node ',i9)
          index = int(order(in)+0.5)
!
! Do we really need to simulate this grid node location?
!
!                  if(sim(index).ne.UNEST) go to 20

!
! Location of the node we are currently working on:
!
        iz = int((index-1)/nxy) + 1
        iy = int((index-(iz-1)*nxy-1)/nx) + 1
        ix = index - (iz-1)*nxy - (iy-1)*nx
        xx = xmn + real(ix-1)*xsiz
        yy = ymn + real(iy-1)*ysiz
        zz = zmn + real(iz-1)*zsiz


        if(sim(index).ne.UNEST)then
              level(index) = 0
              ncloseIndex(index)=0
              tmp(index) = 0.0
        elseif(imbsim.eq.0.and.tmp(index).ne.0.0) then
              id = int(tmp(index)+0.5)
              sim(index) = vr(id,MXCUT)
              level(index) = 0
              ncloseIndex(index)=0
              tmp(index) = 0.0
        else
!
! Now, we'll simulate the point ix,iy,iz.  First, get the close data
! and make sure that there are enough to actually simulate a value,
! we'll only keep the closest "ndmax" data, and look for previously

        icnode(:)=0
        cnodeid(:)=-1
        ncnode=0
! kdtreeigh begin
        nclose = min(nd,ndmax)
        nclosetmp = nclose
        results(:).idx=0
        results(:).dis=99999999.0
        call kdtree2_n_nearest_around_point(real(coord_ISOMAP_trans(1:d_tree,index)),d_tree,tp=tree,idxin=-1,correltime=-1,nn=nclose, results=results)
        tree%REARRANGED_IS_USABLE(useablemap(index)) = .true.  !need to tell the tree that this point is now simulated:
        do i=1,nclosetmp
           if(results(i).dis>radsqd) then
              nclose=nclose-1
           end if
        end do
        ncnode=nclose
        ncnodeIndex(index)=ncnode
            
        !if(threadId.eq.1)then 
        !maxLevel=-1
        !do i=1,nclosetmp
        !   if(results(i).dis<=radsqd) then
        !           if(level(results(i).idx).gt.maxLevel) then
        !                   maxLevel = level(results(i).idx)
        !           end if
        !   end if
        !end do
        !end if


        if(nclose.eq.0)then
                !if(threadId.eq.1)level(index)=0
                ncloseIndex(index)=0
        else
                !if(threadId.eq.1)level(index)=maxLevel+1
                ncloseIndex(index)=nclose
                do i=1,nclosetmp
                        resultsDisIndex(i,index)=results(i).dis
                end do
                do i=1,nclosetmp
                        resultsIdxIndex(i,index)=results(i).idx
                end do
        end if
        !if (maxLevel .gt. numberOfLevels) then
        !        numberOfLevels = maxLevel
        !end if
        sim(index) = 999999.0

        cdfval = real(grnd())
        cdfvalIndex(index)=cdfval

        tmp(index) = 0.0

        end if

!
! END MAIN LOOP OVER NODES:
!
! 20         continue
!
!            tmp(index) = 0.0

      end do

      end if

      end do

      write(*,*)'[FOR](sgs) threadID=',threadId,' finished search'
!$omp barrier


      !if(threadId.eq.1.and.numThreads.gt.1)then
      if(threadId.eq.1)then

          elapsed=secnds(0.0)

          write(*,*)'[FOR](sgs) threadID=1 starting level computation'
          !do in=nloopmax+1,nxyz
          do in=1,nxyz
              index = order(in)
              if(sim(index).eq. 999999.0)then
                  nclosetmp = ncloseIndex(index) 
                  maxLevel=-1
                  
                  do i=1,nclosetmp
                     if(resultsDisIndex(i,index)<=radsqd) then
                             if(level(resultsIdxIndex(i,index))&
                            .gt.maxLevel) then
                                     maxLevel = level(resultsIdxIndex(i,index))
                             end if
                     end if
                  end do
                  if(ncloseIndex(index).eq.0)then
                      level(index)=0
                  else
                      level(index)=maxLevel+1
                  end if
                  !if (maxLevel .gt. numberOfLevels) then
                  !    numberOfLevels = maxLevel
                  !end if
  
              end if
          end do
          t2=secnds(elapsed)
          write(*,*)'[FOR](sgs) threadID=1 finished level computation'
  
      end if
!$omp end parallel

            t2=secnds(t1)
            elapsed=t2
            write(*,'(a,f10.4,a)')  '[FOR](sgs) ',elapsed, &
      's - time to calculate neighbours and levels for simulation'

!stop
            t1=secnds(0.0)

            numberOfLevels=-1
            do i=1,nxyz
                index=order(i)
                !write(*,*)'leveltag=',index,level(index),numberOfLevels
                if (level(index) .gt. numberOfLevels) then
                    numberOfLevels = level(index)
                end if
            end do


!            call system_clock(clock_stop,clock_rate,clock_max)
!            print *,'TIME (make levels1)=',
!     +    real(clock_stop-clock_start,kind=8)/real(clock_rate,kind=8) 

!            call system_clock(clock_start,clock_rate,clock_max)

            levelThreshold = 0
            numberOfLevels = numberOfLevels + 1
            countLev = 0
            lastCount = 0
            allocate(levelCount(numberOfLevels+1))
            allocate(levelStart(numberOfLevels+1))

            do lev = 0,numberOfLevels
               levelStart(lev+1) = countLev + 1
               levelCount(lev+1) = 0 
               do in = 1,nxyz
                  if(level(in).eq.lev)then
                     countLev = countLev +1
                     indexSort(countLev) = in
                     levelCount(lev+1) = levelCount(lev+1) + 1
                  end if
               end do
               levelThreshold = levelThreshold + levelCount(lev+1)

            end do

            levelThreshold = int(0.1*ceiling(real(levelThreshold)/ &
                                real(numberOfLevels-1)))


!            call system_clock(clock_stop,clock_rate,clock_max)
!            print *,'TIME (make levels2)=',
!     +    real(clock_stop-clock_start,kind=8)/real(clock_rate,kind=8) 


!            call system_clock(clock_start,clock_rate,clock_max)

!            do lev=1,numberOfLevels
!               print *,lev,levelCount(lev+1)
!            end do
            !stop

            
!            do i=1,20
!                  cloc(i)=(ciy(i)-1)*nx+cix(i)
!                  write(*,2001) int(cloc(i)),cdif(i),cccdf(i),cccdfo(i)
!            end do
!2001        format(i9,3(f9.3))
!!
!! Write this simulation to the output file:
!!
!            nxysim = 0
!            do ic=1,ncut
!                  ccdf(ic) = 0.0
!            end do
!            do ind=1,nxyz
!                  write(lout,'(f12.4)') sim(ind)
!!
!! Calculate the cdf of the simulated values (for error checking):
!!
!                  if(sim(ind).gt.UNEST) then
!                        nxysim = nxysim + 1
!                        do ic=1,ncut
!                              if(ivtype.eq.0) then
!                                    if(sim(ind).eq.thres(ic)) &
!                                      ccdf(ic)=ccdf(ic)+1.
!                              else
!                                    if(sim(ind).le.thres(ic)) &
!                                      ccdf(ic)=ccdf(ic)+1.
!                              end if
!                        end do
!                  endif
!            end do
!
! Report on the reproduction of the cdf and the number and magnitude
! of order relations violations:
! 
!            write(*,203)    isim
!            write(ldbg,203) isim
!            do icut=1,ncut
!                  ccdf(icut) = ccdf(icut) / max(real(nxysim),1.0)
!                  write(*,204)    icut,cdf(icut),ccdf(icut)
!                  write(ldbg,204) icut,cdf(icut),ccdf(icut)
!            end do
! 203        format(/,' Finished simulation ',i2)
! 204        format('     threshold ',i3,' input cdf = ',f6.4,&
!                       ' realization cdf = ',f6.4)
!            write(*,   300)
!            write(ldbg,300)
! 300        format(/,' Summary of order relations: ')
!            ntot = 0
!            atot = 0.0
!            do icut=1,ncut
!               ntot = ntot + nviol(icut)
!               atot = atot + aviol(icut)
!               aviol(icut) = aviol(icut) / real(max(1,nviol(icut)))
!               write(*,302) icut,nviol(icut),aviol(icut),xviol(icut)
!               write(ldbg,302) icut,nviol(icut),aviol(icut),xviol(icut)
! 302           format('     threshold',i2,' number = ',i6,&
!                      ' average = ',f8.4,' maximum = ',f8.4)
!            end do
!            atot = atot / real(max(1,ntot))
!            btot =(ntot / real(ncut*nxysim)) * 100.0
!            write(ldbg,303) btot,atot
!            write(*,   303) btot,atot
! 303        format(/,' total of ',f18.6,'% with average of ',f8.4)

!
! Initialize the simulation (again):
!

      !Arreglo post memoria simulaciones
             ex_sim_array = -999             
      ! Fin arreglo

            do i=1,nxyz
                  sim(i) = UNEST
                  tmp(i) = 0.0
                  lock(i) = 0
            end do
            lock(data_ind(1:nd)) = 1

            !write(*,*)
            !write(*,*) ' Working on realization number: ',isim
            
!
! Assign the data to the closest grid node:
!
!            TINY = 0.0001
            TINY = 5.1
            do id=1,nd
                  call getindx(nx,xmn,xsiz,dble(x(id)),ix,testind)
                  call getindx(ny,ymn,ysiz,dble(y(id)),iy,testind)
                  call getindx(nz,zmn,zsiz,dble(z(id)),iz,testind)
                  ind  = ix + (iy-1)*nx + (iz-1)*nxy
                  xx   = xmn + real(ix-1)*xsiz
                  yy   = ymn + real(iy-1)*ysiz
                  zz   = zmn + real(iz-1)*zsiz
                  test1 = abs(xx-x(id)) + abs(yy-y(id)) + abs(zz-z(id))
!
! Assign this data to the node (unless there is a closer data):
!
                  atnode(id) = .false.
                  if(sstrat.eq.1)                  atnode(id) = .true.
                  if(sstrat.eq.0.and.test1.le.TINY) atnode(id) = .true.
                  if(atnode(id)) then
                        if(sim(ind).ge.0.0) then
                              id2 = int(sim(ind)+0.5)
                              test2 = abs(xx-x(id2)) + abs(yy-y(id2)) + abs(zz-z(id2))
                              if(test1.le.test2) sim(ind) = real(id)
                              !if(idbg.ge.2) write(ldbg,102) id,id2
                        else
                              sim(ind) = real(id)
                        end if
                  end if
            end do
! 102        format(' WARNING data values ',2i5,' are both assigned to',&
!                 /,'         the same node - taking the closest')
!
! Now, enter the hard data values into the "sim" array and keep the
! data number in the "tmp" array (to be reset when a hard value
! is assigned to that node):
            !write(*,*) ' Enter hard data to sim array: ',isim
            do i=1,nxyz
                  id = int(sim(i)+0.5)
                  if(id.gt.0) then
                        if(id.le.nhd) then
                              sim(i) = vr(id,MXCUT)
                        else
                              tmp(i) = sim(i)
                              sim(i) = UNEST
                        end if
                  end if
            end do
!
! Accumulate the number and magnitude of order relations violations:
!
            !write(*,*) ' Accumulate the number and magnitude of order: ',isim
            nclose = 0
            irepo  = max(1,min((nxyz/10),10000))
            ntviol = 0.0
            atviol = 0.0
            do icut=1,ncut
                  nviol(icut) =  0
                  aviol(icut) =  0.0
                  xviol(icut) = -1.0
            end do
!
! MAIN LOOP OVER ALL THE NODES:
!

            lev=0

            levIni=levelStart(lev+1) 
            levFin=(levelStart(lev+1)+levelCount(lev+1)-1) 
            write(*,*) '[FOR](sgs) Level=',lev,&
            ' Number_per_level=',(levelCount(1))
!            write(*)'level-0:',levelCount(lev+1)
!            write(*)'levelStart-0:',levelStart(lev+1)
            do countLev=levIni,levFin
                  index = indexSort(countLev)
                  !write(*,*)'sim-level-0:',index,sim(index)
                  lock(index)=1
            end do

!      write(*,*) '[FOR](sgs) Level=',lev,' Number_per_level=',(levelCount(1))
!            print *,'level=',lev,'/',numberOfLevels,'ini=',levIni,'fin=',levFin
!
!!!$omp parallel default(firstprivate) shared(x,y,z,vr,level,indexSort,cdfvalIndex,ncnodeIndex,icnodeIndex,cnodeidIndex,cnodexIndex,cnodeyIndex,cnodezIndex,cnodevIndex,cnodetIndex,tmp,sim,order,covtab,beez,lock)
!!            threadId=omp_get_thread_num()+1
!!            numThreads = omp_get_num_threads()
!!            invNumThreads = 1.0/real(numThreads)
!
!!#ifdef _OPENMP
!!            threadId=OMP_get_thread_num()+1
!!#else
!!            threadId=1
!!#endif
!
!!!$omp do schedule(runtime)
!            do count=levIni,levFin
!
!                  index = indexSort(count)
!
!!      !do de RGE
!!            j = 0
!!      !DO DE RGE      
!!            !write(*,*) ' Main loop over all nodes: ',isim
!!            do in=1,nxyz
!!            
!!                  
!!                  if((int(in/irepo)*irepo).eq.in) write(*,104) in
!!! 104              format('   currently on node ',i9)
!!                  index = int(order(in)+0.5)
!                 
!!
!! Do we really need to simulate this grid node location?
!!
!                  if(sim(index).ne.UNEST)then 
!                        lock(index) = 1
!                        go to 200
!                  end if
!                  if(imbsim.eq.0.and.tmp(index).ne.0.0) then
!                        id = int(tmp(index)+0.5)
!                        sim(index) = vr(id,MXCUT)
!                        lock(index) = 1
!                        go to 200
!                  end if
!!
!! Location of the node we are currently working on:
!!
!                  iz = int((index-1)/nxy) + 1
!                  iy = int((index-(iz-1)*nxy-1)/nx) + 1
!                  ix = index - (iz-1)*nxy - (iy-1)*nx
!                  xx = xmn + real(ix-1)*xsiz
!                  yy = ymn + real(iy-1)*ysiz
!                  zz = zmn + real(iz-1)*zsiz
!                  !if(idbg.ge.3)&
!                  !write(ldbg,*) 'Working on grid index:',ix,iy,iz
!                  
!!
!! Now, we'll simulate the point ix,iy,iz.  First, get the close data
!! and make sure that there are enough to actually simulate a value,
!! we'll only keep the closest "ndmax" data, and look for previously
!! simulated grid nodes:
!!
!                  if(sstrat.eq.0) then
!                        !write(*,*) ' srchsupr: ',isim,in
!                        call srchsupr(xx,yy,zz,radsqd,isrot,MAXROT,&
!                             rotmat,nsbtosr,ixsbtosr,iysbtosr,&
!                             izsbtosr,noct,nd,x,y,z,tmpdat,nisb,nxsup,&
!                             xmnsup,xsizsup,nysup,ymnsup,ysizsup,&
!                             nzsup,zmnsup,zsizsup,nclose,close,&
!                             infoct)
!                  !if(in.eq.47362.or.in.eq.72755.or.in.eq.74061.or.in.eq.48441.or.in.eq.77807) then
!                  !      write(*,*) in,(close(i),i=1,nclose)
!                  !      stop
!                  !end if
!                        if(nclose.gt.ndmax) nclose = ndmax
!!                       do i=1,nclose
!!                             iii = int(close(i)+0.5)
!!                             close(i) = real(actloc(iii))
!!                       end do
!                  endif
!
!                  icnode(:)=0
!                  cnodex(:)=0.0
!                  cnodey(:)=0.0
!                  cnodez(:)=0.0
!                  cnodev(:)=0.0
!                  cnodet(:)=0.0
!                  cnodeid(:)=0
!                  ncnode=0
!
! 
!                  !j = j+1
!                  
!                 
!                  !write(*,*) 'coord_ISOMAP_trans(1,1572864,1)=',coord_ISOMAP_trans(1,1572864) 
!                  !write(*,*) ' srchnd: ',isim,in
!                  call srchnd(nloop,index,ix,iy,iz)
!!                  call srchnd(ix,iy,iz,MAXNOD,MAXXYZ,MAXORD,MXYZ,ncnode,
!!     +        maxsec,nctx,ncty,nctz,nlooku,nodmax,nx,ny,nz,UNEST,
!!     +        xmn,ymn,zmn,xsiz,ysiz,zsiz,icnode,ixnode,iynode,iznode,
!!     +        cnodex,cnodey,cnodez,cnodev,cnodet,tmp,sim)
!                  
!                  !write(*,*) ' ex_sim_array: ',isim,in
!                  ex_sim_array(nd+j) = index
!                  
!                  
!                  !if(idbg.ge.3)&
!                  !write(ldbg,*) '  there are ',nclose,' close data', &
!                  !                     ' and ',ncnode,' close nodes'
!!
!! What cdf value are we looking for?
!!
!                  !write(*,*) ' acorni: ',isim,in
!                  zval   = UNEST
!                  cdfval = cdfvalIndex(index)
!!
!! Use the global distribution?
!!
!                  if((nclose+ncnode).le.0) then
!                        !write(*,*) ' beyond: ',isim,in
!                        call beyond(ivtype,ncut,thres,cdf,ng,gcut,gcdf,&
!                                    zmin,zmax,ltail,ltpar,middle,mpar,&
!                                    utail,utpar,zval,cdfval,ierr)
!                  endif
!                  !write(*,*) ' save sim: ',isim,in
!                  sim(index) = zval
!                  lock(index) = 1
!!
!! END MAIN LOOP OVER NODES:
!!
! 200        continue
!            tmp(index) = 0.0
!            
!            end do
!!!$omp end do 
!!!$omp end parallel
!
!
!write(*,*)'MAXDAT-bef-t',MAXDAT

!            invNumThreads = 1.0/real(numThreads)

!
!$omp parallel default(firstprivate) shared(cdf,x,y,z,vr,level,indexSort,cdfvalIndex,ncnodeIndex,ncloseIndex,icnodeIndex,cnodeidIndex,tmp,sim,order,beez,lock,numberOfLevels,levelStart,levelCount,coord_ISOMAP,coord_ISOMAP_trans,resultsDisIndex,resultsIdxIndex,MXYZ,MAXCTX,MAXCTY,MAXCTZ,MAXKR1,MAXROT,MAXDAT,MXCUT,MAXKR2,MAXNOD,MAXCUT,MAXXYZ,MAXNSTvar,nd,nloop,xmn,ymn,zmn,xsiz,ysiz,zsiz,dim,nst,c0,it,cc,aa,atnode,nctx,ncty,nctz,icnode,ixnode,iynode,iznode,nx,ny,nz,NODES_LENGTH,NODES2CAL_LENGTH,d_tree,xyzland,query) private(maxLevel, results, threadId, numThreads, invNumThreads,lev,levIni, levFin,blocknumber, levIniLocal, levFinLocal,countLev,ix,iy,iz,xx,yy,zz,ncnode,nclose, cnodex, cnodey, cnodez,cnodet,cnodev,neighbours,distance,i,j,zval,cdfval,ccdf,ccdfo,ic,aclose,a,r,rr,s,ilock,plock)
!
            threadId=omp_get_thread_num()+1
            numThreads = omp_get_num_threads()
            !threadId=1
            !numThreads =1 
            invNumThreads = 1.0/real(numThreads)


            
            do lev=1,(numberOfLevels)
            !do lev=1,98


            levIni=levelStart(lev+1) 
            levFin=(levelStart(lev+1)+levelCount(lev+1)-1) 
!            write(*,*) '[OMP] lev=',lev,levIni,levFin
            if(threadId.eq.1) write(*,*) '[FOR][OMP](sgs) Level=',lev,&
            ' Number_per_level=',(levelCount(lev+1))

            if(levelCount(lev+1).ge.-1) then

            if(numThreads.gt.1)then
            blocknumber =ceiling(real(levelCount(lev+1)-numThreads+1)* &
                              invNumThreads)
               levIniLocal = levIni+ &
                         blocknumber*(threadId-1)
               levFinLocal = levIni+ &
                         blocknumber*threadId-1
               if(threadId.eq.numThreads) levFinLocal=levFin

            else
               levIniLocal=levIni
               levFinLocal=levFin
            end if

            do countLev=levIniLocal,levFinLocal
            
                  index = indexSort(countLev)
!write(*,*)lev,'count',countLev,'index',index
!
! Location of the node we are currently working on:
!
                  iz = int((index-1)/nxy) + 1
                  iy = int((index-(iz-1)*nxy-1)/nx) + 1
                  ix = index - (iz-1)*nxy - (iy-1)*nx
                  xx = xmn + real(ix-1)*xsiz
                  yy = ymn + real(iy-1)*ysiz
                  zz = zmn + real(iz-1)*zsiz
!                  if(idbg.ge.3)
!     +            write(ldbg,*) 'Working on grid index:',ix,iy,iz

!
! Now, we'll simulate the point ix,iy,iz.  First, get the close data
! and make sure that there are enough to actually simulate a value,
! we'll only keep the closest "ndmax" data, and look for previously
! simulated grid nodes:
!
!                  if(sstrat.eq.0) then
!                        call srchsupr(xx,yy,zz,radsqd,isrot,MAXROT, &
!                            rotmat,nsbtosr,ixsbtosr,iysbtosr, &
!                            izsbtosr,noct,nd,x,y,z,tmpdat,nisb,nxsup, &
!                            xmnsup,xsizsup,nysup,ymnsup,ysizsup, &
!                            nzsup,zmnsup,zsizsup,nclose,close, &
!                            infoct)
!                        if(nclose.gt.ndmax) nclose = ndmax
!!                       do i=1,nclose
!!                             iii = int(close(i)+0.5)
!!                             close(i) = real(actloc(iii))
!!                       end do
!                  endif

! Instead of calling srchnd, we extract neighbour values from 
! arrays ***Index(:), obtained before...

                 ncnode = ncnodeIndex(index)
                 nclose = ncloseIndex(index)
                 !icnode(:) = icnodeIndex(:,index)
                 !cnodeid(:) = cnodeidIndex(:,index)
!if(index.eq.90)then
!write(*,*)'pre:',index,level(index),ncnode,nclose
!end if             

                  !call srchndPopLocal(nloop,index,ix,iy,iz)
                  !nnx=nx
                  !nny=ny
                  !nnz=nz
                  !nnclose=nclose
                  !nMAXORD=MAXORD

                  call srchndPopLocal(nloop,index,ix,iy,iz,nx,ny,&
                  ncnode,nd,nloop,MAXNOD,MAXORD,resultsIdxIndex(:,index), &
                  xmn,ymn,zmn,xsiz,ysiz,zsiz,tmp,cnodex,cnodey,  &
                  cnodez,cnodet)

!              call srchndPop(ix,iy,iz,MAXNOD,MAXXYZ,MAXORD,MXYZ,ncnode,
!     +        maxsec,nctx,ncty,nctz,nlooku,nodmax,nx,ny,nz,UNEST,
!     +        xmn,ymn,zmn,xsiz,ysiz,zsiz,icnode,ixnode,iynode,iznode,
!     +        cnodex,cnodey,cnodez,cnodev,cnodet,tmp,sim,cnodeid)


! The values of cnodev are extracted from sim array using the
! indexes obtained before


            ilock=0
            do while(ilock.eq.0)
               ilock=1
               do i=1,ncnode
                  !ilock=ilock*lock(cnodeid(i))
                  ilock=ilock*lock(resultsIdxIndex(i,index))
               end do
            end do



            !ilock=1
            !plock = lock(resultsIdxIndex(ilock,index)) 
            !do while(ilock.le.nclose)
            !   j=0
            !   do i=1,10000
            !      j=j+2
            !   end do
            !   if(plock.ne.0) then 
            !      ilock=ilock+1
            !      plock = lock(resultsIdxIndex(ilock,index)) 
            !   end if
            !end do

if(index.eq.query)write(*,*)'sim',sim(1140)


            cnodev(:)=0.0
            neighbours(:)=0
            distance(:)=0.0
            do i=1,ncnode
               cnodev(i)=sim(resultsIdxIndex(i,index))
            end do
if(index.eq.query)write(*,*)'resultsIdxIndex',resultsIdxIndex(:,index) 
if(index.eq.query)write(*,*)'cnodev',cnodev(:)
            neighbours=0
            distance=0.0
            do i=1,ncnode
               neighbours(i)=resultsIdxIndex(i,index)
               distance(i)=dble(resultsDisIndex(i,index))
            end do

if(index.eq.query)write(*,*)'cdf',cdf(:)



!
! What cdf value are we looking for?
!
                  zval   = UNEST
                  cdfval = cdfvalIndex(index)
!
! Use the global distribution?
!

!if(index.eq.1157)then
!write(*,*) 'cnode-out' ,index,ix,iy,iz,nclose,ncnode               
!write(*,*) 'cnodex-out',cnodex(:)               
!write(*,*) 'cnodey-out',cnodey(:)               
!write(*,*) 'cnodez-out',cnodez(:)               
!!write(*,*) 'cnodev-out',cnodev(:)               
!end if

                  if((nclose+ncnode).le.0) then
                        call beyond(ivtype,ncut,thres,cdf,ng,gcut,gcdf, &
                                   zmin,zmax,ltail,ltpar,middle,mpar, &
                                   utail,utpar,zval,cdfval,ierr)

                  else
!!
!! Estimate the local distribution by indicator kriging:
!!

                        ccdf(:)=0.0
                        ccdfo(:)=0.0
                        do ic=1,ncut
                              aclose(:)=0
                              a(:)=0.0
                              !call krige(ix,iy,iz,xx,yy,zz,ic,MAXKR1)
!write(*,*)index,close(1),close(2)
!write(*,*)'a-out:',a(:)
!write(*,*)'r-out:',r(:)
!if(index.eq.query)then
!write(*,*)'cdf-pre=',index,ic,cdf(ic)
!write(*,*)'ccdf-pre=',index,ic,ccdf(ic)
!write(*,*)'maxdat-pre=',index,ic,MAXDAT
!                      end if
!                      if(index.eq.697)then
                              !write(*,*)'vr-bef',vr(697,2)
!end if
!write(*,*)'NODES',NODES_LENGTH,xyzland,dim
                              !write(*,*)'neighs-bef',neighbours(:)
                              !write(*,*)'dista-bef',distance(:)
!                              !call krige(ix,iy,iz,xx,yy,zz,ic,cdf(ic),&
                              call krige_1D_v2(index,ix,iy,iz,xx,yy,zz,ic,cdf(ic),&
                              MAXCTX,MAXCTY,MAXCTZ,MAXKR1,MAXDAT,&
                              MXCUT,MAXROT,MAXNST,MAXNOD,MAXCUT,MAXXYZ,MAXKR2,&
                              ccdf(ic),mik,&
                              !nclose,ncnode,ktype,close,distance,&
                              nclose,ncnode,ktype,neighbours,distance,&
                              atnode,vr,aclose,x,y,z,nhd,softdat,&
                              vra,cnodex,cnodey,cnodez,&
                              cnodev,cnodet,ivtype,thres,&
                              nctx,ncty,nctz,icnode,ixnode,iynode,iznode,&
                              nx,ny,nz,xmn,ymn,zmn,xsiz,ysiz,zsiz,&
                              !aa,c0,cc,cmax,& 
                              aa,c0,cc,& 
                              it,nst,& 
                              r,rr,s,a,rotmat,&
                              !rotmat,covtab,&
                              NODES_LENGTH,xyzland,coord_ISOMAP_trans,dim,d_tree,query)
                              
!if(index.eq.query)then
!write(*,*)'cdf-pos=',index,ic,cdf(ic)
!write(*,*)'ccdf-pos=',index,ic,ccdf(ic)
!end if




!     call krige_old(ix,iy,iz,ic,MAXCTX,MAXCTY,MAXCTZ,MAXKR1,          &  
!     MAXROT,MAXDAT,MXCUT,MAXKR2,MAXNOD,MAXCUT,MAXXYZ,MAXNSTvar,idbg, &
!     imbsim,ivtype,ktype,ldbg,mik,nclose,ncnode,nctx,ncty,nctz,   &
!     xx,yy,zz,cdf(ic),ccdf(ic),atnode,softdat,vr,vra,x,y,z,icnode,&
!     iznode,iynode,ixnode,it,nst,aclose,cnodex,cnodey,cnodez,     &
!     cnodev,cnodet,covtab,thres,aa,c0,cc,cmax,close,beez,r,rr,s,a,&
!     rotmat)
                        j=1
                        end do


!
! Correct order relations:
!
if(index.eq.query)then
write(*,*)'ord-pre=',index,ccdf(:)
write(*,*)'ord-pre=',index,ccdfo(:)
end if
                        call ordrel(ivtype,ncut,ccdf,ccdfo,nviol,aviol,&
                                   xviol)
if(index.eq.query)then
write(*,*)'ord-pos=',index,ccdf(:)
write(*,*)'ord-pos=',index,ccdfo(:)
end if

!!
!! Draw from the local distribution:
!!
if(index.eq.query)then
write(*,*)'bey-pre',index,zval,cdfval,ierr,ccdfo(:)
end if
                        call beyond(ivtype,ncut,thres,ccdfo,ng,gcut, &
                                   gcdf,zmin,zmax,ltail,ltpar,middle,&
                                   mpar,utail,utpar,zval,cdfval,ierr)
if(index.eq.query)then
write(*,*)'bey-pos',index,zval,cdfval,ierr
end if

!!
!! Write some debugging information:
!!
!
!!                        if(idbg.ge.3) then
!!                              do ic=1,ncut
!!                              write(ldbg,202) ccdf(ic),ccdfo(ic)
!! 202                          format('  CDF (original and fixed)',2f7.4)
!!                              end do
!!                        endif
                  endif

                  sim(index) = zval
if(index.eq.1140)write(*,*)'sim-pre',sim(1140)

                  !if(level(index).eq.0)then 
                  !write(*,*)'probe1',index,sim(index),level(index)
                  !write(*,*)'probe2',level(resultsIdxIndex(:,index)) 
                  !write(*,*)'probe3',sim(resultsIdxIndex(:,index))  
                  !write(*,*)'probe4',cnodev(:)  
                  !end if

                  !if(index.eq.30020)write(*,*)'zval',index,sim(index)
                  lock(index) = 1
!write(*,*)threadId,index,'passed lock wait'

!
! END MAIN LOOP OVER NODES:
!
            tmp(index) = 0.0 

           end do

           end if

           end do
!$omp end parallel

!write(*,*)'sim(331)=',sim(331)
!write(*,*)'sim(313)=',sim(313)
!write(*,*)'sim(244)=',sim(244)

!
! Write this simulation to the output file:
!
            nxysim = 0
            do ic=1,ncut
                  ccdf(ic) = 0.0
            end do


            do ind=1,nxyz
!
! Calculate the cdf of the simulated values (for error checking):
!

!                  if(sim(ind).gt.UNEST) then
!                        nxysim = nxysim + 1
!                        do ic=1,ncut
!                              if(ivtype.eq.0) then
!                                    if(sim(ind).eq.thres(ic))
!     +                                ccdf(ic)=ccdf(ic)+1.0
!                              else
!                                    if(sim(ind).le.thres(ic))
!     +                                ccdf(ic)=ccdf(ic)+1.0
!                              end if
!                        end do
!                  endif
!

                  write(lout,'(g14.8)') sim(ind)
                  !index = int(order(ind)+0.5)
                  index = ind
                  iz = int((index-1)/nxy) + 1
                  iy = int((index-(iz-1)*nxy-1)/nx) + 1
                  ix = index - (iz-1)*nxy - (iy-1)*nx
                  xx = xmn + real(ix-1)*xsiz
                  yy = ymn + real(iy-1)*ysiz
                  zz = zmn + real(iz-1)*zsiz
!                  write(*,*)'posout',index,int(order(ind)+0.5),xx,yy,zz,sim(ind)


            end do

!            call system_clock(clock_stop,clock_rate,clock_max)
!            print *,'TIME (output)=',
!     +    real(clock_stop-clock_start,kind=8)/real(clock_rate,kind=8) 

!
! END MAIN LOOP OVER SIMULATIONS:
!
      end do
!#endif
      t2=secnds(t1)
      elapsed=t2
      write(*,'(a,f10.4,a)') '[FOR](sgs) ',elapsed, &
      's - time of simulation'

!
! Return to the main program:
!
      return
      end



      subroutine ctable
      !subroutine ctable(MAXNOD,MAXCXY,MAXCTX,MAXCTY,MAXCTZ,MAXROT)
!-----------------------------------------------------------------------
!
!               Establish the Covariance Look up Table
!               **************************************
!
! The idea is to establish a 3-D network that contains the covariance
! value for a range of grid node offsets that should be at as large
! as twice the search radius in each direction.  The reason it has to
! be twice as large as the search radius is because we want to use it
! to compute the data covariance matrix as well as the data-block
! covariance matrix.
!
! Secondly, we want to establish a search for nearby nodes that 
! in order of closeness as defined by the variogram.
!
!
!
! INPUT VARIABLES:
!
!   xsiz,ysiz,zsiz  Definition of the grid being considered
!   MAXCTX,Y,Z      Number of blocks in covariance table
!
!   covariance table parameters
!
!
!
! OUTPUT VARIABLES:  covtab()         Covariance table
!
! EXTERNAL REFERENCES:
!
!   sqdist          Computes 3-D anisotropic squared distance
!   sortem          Sorts multiple arrays in ascending order
!   cova3           Computes the covariance according to a 3-D model
!
!
!
!-----------------------------------------------------------------------
      use geostat
      parameter(TINY=1.0e-10)
      !include  'sisim.inc'
!      include  'sisim_new.inc'
      real*8    sqdist,hsqd
!
! Size of the look-up table:
!
      nctx = min(((MAXCTX-1)/2),(nx-1))
      ncty = min(((MAXCTY-1)/2),(ny-1))
      nctz = min(((MAXCTZ-1)/2),(nz-1))
!
! Initialize the covariance subroutine and cbb at the same time:
!
      call cova3(0.,0.,0.,0.,0.,0.,1,nst,MAXNST,c0,it,cc,aa,1,&
                 MAXROT,rotmat,cmax,cbb)
!
! Now, set up the table and keep track of the node offsets that are
! within the search radius:
!
      ilooku = max((ncut/2),1)
      nlooku = 0
      do icut=1,ncut
      irot = 1 + (icut-1)*MAXNST
      do i=-nctx,nctx
      xx = i * xsiz
      ic = nctx + 1 + i
      do j=-ncty,ncty
      yy = j * ysiz
      jc = ncty + 1 + j
      do k=-nctz,nctz
      zz = k * zsiz
      kc = nctz + 1 + k
if(ic.eq.49.and.jc.eq.26.and.kc.eq.52)then
write(*,*)'covtab-in-1=',xx,yy,zz
      end if
            call cova3(0.,0.,0.,xx,yy,zz,icut,nst,MAXNST,c0,it,cc,aa,&
                       irot,MAXROT,rotmat,cmax,cov)
            covtab(ic,jc,kc,icut) = cov
if(ic.eq.49.and.jc.eq.26.and.kc.eq.52)then
write(*,*)'covtab-in-2=',cov
      end if
            if(icut.eq.ilooku) then
               hsqd = sqdist(0.0,0.0,0.0,xx,yy,zz,isrot,MAXROT,rotmat)
               if(real(hsqd).le.radsqd) then
                  nlooku = nlooku + 1
!
! We subtract the covariance from a large value so that the ascending
! sort subroutine will accomplish the sort we want.  Furthermore, a
! fraction of the distance is also taken off so that we search by
! anisotropic distance once we are beyond the range:
!
                  tmp(nlooku)   =-(covtab(ic,jc,kc,icut)-TINY*hsqd)
                  order(nlooku) =real((kc-1)*MAXCXY+(jc-1)*MAXCTX+ic)
               endif 
            endif
      end do
      end do
      end do
      end do
      ! write(*,*)'nlooku = ',nlooku
!
! Finished setting up the look-up table, now order the nodes such
! that the closest ones, according to variogram distance, are searched
! first. Note: the "loc" array is used because I didn't want to make 
! special allowance for 2 byte integers in the sorting subroutine:
!
      call sortem(1,nlooku,tmp,1,order)
      do il=1,nlooku
            loc = int(order(il))
            iz  = int((loc-1)/MAXCXY) + 1
            iy  = int((loc-(iz-1)*MAXCXY-1)/MAXCTX) + 1
            ix  = loc-(iz-1)*MAXCXY - (iy-1)*MAXCTX
            iznode(il) = iz
            iynode(il) = iy
            ixnode(il) = ix
      end do
      if(nodmax.gt.MAXNOD) then
            write(ldbg,*)
            write(ldbg,*) 'The maximum number of close nodes = ',nodmax
            write(ldbg,*) 'this was reset from your specification due '
            write(ldbg,*) 'to storage limitations.'
            nodmax = MAXNOD
      endif
!
! Debugging output if requested:
!
      if(idbg.le.2) return
      write(ldbg,*)
      write(ldbg,*) 'There are ',nlooku,' nearby nodes that will be '
      write(ldbg,*) 'checked until enough close data are found.'
      write(ldbg,*)
      if(idbg.lt.4) return
      do i=1,nlooku
            xx = (ixnode(i) - nctx - 1) * xsiz
            yy = (iynode(i) - ncty - 1) * ysiz
            zz = (iznode(i) - nctz - 1) * zsiz
            write(ldbg,100) i,xx,yy,zz
      end do
 100  format('Point ',i6,' at ',3f18.6)
!
! All finished:
!
      return
      end



      subroutine srchnd(nloop2,index2,ix,iy,iz)
!-----------------------------------------------------------------------
!
!               Search for nearby Simulated Grid nodes
!               **************************************
!
! The idea is to spiral away from the node being simulated and note all
! the nearby nodes that have been simulated.
!
!
!
! INPUT VARIABLES:
!
!   ix,iy,iz        index of the point currently being simulated
!   sim(nx,ny,nz)   the simulation so far
!   nodmax          the maximum number of nodes that we want
!   nlooku          the number of nodes in the look up table
!   i[x,y,z]node    the relative indices of those nodes.
!   [x,y,z]mn       the origin of the global grid netwrok
!   [x,y,z]siz      the spacing of the grid nodes.
!
!
!
! OUTPUT VARIABLES:
!
!   ncnode          the number of close nodes
!   icnode()        the number in the look up table
!   cnode[x,y,z]()  the location of the nodes
!   cnodev()        the values at the nodes
!
!
!
!-----------------------------------------------------------------------
      use geostat
      use srch
      use global
      !include  'sisim.inc'
!      include  'sisim_new.inc'
      integer, allocatable :: cercanos(:)
!  Creo vecin como variable interna    
      real*8,allocatable    :: vecin(:)
      real temporal
      
      
!
! Edicion Rodrigo Gutierrez
!
      !write(*,*)'srchnd: dim=',dim,' nodmax=',nodmax 
      allocate(vecin(dim))
      
      allocate(cercanos(nodmax))
      

      cercanos = 0
      vecin = 0
      
      
! Coordenadas del índice en cuestion.
!  d_tree = dimensiones máximas a usar.

      !write(*,*)'srchnd: cargando coord_ISOMAP_trans' 
      do iii=1,dim
          !write(*,*)'srchnd: coord_ISOMAP_trans(',iii,',',index2,')' 
          !vecin(iii)=coord_ISOMAP(index2,iii)
          vecin(iii)=coord_ISOMAP_trans(iii,index2)
      end do
         
! Busqueda de vecinos del índice en cuestion.

      !write(*,*)'ix = ',ix
      !write(*,*)'iy = ',iy
      !write(*,*)'iz = ',iz
      !write(*,*)'vecin = ',vecin
      !write(*,*)'dim = ',dim 
      !write(*,*)'nodmax = ',nodmax
      
      !write(*,*)'ex_results_d = ',ex_results_d
      
      !write(*,*)'ex_sim_array = ',ex_sim_array
      !write(*,*)'sim = ',sim
      !write(*,*)'srchnd: entrando a exhaustive_srch' 
      call exhaustive_srch(vecin,dim,nodmax,cercanos,&
                              ex_results_d,nloop2)

      !write(*,*)'cercanos = ',cercanos
      
       

       
      !              k = int(index/(nx*ny)) + 1
      !"      !if(k.lt.1.or.k.gt.nz) go to 1
      !      j = (int((index-(k-1)*nx*ny)/nx)) + 1
      !      !if(j.lt.1.or.j.gt.ny) go to 1
      !      i = index - (k-1)*nx*ny - (j-1)*nx
      
!
! Fin Edición Rodrigo Gutierrez
!      
!
! Consider all the nearby nodes until enough have been found:
!
      !write(*,*)'srchnd: continuar srchnd' 
      ncnode = 0
      ncsec  = 0
      do il=1,nlooku
            if(ncnode.eq.nodmax) return
            i = ix + (int(ixnode(il))-nctx-1)
            if(i.lt.1.or.i.gt.nx) go to 1
            j = iy + (int(iynode(il))-ncty-1)
            if(j.lt.1.or.j.gt.ny) go to 1
            k = iz + (int(iznode(il))-nctz-1)
            if(k.lt.1.or.k.gt.nz) go to 1
!
! Check this potentially informed grid node:
!
            index = (k-1)*nx*ny + (j-1)*nx + i
            !write(*,*)'index lookup = ',index
            do i=1,nodmax
               
               if(cercanos(i).eq.index) then
                  
                  index = cercanos(i)
                  
                  
                  go to 30

               end if
               
               if(i.eq.nodmax) then 
                  go to 1
               end if
               
            end do
            
 30         continue     
            
            if(sim(index).gt.UNEST.or.tmp(index).gt.0.5) then
                  if(sim(index).le.UNEST.and.tmp(index).gt.0.5) then
                        ncsec  = ncsec + 1
                        if(ncsec.gt.maxsec) go to 1
                  end if
                  ncnode         = ncnode + 1
                  !write(*,*)'zona peligrosa, index = ',index
                  icnode(ncnode) = il
                  cnodex(ncnode) = xmn + real(i-1)*xsiz
                  cnodey(ncnode) = ymn + real(j-1)*ysiz
                  cnodez(ncnode) = zmn + real(k-1)*zsiz
                  cnodev(ncnode) = sim(index)
                  cnodet(ncnode) = tmp(index)
            endif
            
            
 1          continue
      end do
!
! Return to calling program:
!
      return
      end

      !subroutine srchndPushLocal(nloop2,index2,ix,iy,iz,ncnode,icnode,cnodeid)
      subroutine srchndPushLocal(nloop2,index2,ix,iy,iz)
!-----------------------------------------------------------------------
!
!               Search for nearby Simulated Grid nodes
!               **************************************
!
! The idea is to spiral away from the node being simulated and note all
! the nearby nodes that have been simulated.
!
!
!
! INPUT VARIABLES:
!
!   ix,iy,iz        index of the point currently being simulated
!   sim(nx,ny,nz)   the simulation so far
!   nodmax          the maximum number of nodes that we want
!   nlooku          the number of nodes in the look up table
!   i[x,y,z]node    the relative indices of those nodes.
!   [x,y,z]mn       the origin of the global grid netwrok
!   [x,y,z]siz      the spacing of the grid nodes.
!
!
!
! OUTPUT VARIABLES:
!
!   ncnode          the number of close nodes
!   icnode()        the number in the look up table
!   cnode[x,y,z]()  the location of the nodes
!   cnodev()        the values at the nodes
!
!
!
!-----------------------------------------------------------------------
      use geostat
      use srch
      use global
      !include  'sisim.inc'
!      include  'sisim_new.inc'

      !integer ncnode, icnode
      !integer icnode(MAXNOD)
      !integer cnodeid(MAXNOD)
      !integer i,j,k,il,indexloc,ncsec


      integer, allocatable :: cercanos(:)
!  Creo vecin como variable interna    
      real*8,allocatable    :: vecin(:)
      real temporal
      
      
!
! Edicion Rodrigo Gutierrez
!
      !write(*,*)'srchnd: dim=',dim,' nodmax=',nodmax 
      allocate(vecin(dim))
      
      allocate(cercanos(nodmax))
      

      cercanos = 0
      vecin = 0
      
      
! Coordenadas del índice en cuestion.
!  d_tree = dimensiones máximas a usar.

      !write(*,*)'srchnd: cargando coord_ISOMAP_trans' 
      do iii=1,dim
          !write(*,*)'srchnd: coord_ISOMAP_trans(',iii,',',index2,')' 
          !vecin(iii)=coord_ISOMAP(index2,iii)
          vecin(iii)=coord_ISOMAP_trans(iii,index2)
      end do
         
! Busqueda de vecinos del índice en cuestion.

      !write(*,*)'ix = ',ix
      !write(*,*)'iy = ',iy
      !write(*,*)'iz = ',iz
      !write(*,*)'vecin = ',vecin
      !write(*,*)'dim = ',dim 
      !write(*,*)'nodmax = ',nodmax
      
      !write(*,*)'ex_results_d = ',ex_results_d
      
      !write(*,*)'ex_sim_array = ',ex_sim_array
      !write(*,*)'sim = ',sim
      !write(*,*)'srchnd: entrando a exhaustive_srch' 
      call exhaustive_srch(vecin,dim,nodmax,cercanos,&
                              ex_results_d,nloop2)

      !write(*,*)'cercanos = ',cercanos
      
       

       
      !              k = int(index/(nx*ny)) + 1
      !"      !if(k.lt.1.or.k.gt.nz) go to 1
      !      j = (int((index-(k-1)*nx*ny)/nx)) + 1
      !      !if(j.lt.1.or.j.gt.ny) go to 1
      !      i = index - (k-1)*nx*ny - (j-1)*nx
      
!
! Fin Edición Rodrigo Gutierrez
!      
!
! Consider all the nearby nodes until enough have been found:
!
      !write(*,*)'srchnd: continuar srchnd' 
      ncnode = 0
      ncsec  = 0
      do il=1,nlooku
            if(ncnode.eq.nodmax) return
            i = ix + (int(ixnode(il))-nctx-1)
            if(i.lt.1.or.i.gt.nx) go to 1
            j = iy + (int(iynode(il))-ncty-1)
            if(j.lt.1.or.j.gt.ny) go to 1
            k = iz + (int(iznode(il))-nctz-1)
            if(k.lt.1.or.k.gt.nz) go to 1
!
! Check this potentially informed grid node:
!
            indexloc = (k-1)*nx*ny + (j-1)*nx + i
            !write(*,*)'index lookup = ',index
            do i=1,nodmax
               
               if(cercanos(i).eq.indexloc) then
                  
                  indexloc = cercanos(i)
                  
                  
                  go to 30

               end if
               
               if(i.eq.nodmax) then 
                  go to 1
               end if
               
            end do
            
 30         continue     
            
            if(sim(indexloc).gt.UNEST.or.tmp(indexloc).gt.0.5) then
                  if(sim(indexloc).le.UNEST.and.tmp(indexloc).gt.0.5) then
                        ncsec  = ncsec + 1
                        if(ncsec.gt.maxsec) go to 1
                  end if
                  ncnode         = ncnode + 1
                  !write(*,*)'zona peligrosa, index = ',index
                  !icnode(ncnode) = il
                  !cnodex(ncnode) = xmn + real(i-1)*xsiz
                  !cnodey(ncnode) = ymn + real(j-1)*ysiz
                  !cnodez(ncnode) = zmn + real(k-1)*zsiz
                  !cnodev(ncnode) = sim(index)
                  !cnodet(ncnode) = tmp(index)
                  cnodeid(ncnode) = indexloc
                  icnode(ncnode) = il
            endif
            
            
 1          continue
      end do
!
! Return to calling program:
!
      return
      end

      subroutine srchndPopLocal(nloop2,index2,ix,iy,iz,nx,ny,ncnode,nd,nloop,MAXNOD,MAXORD,results,xmn,ymn,zmn,xsiz,ysiz,zsiz,tmp,cnodex,cnodey,cnodez,cnodet)
!-----------------------------------------------------------------------
      !use geostat
      !use srch
      !use global
      !!include  'sisim.inc'
      !include  'sisim_new.inc'

      implicit none


      integer,intent(in) :: nloop2,index2,ix,iy,iz,nx,ny,ncnode,nd,nloop,MAXNOD,MAXORD
      !integer,intent(in) :: resultsIdxIndex(nd,nloop)
      integer,intent(in) :: results(nd)
      !integer,dimension(:,:),intent(in) :: resultsIdxIndex
      real*8,intent(in) :: xmn,ymn,zmn,xsiz,ysiz,zsiz
      real,intent(in) :: tmp(MAXORD)
      real,intent(inout) :: cnodex(MAXNOD),cnodey(MAXNOD),cnodez(MAXNOD),cnodet(MAXNOD)


      integer nxy,il,indexloc,i,j,k

!if(index2.eq.32256)then                  
!write(*,*)'neig-pos2: ',index2,results(:)
!end if

!if(index2.eq.1157)then
!write(*,*)'pos:',index2,ncnode,xmn,ymn,zmn,xsiz,ysiz,zsiz
!end if             
!write(*,*)'resultsIdxIndex: ',resultsIdxIndex(:,index2) 
      nxy=nx*ny
      do il=1,ncnode

            !indexloc = cnodeid(il)
            !indexloc = resultsIdxIndex(il,index2)
            indexloc = results(il)
            k = int((indexloc-1)/nxy) + 1
            j = int((indexloc-(k-1)*nxy-1)/nx) + 1
            i = indexloc - (k-1)*nxy - (j-1)*nx

!if(index2.eq.1157)then
!!write(*,*)'pos2:',ncnode,il,indexloc,k,j,i
!write(*,*)'cnodex-loop:',ncnode,il,indexloc,xmn,real(i-1),xsiz
!write(*,*)'cnodey-loop:',ncnode,il,indexloc,ymn,real(j-1),ysiz
!write(*,*)'cnodez-loop:',ncnode,il,indexloc,zmn,real(k-1),zsiz
!end if            
!            ncnode         = ncnode + 1
!            cnodeid(il) = indexloc
!            icnode(il) = il
!write(*,*)'il:',index2,il,MAXNOD,indexloc
!write(*,*)'il:',index2,il,MAXNOD,indexloc,tmp(indexloc)
            cnodex(il) = xmn + real(i-1)*xsiz
            cnodey(il) = ymn + real(j-1)*ysiz
            cnodez(il) = zmn + real(k-1)*zsiz
!!            cnodev(il) = sim(indexloc)
            cnodet(il) = tmp(indexloc)
!!if(index2.eq.1157)then
!!write(*,*)'cnodex-post:',cnodex(il)
!!write(*,*)'cnodey-post:',cnodey(il)
!!write(*,*)'cnodez-post:',cnodez(il)
!!end if
      end do
!if(index2.eq.1157)then
!write(*,*) 'cnode-in' ,index2,ix,iy,iz,ncnode               
!write(*,*) 'cnodex-in',cnodex(:)               
!write(*,*) 'cnodey-in',cnodey(:)               
!write(*,*) 'cnodez-in',cnodez(:)               
!!write(*,*) 'cnodev-in',cnodev(:)               
!end if



!
! Return to calling program:
!
      return
      end



!      subroutine srchndPush(ix,iy,iz,MAXNOD,MAXXYZ,MAXORD,MXYZ,ncnode,
!     +        maxsec,nctx,ncty,nctz,nlooku,nodmax,nx,ny,nz,UNEST,
!     +        xmn,ymn,zmn,xsiz,ysiz,zsiz,icnode,ixnode,iynode,iznode,
!     +        cnodex,cnodey,cnodez,cnodev,cnodet,tmp,sim,cnodeid)
!
!
!      implicit none
!
!      integer ix,iy,iz,MAXNOD,MAXXYZ,MAXORD,MXYZ,ncnode,maxsec,
!     +        nctx,ncty,nctz,nlooku,nodmax,nx,ny,nz 
!      real    UNEST,xmn,ymn,zmn,xsiz,ysiz,zsiz
!      integer icnode(MAXNOD),ixnode(MAXXYZ),iynode(MAXXYZ),
!     +        iznode(MAXXYZ)
!      real    cnodex(MAXNOD),cnodey(MAXNOD),cnodez(MAXNOD),
!     + cnodev(MAXNOD),cnodet(MAXNOD),tmp(MAXORD),sim(MXYZ)
!      integer cnodeid(MAXNOD)
!
!      integer i,j,k,il,indexloc,ncsec
!c
!c Consider all the nearby nodes until enough have been found:
!c
!      ncnode = 0
!      ncsec  = 0
!
!
!      do 2 il=1,nlooku
!
!            if(ncnode.eq.nodmax) return
!            i = ix + (int(ixnode(il))-nctx-1)
!            if(i.lt.1.or.i.gt.nx) go to 2
!            j = iy + (int(iynode(il))-ncty-1)
!            if(j.lt.1.or.j.gt.ny) go to 2
!            k = iz + (int(iznode(il))-nctz-1)
!            if(k.lt.1.or.k.gt.nz) go to 2
!
!            indexloc = (k-1)*nx*ny + (j-1)*nx + i
!
!            if(sim(indexloc).gt.UNEST) then
!                  ncnode         = ncnode + 1
!                  cnodeid(ncnode) = indexloc
!                  icnode(ncnode) = il
!            endif
! 2    continue
!
!!
!! Return to calling program:
!!
!      return
!      end

!      subroutine srchndPop(ix,iy,iz,MAXNOD,MAXXYZ,MAXORD,MXYZ,ncnode,
!     +        maxsec,nctx,ncty,nctz,nlooku,nodmax,nx,ny,nz,UNEST,
!     +        xmn,ymn,zmn,xsiz,ysiz,zsiz,icnode,ixnode,iynode,iznode,
!     +        cnodex,cnodey,cnodez,cnodev,cnodet,tmp,sim,cnodeid)
!
!
!      implicit none
!
!      integer ix,iy,iz,MAXNOD,MAXXYZ,MAXORD,MXYZ,ncnode,maxsec,
!     +        nctx,ncty,nctz,nlooku,nodmax,nx,ny,nz 
!      real    UNEST,xmn,ymn,zmn,xsiz,ysiz,zsiz
!      integer icnode(MAXNOD),ixnode(MAXXYZ),iynode(MAXXYZ),
!     +        iznode(MAXXYZ)
!      real    cnodex(MAXNOD),cnodey(MAXNOD),cnodez(MAXNOD),
!     + cnodev(MAXNOD),cnodet(MAXNOD),tmp(MAXORD),sim(MXYZ)
!      integer cnodeid(MAXNOD)
!
!      integer i,j,k,il,indexloc,ncsec,nxy
!c
!c Consider all the nearby nodes until enough have been found:
!c
!      nxy=nx*ny
!      do il=1,ncnode
!
!            indexloc = cnodeid(il)
!            k = int((indexloc-1)/nxy) + 1
!            j = int((indexloc-(k-1)*nxy-1)/nx) + 1
!            i = indexloc - (k-1)*nxy - (j-1)*nx
!
!c            ncnode         = ncnode + 1
!c            cnodeid(il) = indexloc
!c            icnode(il) = il
!            cnodex(il) = xmn + real(i-1)*xsiz
!            cnodey(il) = ymn + real(j-1)*ysiz
!            cnodez(il) = zmn + real(k-1)*zsiz
!c            cnodev(il) = sim(indexloc)
!            cnodet(il) = tmp(indexloc)
!
!      end do
!c
!c Return to calling program:
!c
!      return
!      end
!
!

      subroutine krige_1D(index2,ix,iy,iz,xx,yy,zz,icut,gmean,&
                       MAXCTX,MAXCTY,MAXCTZ,MAXKR1,MAXDAT,&
                       MXCUT,MAXROT,MAXNST,MAXNOD,MAXCUT,MAXXYZ,MAXKR2,&
                       cmean,mik,&
                       !nclose,ncnode,ktype,close,distance,&
                       nclose,ncnode,ktype,neighbours,distance,&
                       atnode,vr,aclose,x,y,z,nhd,softdat,&
                       vra,cnodex,cnodey,cnodez,&
                       cnodev,cnodet,ivtype,thres,&
                       nctx,ncty,nctz,icnode,ixnode,iynode,iznode,&
                       nx,ny,nz,xmn,ymn,zmn,xsiz,ysiz,zsiz,&
                       !aa,c0,cc,cmax,& 
                       aa,c0,cc,& 
                       it,nst,& 
                       r,rr,s,a,rotmat,covtab,&
                       !rotmat,covtab,&
                       NODES,xyzland,coord_ISOMAP_trans,dim,d_tree)
                       !r,rr,s,a,rotmat)


!-----------------------------------------------------------------------
!
!            Builds and Solves the SK or OK Kriging System
!            *********************************************
!
! INPUT VARIABLES:
!
!   ix,iy,iz        index of the point currently being simulated
!   xx,yy,zz        location of the point currently being simulated
!   icut            cutoff number to use for either the covariance look
!                     up table or the covariance calculation
!
!
!
! OUTPUT VARIABLES:
!
!   cmean           kriged estimate
!
!
! 
! EXTERNAL REFERENCES: ksol   Gaussian elimination system solution
!
!
! NOTE: 1. the array "aclose" is used to flag those samples which exist
!          at the cutoff currently being kriged.
!
!
!-----------------------------------------------------------------------
      implicit none
      integer,intent(in) ::index2,ix,iy,iz,icut,MAXCTX,MAXCTY,MAXCTZ,MAXKR1
      integer,intent(in)::MAXROT,MAXDAT,MXCUT,MAXNST,MAXNOD,MAXCUT,MAXXYZ,MAXKR2
      integer,intent(in) ::mik,nclose,ncnode,ktype,nhd,ivtype
      integer,intent(in) ::nctx,ncty,nctz
      real*8,intent(in) :: xx,yy,zz
      real,intent(in) :: gmean
      integer,intent(in) :: neighbours(MAXKR1)
      real*8,intent(in) :: distance(MAXKR1) 
      logical,intent(in) :: atnode(MAXDAT)
      real,intent(in) :: vr(MAXDAT,MXCUT),x(MAXDAT),y(MAXDAT),z(MAXDAT)
      integer,intent(inout) :: aclose(MAXKR1),softdat(MAXKR1)
      real,intent(inout) :: vra(MAXKR1)
      real,intent(in) :: cnodex(MAXNOD),cnodey(MAXNOD),cnodez(MAXNOD)
      real,intent(in) :: cnodev(MAXNOD),cnodet(MAXNOD),thres(MAXCUT) 
      integer,intent(in) :: icnode(MAXNOD)
      integer,intent(in) :: iznode(MAXXYZ),iynode(MAXXYZ),ixnode(MAXXYZ)
      integer,intent(in) :: nx,ny,nz
      real*8,intent(in) :: xmn,ymn,zmn,xsiz,ysiz,zsiz
      real, intent(in) :: aa(MAXCUT*MAXNST),c0(MAXCUT),cc(MAXCUT*MAXNST) !,cmax(MAXCUT) 
      real,intent(in) :: it(MAXCUT*MAXNST),nst(MAXCUT),covtab(MAXCTX,MAXCTY,MAXCTZ,MAXCUT)
      real*8,intent(inout) ::  r(MAXKR1),rr(MAXKR1),s(MAXKR1),a(MAXKR2),rotmat(MAXROT,3,3) 
      !real*8,intent(inout) ::  rotmat(MAXROT,3,3) 
      !real*8,intent(inout) :: a(MAXKR2),rotmat(MAXROT,3,3) 
     
      !integer,intent(in) :: NCLOSEvar,nloop 
      !integer,intent(in):: resultsIdxIndex(NCLOSEvar,nloop)
      integer,intent(in):: NODES,xyzland
      !real*8,intent(in) :: coord_ISOMAP_trans(NODES,xyzland) 
      real*8,intent(in) :: coord_ISOMAP_trans(xyzland,NODES) 
      integer, intent(in):: dim,d_tree

      real,intent(inout) :: cmean

      !use      geostat
      !include 'sisim.inc'
      !include 'sisim_new.inc'
      !integer aclose(nclose)
      !real*8, allocatable ::  r(:),rr(:),s(:),a(:)
      logical krig,somesoft,bothsoft,testind
      integer mclose,index,na,neq,i,j,irot,in,j1,iii,ind,ising,k
      integer ix1,iy1,iz1,ix2,iy2,iz2,ix11,iy11,iz11,ix22,iy22,iz22,ii,jj,kk,ppp
      real x1,y1,z1,x2,y2,z2,sumwt,cov
      real*8 dist,covd,sumdiff,sumdiff2
      integer ind1,ind2,notnull
      real cmax
!
! Size of the kriging system:  Some of the data values may be missing
! which would correspond to a constraint interval.  Note that there
! should not be any missing values if the median approximation is being
! considered.  The variable ``krig'' is used
! to flag whether kriging is to be done or if the previous weights are
! to be used.
!


      !allocate(r(MAXKR1))
      !allocate(rr(MAXKR1))
      !allocate(s(MAXKR1))
      !allocate(a(MAXKR2))
      r(:)=0.0
      rr(:)=0.0
      s(:)=0.0
      a(:)=0.0

!      write(*,*)'Inside krige',MAXKR1
!      write(*,*)'neighs:',neighbours(:)
!      write(*,*)'dista:',distance(:)
      !allocate(aclose(MAXKR1))
      !write(*,*)'passed'
      somesoft = .false.
      krig     = .true.
      if(mik.eq.1.and.icut.gt.1) krig = .false.
      if(krig) then
            mclose = 0
!write(*,*)'krige:',index2,nclose
            do i=1,nclose
                  index     =  neighbours(i)
                  write(*,*)'index',i,index,MAXDAT
                  !write(*,*)'neig',index,i,neighbours(i)
                  !write(*,*)'atnode',index,i,atnode(index)
                  !write(*,*)'vr-size:',shape(vr),icut
                  !write(*,*)'vr',index,i,icut,vr(index,icut)
      !write(*,*)'Processing:',i,index
                  if(.not.atnode(index)) then
                  if(vr(index,icut).ge.0.0) then
      write(*,*)'Processing-inside:',i,index,mclose
                        mclose = mclose + 1
                        !write(*,*)'aclose:',mclose,MAXKR1,index
                        aclose(mclose) = index
                  endif
                  endif
            end do
            na  = mclose + ncnode
            neq = na + ktype
      endif
!
! There are no data yet:
!
      irot   = 1 + (icut-1)*MAXNST
      !write(*,*)'inside',na,neq


!! From SGS_LVA_LEVELS
!      !a(1) = 0.0
!      in = 0
!      do j=1,na
!         !ind1=resultsIdxIndex(i,index)
!         !ind1=0
!         ind1=neighbours(j)
!!         if(j.le.mclose) then
!!             index  = aclose(j)
!!             vra(j) = vr(index,icut)
!!             x1     = x(index)
!!             y1     = y(index)
!!             z1     = z(index)
!!             if(index.gt.nhd) softdat(j) = .true.
!
!         do i=1,j
!            !ind2=resultsIdxIndex(j,index)
!            !ind2=0
!            ind2=neighbours(i)
!            sumdiff=0.0
!            do k=1,dim
!!                write(*,*)'ent-cova3-ind1',index2,k,ind1,coord_ISOMAP_trans(k,ind1) 
!!                write(*,*)'ent-cova3-ind2',index2,k,ind2,coord_ISOMAP_trans(k,ind2) 
!                sumdiff=sumdiff + (coord_ISOMAP_trans(k,ind1)-coord_ISOMAP_trans(k,ind2))*(coord_ISOMAP_trans(k,ind1)-coord_ISOMAP_trans(k,ind2)) 
!            end do
!!            write(*,*)'ent-cova3-0',sumdiff
!            dist=dsqrt (sumdiff)
!            !dist=dsqrt ( sum( ( coord_ISOMAP_trans(1:dim,ind1)-coord_ISOMAP_trans(1:dim,ind2) ) **2    )    )
!            !write(*,*)'ent-cova3-0', coord_ISOMAP_trans(1:dim,ind1),coord_ISOMAP_trans(1:dim,ind2)
!            call cova3_1D_sp(real(dist),icut,nst,MAXNST,MAXCUT,c0,it,cc,aa,cmax,cov) 
!!            write(*,*)'ent-cova3-1',dist,dble(cov)
!            in = in + 1
!            write(*,*)'debug',index2,in,MAXKR2
!            a(in) = dble(cov)
!         end do
!!write(*,*)'kkk',j,i,ind1
!         vra(j)=cnodev(j) 
!         if(d_tree/=dim) then
!             dist=distance(j) + sum( ( coord_ISOMAP_trans(d_tree+1:dim,ind1)-coord_ISOMAP_trans(d_tree+1:dim,index) ) **2    ) !dist in tree is only for the first few dims
!         else
!             dist=distance(j)
!         end if
!         dist=dsqrt(dist)
!         call cova3_1D_sp(real(dist),icut,nst,MAXNST,MAXCUT,c0,it,cc,aa,cmax,cov)
!!         write(*,*)'ent-cova3-2',ist,dble(cov)
!         r(j)=dble(cov) 
!
!      end do
!      if(index2.eq.265)then
!write(*,*)'nclose',nclose
!write(*,*)'na',na
!write(*,*)'in',in
!write(*,*)'a',a(:)
!write(*,*)'r',r(:)
!      end if
!!     !do i=1,na
!!     !   a(neq*(i-1)+na+1) = dble(cmax)
!!     !   a(neq*na+i)       = dble(cmax)
!!     !end do
!!
!!
!!! END From SGS_LVA_LEVELS


!      a=0.0
!
! Set up kriging matrices:
!
      in = 0
      j1 = 0
      do 1 j=1,na
            softdat(j) = .false.
!
! Sort out the actual location of point "j"
!
            if(j.le.mclose) then
                  index  = aclose(j)
                  vra(j) = vr(index,icut)
                  x1     = x(index)
                  y1     = y(index)
                  z1     = z(index)
                  if(index.gt.nhd) softdat(j) = .true.
            else
!
! It is a previously simulated node (keep index for table look-up):
!
                  index  = j-mclose
                  x1     = cnodex(index)
                  y1     = cnodey(index)
                  z1     = cnodez(index)
!
! Is this node informed by a hard datum or a soft datum?
!
                  if(cnodet(index).le.0.5) then
                        if(ivtype.eq.0) then
                           vra(j) = 0.0
                           if(int(cnodev(index)+0.5).eq.int(thres(icut)+0.5)) vra(j) = 1.0
                        else
                           vra(j) = 1.0
                           if(cnodev(index).gt.thres(icut)) vra(j) = 0.0
                        end if
                  else
                        iii    = int(cnodet(index)+0.5)
                        vra(j) = vr(iii,icut)
                        softdat(j) = .true.
                  end if
                  call getindx(nx,xmn,xsiz,dble(cnodex(index)),ix11,testind)
                  call getindx(ny,ymn,ysiz,dble(cnodey(index)),iy11,testind)
                  call getindx(nz,zmn,zsiz,dble(cnodez(index)),iz11,testind)
!if(ix.eq.25.and.iy.eq.11.and.iz.eq.28)then
!write(*,*)'node-n:',ix,iy,iz,nx,ny,nz
!write(*,*)'node-mn:',ix,iy,iz,xmn,ymn,zmn
!write(*,*)'node-siz:',ix,iy,iz,xsiz,ysiz,zsiz
!write(*,*)'node-cnode:',ix,iy,iz,dble(cnodex(index)),dble(cnodey(index)),dble(cnodez(index))
!write(*,*)'node-iout:',ix,iy,iz,ix11,iy11,iz11
!          end if
                  !ind    = icnode(index)
                  !ix1    = ix + (int(ixnode(ind))-nctx-1)
                  !iy1    = iy + (int(iynode(ind))-ncty-1)
                  !iz1    = iz + (int(iznode(ind))-nctz-1)

                  ix1    = ix + (int(ix11)-nctx-1)
                  iy1    = iy + (int(iy11)-ncty-1)
                  iz1    = iz + (int(iz11)-nctz-1)
            end if
!
! Only set up the matrix and the RHS if kriging:
!
            if(krig) then
               do 2 i=1,j
!
! Sort out the actual location of point "i"
!
                  if(i.le.mclose) then
                        index  = aclose(i)
                        x2     = x(index)
                        y2     = y(index)
                        z2     = z(index)
                  else
!
! It is a previously simulated node (keep index for table look-up):
!
                        index  = i-mclose
                        x2     = cnodex(index)
                        y2     = cnodey(index)
                        z2     = cnodez(index)
!                        ind    = icnode(index)
!                        ix2    = ix + (int(ixnode(ind))-nctx-1)
!                        iy2    = iy + (int(iynode(ind))-ncty-1)
!                        iz2    = iz + (int(iznode(ind))-nctz-1)
                        call getindx(nx,xmn,xsiz,dble(cnodex(index)),ix22,testind)
                        call getindx(ny,ymn,ysiz,dble(cnodey(index)),iy22,testind)
                        call getindx(nz,zmn,zsiz,dble(cnodez(index)),iz22,testind)
                        ix2    = ix + (int(ix22)-nctx-1)
                        iy2    = iy + (int(iy22)-ncty-1)
                        iz2    = iz + (int(iz22)-nctz-1)
                  endif
!
! Now, get the covariance value:
!
                  in = in + 1
!
! Decide whether or not to use the covariance look-up table:
!
!if(ix.eq.25.and.iy.eq.11.and.iz.eq.28)then
!write(*,*)'node: points',ix,iy,iz
!write(*,*)'node: coord1',x1,y1,z1
!write(*,*)'node: coord2',x2,y2,z2
!          end if
                  if(j.le.mclose.or.i.le.mclose) then
                  ppp=0
                        call cova3(x1,y1,z1,x2,y2,z2,icut,nst,MAXNST,&
                             c0,it,cc,aa,irot,MAXROT,rotmat,cmax,cov)
                        a(in) = dble(cov)
!if(ix.eq.25.and.iy.eq.11.and.iz.eq.28)write(*,*)'a2: use cova3',ix,iy,iz 
                  else
!
! Try to use the covariance look-up (if the distance is in range):
!
                        ii = nctx + 1 + (ix1 - ix2)
                        jj = ncty + 1 + (iy1 - iy2)
                        kk = nctz + 1 + (iz1 - iz2)
!if(ix.eq.25.and.iy.eq.11.and.iz.eq.28)write(*,*)'box1',ix,iy,iz,ii,jj,kk
!if(ix.eq.25.and.iy.eq.11.and.iz.eq.28)write(*,*)'box2',ix,iy,iz,nctx,ncty,nctz
!if(ix.eq.25.and.iy.eq.11.and.iz.eq.28)write(*,*)'box3',ix,iy,iz,MAXCTX,MAXCTY,MAXCTZ
                        if(ii.lt.1.or.ii.gt.MAXCTX.or.jj.lt.1.or.jj.gt.MAXCTY.or.kk.lt.1.or.kk.gt.MAXCTZ) then
                  ppp=0
                              call cova3(x1,y1,z1,x2,y2,z2,icut,nst,&
                                         MAXNST,c0,it,cc,aa,irot,MAXROT,&
                                         rotmat,cmax,cov)
                              a(in) = dble(cov)
!if(ix.eq.25.and.iy.eq.11.and.iz.eq.28)write(*,*)'a: use cova3',ix,iy,iz 
                        else
                  ppp=0
                              a(in) = dble(covtab(ii,jj,kk,icut))
!if(ix.eq.25.and.iy.eq.11.and.iz.eq.28)write(*,*)'a: use covtab',ix,iy,iz 
                        endif
                  endif
!if(ix.eq.25.and.iy.eq.11.and.iz.eq.28)then
!write(*,*)'node:',ix,iy,iz,'a',a(in)
!end if
 2          continue
!
! Get the RHS value (possibly with covariance look-up table):
!
            if(j.le.mclose) then
                  ppp=0
                  call cova3(xx,yy,zz,x1,y1,z1,icut,nst,MAXNST,&
                             c0,it,cc,aa,irot,MAXROT,rotmat,cmax,cov)
                  r(j) = dble(cov)
!if(ix.eq.25.and.iy.eq.11.and.iz.eq.28)write(*,*)'r: use cova3',ix,iy,iz 
            else
!
! Try to use the covariance look-up (if the distance is in range):
!
                  ii = nctx + 1 + (ix - ix1)
                  jj = ncty + 1 + (iy - iy1)
                  kk = nctz + 1 + (iz - iz1)
                  if(ii.lt.1.or.ii.gt.MAXCTX.or.jj.lt.1.or.jj.gt.MAXCTY.or.kk.lt.1.or.kk.gt.MAXCTZ) then
                  ppp=0
                        call cova3(xx,yy,zz,x1,y1,z1,icut,nst,MAXNST,&
                             c0,it,cc,aa,irot,MAXROT,rotmat,cmax,cov)
                        r(j) = dble(cov)
!if(ix.eq.25.and.iy.eq.11.and.iz.eq.28)write(*,*)'r2: use covtab',ix,iy,iz 
                  else
                  ppp=0
                        r(j) = dble(covtab(ii,jj,kk,icut))
!if(ix.eq.25.and.iy.eq.11.and.iz.eq.28)write(*,*)'r: use covtab',ix,iy,iz 
                  endif
            endif
            rr(j) = r(j)
!if(ix.eq.25.and.iy.eq.11.and.iz.eq.28)then
!write(*,*)'node:',ix,iy,iz,'r',r(j)
!end if
!
! End ``if'' block (true if kriging)
!
         end if
!
! End loop over all of the nearby data
!
      if(softdat(j)) somesoft = .true.
 1    continue


!!
!! If we are doing Markov-Bayes are there are soft data we need to
!! correct some of the covariance values in the kriging matrix:
!!
!      if(imbsim.eq.1.and.somesoft) then
!            in = 0
!            do j=1,na
!                  do i=1,j
!                        in = in + 1
!                        bothsoft = .false.
!                        if(softdat(j).and.softdat(i)) bothsoft = .true.
!!
!! Correct for soft-soft covariance or soft-hard covariance:
!!
!                        if(bothsoft) then
!                              a(in) = a(in)*dble(beez(icut))
!                              if(i.ne.j) a(in) = a(in)*dble(beez(icut))
!                        else
!                              if(softdat(j).or.softdat(i)) a(in) = a(in)*dble(beez(icut))
!                        end if
!                  end do
!!
!! Correct the right hand side for soft-hard covariance:
!!
!                  if(softdat(j)) then
!                        r(j)  = r(j)*dble(beez(icut))
!                        rr(j) = r(j)
!                  end if
!            end do
!      end if

!
! Addition of OK constraint:
!
      if(krig.and.ktype.eq.1) then
            do i=1,na
                  in    = in + 1
                  a(in) = 1.0
            end do
            in      = in + 1
            a(in)   = 0.0
            r(neq)  = 1.0
            rr(neq) = 1.0
      endif
!
! Write out the kriging Matrix if Seriously Debugging:
!
!      if(krig.and.idbg.ge.4) then
!            write(ldbg,101) ix,iy,iz
!            is = 1
!            do i=1,neq
!                  ie = is + i - 1
!                  !write(ldbg,102) i,r(i),(a(j),j=is,ie)
!                  is = is + i
!            end do
! 101        format(/,'Kriging Matrices for Node: ',3i4,&
!                     ' RHS first')
!! 102        format('    r(',i2,') =',f7.4,'  a= ',9(10f7.4))
!      endif

!
! Solve the Kriging System:
!
      if(krig) then
            if(neq.eq.1.and.ktype.eq.0) then
                  s(1)  = r(1) / a(1)
                  ising = 0
            else
!write(*,*)'a:',a(:)
!if(index2.eq.22676)then
!write(*,*)'a:',a(1:in)
!write(*,*)'r:',r(1:na)
!            end if
                  call ksol(1,neq,1,a,r,s,ising)
!if(index2.eq.22676)then
!write(*,*)'s:',s(1:na)
!            end if
            endif
      endif
!write(*,*)ising,s(:)
!
! Write a warning if the matrix is singular:
!
      if(ising.ne.0) then
!if(ix.eq.25.and.iy.eq.11.and.iz.eq.28)then
write(*,*)'ising.ne.0'
!      end if
!            if(idbg.ge.1) then
!                  write(ldbg,*) 'WARNING SISIM: singular matrix'
!                  write(ldbg,*) '              for block',ix,iy,iz
!            endif
            cmean  = 0.0
      !deallocate(r)
      !deallocate(rr)
      !deallocate(s)
      !deallocate(a)
            return
      endif
!
! Write out the kriging Matrix if Seriously Debugging:
!
!      if(krig.and.idbg.ge.4) then
!            do i=1,na
!                  write(ldbg,140) i,s(i)
!            end do
! 140        format(' Kriging weight for data: ',i4,' = ',f8.4)
!      endif

!
! Compute the estimate, the sum of weights, correct for SK, and return:
!
      cmean = 0.0
      sumwt = 0.0
      notnull = 0
!if(index2.eq.22676)then
!write(*,*)'vra:',vra(1:na)
!            end if
      do i=1,na
      if(vra(i).ne.-999.0)then
            cmean = cmean + real(s(i)) * vra(i)
            sumwt = sumwt + real(s(i))
            notnull = notnull + 1
      else
            notnull = notnull - 1 
      end if
      end do
!      if(notnull.lt.1)write(*,*)'WARNING: less than 1 not nulls in vra'
!write(*,*)'ktype',ktype
      if(ktype.eq.0) cmean = cmean + (1.0-sumwt)*gmean
      !deallocate(r)
      !deallocate(rr)
      !deallocate(s)
      !deallocate(a)
      return
      end

      subroutine krige_1D_v2(index2,ix,iy,iz,xx,yy,zz,icut,gmean,&
                       MAXCTX,MAXCTY,MAXCTZ,MAXKR1,MAXDAT,&
                       MXCUT,MAXROT,MAXNST,MAXNOD,MAXCUT,MAXXYZ,MAXKR2,&
                       cmean,mik,&
                       !nclose,ncnode,ktype,close,distance,&
                       nclose,ncnode,ktype,neighbours,distance,&
                       atnode,vr,aclose,x,y,z,nhd,softdat,&
                       vra,cnodex,cnodey,cnodez,&
                       cnodev,cnodet,ivtype,thres,&
                       nctx,ncty,nctz,icnode,ixnode,iynode,iznode,&
                       nx,ny,nz,xmn,ymn,zmn,xsiz,ysiz,zsiz,&
                       !aa,c0,cc,cmax,& 
                       aa,c0,cc,& 
                       it,nst,& 
                       r,rr,s,a,rotmat,&
                       !rotmat,covtab,&
                       NODES,xyzland,coord_ISOMAP_trans,dim,d_tree,query)
                       !r,rr,s,a,rotmat)


!-----------------------------------------------------------------------
!
!            Builds and Solves the SK or OK Kriging System
!            *********************************************
!
! INPUT VARIABLES:
!
!   ix,iy,iz        index of the point currently being simulated
!   xx,yy,zz        location of the point currently being simulated
!   icut            cutoff number to use for either the covariance look
!                     up table or the covariance calculation
!
!
!
! OUTPUT VARIABLES:
!
!   cmean           kriged estimate
!
!
! 
! EXTERNAL REFERENCES: ksol   Gaussian elimination system solution
!
!
! NOTE: 1. the array "aclose" is used to flag those samples which exist
!          at the cutoff currently being kriged.
!
!
!-----------------------------------------------------------------------
      implicit none
      integer,intent(in) ::index2,ix,iy,iz,icut,MAXCTX,MAXCTY,MAXCTZ,MAXKR1
      integer,intent(in)::MAXROT,MAXDAT,MXCUT,MAXNST,MAXNOD,MAXCUT,MAXXYZ,MAXKR2
      integer,intent(in) ::mik,nclose,ncnode,ktype,nhd,ivtype
      integer,intent(in) ::nctx,ncty,nctz
      real*8,intent(in) :: xx,yy,zz
      real,intent(in) :: gmean
      integer,intent(in) :: neighbours(MAXKR1)
      real*8,intent(in) :: distance(MAXKR1) 
      logical,intent(in) :: atnode(MAXDAT)
      real,intent(in) :: vr(MAXDAT,MXCUT),x(MAXDAT),y(MAXDAT),z(MAXDAT)
      integer,intent(inout) :: aclose(MAXKR1),softdat(MAXKR1)
      real,intent(inout) :: vra(MAXKR1)
      real,intent(in) :: cnodex(MAXNOD),cnodey(MAXNOD),cnodez(MAXNOD)
      real,intent(in) :: cnodev(MAXNOD),cnodet(MAXNOD),thres(MAXCUT) 
      integer,intent(in) :: icnode(MAXNOD)
      integer,intent(in) :: iznode(MAXXYZ),iynode(MAXXYZ),ixnode(MAXXYZ)
      integer,intent(in) :: nx,ny,nz
      real*8,intent(in) :: xmn,ymn,zmn,xsiz,ysiz,zsiz
      real, intent(in) :: aa(MAXCUT*MAXNST),c0(MAXCUT),cc(MAXCUT*MAXNST)!,cmax(MAXCUT) 
      integer,intent(in) :: it(MAXCUT*MAXNST),nst(MAXCUT)
      real*8,intent(inout) ::  r(MAXKR1),rr(MAXKR1),s(MAXKR1),a(MAXKR2),rotmat(MAXROT,3,3) 
      !real*8,intent(inout) ::  rotmat(MAXROT,3,3) 
      !real*8,intent(inout) :: a(MAXKR2),rotmat(MAXROT,3,3) 
     
      !integer,intent(in) :: NCLOSEvar,nloop 
      !integer,intent(in):: resultsIdxIndex(NCLOSEvar,nloop)
      integer,intent(in):: NODES,xyzland
      !real*8,intent(in) :: coord_ISOMAP_trans(NODES,xyzland) 
      real*8,intent(in) :: coord_ISOMAP_trans(xyzland,NODES) 
      integer, intent(in):: dim,d_tree,query

      real,intent(inout) :: cmean

      !use      geostat
      !include 'sisim.inc'
      !include 'sisim_new.inc'
      !integer aclose(nclose)
      !real*8, allocatable ::  r(:),rr(:),s(:),a(:)
      logical krig,somesoft,bothsoft,testind
      integer mclose,index,na,neq,i,j,irot,in,j1,iii,ind,ising,k
      integer ix1,iy1,iz1,ix2,iy2,iz2,ix11,iy11,iz11,ix22,iy22,iz22,ii,jj,kk,ppp
      real x1,y1,z1,x2,y2,z2,sumwt,cov
      real*8 dist,covd,sumdiff,sumdiff2
      integer ind1,ind2,notnull
      real cmax
!
! Size of the kriging system:  Some of the data values may be missing
! which would correspond to a constraint interval.  Note that there
! should not be any missing values if the median approximation is being
! considered.  The variable ``krig'' is used
! to flag whether kriging is to be done or if the previous weights are
! to be used.
!


      cmax=0.0

      !write(*,*)'cmax',cmax,c0(icut)
      !allocate(r(MAXKR1))
      !allocate(rr(MAXKR1))
      !allocate(s(MAXKR1))
      !allocate(a(MAXKR2))
      r(:)=0.0
      rr(:)=0.0
      s(:)=0.0
      a(:)=0.0

!      write(*,*)'Inside krige',MAXKR1
!      write(*,*)'neighs:',neighbours(:)
!      write(*,*)'dista:',distance(:)
      !allocate(aclose(MAXKR1))
      !write(*,*)'passed'
      somesoft = .false.
      krig     = .true.
      if(mik.eq.1.and.icut.gt.1) krig = .false.
      if(krig) then
            mclose = 0
!write(*,*)'krige:',index2,nclose
            do i=1,nclose
                  index     =  neighbours(i)
                  !write(*,*)'index',i,index,MAXDAT
                  !write(*,*)'neig',index,i,neighbours(i)
                  !write(*,*)'atnode',index,i,atnode(index)
                  !write(*,*)'vr-size:',shape(vr),icut
                  !write(*,*)'vr',index,i,icut,vr(index,icut)
                  !write(*,*)'Processing:',i,index
!                  if(.not.atnode(index)) then
!                  if(vr(index,icut).ge.0.0) then
!      write(*,*)'Processing-inside:',i,index,mclose
!                        mclose = mclose + 1
!                        !write(*,*)'aclose:',mclose,MAXKR1,index
!                        aclose(mclose) = index
!                        vra(mclose) = vr(index,icut)
!                  endif
!                  endif
            end do
            na  = mclose + ncnode
            na  = ncnode
            neq = na + ktype
      endif
!
! There are no data yet:
!
      irot   = 1 + (icut-1)*MAXNST
      !write(*,*)'inside',na,neq


!! From SGS_LVA_LEVELS
!      !a(1) = 0.0
!      in = 0
!      do j=1,na
!         !ind1=resultsIdxIndex(i,index)
!         !ind1=0
!         ind1=neighbours(j)
!!         if(j.le.mclose) then
!!             index  = aclose(j)
!!             vra(j) = vr(index,icut)
!!             x1     = x(index)
!!             y1     = y(index)
!!             z1     = z(index)
!!             if(index.gt.nhd) softdat(j) = .true.
!
!         do i=1,j
!            !ind2=resultsIdxIndex(j,index)
!            !ind2=0
!            ind2=neighbours(i)
!            sumdiff=0.0
!            do k=1,dim
!!                write(*,*)'ent-cova3-ind1',index2,k,ind1,coord_ISOMAP_trans(k,ind1) 
!!                write(*,*)'ent-cova3-ind2',index2,k,ind2,coord_ISOMAP_trans(k,ind2) 
!                sumdiff=sumdiff + (coord_ISOMAP_trans(k,ind1)-coord_ISOMAP_trans(k,ind2))*(coord_ISOMAP_trans(k,ind1)-coord_ISOMAP_trans(k,ind2)) 
!            end do
!!            write(*,*)'ent-cova3-0',sumdiff
!            dist=dsqrt (sumdiff)
!            !dist=dsqrt ( sum( ( coord_ISOMAP_trans(1:dim,ind1)-coord_ISOMAP_trans(1:dim,ind2) ) **2    )    )
!            !write(*,*)'ent-cova3-0', coord_ISOMAP_trans(1:dim,ind1),coord_ISOMAP_trans(1:dim,ind2)
!            call cova3_1D_sp(real(dist),icut,nst,MAXNST,MAXCUT,c0,it,cc,aa,cmax,cov) 
!!            write(*,*)'ent-cova3-1',dist,dble(cov)
!            in = in + 1
!            write(*,*)'debug',index2,in,MAXKR2
!            a(in) = dble(cov)
!         end do
!!write(*,*)'kkk',j,i,ind1
!         vra(j)=cnodev(j) 
!         if(d_tree/=dim) then
!             dist=distance(j) + sum( ( coord_ISOMAP_trans(d_tree+1:dim,ind1)-coord_ISOMAP_trans(d_tree+1:dim,index) ) **2    ) !dist in tree is only for the first few dims
!         else
!             dist=distance(j)
!         end if
!         dist=dsqrt(dist)
!         call cova3_1D_sp(real(dist),icut,nst,MAXNST,MAXCUT,c0,it,cc,aa,cmax,cov)
!!         write(*,*)'ent-cova3-2',ist,dble(cov)
!         r(j)=dble(cov) 
!
!      end do
!      if(index2.eq.265)then
!write(*,*)'nclose',nclose
!write(*,*)'na',na
!write(*,*)'in',in
!write(*,*)'a',a(:)
!write(*,*)'r',r(:)
!      end if
!!     !do i=1,na
!!     !   a(neq*(i-1)+na+1) = dble(cmax)
!!     !   a(neq*na+i)       = dble(cmax)
!!     !end do
!!
!!
!!! END From SGS_LVA_LEVELS


!      a=0.0
!
! Set up kriging matrices:
!
      in = 0
      j1 = 0
      do 1 j=1,na
            softdat(j) = .false.
!
! Sort out the actual location of point "j"
!
            !if(j.le.mclose) then
            !      index  = aclose(j)
            !      vra(j) = vr(index,icut)
            !      x1     = x(index)
            !      y1     = y(index)
            !      z1     = z(index)
            !      if(index.gt.nhd) softdat(j) = .true.
            !      ind1 = index
            !      !ind1 = neighbours(j-mclose)
            !      if(ind1.eq.0)write(*,*)'WARN ind1 mclose',ind1
            !      write(*,*)'j',j,'le mclose'
            !else
!
! It is a previously simulated node (keep index for table look-up):
!
                  index  = j-mclose
                  x1     = cnodex(index)
                  y1     = cnodey(index)
                  z1     = cnodez(index)
!
! Is this node informed by a hard datum or a soft datum?
!
                  if(cnodet(index).le.0.5) then
                        if(ivtype.eq.0) then
                           vra(j) = 0.0
                           if(int(cnodev(index)+0.5).eq.int(thres(icut)+0.5)) vra(j) = 1.0
                        else
                           vra(j) = 1.0
                           if(cnodev(index).gt.thres(icut)) vra(j) = 0.0
                        end if
                  else
                        iii    = int(cnodet(index)+0.5)
                        vra(j) = vr(iii,icut)
                        softdat(j) = .true.
                  end if
                  call getindx(nx,xmn,xsiz,dble(cnodex(index)),ix11,testind)
                  call getindx(ny,ymn,ysiz,dble(cnodey(index)),iy11,testind)
                  call getindx(nz,zmn,zsiz,dble(cnodez(index)),iz11,testind)
!if(ix.eq.25.and.iy.eq.11.and.iz.eq.28)then
!write(*,*)'node-n:',ix,iy,iz,nx,ny,nz
!write(*,*)'node-mn:',ix,iy,iz,xmn,ymn,zmn
!write(*,*)'node-siz:',ix,iy,iz,xsiz,ysiz,zsiz
!write(*,*)'node-cnode:',ix,iy,iz,dble(cnodex(index)),dble(cnodey(index)),dble(cnodez(index))
!write(*,*)'node-iout:',ix,iy,iz,ix11,iy11,iz11
!          end if
                  !ind    = icnode(index)
                  !ix1    = ix + (int(ixnode(ind))-nctx-1)
                  !iy1    = iy + (int(iynode(ind))-ncty-1)
                  !iz1    = iz + (int(iznode(ind))-nctz-1)

                  ix1    = ix + (int(ix11)-nctx-1)
                  iy1    = iy + (int(iy11)-ncty-1)
                  iz1    = iz + (int(iz11)-nctz-1)
                  ind1 = neighbours(j-mclose)
                  if(ind1.eq.0)write(*,*)'WARN ind1 else',ind1
                  !write(*,*)'j',j,'gt mclose'
            !end if
!
! Only set up the matrix and the RHS if kriging:
!
            if(krig) then
               do 2 i=1,j
!
! Sort out the actual location of point "i"
!
                  !if(i.le.mclose) then
                  !      index  = aclose(i)
                  !      x2     = x(index)
                  !      y2     = y(index)
                  !      z2     = z(index)
                  !      ind2 = index
                  !      !ind2 = neighbours(i-mclose)
                  !      if(ind2.eq.0)write(*,*)'WARN ind2 mclose',ind2
                  !      write(*,*)'i',i,'le mclose'
                  !else
!
! It is a previously simulated node (keep index for table look-up):
!
                        index  = i-mclose
                        x2     = cnodex(index)
                        y2     = cnodey(index)
                        z2     = cnodez(index)
!                        ind    = icnode(index)
!                        ix2    = ix + (int(ixnode(ind))-nctx-1)
!                        iy2    = iy + (int(iynode(ind))-ncty-1)
!                        iz2    = iz + (int(iznode(ind))-nctz-1)
                        call getindx(nx,xmn,xsiz,dble(cnodex(index)),ix22,testind)
                        call getindx(ny,ymn,ysiz,dble(cnodey(index)),iy22,testind)
                        call getindx(nz,zmn,zsiz,dble(cnodez(index)),iz22,testind)
                        ix2    = ix + (int(ix22)-nctx-1)
                        iy2    = iy + (int(iy22)-ncty-1)
                        iz2    = iz + (int(iz22)-nctz-1)
                        
                        ind2 = neighbours(i-mclose)
                        if(ind2.eq.0)write(*,*)'WARN ind2 else',ind2
                        !write(*,*)'i',i,'gt mclose'

                  !endif
!
! Now, get the covariance value:
!
                  in = in + 1
!
! Decide whether or not to use the covariance look-up table:
!
                  !if(j.le.mclose.or.i.le.mclose) then
                  ppp=0
                        !call cova3(x1,y1,z1,x2,y2,z2,icut,nst,MAXNST,&
                        !     c0,it,cc,aa,irot,MAXROT,rotmat,cmax,cov)
                        !a(in) = dble(cov)


                  sumdiff=0.0
                  do k=1,dim
                      sumdiff=sumdiff + (coord_ISOMAP_trans(k,ind1)-coord_ISOMAP_trans(k,ind2))*(coord_ISOMAP_trans(k,ind1)-coord_ISOMAP_trans(k,ind2)) 
                  end do
                  dist = dsqrt (sumdiff)
                  call cova3_1D_sp(real(dist),icut,nst,MAXNST,MAXCUT,c0,it,cc,aa,cmax,cov) 
                  a(in) = dble(cov)

                  !write(*,*)'deb',ind1,ind2,dist,sumdiff,a(in)
                  !write(*,*)'icut',icut
                  !write(*,*)'nst',nst
                  !write(*,*)'MAXNST',MAXNST
                  !write(*,*)'MAXCUT',MAXCUT
                  !write(*,*)'c0',c0
                  !write(*,*)'it',it
                  !write(*,*)'cc',cc
                  !write(*,*)'aa',aa
                  !write(*,*)'cmax',cmax


                 !else
!!
!! Try to use the covariance look-up (if the distance is in range):
!!
!                        ii = nctx + 1 + (ix1 - ix2)
!                        jj = ncty + 1 + (iy1 - iy2)
!                        kk = nctz + 1 + (iz1 - iz2)
!!if(ix.eq.25.and.iy.eq.11.and.iz.eq.28)write(*,*)'box1',ix,iy,iz,ii,jj,kk
!!if(ix.eq.25.and.iy.eq.11.and.iz.eq.28)write(*,*)'box2',ix,iy,iz,nctx,ncty,nctz
!!if(ix.eq.25.and.iy.eq.11.and.iz.eq.28)write(*,*)'box3',ix,iy,iz,MAXCTX,MAXCTY,MAXCTZ
!                        if(ii.lt.1.or.ii.gt.MAXCTX.or.jj.lt.1.or.jj.gt.MAXCTY.or.kk.lt.1.or.kk.gt.MAXCTZ) then
!                  ppp=0
!                              call cova3(x1,y1,z1,x2,y2,z2,icut,nst,&
!                                         MAXNST,c0,it,cc,aa,irot,MAXROT,&
!                                         rotmat,cmax,cov)
!                              a(in) = dble(cov)
!!if(ix.eq.25.and.iy.eq.11.and.iz.eq.28)write(*,*)'a: use cova3',ix,iy,iz 
!                        else
!                  ppp=0
!                              a(in) = dble(covtab(ii,jj,kk,icut))
!!if(ix.eq.25.and.iy.eq.11.and.iz.eq.28)write(*,*)'a: use covtab',ix,iy,iz 
!                        endif
!                  endif
!if(ix.eq.25.and.iy.eq.11.and.iz.eq.28)then
!write(*,*)'node:',ix,iy,iz,'a',a(in)
!end if
 2          continue

!      do i=1,na
!         a(neq*(i-1)+na+1) = dble(cmax)
!         a(neq*na+i)       = dble(cmax)
!      end do






!
! Get the RHS value (possibly with covariance look-up table):
!
            !if(j.le.mclose) then
            !      ppp=0
            !      call cova3(xx,yy,zz,x1,y1,z1,icut,nst,MAXNST,&
            !                 c0,it,cc,aa,irot,MAXROT,rotmat,cmax,cov)
            !      r(j) = dble(cov)
!if(ix.eq.25!.and.iy.eq.11.and.iz.eq.28)write(*,*)'r: use cova3',ix,iy,iz 
            !else
!
! Try to use! the covariance look-up (if the distance is in range):
!
            !      ii = nctx + 1 + (ix - ix1)
            !      jj = ncty + 1 + (iy - iy1)
            !      kk = nctz + 1 + (iz - iz1)
            !      if(ii.lt.1.or.ii.gt.MAXCTX.or.jj.lt.1.or.jj.gt.MAXCTY.or.kk.lt.1.or.kk.gt.MAXCTZ) then
            !      ppp=0
            !            call cova3(xx,yy,zz,x1,y1,z1,icut,nst,MAXNST,&
            !                 c0,it,cc,aa,irot,MAXROT,rotmat,cmax,cov)
            !            r(j) = dble(cov)
!if(ix.eq.25!.and.iy.eq.11.and.iz.eq.28)write(*,*)'r2: use covtab',ix,iy,iz 
            !      else
            !      ppp=0
            !            r(j) = dble(covtab(ii,jj,kk,icut))
!if(ix.eq.25!.and.iy.eq.11.and.iz.eq.28)write(*,*)'r: use covtab',ix,iy,iz 
            !      endif
            !endif


            if(d_tree/=dim) then
                dist=distance(j) + sum( (coord_ISOMAP_trans(d_tree+1:dim,ind1)-coord_ISOMAP_trans(d_tree+1:dim,index2) ) **2    ) !dist in tree is only for the first few dims
            else
                dist=distance(j)
            end if
            dist=dsqrt(dist)
            call cova3_1D_sp(real(dist),icut,nst,MAXNST,MAXCUT,c0,it,cc,aa,cmax,cov)
            r(j)=dble(cov) 

                  !write(*,*)'r-deb',ind1,ind2,dist,r(j)
                  !write(*,*)'r-icut',icut
                  !write(*,*)'r-nst',nst
                  !write(*,*)'r-MAXNST',MAXNST
                  !write(*,*)'r-MAXCUT',MAXCUT
                  !write(*,*)'r-c0',c0
                  !write(*,*)'r-it',it
                  !write(*,*)'r-cc',cc
                  !write(*,*)'r-aa',aa
                  !write(*,*)'r-cmax',cmax



            rr(j) = r(j)
!if(ix.eq.25.and.iy.eq.11.and.iz.eq.28)then
!write(*,*)'node:',ix,iy,iz,'r',r(j)
!end if
!
! End ``if'' block (true if kriging)
!
         end if
!
! End loop over all of the nearby data
!
      if(softdat(j)) somesoft = .true.
 1    continue


!!
!! If we are doing Markov-Bayes are there are soft data we need to
!! correct some of the covariance values in the kriging matrix:
!!
!      if(imbsim.eq.1.and.somesoft) then
!            in = 0
!            do j=1,na
!                  do i=1,j
!                        in = in + 1
!                        bothsoft = .false.
!                        if(softdat(j).and.softdat(i)) bothsoft = .true.
!!
!! Correct for soft-soft covariance or soft-hard covariance:
!!
!                        if(bothsoft) then
!                              a(in) = a(in)*dble(beez(icut))
!                              if(i.ne.j) a(in) = a(in)*dble(beez(icut))
!                        else
!                              if(softdat(j).or.softdat(i)) a(in) = a(in)*dble(beez(icut))
!                        end if
!                  end do
!!
!! Correct the right hand side for soft-hard covariance:
!!
!                  if(softdat(j)) then
!                        r(j)  = r(j)*dble(beez(icut))
!                        rr(j) = r(j)
!                  end if
!            end do
!      end if

!
! Addition of OK constraint:
!
      if(krig.and.ktype.eq.1) then
            do i=1,na
                  in    = in + 1
                  a(in) = 1.0
            end do
            in      = in + 1
            a(in)   = 0.0
            r(neq)  = 1.0
            rr(neq) = 1.0
      endif
!
! Write out the kriging Matrix if Seriously Debugging:
!
!      if(krig.and.idbg.ge.4) then
!            write(ldbg,101) ix,iy,iz
!            is = 1
!            do i=1,neq
!                  ie = is + i - 1
!                  !write(ldbg,102) i,r(i),(a(j),j=is,ie)
!                  is = is + i
!            end do
! 101        format(/,'Kriging Matrices for Node: ',3i4,&
!                     ' RHS first')
!! 102        format('    r(',i2,') =',f7.4,'  a= ',9(10f7.4))
!      endif

!
! Solve the Kriging System:
!
      if(krig) then
            if(neq.eq.1.and.ktype.eq.0) then
                  s(1)  = r(1) / a(1)
                  ising = 0
            else
!write(*,*)'a:',a(:)
!if(index2.eq.22676)then
!write(*,*)'a:',a(1:in)
!write(*,*)'r:',r(1:na)
!            end if

!write(*,*)'neq',neq
if(index2.eq.query)then
write(*,*)'a',a(:)
write(*,*)'r',r(:)
end if
                    call ksol(1,neq,1,MAXKR1,MAXKR2,a,r,s,ising)
                    !call ksolopt(1,neq,1,MAXKR1,MAXKR2,a,r,s,ising)
if(index2.eq.query)then
write(*,*)'s',s(:)
end if
!if(index2.eq.22676)then
!write(*,*)'s:',s(1:na)
!            end if
            endif
      endif
!write(*,*)ising,s(:)
!
! Write a warning if the matrix is singular:
!
      if(ising.ne.0) then
!if(ix.eq.25.and.iy.eq.11.and.iz.eq.28)then
write(*,*)'ising.ne.0'
!      end if
!            if(idbg.ge.1) then
!                  write(ldbg,*) 'WARNING SISIM: singular matrix'
!                  write(ldbg,*) '              for block',ix,iy,iz
!            endif
            cmean  = 0.0
      !deallocate(r)
      !deallocate(rr)
      !deallocate(s)
      !deallocate(a)
            return
      endif
!
! Write out the kriging Matrix if Seriously Debugging:
!
!      if(krig.and.idbg.ge.4) then
!            do i=1,na
!                  write(ldbg,140) i,s(i)
!            end do
! 140        format(' Kriging weight for data: ',i4,' = ',f8.4)
!      endif

!
! Compute the estimate, the sum of weights, correct for SK, and return:
!
      cmean = 0.0
      sumwt = 0.0
      notnull = 0
!if(index2.eq.22676)then
!write(*,*)'vra:',vra(1:na)
!            end if
if(index2.eq.query)then
write(*,*)'vra',vra(:)
end if
      do i=1,na
!write(*,*)'vra',vra(:)
!write(*,*)'s',s(:)
      if(vra(i).ne.-999.0)then
            cmean = cmean + real(s(i)) * vra(i)
if(index2.eq.query)write(*,*)'cmean',cmean,real(s(i)),vra(i)
            sumwt = sumwt + real(s(i))
            notnull = notnull + 1
      else
            notnull = notnull - 1 
      end if
      end do
!      if(notnull.lt.1)write(*,*)'WARNING: less than 1 not nulls in vra'
!write(*,*)'ktype',ktype
      if(ktype.eq.0)then
          cmean = cmean + (1.0-sumwt)*gmean
if(index2.eq.query)write(*,*)'cmean',cmean,(1.0-sumwt),sumwt,gmean
      end if
      !deallocate(r)
      !deallocate(rr)
      !deallocate(s)
      !deallocate(a)
      return
      end





      !subroutine krige(ix,iy,iz,xx,yy,zz,icut,MAXKR1)
      subroutine krige(ix,iy,iz,xx,yy,zz,icut,gmean,&
                       MAXCTX,MAXCTY,MAXCTZ,MAXKR1,MAXDAT,&
                       MXCUT,MAXROT,MAXNST,MAXNOD,MAXCUT,MAXXYZ,MAXKR2,&
                       cmean,mik,&
                       nclose,ncnode,ktype,close,&
                       atnode,vr,aclose,x,y,z,nhd,softdat,&
                       vra,cnodex,cnodey,cnodez,&
                       cnodev,cnodet,ivtype,thres,&
                       nctx,ncty,nctz,icnode,ixnode,iynode,iznode,&
                       nx,ny,nz,xmn,ymn,zmn,xsiz,ysiz,zsiz,&
                       aa,c0,cc,cmax,& 
                       it,nst,& 
                       r,rr,s,a,rotmat,covtab)
                       !r,rr,s,a,rotmat)


!-----------------------------------------------------------------------
!
!            Builds and Solves the SK or OK Kriging System
!            *********************************************
!
! INPUT VARIABLES:
!
!   ix,iy,iz        index of the point currently being simulated
!   xx,yy,zz        location of the point currently being simulated
!   icut            cutoff number to use for either the covariance look
!                     up table or the covariance calculation
!
!
!
! OUTPUT VARIABLES:
!
!   cmean           kriged estimate
!
!
! 
! EXTERNAL REFERENCES: ksol   Gaussian elimination system solution
!
!
! NOTE: 1. the array "aclose" is used to flag those samples which exist
!          at the cutoff currently being kriged.
!
!
!-----------------------------------------------------------------------
      implicit none
      integer,intent(in) ::ix,iy,iz,icut,MAXCTX,MAXCTY,MAXCTZ,MAXKR1
      integer,intent(in) ::MAXROT,MAXDAT,MXCUT,MAXNST,MAXNOD,MAXCUT,MAXXYZ,MAXKR2
      integer,intent(in) ::mik,nclose,ncnode,ktype,nhd,ivtype
      integer,intent(in) ::nctx,ncty,nctz
      real*8,intent(in) :: xx,yy,zz
      real,intent(in) :: gmean
      real,intent(in) :: close(MAXDAT) 
      logical,intent(in) :: atnode(MAXDAT)
      real,intent(in) :: vr(MAXDAT,MXCUT),x(MAXDAT),y(MAXDAT),z(MAXDAT)
      integer,intent(inout) :: aclose(MAXKR1),softdat(MAXKR1)
      real,intent(inout) :: vra(MAXKR1)
      real,intent(in) :: cnodex(MAXNOD),cnodey(MAXNOD),cnodez(MAXNOD)
      real,intent(in) :: cnodev(MAXNOD),cnodet(MAXNOD),thres(MAXCUT) 
      integer,intent(in) :: icnode(MAXNOD)
      integer,intent(in) :: iznode(MAXXYZ),iynode(MAXXYZ),ixnode(MAXXYZ)
      integer,intent(in) :: nx,ny,nz
      real*8,intent(in) :: xmn,ymn,zmn,xsiz,ysiz,zsiz
      real, intent(in) :: aa(MAXCUT*MAXNST),c0(MAXCUT),cc(MAXCUT*MAXNST),cmax(MAXCUT) 
      real,intent(in) ::  it(MAXCUT*MAXNST),nst(MAXCUT),covtab(MAXCTX,MAXCTY,MAXCTZ,MAXCUT) 
      real*8,intent(inout) ::  r(MAXKR1),rr(MAXKR1),s(MAXKR1),a(MAXKR2),rotmat(MAXROT,3,3) 
      !real*8,intent(inout) :: a(MAXKR2),rotmat(MAXROT,3,3) 


      real,intent(inout) :: cmean

      !use      geostat
      !include 'sisim.inc'
      !include 'sisim_new.inc'
      !integer aclose(nclose)
      logical krig,somesoft,bothsoft,testind
      integer mclose,index,na,neq,i,j,irot,in,j1,iii,ind,ising
      integer ix1,iy1,iz1,ix2,iy2,iz2,ix11,iy11,iz11,ix22,iy22,iz22,ii,jj,kk,ppp
      real x1,y1,z1,x2,y2,z2,sumwt,cov 
!
! Size of the kriging system:  Some of the data values may be missing
! which would correspond to a constraint interval.  Note that there
! should not be any missing values if the median approximation is being
! considered.  The variable ``krig'' is used
! to flag whether kriging is to be done or if the previous weights are
! to be used.
!

      !write(*,*)'Inside krige',MAXKR1
      !allocate(aclose(MAXKR1))
      !write(*,*)'passed'
      somesoft = .false.
      krig     = .true.
      if(mik.eq.1.and.icut.gt.1) krig = .false.
      if(krig) then
            mclose = 0
            do i=1,nclose
                  index     =  int(close(i))
      !write(*,*)'Processing:',i,index
                  if(.not.atnode(index).and.vr(index,icut).ge.0.0) then
      !write(*,*)'Processing-inside:',i,index,mclose
                        mclose = mclose + 1
                        !write(*,*)'aclose:',mclose,MAXKR1,index
                        aclose(mclose) = index
                  endif
            end do
            na  = mclose + ncnode
            neq = na + ktype
      endif
!
! There are no data yet:
!
      irot   = 1 + (icut-1)*MAXNST
!
! Set up kriging matrices:
!
      in = 0
      j1 = 0
      do 1 j=1,na
            softdat(j) = .false.
!
! Sort out the actual location of point "j"
!
            if(j.le.mclose) then
                  index  = aclose(j)
                  vra(j) = vr(index,icut)
                  x1     = x(index)
                  y1     = y(index)
                  z1     = z(index)
                  if(index.gt.nhd) softdat(j) = .true.
            else
!
! It is a previously simulated node (keep index for table look-up):
!
                  index  = j-mclose
                  x1     = cnodex(index)
                  y1     = cnodey(index)
                  z1     = cnodez(index)
!
! Is this node informed by a hard datum or a soft datum?
!
                  if(cnodet(index).le.0.5) then
                        if(ivtype.eq.0) then
                           vra(j) = 0.0
                           if(int(cnodev(index)+0.5).eq.int(thres(icut)+0.5)) vra(j) = 1.0
                        else
                           vra(j) = 1.0
                           if(cnodev(index).gt.thres(icut)) vra(j) = 0.0
                        end if
                  else
                        iii    = int(cnodet(index)+0.5)
                        vra(j) = vr(iii,icut)
                        softdat(j) = .true.
                  end if
                  call getindx(nx,xmn,xsiz,dble(cnodex(index)),ix11,testind)
                  call getindx(ny,ymn,ysiz,dble(cnodey(index)),iy11,testind)
                  call getindx(nz,zmn,zsiz,dble(cnodez(index)),iz11,testind)
if(ix.eq.25.and.iy.eq.11.and.iz.eq.28)then
write(*,*)'node-n:',ix,iy,iz,nx,ny,nz
write(*,*)'node-mn:',ix,iy,iz,xmn,ymn,zmn
write(*,*)'node-siz:',ix,iy,iz,xsiz,ysiz,zsiz
write(*,*)'node-cnode:',ix,iy,iz,dble(cnodex(index)),dble(cnodey(index)),dble(cnodez(index))
write(*,*)'node-iout:',ix,iy,iz,ix11,iy11,iz11
          end if
                  !ind    = icnode(index)
                  !ix1    = ix + (int(ixnode(ind))-nctx-1)
                  !iy1    = iy + (int(iynode(ind))-ncty-1)
                  !iz1    = iz + (int(iznode(ind))-nctz-1)

                  ix1    = ix + (int(ix11)-nctx-1)
                  iy1    = iy + (int(iy11)-ncty-1)
                  iz1    = iz + (int(iz11)-nctz-1)
            end if
!
! Only set up the matrix and the RHS if kriging:
!
            if(krig) then
               do 2 i=1,j
!
! Sort out the actual location of point "i"
!
                  if(i.le.mclose) then
                        index  = aclose(i)
                        x2     = x(index)
                        y2     = y(index)
                        z2     = z(index)
                  else
!
! It is a previously simulated node (keep index for table look-up):
!
                        index  = i-mclose
                        x2     = cnodex(index)
                        y2     = cnodey(index)
                        z2     = cnodez(index)
!                        ind    = icnode(index)
!                        ix2    = ix + (int(ixnode(ind))-nctx-1)
!                        iy2    = iy + (int(iynode(ind))-ncty-1)
!                        iz2    = iz + (int(iznode(ind))-nctz-1)
                        call getindx(nx,xmn,xsiz,dble(cnodex(index)),ix22,testind)
                        call getindx(ny,ymn,ysiz,dble(cnodey(index)),iy22,testind)
                        call getindx(nz,zmn,zsiz,dble(cnodez(index)),iz22,testind)
                        ix2    = ix + (int(ix22)-nctx-1)
                        iy2    = iy + (int(iy22)-ncty-1)
                        iz2    = iz + (int(iz22)-nctz-1)
                  endif
!
! Now, get the covariance value:
!
                  in = in + 1
!
! Decide whether or not to use the covariance look-up table:
!
if(ix.eq.25.and.iy.eq.11.and.iz.eq.28)then
write(*,*)'node: points',ix,iy,iz
write(*,*)'node: coord1',x1,y1,z1
write(*,*)'node: coord2',x2,y2,z2
          end if
                  if(j.le.mclose.or.i.le.mclose) then
                  ppp=0
                        call cova3(x1,y1,z1,x2,y2,z2,icut,nst,MAXNST,&
                             c0,it,cc,aa,irot,MAXROT,rotmat,cmax,cov)
                        a(in) = dble(cov)
if(ix.eq.25.and.iy.eq.11.and.iz.eq.28)write(*,*)'a2: use cova3',ix,iy,iz 
                  else
!
! Try to use the covariance look-up (if the distance is in range):
!
                        ii = nctx + 1 + (ix1 - ix2)
                        jj = ncty + 1 + (iy1 - iy2)
                        kk = nctz + 1 + (iz1 - iz2)
if(ix.eq.25.and.iy.eq.11.and.iz.eq.28)write(*,*)'box1',ix,iy,iz,ii,jj,kk
if(ix.eq.25.and.iy.eq.11.and.iz.eq.28)write(*,*)'box2',ix,iy,iz,nctx,ncty,nctz
if(ix.eq.25.and.iy.eq.11.and.iz.eq.28)write(*,*)'box3',ix,iy,iz,MAXCTX,MAXCTY,MAXCTZ
                        if(ii.lt.1.or.ii.gt.MAXCTX.or.jj.lt.1.or.jj.gt.MAXCTY.or.kk.lt.1.or.kk.gt.MAXCTZ) then
                  ppp=0
                              call cova3(x1,y1,z1,x2,y2,z2,icut,nst,&
                                         MAXNST,c0,it,cc,aa,irot,MAXROT,&
                                         rotmat,cmax,cov)
                              a(in) = dble(cov)
if(ix.eq.25.and.iy.eq.11.and.iz.eq.28)write(*,*)'a: use cova3',ix,iy,iz 
                        else
                  ppp=0
                              a(in) = dble(covtab(ii,jj,kk,icut))
if(ix.eq.25.and.iy.eq.11.and.iz.eq.28)write(*,*)'a: use covtab',ix,iy,iz 
                        endif
                  endif
if(ix.eq.25.and.iy.eq.11.and.iz.eq.28)then
write(*,*)'node:',ix,iy,iz,'a',a(in)
end if
 2          continue
!
! Get the RHS value (possibly with covariance look-up table):
!
            if(j.le.mclose) then
                  ppp=0
                  call cova3(xx,yy,zz,x1,y1,z1,icut,nst,MAXNST,&
                             c0,it,cc,aa,irot,MAXROT,rotmat,cmax,cov)
                  r(j) = dble(cov)
if(ix.eq.25.and.iy.eq.11.and.iz.eq.28)write(*,*)'r: use cova3',ix,iy,iz 
            else
!
! Try to use the covariance look-up (if the distance is in range):
!
                  ii = nctx + 1 + (ix - ix1)
                  jj = ncty + 1 + (iy - iy1)
                  kk = nctz + 1 + (iz - iz1)
                  if(ii.lt.1.or.ii.gt.MAXCTX.or.jj.lt.1.or.jj.gt.MAXCTY.or.kk.lt.1.or.kk.gt.MAXCTZ) then
                  ppp=0
                        call cova3(xx,yy,zz,x1,y1,z1,icut,nst,MAXNST,&
                             c0,it,cc,aa,irot,MAXROT,rotmat,cmax,cov)
                        r(j) = dble(cov)
if(ix.eq.25.and.iy.eq.11.and.iz.eq.28)write(*,*)'r2: use covtab',ix,iy,iz 
                  else
                  ppp=0
                        r(j) = dble(covtab(ii,jj,kk,icut))
if(ix.eq.25.and.iy.eq.11.and.iz.eq.28)write(*,*)'r: use covtab',ix,iy,iz 
                  endif
            endif
            rr(j) = r(j)
if(ix.eq.25.and.iy.eq.11.and.iz.eq.28)then
write(*,*)'node:',ix,iy,iz,'r',r(j)
end if
!
! End ``if'' block (true if kriging)
!
         end if
!
! End loop over all of the nearby data
!
      if(softdat(j)) somesoft = .true.
 1    continue
!!
!! If we are doing Markov-Bayes are there are soft data we need to
!! correct some of the covariance values in the kriging matrix:
!!
!      if(imbsim.eq.1.and.somesoft) then
!            in = 0
!            do j=1,na
!                  do i=1,j
!                        in = in + 1
!                        bothsoft = .false.
!                        if(softdat(j).and.softdat(i)) bothsoft = .true.
!!
!! Correct for soft-soft covariance or soft-hard covariance:
!!
!                        if(bothsoft) then
!                              a(in) = a(in)*dble(beez(icut))
!                              if(i.ne.j) a(in) = a(in)*dble(beez(icut))
!                        else
!                              if(softdat(j).or.softdat(i)) a(in) = a(in)*dble(beez(icut))
!                        end if
!                  end do
!!
!! Correct the right hand side for soft-hard covariance:
!!
!                  if(softdat(j)) then
!                        r(j)  = r(j)*dble(beez(icut))
!                        rr(j) = r(j)
!                  end if
!            end do
!      end if

!
! Addition of OK constraint:
!
      if(krig.and.ktype.eq.1) then
            do i=1,na
                  in    = in + 1
                  a(in) = 1.0
            end do
            in      = in + 1
            a(in)   = 0.0
            r(neq)  = 1.0
            rr(neq) = 1.0
      endif
!
! Write out the kriging Matrix if Seriously Debugging:
!
!      if(krig.and.idbg.ge.4) then
!            write(ldbg,101) ix,iy,iz
!            is = 1
!            do i=1,neq
!                  ie = is + i - 1
!                  !write(ldbg,102) i,r(i),(a(j),j=is,ie)
!                  is = is + i
!            end do
! 101        format(/,'Kriging Matrices for Node: ',3i4,&
!                     ' RHS first')
!! 102        format('    r(',i2,') =',f7.4,'  a= ',9(10f7.4))
!      endif

!
! Solve the Kriging System:
!
      if(krig) then
            if(neq.eq.1.and.ktype.eq.0) then
                  s(1)  = r(1) / a(1)
                  ising = 0
            else
!write(*,*)'a:',a(:)
!write(*,*)'r:',r(:)
                  call ksol(1,neq,1,a,r,s,ising)
            endif
      endif
!write(*,*)ising,s(:)
!
! Write a warning if the matrix is singular:
!
      if(ising.ne.0) then
if(ix.eq.25.and.iy.eq.11.and.iz.eq.28)then
write(*,*)'ising.ne.0'
      end if
!            if(idbg.ge.1) then
!                  write(ldbg,*) 'WARNING SISIM: singular matrix'
!                  write(ldbg,*) '              for block',ix,iy,iz
!            endif
            cmean  = 0.0
            return
      endif
!
! Write out the kriging Matrix if Seriously Debugging:
!
!      if(krig.and.idbg.ge.4) then
!            do i=1,na
!                  write(ldbg,140) i,s(i)
!            end do
! 140        format(' Kriging weight for data: ',i4,' = ',f8.4)
!      endif

!
! Compute the estimate, the sum of weights, correct for SK, and return:
!
      cmean = 0.0
      sumwt = 0.0
      do i=1,na
            cmean = cmean + real(s(i)) * vra(i)
            sumwt = sumwt + real(s(i))
      end do
      if(ktype.eq.0) cmean = cmean + (1.0-sumwt)*gmean
      return
      end

! Editado por RGE











      subroutine fix_angles ( angles )
        parameter(DEG2RAD=3.141592654/180.0,EPSLON=1.e-20)
    
        real(kind=8), dimension ( 3 ), intent (inout) :: angles
        !1st - angle between the major axis of anisotropy and the
        !      E-W axis. Note: Counter clockwise is positive.
        !2nd - angle between major axis and the horizontal plane.
        !      (The dip of the ellipsoid measured positive down)
        !3rd - Angle of rotation of minor axis about the major axis
        !      of the ellipsoid.
        if(angles(1) >= 0.0 .and. angles(1) < 270.0) then
            angles(1) = ( 90.0 - angles(1)) * DEG2RAD ; else
            angles(1) = (450.0 - angles(1)) * DEG2RAD ; endif
        angles(2) = -1.0 * angles(2) * DEG2RAD
        angles(3) =  1.0 * angles(3) * DEG2RAD
        Return
        end subroutine fix_angles
        
!
! Fin Editado por RGE
!


      subroutine makepar
!-----------------------------------------------------------------------
!
!                      Write a Parameter File
!                      **********************
!
!
!
!-----------------------------------------------------------------------
      lun = 99
      open(lun,file='sisim_lva.par',status='UNKNOWN')
      write(lun,10)
 10   format('                  Parameters for SISIM_LVA',/, &
             '                  ********************',/,/,&
             'START OF PARAMETERS:')

      write(lun,11)
 11   format('0                             ', &
             '-1=continuous(cdf), 0=categorical(pdf)')
      write(lun,12)
 12   format('4                             ', &
             '-number thresholds/categories')
      write(lun,13)
 13   format('1   2   3   4   ', &
             '-   thresholds / categories')
      write(lun,14)
 14   format('0.009  0.009  0.012  0.97  ', &
             '-   global cdf / pdf')
      write(lun,15)
 15   format('jeronimo_corto.txt           ', &
             '-file with data')
      write(lun,16)
 16   format('1   2   3   4                 ', &
             '-   columns for X,Y,Z, and variable')
      write(lun,17)
 17   format('direct.ik                     ', &
             '-file with soft indicator input')
      write(lun,18)
 18   format('1   2   0   3 4 5 6 7         ', &
             '-   columns for X,Y,Z, and indicators')
      write(lun,19)
 19   format('0                             ', &
             '-   Markov-Bayes simulation (0=no,1=yes)')
      write(lun,20)
 20   format('0.61  0.54  0.56  0.53  0.29  ', &
             '-      calibration B(z) values')
      write(lun,21)
 21   format('-1.0e21    1.0e21             ', &
             '-trimming limits')
      write(lun,22)
 22   format('0.0   30.0                    ', &
             '-minimum and maximum data value')
      write(lun,23)
 23   format('1      0.0                    ', &
             '-   lower tail option and parameter')
      write(lun,24)
 24   format('1      1.0                    ', &
             '-   middle     option and parameter')
      write(lun,25)
 25   format('1     30.0                    ', &
             '-   upper tail option and parameter')
      write(lun,26)
 26   format('cluster.dat                   ', &
             '-   file with tabulated values')
      write(lun,27)
 27   format('3   0                         ', &
             '-      columns for variable, weight')
      write(lun,28)
 28   format('0                             ', &
             '-debugging level: 0,1,2,3')
      write(lun,29)
 29   format('sisim.dbg                     ', &
             '-file for debugging output')
      write(lun,30)
 30   format('sisim_lva.out                     ', &
             '-file for simulation output')
      write(lun,31)
 31   format('1                             ', &
             '-number of realizations')
      write(lun,32)
 32   format('50   0.5    1.0               ', &
             '-nx,xmn,xsiz')
      write(lun,33)
 33   format('50   0.5    1.0               ', &
             '-ny,ymn,ysiz')
      write(lun,34)
 34   format('1    1.0   10.0               ', &
             '-nz,zmn,zsiz')
      write(lun,35)
 35   format('69069                         ', &
             '-random number seed')
      write(lun,36)
 36   format('12                            ', &
             '-maximum original data  for each kriging')
      write(lun,37)
 37   format('12                            ', &
             '-maximum previous nodes for each kriging')
      write(lun,38)
 38   format('1                             ', &
             '-maximum soft indicator nodes for kriging')
      write(lun,39)
 39   format('1                             ', &
             '-assign data to nodes? (0=no,1=yes)')
      write(lun,40)
 40   format('0     3                       ', &
             '-multiple grid search? (0=no,1=yes),num')
      write(lun,41)
 41   format('0                             ', &
             '-maximum per octant    (0=not used)')
      write(lun,42)
 42   format('20.0  20.0  20.0              ', &
             '-maximum search radii')
      write(lun,43)
 43   format(' 0.0   0.0   0.0              ', &
             '-angles for search ellipsoid')
      write(lun,44)
 44   format('51    51    11                ', &
             '-size of covariance lookup table')
      write(lun,47)
 47   format('0    2.5                      ', &
             '-0=full IK, 1=median approx. (cutoff)')
      write(lun,48)
 48   format('0                             ', &
             '-0=SK, 1=OK')
      write(lun,49)
 49   format('1    0.15                     ', &
             '-One   nst, nugget effect')
      write(lun,50)
 50   format('1    0.85 0.0   0.0   0.0     ', &
             '-      it,cc,ang1,ang2,ang3')
      write(lun,51)
 51   format('         10.0  10.0  10.0     ', &
             '-      a_hmax, a_hmin, a_vert')
      write(lun,52)
 52   format('1    0.10                     ', &
             '-Two   nst, nugget effect')
      write(lun,53)
 53   format('1    0.90 0.0   0.0   0.0     ', &
             '-      it,cc,ang1,ang2,ang3')
      write(lun,54)
 54   format('         10.0  10.0  10.0     ', &
             '-      a_hmax, a_hmin, a_vert')
      write(lun,55)
 55   format('1    0.10                     ', &
             '-Three nst, nugget effect')
      write(lun,56)
 56   format('1    0.90 0.0   0.0   0.0     ', &
             '-      it,cc,ang1,ang2,ang3')
      write(lun,57)
 57   format('         10.0  10.0  10.0     ', &
             '-      a_hmax, a_hmin, a_vert')
      write(lun,58)
 58   format('1    0.10                     ', &
             '-Four  nst, nugget effect')
      write(lun,59)
 59   format('1    0.90 0.0   0.0   0.0     ', &
             '-      it,cc,ang1,ang2,ang3')
      write(lun,60)
 60   format('         10.0  10.0  10.0     ', &
             '-      a_hmax, a_hmin, a_vert')
      write(lun,61)
 61   format('LVA_field_25_menos22.out     ', &
             '-file containing the LVA grid')
      write(lun,62)
 62   format('1 2 3 4 5     ', &
             '                  -LVA grid columns')
      write(lun,63)
 63   format('68.0  460301  25.0     ', &
             '             -nx,xmn,xsiz (LVA GRID)')
      write(lun,64)
 64   format('32.0  7069512  25.0     ', &
             '             -ny,ymn,ysiz (LVA GRID)')
      write(lun,65)
 65   format('4.0  3065  25.0     ', &
             '             -nz,zmn,zsiz (LVA GRID)')
      write(lun,66)
 66   format('1     ', &
             '                  -noffsets for graph')
      write(lun,67)
 67   format('2     ', &
             '                  - MDS? 2=L-ISOMAP 3=read dist')
      write(lun,68)
 68   format('2  2  2     ', &
             '                  -number of landmark points in x,y,z')
      write(lun,69)
 69   format('-1     ', &
             '                  -max nº of dim (set -1 to use max)')
      write(lun,70)
 70   format('200     ', &
             '                  -maximum search radii')
      write(lun,71)
 71   format('-1     ', &
             '                 -max nº of dimensions to use in search')

      close(lun)
      return
      end
      subroutine srchsupr(xloc,yloc,zloc,radsqd,irot,MAXROT,rotmat,&
                          nsbtosr,ixsbtosr,iysbtosr,izsbtosr,noct,nd,&
                          x,y,z,tmp,nisb,nxsup,xmnsup,xsizsup,&
                          nysup,ymnsup,ysizsup,nzsup,zmnsup,zsizsup,&
                          nclose,close,infoct)
!-----------------------------------------------------------------------
!
!              Search Within Super Block Search Limits
!              ***************************************
!
!
! This subroutine searches through all the data that have been tagged in
! the super block subroutine.  The close data are passed back in the
! index array "close".  An octant search is allowed.
!
!
!
! INPUT VARIABLES:
!
!   xloc,yloc,zloc   location of point being estimated/simulated
!   radsqd           squared search radius
!   irot             index of the rotation matrix for searching
!   MAXROT           size of rotation matrix arrays
!   rotmat           rotation matrices
!   nsbtosr          Number of super blocks to search
!   ixsbtosr         X offsets for super blocks to search
!   iysbtosr         Y offsets for super blocks to search
!   izsbtosr         Z offsets for super blocks to search
!   noct             If >0 then data will be partitioned into octants
!   nd               Number of data
!   x(nd)            X coordinates of the data
!   y(nd)            Y coordinates of the data
!   z(nd)            Z coordinates of the data
!   tmp(nd)          Temporary storage to keep track of the squared
!                      distance associated with each data
!   nisb()                Array with cumulative number of data in each
!                           super block.
!   nxsup,xmnsup,xsizsup  Definition of the X super block grid
!   nysup,ymnsup,ysizsup  Definition of the X super block grid
!   nzsup,zmnsup,zsizsup  Definition of the X super block grid
!
!
!
! OUTPUT VARIABLES:
!
!   nclose           Number of close data
!   close()          Index of close data
!   infoct           Number of informed octants (only computes if
!                      performing an octant search)
!
!
!
! EXTERNAL REFERENCES:
!
!   sqdist           Computes anisotropic squared distance
!   sortem           Sorts multiple arrays in ascending order
!
!
!
!-----------------------------------------------------------------------
      real    x(*),y(*),z(*),tmp(*),close(*)
      real*8  rotmat(MAXROT,3,3),hsqd,sqdist
      integer nisb(*),inoct(8)
      integer ixsbtosr(*),iysbtosr(*),izsbtosr(*)
      logical inflag
!
! Determine the super block location of point being estimated:
!
      call getindx(nxsup,xmnsup,xsizsup,xloc,ix,inflag)
      call getindx(nysup,ymnsup,ysizsup,yloc,iy,inflag)
      call getindx(nzsup,zmnsup,zsizsup,zloc,iz,inflag)
!
! Loop over all the possible Super Blocks:
!
      nclose = 0
      do 1 isup=1,nsbtosr
!
! Is this super block within the grid system:
!
            ixsup = ix + ixsbtosr(isup)
            iysup = iy + iysbtosr(isup)
            izsup = iz + izsbtosr(isup)
            if(ixsup.le.0.or.ixsup.gt.nxsup.or.iysup.le.0.or.iysup.gt.nysup.or.izsup.le.0.or.izsup.gt.nzsup) go to 1
!
! Figure out how many samples in this super block:
!
            ii = ixsup + (iysup-1)*nxsup + (izsup-1)*nxsup*nysup
            if(ii.eq.1) then
                  nums = nisb(ii)
                  i    = 0
            else
                  nums = nisb(ii) - nisb(ii-1)
                  i    = nisb(ii-1)
            endif
!
! Loop over all the data in this super block:
!
            do 2 ii=1,nums
                  i = i + 1
!
! Check squared distance:
!
                  hsqd = sqdist(xloc,yloc,zloc,x(i),y(i),z(i),irot,MAXROT,rotmat)
                  if(real(hsqd).gt.radsqd) go to 2
!
! Accept this sample:
!
                  nclose = nclose + 1
                  close(nclose) = real(i)
                  tmp(nclose)  = real(hsqd)
                  if (xloc.eq.181.5.and.yloc.eq.41.5)then
                        write(*,*) x(i),y(i)
                  end if
 2          continue
 1    continue
!
! Sort the nearby samples by distance to point being estimated:
!
      call sortem(1,nclose,tmp,1,close)
!
! If we aren't doing an octant search then just return:
!
      if(noct.le.0) return
!
! PARTITION THE DATA INTO OCTANTS:
!
      do i=1,8
            inoct(i) = 0
      end do
!
! Now pick up the closest samples in each octant:
!
      nt = 8*noct
      na = 0
      do j=1,nclose
            i  = int(close(j))
            h  = tmp(j)
            dx = x(i) - xloc
            dy = y(i) - yloc
            dz = z(i) - zloc
            if(dz.lt.0.) go to 5
            iq=4
            if(dx.le.0.0 .and. dy.gt.0.0) iq=1
            if(dx.gt.0.0 .and. dy.ge.0.0) iq=2
            if(dx.lt.0.0 .and. dy.le.0.0) iq=3
            go to 6
 5          iq=8
            if(dx.le.0.0 .and. dy.gt.0.0) iq=5
            if(dx.gt.0.0 .and. dy.ge.0.0) iq=6
            if(dx.lt.0.0 .and. dy.le.0.0) iq=7
 6          continue
            inoct(iq) = inoct(iq) + 1
!
! Keep this sample if the maximum has not been exceeded:
!
            if(inoct(iq).le.noct) then
                  na = na + 1
                  close(na) = i
                  tmp(na)   = h
                  if(na.eq.nt) go to 7
            endif
      end do
!
! End of data selection. Compute number of informed octants and return:
!
 7    nclose = na
      infoct = 0
      do i=1,8
            if(inoct(i).gt.0) infoct = infoct + 1
      end do
!
! Finished:
!
      return
      end
