module global
!just contains some global variables, mostly for trouble shooting
    implicit none
    
    real*8, public :: F_squared
    real*8, allocatable, dimension (:,:,:,:,:,:) :: template_D !template to speed up dist cals
    real*8, allocatable, dimension (:,:,:) :: norm_D !template of lengths of offsets
    real*8, allocatable, dimension (:,:) :: data_cov ! store the distance between all points
    real*8, allocatable, dimension (:,:) :: a_JM,d2lanmark,vec,vec1,val,subB2,subB22,vec_test
    real*8, allocatable, target, dimension (:,:) :: coord_ISOMAP,coord_ISOMAP_trans
    real*4, allocatable, target, dimension (:,:) :: coord_ISOMAP_trans_tmp 
    real*4, allocatable, target, dimension (:,:) :: coord_ISOMAP_trans_rearranged 
    real*8, allocatable, dimension (:,:) :: temp_coord,coord_LLE,anticline_ex_dij_orij
    real*8, allocatable, dimension (:,:) :: subB
    real*8, allocatable :: evalues(:), vectors(:,:),x_loc(:),y_loc(:),z_loc(:)
    integer xland,yland,zland,xspace,yspace,zspace,add_paths_data,add_paths_land

    integer, allocatable :: ex_results(:),ex_sim_array(:)
    real*8, allocatable :: est(:),ex_results_d(:)
    
    real*8, allocatable, dimension (:) :: r_JM, r_dist,dist_temp
    real*8, allocatable, dimension (:) :: x_JM,r_tmp
    integer, allocatable, dimension (:) :: bad_data,IPIV,data_ind
    integer, allocatable :: landpts(:)
    real(kind=8) :: bbb = 0.5 !for John Manchuk's robust solver    
    real*8, pointer, dimension (:) :: hash_dist ! hash table to store the distances
    integer, pointer, dimension (:,:) :: keys ! hash table to store the keys
    integer ndata,cal_stress,cur_node
    real*8 post_MDS_dist,sum_dist
    integer :: nbad_matrix !count the number of bad matricies
    integer :: ncol !count the number of colisions for the hash tabel
    integer :: n_seg ! number of times optimization segments the line
    integer :: max_seg ! if greater than -1, number of times to segment the line
    integer :: n_hash_table !count number of entries in the hash table
    integer T_size ! the size of the hash table
    integer max_hash !most number of entries in the hash table
    real*8  ::  min_dist_tol,stress !for timers
    real*4 :: elapsed, total, start_time,sum_time
    integer, parameter :: time_out= 11,time_out2= 13, T_eigs=145, T_matrix=146, T_dijk=147,T_search=148
    logical err_beta,inv_dist
    character*40 error_message
    real*8 min_eigen
    integer n_indef,na_total,na_cnt,good,graph_opt,ldistout,all_dist
    integer  indef_mat,c1,c2,c3,assign_nodes,xind,yind,zind,LLE_opt
    integer xyzland
    real*8 :: DIAG1,dist,pp1(20),pp2(20)
    integer dim,d_tree
    real*8 edge_dist
    integer edge_NODE,seedd,nreal,MDS_opt
    real*8 power,rad_inv_dis
    real*8 pp
    real*4 exhaustive_srch_time,time_temp,kd_time,xpp
    
    logical use_kd_tree
    
    !integer nd   (este parámetro está en sisim.inc)
    integer n_searched
    
    logical,allocatable,target :: is_usable(:),is_usable_rearranged(:)
        
    !real(kind(0.0)), target, allocatable :: input_data(:,:)
    !real*8, target, allocatable :: input_data(:,:)
    real*8, target, allocatable :: input_data(:,:)
   
    logical, target, allocatable :: input_usable(:)

    integer ivtype,ncut,test 
    integer, parameter :: MAXNST=4    
    integer, parameter :: MAXDT=9    
    real*8 VERSION
    integer MAXNOD,MAXCXY,MAXCTX,MAXCTY,MAXCTZ,MAXKR1,MAXSBX, &
                    MAXSBY,MAXSBZ,MAXCUT,MXCUT,MAXROT
    integer MAXDAT, MAXTAB
    integer MAXX,MAXY,MAXZ,MXYZ,MAXXYZ,MAXSAM,MAXKR2,MAXSB,MXSXY,MXSX,MAXORD
    real av,ss
    integer lin,lout,ldbg
    integer ixl,iyl,izl,ivrl
    integer ixs,iys,izs,imbsim
    real tmin,tmax
    real*8 zmin,zmax
    integer ltail,utail,middle
    real utpar,ltpar,mpar
    integer itabvr,itabwt,idbg
    integer nsim,nx,ny,nz,nxy,nxyz
    real*8 xmn,ymn,zmn,xsiz,ysiz,zsiz,radius,radius1,radius2,skmean,unbias,radsqd
    real*8 radius_lva,radsqd_lva
    real*8 sanis1,sanis2,sang1,sang2,sang3
    real*8 sanis1_lva,sanis2_lva,sang1_lva,sang2_lva,sang3_lva
    integer ixv(1),idum
    integer ndmax,nodmax
    integer sstrat,mults,nmult
    integer noct
    real*8,parameter :: EPSLON=0.000001
    real*8,parameter :: UNEST=-999.0
    integer mxctx,mxcty,mxctz
    integer ktype,istart,mik,maxsec
    real cutmik
    real*8 aa1,aa2
    integer ia1g,ia2g,ia3g,ir1,ir2
    integer inflag
    integer locx,locy,locz
    integer nd,ng,nhd,nsec
    real vrt,xd,tcdf,cp,oldcp

!    real xmn,ymn,zmn,xsiz,ysiz,zsiz
!    integer nx,ny,nz
!    integer,public :: MAXNOD,MAXCXY,MAXCTX,MAXCTY,MAXCTZ,MAXKR1,MAXSBX, &
!                    MAXSBY,MAXSBZ,MAXCUT,MXCUT,MAXROT
!    integer,parameter,public :: MAXNST = kind(9)

end module global
