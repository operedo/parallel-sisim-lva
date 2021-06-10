module srch
  use global
!  implicit none

contains 

!kdtree2_n_nearest_around_point(real   (coord_ISOMAP(index,1:d_tree)),d_tree,nn=nclose, results=results)
  subroutine exhaustive_srch           (vecin,dim,nn,results,results_d,n2srch)
    ! Find the 'nn' vectors nearest to  'vecin',
    ! returing results in results(:)
    !n2srch has the number of previously simulated nodes to search in the global array ex_sim_array
    
    integer                        :: idxin, correltime, nn, dim,results(nn),n2srch,i,j,k,ind,tyty(1),largest
    real*8                           :: vecin(dim),d,largest_d,results_d(nn),results_temp(nn)
    logical        ::          skip

    !do an exhaustive search

    results_d=1e21
    !get the largest dist:
    largest   = 1
    largest_d = 1e21
    results_temp=0
    
    do i=1,n2srch
        !get the distance to this node
        
        ind = ex_sim_array(i)
        
        if(ind == -999) exit !(we are done now)
    
        !write(*,*)'Entering coord_ISOMAP_trans(1:',dim,',',ind,')' 
        !write(*,*)'coord_ISOMAP_trans(1:',dim,',',ind,')=',coord_ISOMAP_trans(1:dim,ind) 
        !d= sum( ( coord_ISOMAP(ind,1:dim)-vecin ) **2    )
        d= sum( ( coord_ISOMAP_trans(1:dim,ind)-vecin ) **2    )
        !write(*,*)'Finishing coord_ISOMAP_trans(1:',dim,',',ind,')' 
        
        if(d<largest_d) then  !if d==0.0 then we are at a simualtion location and skip

            !do we have this location already?
            skip=.false.
            do k=1,nn
                if(results_temp(k) ==real(ind)) then 
                    skip = .true.
                    exit
                end if
            end do

            if(not(skip) ) then 
                results_d(largest) = d
                results_temp(largest) = real(ind)
                
                !get the new largest value
                tyty = MAXLOC(results_d)
                largest=tyty(1)
                largest_d=results_d(largest)
            end if
        end if
    end do

    !call sortem(1,nn,results_d,1,results_temp) 
    call sortem2(1,nn,results_d,1,results_temp)

    if(minval(results_temp) == -999) then !need to get new nclose
        do i=1,nn
            if(results_temp(i) == -999) then
                nn=i-1
                exit
            end if
        end do
    end if


    results=int(results_temp)

end subroutine exhaustive_srch


end module srch
