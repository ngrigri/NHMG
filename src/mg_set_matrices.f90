module mg_set_matrices

  use mg_mpi
  use mg_tictoc
  use mg_namelist
  use mg_grids
  use mg_mpi_exchange
  use mg_gather
  use mg_netcdf_out

  implicit none

contains

  !-----------------------------------------------------------------------------------
  subroutine set_matrices()

    ! Define matrix coefficients cA
    ! Coefficients are stored in order of diagonals
    ! cA(1,:,:,:)      -> p(k,j,i)
    ! cA(2,:,:,:)      -> p(k-1,j,i)
    ! cA(3,:,:,:)      -> p(k+1,j-1,i)
    ! cA(4,:,:,:)      -> p(k,j-1,i)
    ! cA(5,:,:,:)      -> p(k-1,j-1,i)
    ! cA(6,:,:,:)      -> p(k+1,j,i-1)
    ! cA(7,:,:,:)      -> p(k,j,i-1)
    ! cA(8,:,:,:)      -> p(k-1,j,i-1)

    integer(kind=ip) :: lev
    integer(kind=ip) :: k,j,i
    integer(kind=ip) :: nx,ny,nz

    real(kind=rp), dimension(:,:),     pointer :: dx,dy
    real(kind=rp), dimension(:,:,:),   pointer :: zr,zw
    real(kind=rp), dimension(:,:,:),   pointer :: dzw
    real(kind=rp), dimension(:,:,:),   pointer :: Arx,Ary
    real(kind=rp), dimension(:,:,:),   pointer :: zxdy,zydx
    real(kind=rp), dimension(:,:,:),   pointer :: cw
    real(kind=rp), dimension(:,:,:,:), pointer :: cA

    real(kind=rp), parameter :: one  = 1._rp
    real(kind=rp), parameter :: eigh = one/8._rp
    real(kind=rp), parameter :: qrt  = 0.25_rp
    real(kind=rp), parameter :: hlf  = 0.5_rp

    integer(kind=ip), save :: iter_mat=-1
    iter_mat = iter_mat + 1

    if (myrank==0) write(*,*)'- set matrices:',iter_mat

    do lev = 1, nlevs

       if (myrank==0) write(*,*)'   lev=',lev

       nx = grid(lev)%nx
       ny = grid(lev)%ny
       nz = grid(lev)%nz

       dx => grid(lev)%dx
       dy => grid(lev)%dy
       zr => grid(lev)%zr
       zw => grid(lev)%zw
       dzw => grid(lev)%dzw
       Arx => grid(lev)%Arx
       Ary => grid(lev)%Ary
       cw  => grid(lev)%cw
       zxdy => grid(lev)%zxdy
       zydx => grid(lev)%zydx
       cA => grid(lev)%cA 

       !! interaction coeff with neighbours
   
       ! lower level      
       k = 1
       do i = 1,nx
          do j = 1,ny+1
             cA(3,k,j,i) = qrt * ( &
                  ( hlf * (zr(k+1,j+1,i  )-zr(k+1,j-1,i  )) / dy(j,i) ) * dx(j,i) + &
                  ( hlf * (zr(k,j,i  )-zr(k,j-2,i  )) / dy(j-1,i) ) * dx(j-1,i) )
             cA(4,k,j,i) =                                                                &
                                ! couples with j-1
                  ( qrt                                                                 * &
                  ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j-1,i) - zw(k,j-1,i) )             * &
                  (dx(j,i)+dx(j-1,i)) )/ (hlf * (dy(j,i)+dy(j-1,i)))                      & 
                                ! topo terms 
                  - ( (( hlf * (zr(k,j+1,i)-zr(k,j-1,i)) / dy(j,i) ) * dx(j,i))**2 / &
                  ( cw(k,j  ,i) + cw(k+1,j  ,i) )     &
                  +   (( hlf * (zr(k,j,i)-zr(k,j-2,i)) / dy(j-1,i) ) * dx(j-1,i))**2 / &
                  ( cw(k,j-1,i) + cw(k+1,j-1,i) )   ) &
                                ! from j,k cross terms
                  - qrt * ( &
                  (( hlf * (zr(k,j,i  )-zr(k,j-2,i  )) / dy(j-1,i) ) * dx(j-1,i)) -  &
                  ( hlf * (zr(k,j+1,i  )-zr(k,j-1,i  )) / dy(j,i) ) * dx(j,i))  
          enddo
       enddo
       do i = 1,nx+1
          do j = 1,ny
             cA(6,k,j,i) = qrt * ( &
                  ( hlf * (zr(k+1,j  ,i+1)-zr(k+1,j  ,i-1)) / dx(j,i) ) * dy(j,i) + &
                  ( hlf * (zr(k,j  ,i)-zr(k,j  ,i-2)) / dx(j,i-1) ) * dy(j,i-1) )  ! couples with k+1 i-1
             cA(7,k,j,i) =                                                                &
                                ! couples with i-1
                  ( qrt                                                                 * &
                  ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j,i-1) - zw(k,j,i-1) )             * &
                  (dy(j,i)+dy(j,i-1)) )                                                 / &
                  ( hlf * (dx(j,i)+dx(j,i-1)) )                                           &  
                                ! topo terms
                  - ( (( hlf * (zr(k,j  ,i+1)-zr(k,j  ,i-1)) / dx(j,i) ) * dy(j,i))**2 / &
                  ( cw(k,j,i  ) + cw(k+1,j,i  ) )     &
                  +   (( hlf * (zr(k,j  ,i)-zr(k,j  ,i-2)) / dx(j,i-1) ) * dy(j,i-1))**2 / &
                  ( cw(k,j,i-1) + cw(k+1,j,i-1) )   ) &
                                ! from i,k cross terms
                  - qrt * ( &
                  ( hlf * (zr(k,j  ,i  )-zr(k,j  ,i-2)) / dx(j,i-1) ) * dy(j,i-1) - &
                  ( hlf * (zr(k,j  ,i+1)-zr(k,j  ,i-1)) / dx(j,i  ) ) * dy(j,i) )
          enddo
       enddo
       do i = 1,nx+1
          do j = 0,ny
             ! only for k==1, couples with j+1,i-1
             cA(5,k,j,i) = &
                  + hlf * &
                  (( hlf * (zr(k,j+1,i+1)-zr(k,j+1,i-1)) / dx(j+1,i) ) * dy(j+1,i)) * &
                  (( hlf * (zr(k,j+2,i  )-zr(k,j  ,i  )) / dy(j+1,i) ) * dx(j+1,i)) / &
                  ( cw(k,j+1,i  ) + cw(k+1,j+1,i  ))  &
                  + hlf * &
                  (( hlf * (zr(k,j  ,i)-zr(k,j  ,i-2)) / dx(j,i-1) ) * dy(j,i-1)) * &
                  (( hlf * (zr(k,j+1,i-1)-zr(k,j-1,i-1)) / dy(j,i-1) ) * dx(j,i-1)) / &
                  ( cw(k,j  ,i-1) + cw(k+1,j  ,i-1))
          enddo
       enddo
       do i = 1,nx+1
          do j = 1,ny+1
             ! only for k==1, couples with j-1,i-1
             cA(8,k,j,i) = &
                  - hlf * &
                  (( hlf * (zr(k,j-1,i+1)-zr(k,j-1,i-1)) / dx(j-1,i) ) * dy(j-1,i)) * &
                  (( hlf * (zr(k,j,i  )-zr(k,j-2,i  )) / dy(j-1,i) ) * dx(j-1,i)) / &
                  (cw(k,j-1,i  ) + cw(k+1,j-1,i  )) & 
                  - hlf * &
                  (( hlf * (zr(k,j,i)-zr(k,j,i-2)) / dx(j,i-1) ) * dy(j,i-1)) * &
                  (( hlf * (zr(k,j+1,i-1)-zr(k,j-1,i-1)) / dy(j,i-1) ) * dx(j,i-1)) / &
                  (cw(k,j  ,i-1) + cw(k+1,j  ,i-1)) 
          enddo
       enddo

       ! interior levels
       do i = 1,nx
          do j = 1,ny
             do k = 2,nz-1 
                cA(2,k,j,i) =  cw(k,j,i)                                  ! couples with k-1
             enddo
          enddo
       enddo
       do i = 1,nx
          do j = 1,ny+1
             do k = 2,nz-1 
                cA(3,k,j,i) =  qrt * ( &
                     ( hlf * (zr(k+1,j+1,i  )-zr(k+1,j-1,i  )) / dy(j,i) ) * dx(j,i) + &
                     ( hlf * (zr(k,j,i  )-zr(k,j-2,i  )) / dy(j-1,i) ) * dx(j-1,i) )     ! couples with k+1 j-1
                cA(4,k,j,i) =  ( qrt * &
                     ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j-1,i) - zw(k,j-1,i) ) * &
                     (dx(j,i)+dx(j-1,i)) ) / ( hlf * (dy(j,i)+dy(j-1,i)) ) ! couples with j-1
                cA(5,k,j,i) =- qrt * ( &
                     (( hlf * (zr(k-1,j+1,i  )-zr(k-1,j-1,i  )) / dy(j,i) ) * dx(j,i)) + &
                     (( hlf * (zr(k,j,i  )-zr(k,j-2,i  )) / dy(j-1,i) ) * dx(j-1,i)) )     ! couples with k-1 j-1
             enddo
          enddo
       enddo
       do i = 1,nx+1
          do j = 1,ny 
             do k = 2,nz-1 
                cA(6,k,j,i) =  qrt * ( &
                     (( hlf * (zr(k+1,j  ,i+1)-zr(k+1,j  ,i-1)) / dx(j,i) ) * dy(j,i)) + &
                     (( hlf * (zr(k,j  ,i)-zr(k,j  ,i-2)) / dx(j,i-1) ) * dy(j,i-1)) )     ! Couples with k+1 i-1
                cA(7,k,j,i) =   (qrt * &
                     ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j,i-1) - zw(k,j,i-1) ) * &
                     (dy(j,i)+dy(j,i-1)) ) / ( hlf * (dx(j,i)+dx(j,i-1)) ) ! Couples with i-1
                cA(8,k,j,i) =- qrt * ( &
                     (( hlf * (zr(k-1,j  ,i+1)-zr(k-1,j  ,i-1)) / dx(j,i) ) * dy(j,i)) + &
                     (( hlf * (zr(k,j  ,i)-zr(k,j  ,i-2)) / dx(j,i-1) ) * dy(j,i-1))  )     ! Couples with k-1 i-1
             enddo
          enddo
       enddo

       ! upper level
       k = nz
       do i = 1,nx    
          do j = 1,ny 
             cA(2,k,j,i) = cw(k,j,i)                                    ! couples with k-1
          enddo
       enddo
       do i = 1,nx
          do j = 1,ny+1
             cA(4,k,j,i) = ( qrt * &
                  ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j-1,i) - zw(k,j-1,i) ) * &
                  (dx(j,i)+dx(j-1,i)) ) / ( hlf * (dy(j,i)+dy(j-1,i)) ) & ! couples with j-1
                  + qrt * ( &
                  - (( hlf * (zr(k,j  ,i)-zr(k,j-2,i)) / dy(j-1,i) ) * dx(j-1,i)) &
                  + (( hlf * (zr(k,j+1,i)-zr(k,j-1,i)) / dy(j  ,i) ) * dx(j  ,i)) )
             cA(5,k,j,i) =- qrt * ( &
                  (( hlf * (zr(k-1,j+1,i  )-zr(k-1,j-1,i  )) / dy(j,i) ) * dx(j,i)) + &
                  (( hlf * (zr(k,j,i  )-zr(k,j-2,i  )) / dy(j-1,i) ) * dx(j-1,i)) )     ! couples with k-1 j-1
          enddo
       enddo
       do i = 1,nx+1
          do j = 1,ny 
             cA(7,k,j,i) = (qrt * &
                  ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j,i-1) - zw(k,j,i-1) ) * &
                  (dy(j,i)+dy(j,i-1)) ) / ( hlf * (dx(j,i)+dx(j,i-1)) ) & ! Couples with i-1
                  + qrt * ( &
                  - (( hlf * (zr(k,j,i  )-zr(k,j,i-2)) / dx(j,i-1)) * dy(j,i-1)) &
                  + (( hlf * (zr(k,j,i+1)-zr(k,j,i-1)) / dx(j,i  )) * dy(j,i  )) )
             cA(8,k,j,i) =- qrt * ( &
                  (( hlf * (zr(k-1,j  ,i+1)-zr(k-1,j  ,i-1)) / dx(j,i) ) * dy(j,i)) + &
                  (( hlf * (zr(k,j  ,i)-zr(k,j  ,i-2)) / dx(j,i-1) ) * dy(j,i-1)) )     ! Couples with k-1 i-1
          enddo
       enddo

       call fill_halo(lev,cA)

       !! interaction coeff with itself
       do i = 1,nx
          do j = 1,ny

             k = 1 !lower level
             cA(1,k,j,i) =                     &
                  -cA(2,k+1,j,i)               &
                  -cA(4,k,j,i)-cA(4,k,j+1,i)   &
                  -cA(7,k,j,i)-cA(7,k,j,i+1)   &
                  -cA(6,k,j,i)-cA(8,k+1,j,i+1) &
                  -cA(3,k,j,i)-cA(5,k+1,j+1,i) &
                  -cA(5,k,j,i)-cA(5,k,j-1,i+1) &
                  -cA(8,k,j,i)-cA(8,k,j+1,i+1)

             do k = 2,nz-1 !interior levels
                cA(1,k,j,i) = &
                     -cA(2,k,j,i)-cA(2,k+1,j,i)   &
                     -cA(4,k,j,i)-cA(4,k,j+1,i)   &
                     -cA(7,k,j,i)-cA(7,k,j,i+1)   &
                     -cA(6,k,j,i)-cA(6,k-1,j,i+1) &
                     -cA(8,k,j,i)-cA(8,k+1,j,i+1) & 
                     -cA(3,k,j,i)-cA(3,k-1,j+1,i) &
                     -cA(5,k,j,i)-cA(5,k+1,j+1,i)

             enddo

             k = nz !upper level
             cA(1,k,j,i) =                    &
                  - cA(2,k,j,i)               &
                  - cw(k+1,j,i)               &
                  + hlf * (hlf * (zr(k,j  ,i+2) - zr(k,j  ,i  )) / dx(j  ,i+1)) * dy(j  ,i+1) & ! + 0.25*zxdy(nz,j,i+1) + 0.25*zxdy(nz,j,i+1)
                  - hlf * (hlf * (zr(k,j  ,i  ) - zr(k,j  ,i-2)) / dx(j  ,i-1)) * dy(j  ,i-1) & ! - 0.25*zxdy(nz,j,i-1) - 0.25*zxdy(nz,j,i-1)
                  + hlf * (hlf * (zr(k,j+2,i  ) - zr(k,j  ,i  )) / dy(j+1,i  )) * dx(j+1,i  ) & ! + 0.25*zydx(nz,j+1,i) + 0.25*zydx(nz,j+1,i)
                  - hlf * (hlf * (zr(k,j  ,i  ) - zr(k,j-2,i  )) / dy(j-1,i  )) * dx(j-1,i  ) & ! - 0.25*zydx(nz,j-1,i) - 0.25*zydx(nz,j-1,i)
                  - cA(4,k,j,i)-cA(4,k,j+1,i) &
                  - cA(7,k,j,i)-cA(7,k,j,i+1) &
                  - cA(6,k-1,j,i+1)           &
                  - cA(8,k,j,i)               &
                  - cA(3,k-1,j+1,i)           &
                  - cA(5,k,j,i)
          enddo ! j
       enddo ! i

       if (netcdf_output) then
          if (myrank==0) write(*,*)'       write cA in a netcdf file'
          call write_netcdf(grid(lev)%cA,vname='ca',netcdf_file_name='cA.nc',rank=myrank,iter=lev)
       endif

    enddo

  end subroutine set_matrices

end module mg_set_matrices
