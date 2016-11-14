module mg_define_matrices

  use mg_mpi
  use mg_tictoc
  use mg_namelist
  use mg_grids
  use mg_zr_zw
  use mg_mpi_exchange
  use mg_gather
  use mg_netcdf_out

  implicit none

  interface define_matrices
     module procedure           &
          define_matrices_topo
  end interface define_matrices

  !NG comment: constants in a mg_cst.f90 file ?
  real(kind=rp), parameter :: one  = 1._rp
  real(kind=rp), parameter :: eigh = one/8._rp
  real(kind=rp), parameter :: qrt  = 0.25_rp
  real(kind=rp), parameter :: hlf  = 0.5_rp

contains

  !-------------------------------------------------------------------------  
  subroutine define_matrices_topo(zeta, h)

    real(kind=rp), dimension(:,:), intent(in) :: zeta
    real(kind=rp), dimension(:,:), intent(in) :: h

    integer(kind=ip)::  lev

    real(kind=rp), dimension(:,:), pointer :: zetaf
    real(kind=rp), dimension(:,:), pointer :: zetac

    real(kind=rp), dimension(:,:), pointer :: hf
    real(kind=rp), dimension(:,:), pointer :: hc

    real(kind=rp), dimension(:,:,:), pointer :: zrc
    real(kind=rp), dimension(:,:,:), pointer :: zwc

    integer(kind=ip) :: ny,nx
    integer(kind=ip) :: nyf,nxf
    integer(kind=ip) :: nyc,nxc

    if (myrank==0) write(*,*)'- define matrices from topography (h):'

    do lev = 1, nlevs

       if (myrank==0) write(*,*)'   lev=',lev

       nx=grid(lev)%nx
       ny=grid(lev)%ny

       if (lev == 1) then
     
          grid(lev)%zeta(0:ny+1,0:nx+1) = zeta
          grid(lev)%h (0:ny+1,0:nx+1) =  h

       else

          nxf =grid(lev-1)%nx
          nyf =grid(lev-1)%ny

          zetaf  => grid(lev-1)%zeta
          hf  => grid(lev-1)%h

          if (grid(lev)%gather == 1) then
             nxc= nx/grid(lev)%ngx
             nyc= ny/grid(lev)%ngy
             allocate(zetac(0:nyc+1,0:nxc+1))
             allocate( hc(0:nyc+1,0:nxc+1))
          else
             nxc = nx
             nyc = ny        
             zetac  => grid(lev)%zeta
             hc  => grid(lev)%h
          endif

          ! Coarsen zeta and h

          zetac(1:nyc,1:nxc)  = qrt      * ( &
               zetaf(1:nyf  :2,1:nxf  :2)  + &
               zetaf(2:nyf+1:2,1:nxf  :2)  + &
               zetaf(1:nyf  :2,2:nxf+1:2)  + &
               zetaf(2:nyf+1:2,2:nxf+1:2)  )

          hc(1:nyc,1:nxc)  = qrt      * ( &
               hf(1:nyf  :2,1:nxf  :2)  + &
               hf(2:nyf+1:2,1:nxf  :2)  + &
               hf(1:nyf  :2,2:nxf+1:2)  + &
               hf(2:nyf+1:2,2:nxf+1:2)  )

          if (grid(lev)%gather == 1) then

             call gather(lev,zetac,grid(lev)%zeta)
             call gather(lev, hc,grid(lev)%h)

             deallocate(zetac)
             deallocate( hc)
          endif

       endif ! lev == 1

       call fill_halo(lev, grid(lev)%zeta)
       call fill_halo(lev, grid(lev)%h)

       zrc => grid(lev)%zr
       zwc => grid(lev)%zw

       ! Compute zr and zw
       call setup_zr_zw                    (  & 
            nhhc,nhtheta_b,nhtheta_s,     &
            grid(lev)%zeta,grid(lev)%h,   &  ! input args
            grid(lev)%zr, grid(lev)%zw,   &  ! output args
            coord_type='new_s_coord'      )    ! optional

       call fill_halo(lev,grid(lev)%zr) ! Special fill_halo nh = 2
       call fill_halo(lev,grid(lev)%zw) ! Special fill_halo nh = 2

       if (netcdf_output) then
          call write_netcdf(grid(lev)%zeta,vname='zeta',netcdf_file_name='zeta.nc',rank=myrank,iter=lev)
          call write_netcdf(grid(lev)%h ,vname='h' ,netcdf_file_name='h.nc' ,rank=myrank,iter=lev)
          call write_netcdf(grid(lev)%zr,vname='zr',netcdf_file_name='zr.nc',rank=myrank,iter=lev)
          call write_netcdf(grid(lev)%zw,vname='zw',netcdf_file_name='zw.nc',rank=myrank,iter=lev)
       endif

       ! Define matrix coefficients from zr and zw coarsened
       call define_matrix(lev, grid(lev)%zr, grid(lev)%zw)

    enddo ! lev

  end subroutine define_matrices_topo

  !-----------------------------------------------------------------------------------
  subroutine define_matrix(lev, zr, zw)

    integer(kind=ip),intent(in):: lev
    real(kind=rp), dimension(:,:,:), pointer, intent(in) :: zr, zw

    ! Define matrix coefficients cA
    ! Coefficients are stored in order of diagonals
    ! cA(1,:,:,:)      -> p(k,j,i)
    ! cA(2,:,:,:)      -> p(k-1,j,i)
    ! cA(3,:,:,:)      -> p(k+1,j-1,i)
    ! cA(4,:,:,:)      -> p(k,j-1,i)
    ! cA(5,:,:,:)      -> p(k-1,j-1,i)
    ! cA(6,:,:,:)      -> p(k+1,j,i-1)
    ! cA(7,:,:,:)      -> p(k,j,i-1)xd
    ! cA(8,:,:,:)      -> p(k-1,j,i-1)

    integer(kind=ip):: k, j, i
    integer(kind=ip):: nx, ny, nz

    real(kind=rp) :: Arz
    real(kind=rp), dimension(:,:),     pointer :: dx,dy
    real(kind=rp), dimension(:,:,:),   pointer :: dzw
    real(kind=rp), dimension(:,:,:),   pointer :: zydx,zxdy
    real(kind=rp), dimension(:,:,:),   pointer :: cw
    real(kind=rp), dimension(:,:,:,:), pointer :: cA

    nx = grid(lev)%nx
    ny = grid(lev)%ny
    nz = grid(lev)%nz

    dx => grid(lev)%dx 
    dy => grid(lev)%dy 

    cA => grid(lev)%cA 

    !- Used in compute_rhs -!
    if (lev == 1) then
       dzw => grid(lev)%dzw
       do i = 0,nx+1
          do j = 0,ny+1
             dzw(1,j,i) = zr(1,j,i) - zw(1,j,i) !!
             do k = 2,nz
                dzw(k,j,i) = zr(k,j,i) - zr(k-1,j,i) !!  cell height at w-points
             enddo
             dzw(nz+1,j,i) = zw(nz+1,j,i) - zr(nz,j,i) !!
          enddo
       enddo

       !! Slopes in x- and y-direction defined at rho-points
       zxdy => grid(lev)%zxdy
       zydx => grid(lev)%zydx
       do i = 0,nx+1
          do j = 0,ny+1
             do k = 1, nz
                zydx(k,j,i) = hlf * (( zr(k,j+1,i  ) - zr(k,j-1,i  ) ) / dy(j,i) ) * dx(j,i)
                zxdy(k,j,i) = hlf * (( zr(k,j  ,i+1) - zr(k,j  ,i-1) ) / dx(j,i) ) * dy(j,i)
             enddo
          enddo
       enddo
    endif

    !- Used also in compute_rhs -!
    cw => grid(lev)%cw
    do i = 0,nx+1
       do j = 0,ny+1

          Arz = dx(j,i)*dy(j,i)

          k=1
          cw(k,j,i) = ( Arz / (zr(k,j,i)-zw(k,j,i)) ) * &
               (one + &
               ( hlf * (zw(k,j  ,i+1)-zw(k,j  ,i-1)) / dx(j,i) )**2 + &
               ( hlf * (zw(k,j+1,i  )-zw(k,j-1,i  )) / dy(j,i) )**2 )

          do k = 2,nz
             cw(k,j,i) = ( Arz / (zr(k,j,i)-zr(k-1,j,i)) ) * &
                  (one + &
                  ( hlf * (zw(k,j  ,i+1)-zw(k,j  ,i-1)) / dx(j,i) )**2 + &
                  ( hlf * (zw(k,j+1,i  )-zw(k,j-1,i  )) / dy(j,i) )**2 )

          enddo

          k=nz+1
          cw(k,j,i) = ( Arz / (zw(k,j,i)-zr(k-1,j,i)) ) * &
               (one + &
               ( hlf * (zw(k,j  ,i+1)-zw(k,j  ,i-1)) / dx(j,i) )**2 + &
               ( hlf * (zw(k,j+1,i  )-zw(k,j-1,i  )) / dy(j,i) )**2 )

       enddo
    enddo

!! interaction coeff with neighbours
!XX
          !---------!
          !- K = 1 -! lower level
          !---------!
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

!XX
          !---------------!
          !- K = 2, nz-1 -! interior levels
          !---------------!

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

!XX
          !----------!
          !- K = nz -! upper level
          !----------!
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

  end subroutine define_matrix

end module mg_define_matrices
