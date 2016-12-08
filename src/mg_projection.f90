module mg_projection

  use mg_cst
  use mg_mpi
  use mg_tictoc
  use mg_namelist
  use mg_grids
  use mg_mpi_exchange
  use mg_gather
  use mg_netcdf_out

  implicit none

contains
  !-------------------------------------------------------------------------     
  subroutine set_rhs(u,v,w)

    real(kind=rp), dimension(:,:,:), pointer, intent(in) :: u,v,w

    integer(kind=ip):: k,j,i
    integer(kind=ip):: nx,ny,nz

    real(kind=rp) :: gamma
    real(kind=rp), dimension(:,:)  , pointer :: dx,dy
    real(kind=rp), dimension(:,:),   pointer :: dxu,dyv
    real(kind=rp), dimension(:,:,:), pointer :: dzw
    real(kind=rp), dimension(:,:,:), pointer :: Arx,Ary
    real(kind=rp), dimension(:,:),   pointer :: Arz
    real(kind=rp), dimension(:,:,:), pointer :: zxdy,zydx
    real(kind=rp), dimension(:,:,:), pointer :: alpha
    real(kind=rp), dimension(:,:)  , pointer :: beta
    real(kind=rp), dimension(:,:,:), pointer :: uf,vf,wf
    real(kind=rp), dimension(:,:,:), pointer :: rhs

    if (myrank==0) write(*,*)'   - set rhs:'

    nx = grid(1)%nx
    ny = grid(1)%ny
    nz = grid(1)%nz

    dx    => grid(1)%dx    !
    dy    => grid(1)%dy    !
    dxu   => grid(1)%dxu   !
    dyv   => grid(1)%dyv   !
    dzw   => grid(1)%dzw   !
    Arx   => grid(1)%Arx   !
    Ary   => grid(1)%Ary   !
    Arz   => grid(1)%Arz   !
    alpha => grid(1)%alpha !
    beta  => grid(1)%beta  !
    zxdy  => grid(1)%zxdy  !
    zydx  => grid(1)%zydx  !

    rhs => grid(1)%b
    rhs(:,:,:) = zero

    !- UF -!

    uf => grid(1)%dummy3Dnz
    uf(:,:,:) = zero

    ! lower level
    k = 1
    do i = 1,nx+1  
       do j = 1,ny 

          gamma = one - qrt * ( &
               (zxdy(k,j,i  )/dy(j,i  ))**2/alpha(k,j,i  ) + &
               (zxdy(k,j,i-1)/dy(j,i-1))**2/alpha(k,j,i-1) )

          uf(k,j,i) = gamma * Arx(k,j,i) * u(k,j,i) &
               - qrt * ( &
               + zxdy(k,j,i  ) * dzw(k+1,j,i  ) * w(k+1,j,i  ) &
               + zxdy(k,j,i-1) * dzw(k+1,j,i-1) * w(k+1,j,i-1) )  &
               - beta(j,i-1)   * dyv(j  ,i-1)   * v(k,j  ,i-1) &
               - beta(j,i-1)   * dyv(j+1,i-1)   * v(k,j+1,i-1) &
               - beta(j,i  )   * dyv(j  ,i  )   * v(k,j  ,i  ) &
               - beta(j,i  )   * dyv(j+1,i  )   * v(k,j+1,i  )

       enddo
    enddo

    !interior levels
    do i = 1,nx+1  
       do j = 1,ny 
          do k = 2,nz-1 

             uf(k,j,i) = Arx(k,j,i) * u(k,j,i) &
                  - qrt * ( &
                  + zxdy(k,j,i  ) * dzw(k  ,j,i  ) * w(k  ,j,i  ) &
                  + zxdy(k,j,i  ) * dzw(k+1,j,i  ) * w(k+1,j,i  ) &
                  + zxdy(k,j,i-1) * dzw(k  ,j,i-1) * w(k  ,j,i-1) &
                  + zxdy(k,j,i-1) * dzw(k+1,j,i-1) * w(k+1,j,i-1) )
          enddo
       enddo
    enddo

    k = nz !upper level
    do i = 1,nx+1  
       do j = 1,ny 

          uf(k,j,i) = Arx(k,j,i) * u(k,j,i) &
               - qrt * ( &
               + zxdy(k,j,i  ) *       dzw(k  ,j,i  ) * w(k  ,j,i  ) &
               + zxdy(k,j,i  ) * two * dzw(k+1,j,i  ) * w(k+1,j,i  ) &
               + zxdy(k,j,i-1) *       dzw(k  ,j,i-1) * w(k  ,j,i-1) &
               + zxdy(k,j,i-1) * two * dzw(k+1,j,i-1) * w(k+1,j,i-1) )
       enddo
    enddo

    call fill_halo(1,uf,lbc_null='u')

    do i = 1,nx
       do j = 1,ny 
          do k = 1,nz
             rhs(k,j,i) = uf(k,j,i+1) - uf(k,j,i) 
          enddo
       enddo
    enddo

    uf => null()

    !- VF -!

    vf => grid(1)%dummy3Dnz
    vf(:,:,:) = zero

    !lower level
    k = 1
    do i = 1,nx
       do j = 1,ny+1

          gamma = one - qrt * (  &
               (zydx(k,j  ,i)/dx(j  ,i))**2/alpha(k,j  ,i  ) + &
               (zydx(k,j-1,i)/dx(j-1,i))**2/alpha(k,j-1,i) )

          vf(k,j,i) = gamma * Ary(k,j,i) * v(k,j,i) &
               - qrt * ( &
               + zydx(k,j  ,i) * dzw(k+1,j  ,i) * w(k+1,j  ,i) &
               + zydx(k,j-1,i) * dzw(k+1,j-1,i) * w(k+1,j-1,i) ) &
               - beta(j-1,i)   * dxu(j-1,i  )   * u(k,j-1,i  ) &
               - beta(j-1,i)   * dxu(j-1,i+1)   * u(k,j-1,i+1) &
               - beta(j  ,i)   * dxu(j  ,i  )   * u(k,j  ,i  ) &
               - beta(j  ,i)   * dxu(j  ,i+1)   * u(k,j  ,i+1)
       enddo
    enddo

    !interior levels
    do i = 1,nx
       do j = 1,ny+1
          do k = 2,nz-1

             vf(k,j,i) = Ary(k,j,i) * v(k,j,i) &
                  - qrt * ( &
                  + zydx(k,j  ,i) * dzw(k  ,j  ,i) * w(k  ,j  ,i) &
                  + zydx(k,j  ,i) * dzw(k+1,j  ,i) * w(k+1,j  ,i) &
                  + zydx(k,j-1,i) * dzw(k  ,j-1,i) * w(k  ,j-1,i) &
                  + zydx(k,j-1,i) * dzw(k+1,j-1,i) * w(k+1,j-1,i) )
          enddo
       enddo
    enddo

    k = nz !upper level
    do i = 1,nx
       do j = 1,ny+1

          vf(k,j,i) = Ary(k,j,i) * v(k,j,i) &
               - qrt * ( &
               + zydx(k,j  ,i)       * dzw(k  ,j  ,i) * w(k  ,j  ,i) &
               + zydx(k,j  ,i) * two * dzw(k+1,j  ,i) * w(k+1,j  ,i) &
               + zydx(k,j-1,i)       * dzw(k  ,j-1,i) * w(k  ,j-1,i) &
               + zydx(k,j-1,i) * two * dzw(k+1,j-1,i) * w(k+1,j-1,i) ) 
       enddo
    enddo

    call fill_halo(1,vf,lbc_null='v')

    do i = 1,nx
       do j = 1,ny 
          do k = 1,nz
             rhs(k,j,i) =  rhs(k,j,i) + vf(k,j+1,i) - vf(k,j,i)
          enddo
       enddo
    enddo

    vf => null()

    !- WF -!

    wf => grid(1)%dummy3Dnzp
    wf(:,:,:) = zero

    !bottom
    k = 1
    do i = 1,nx
       do j = 1,ny
          wf(k,j,i) = zero
       enddo
    enddo

    !interior levels
    do i = 1,nx
       do j = 1,ny
          do k = 2,nz

             wf(k,j,i) = hlf * (alpha(k-1,j,i) + alpha(k,j,i)) * Arz(j,i) * w(k,j,i) &
                  - qrt * ( &
                  + zxdy(k  ,j,i) * dxu(j,i  ) * u(k  ,j,i  ) &
                  + zxdy(k  ,j,i) * dxu(j,i+1) * u(k  ,j,i+1) &
                  + zxdy(k-1,j,i) * dxu(j,i  ) * u(k-1,j,i  ) &
                  + zxdy(k-1,j,i) * dxu(j,i+1) * u(k-1,j,i+1) ) &
                  - qrt * ( &
                  + zydx(k  ,j,i) * dyv(j  ,i) * v(k  ,j  ,i) &
                  + zydx(k  ,j,i) * dyv(j+1,i) * v(k  ,j+1,i) &
                  + zydx(k-1,j,i) * dyv(j  ,i) * v(k-1,j  ,i) &
                  + zydx(k-1,j,i) * dyv(j+1,i) * v(k-1,j+1,i) )
          enddo
       enddo
    enddo

    !surface
    k = nz+1 
    do i = 1,nx
       do j = 1,ny

          wf(k,j,i) = alpha(k-1,j,i) * Arz(j,i) * w(k,j,i) &
               - hlf * ( &
               + zxdy(k-1,j,i) * dxu(j,i  ) * u(k-1,j,i  ) &
               + zxdy(k-1,j,i) * dxu(j,i+1) * u(k-1,j,i+1) ) &
               - hlf * ( &
               + zydx(k-1,j,i) * dyv(j  ,i) * v(k-1,j  ,i) &
               + zydx(k-1,j,i) * dyv(j+1,i) * v(k-1,j+1,i) )

       enddo
    enddo

    do i = 1,nx
       do j = 1,ny 
          do k = 1,nz
             rhs(k,j,i) = rhs(k,j,i) + wf(k+1,j,i) - wf(k,j,i)
          enddo
       enddo
    enddo

    wf => null()

  end subroutine set_rhs

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
    real(kind=rp), dimension(:,:),     pointer :: dxu,dyv
    real(kind=rp), dimension(:,:,:),   pointer :: dzw
    real(kind=rp), dimension(:,:,:),   pointer :: Arx,Ary
    real(kind=rp), dimension(:,:)  ,   pointer :: Arz
    real(kind=rp), dimension(:,:,:),   pointer :: zxdy,zydx
    real(kind=rp), dimension(:,:,:),   pointer :: alpha
    real(kind=rp), dimension(:,:),     pointer :: beta
    real(kind=rp), dimension(:,:),     pointer :: gamu
    real(kind=rp), dimension(:,:),     pointer :: gamv
    real(kind=rp), dimension(:,:,:,:), pointer :: cA

    if (myrank==0) write(*,*)'   - set matrices:'

    do lev = 1, nlevs

       nx = grid(lev)%nx
       ny = grid(lev)%ny
       nz = grid(lev)%nz

       dx    => grid(lev)%dx    !
       dy    => grid(lev)%dy    !
       dxu   => grid(lev)%dxu   !
       dyv   => grid(lev)%dyv   !
       dzw   => grid(lev)%dzw   !
       Arx   => grid(lev)%Arx   !
       Ary   => grid(lev)%Ary   !
       Arz   => grid(lev)%Arz   !
       alpha => grid(lev)%alpha !
       beta  => grid(lev)%beta  !
       gamu  => grid(lev)%gamu  !
       gamv  => grid(lev)%gamv  !
       zxdy  => grid(lev)%zxdy  !
       zydx  => grid(lev)%zydx  !
       cA    => grid(lev)%cA    !

       !! interaction coeff with neighbours !!

       !---------------!
       !- lower level -!
       !---------------!
       k = 1
       do i = 1,nx
          do j = 1,ny+1
             ! couples with k+1,j-1
             cA(3,k,j,i) = qrt * ( zydx(k+1,j,i) + zydx(k,j-1,i) )
             ! couples with j-1
             cA(4,k,j,i) = hlf * (gamv(j,i) + gamv(j-1,i)) * Ary(k,j,i) / dyv(j,i) &
                  - qrt * ( zydx(k,j-1,i) - zydx(k,j,i) )
          enddo
       enddo

       do i = 1,nx+1
          do j = 1,ny
             ! couples with k+1,i-1
             cA(6,k,j,i) = qrt * ( zxdy(k+1,j,i) + zxdy(k,j,i-1) )
             ! couples with i-1
             cA(7,k,j,i) = hlf * (gamu(j,i) + gamu(j,i-1))  * Arx(k,j,i) / dxu(j,i) &
                  - qrt * ( zxdy(k,j,i-1) - zxdy(k,j,i) )
          enddo
       enddo

       do i = 1,nx+1
          do j = 0,ny
             ! only for k==1, couples with j+1,i-1
             cA(5,k,j,i) = beta(j,i-1) + beta(j+1,i)
          enddo
       enddo

       do i = 1,nx+1
          do j = 1,ny+1
             ! only for k==1, couples with j-1,i-1
             cA(8,k,j,i) = - beta(j,i-1) - beta(j-1,i)
          enddo
       enddo

       !-------------------!
       !- interior levels -!
       !-------------------!
       do i = 1,nx
          do j = 1,ny
             do k = 2,nz-1 
                ! couples with k-1
                cA(2,k,j,i) = &
                     ( Arz(j,i) / dzw(k,j,i) ) * &
                     hlf * (alpha(k,j,i)+alpha(k-1,j,i))
             enddo
          enddo
       enddo

       do i = 1,nx
          do j = 1,ny+1
             do k = 2,nz-1 
                ! couples with k+1,j-1
                cA(3,k,j,i) = qrt * ( zydx(k+1,j,i) + zydx(k,j-1,i) )
                ! couples with j-1
                cA(4,k,j,i) =  Ary(k,j,i) / dyv(j,i) 
                ! couples with k-1,j-1
                cA(5,k,j,i) = - qrt * ( zydx(k-1,j,i) + zydx(k,j-1,i) )
             enddo
          enddo
       enddo

       do i = 1,nx+1
          do j = 1,ny 
             do k = 2,nz-1 
                ! couples with k+1,i-1
                cA(6,k,j,i) = qrt * ( zxdy(k+1,j,i) + zxdy(k,j,i-1) )
                ! couples with i-1
                cA(7,k,j,i) = Arx(k,j,i) / dxu(j,i) 
                ! couples with k-1,i-1
                cA(8,k,j,i) = - qrt * ( zxdy(k-1,j,i) + zxdy(k,j,i-1) )
             enddo
          enddo
       enddo

       !---------------!
       !- upper level -!
       !---------------!
       k = nz
       do i = 1,nx    
          do j = 1,ny 
             ! couples with k-1
             cA(2,k,j,i) = (Arz(j,i)/dzw(k,j,i)) * hlf * ( alpha(k,j,i) + alpha(k-1,j,i) )
          enddo
       enddo

       do i = 1,nx
          do j = 1,ny+1
             ! couples with j-1
             cA(4,k,j,i) = Ary(k,j,i) / dyv(j,i) &
                  + qrt * ( - zydx(k,j-1,i) + zydx(k,j,i) )
             ! couples with k-1,j-1
             cA(5,k,j,i) = - qrt * ( zydx(k-1,j,i) + zydx(k,j-1,i) )
          enddo
       enddo

       do i = 1,nx+1
          do j = 1,ny 
             ! couples with i-1
             cA(7,k,j,i) = Arx(k,j,i) / dxu(j,i) &
                  + qrt * ( -zxdy(k,j,i-1) + zxdy(k,j,i) )
             ! couples with k-1,i-1
             cA(8,k,j,i) = - qrt * ( zxdy(k-1,j,i) + zxdy(k,j,i-1) )
          enddo
       enddo

       call fill_halo(lev,cA)

       !! self-interaction coeff !!

       do i = 1,nx
          do j = 1,ny

             k = 1 !lower level
             cA(1,k,j,i) = &
                  - Arz(j,i) / dzw(k+1,j,i) * hlf * (alpha(k+1,j,i) + alpha(k,j,i)) &
                  - Arx(k,j,i  )/dxu(j,i  ) * hlf * (gamu(j,i) + gamu(j  ,i-1)) &
                  - Arx(k,j,i+1)/dxu(j,i+1) * hlf * (gamu(j,i) + gamu(j  ,i+1)) &
                  - Ary(k,j  ,i)/dyv(j  ,i) * hlf * (gamv(j,i) + gamv(j-1,i  )) &
                  - Ary(k,j+1,i)/dyv(j+1,i) * hlf * (gamv(j,i) + gamv(j+1,i  ))

             do k = 2,nz-1 !interior levels
                cA(1,k,j,i) = &
                     - Arz(j,i) / dzw(k+1,j,i) * hlf * (alpha(k+1,j,i) + alpha(k,j,i)) &
                     - Arz(j,i) / dzw(k  ,j,i) * hlf * (alpha(k-1,j,i) + alpha(k,j,i)) &
                     - Arx(k,j,i  )/dxu(j,i  )  &
                     - Arx(k,j,i+1)/dxu(j,i+1)  &
                     - Ary(k,j  ,i)/dyv(j  ,i)  &
                     - Ary(k,j+1,i)/dyv(j+1,i)
             enddo

             k=nz ! upper level
             cA(1,k,j,i) = &
                  - Arz(j,i) / dzw(k+1,j,i) * alpha(k,j,i) &
                  - Arz(j,i) / dzw(k  ,j,i) * hlf * (alpha(k-1,j,i) + alpha(k,j,i)) &
                  - Arx(k,j,i  )/dxu(j,i  )  &
                  - Arx(k,j,i+1)/dxu(j,i+1)  &
                  - Ary(k,j  ,i)/dyv(j  ,i)  &
                  - Ary(k,j+1,i)/dyv(j+1,i)

          enddo
       enddo

       if (netcdf_output) then
          if (myrank==0) write(*,*)'       write cA in a netcdf file'
          call write_netcdf(grid(lev)%cA,vname='ca',netcdf_file_name='cA.nc',rank=myrank,iter=lev)
       endif

    enddo

  end subroutine set_matrices

  !-------------------------------------------------------------------------     
  subroutine correct_uvw(u,v,w)

    real(kind=rp), dimension(:,:,:), pointer, intent(inout) :: u,v,w

    integer(kind=ip):: k, j, i
    integer(kind=ip):: nx, ny, nz

    real(kind=rp), dimension(:,:)  , pointer :: dxu,dyv
    real(kind=rp), dimension(:,:,:), pointer :: dzw
    real(kind=rp), dimension(:,:,:), pointer :: p

    nx = grid(1)%nx
    ny = grid(1)%ny
    nz = grid(1)%nz

    dxu => grid(1)%dxu
    dyv => grid(1)%dyv
    dzw => grid(1)%dzw

    if (myrank==0) write(*,*)'   - correct u,v,w:'

    !! Correct
    p => grid(1)%p

    do i = 1,nx+1
!    do i = 2,nx
!       do j = 0,ny+1
        do j = 1,ny
          do k = 1,nz
             u(k,j,i) = u(k,j,i) - one / dxu(j,i) * (p(k,j,i)-p(k,j,i-1))
          enddo
       enddo
    enddo

!    do i = 0,nx+1
    do i = 1,nx
       do j = 1,ny+1 
!       do j = 2,ny 
          do k = 1,nz
          v(k,j,i) = v(k,j,i) - one / dyv(j,i) * (p(k,j,i)-p(k,j-1,i  ))
          enddo
       enddo
    enddo

!    do i = 0,nx+1
    do i = 1,nx
!       do j = 0,ny+1
       do j = 1,ny

          k = 1 !bottom

          do k = 2,nz !interior levels
             w(k,j,i) = w(k,j,i) - one / dzw(k,j,i) * (p(k,j,i)-p(k-1,j,i))
          enddo

          k = nz+1 !surface
          w(k,j,i) = w(k,j,i) - one / dzw(k,j,i) * (-p(k-1,j,i))

       enddo
    enddo

    call fill_halo(1,w)

  end subroutine correct_uvw

end module mg_projection
