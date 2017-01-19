module mg_fluxes

  use mg_cst
  use mg_mpi
  use mg_tictoc
  use mg_namelist
  use mg_grids
  use mg_mpi_exchange
  use mg_netcdf_out

  implicit none

contains
  !-------------------------------------------------------------------------     
  subroutine set_fluxes(u,v,w,uf,vf)

    real(kind=rp), dimension(:,:,:), pointer, intent(in)  :: u,v,w
    real(kind=rp), dimension(:,:,:), pointer, intent(out) :: uf,vf

    integer(kind=ip) :: k,j,i
    integer(kind=ip) :: nx,ny,nz
    real(kind=rp)    :: gamma 

    real(kind=rp), dimension(:,:)  , pointer :: dx,dy
    real(kind=rp), dimension(:,:)  , pointer :: dxu,dyv
    real(kind=rp), dimension(:,:,:), pointer :: dzw
    real(kind=rp), dimension(:,:,:), pointer :: Arx,Ary
    real(kind=rp), dimension(:,:,:), pointer :: zxdy,zydx
    real(kind=rp), dimension(:,:,:), pointer :: alpha
    real(kind=rp), dimension(:,:)  , pointer :: beta

    nx = grid(1)%nx
    ny = grid(1)%ny
    nz = grid(1)%nz

    dx    => grid(1)%dx
    dy    => grid(1)%dy
    dxu   => grid(1)%dxu
    dyv   => grid(1)%dyv
    dzw   => grid(1)%dzw
    Arx   => grid(1)%Arx
    Ary   => grid(1)%Ary
    alpha => grid(1)%alpha
    beta  => grid(1)%beta
    zxdy  => grid(1)%zxdy
    zydx  => grid(1)%zydx

    !- uf -!

    do i = 1,nx+1  
       do j = 0,ny+1

          k = 1
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

          do k = 2,nz-1 
             uf(k,j,i) = Arx(k,j,i) * u(k,j,i) &
                  - qrt * ( &
                  + zxdy(k,j,i  ) * dzw(k  ,j,i  ) * w(k  ,j,i  ) &
                  + zxdy(k,j,i  ) * dzw(k+1,j,i  ) * w(k+1,j,i  ) &
                  + zxdy(k,j,i-1) * dzw(k  ,j,i-1) * w(k  ,j,i-1) &
                  + zxdy(k,j,i-1) * dzw(k+1,j,i-1) * w(k+1,j,i-1) )
          enddo

          k = nz
          uf(k,j,i) = Arx(k,j,i) * u(k,j,i) &
               - qrt * ( &
               + zxdy(k,j,i  ) *       dzw(k  ,j,i  ) * w(k  ,j,i  ) &
               + zxdy(k,j,i  ) * two * dzw(k+1,j,i  ) * w(k+1,j,i  ) &
               + zxdy(k,j,i-1) *       dzw(k  ,j,i-1) * w(k  ,j,i-1) &
               + zxdy(k,j,i-1) * two * dzw(k+1,j,i-1) * w(k+1,j,i-1) )
       enddo
    enddo

    !- vf -!

    do i = 0,nx+1
       do j = 1,ny+1

          k = 1
          gamma = one - qrt * ( &
               (zydx(k,j  ,i)/dx(j  ,i))**2/alpha(k,j,i  ) + &
               (zydx(k,j-1,i)/dx(j-1,i))**2/alpha(k,j-1,i) )
          vf(k,j,i) = gamma * Ary(k,j,i) * v(k,j,i) &
               - qrt * ( &
               + zydx(k,j  ,i) * dzw(k+1,j  ,i) * w(k+1,j  ,i  ) &
               + zydx(k,j-1,i) * dzw(k+1,j-1,i) * w(k+1,j-1,i  ) ) &
               - beta(j-1,i)   * dxu(j-1,i  )   * u(k  ,j-1,i  ) &
               - beta(j-1,i)   * dxu(j-1,i+1)   * u(k  ,j-1,i+1) &
               - beta(j  ,i)   * dxu(j  ,i  )   * u(k  ,j  ,i  ) &
               - beta(j  ,i)   * dxu(j  ,i+1)   * u(k  ,j  ,i+1)

          do k = 2,nz-1
             vf(k,j,i) = Ary(k,j,i) * v(k,j,i) &
                  - qrt * ( &
                  + zydx(k,j  ,i) * dzw(k  ,j  ,i) * w(k  ,j  ,i) &
                  + zydx(k,j  ,i) * dzw(k+1,j  ,i) * w(k+1,j  ,i) &
                  + zydx(k,j-1,i) * dzw(k  ,j-1,i) * w(k  ,j-1,i) &
                  + zydx(k,j-1,i) * dzw(k+1,j-1,i) * w(k+1,j-1,i) )
          enddo

          k = nz
          vf(k,j,i) = Ary(k,j,i) * v(k,j,i) &
               - qrt * ( &
               + zydx(k,j  ,i)       * dzw(k  ,j  ,i) * w(k  ,j  ,i) &
               + zydx(k,j  ,i) * two * dzw(k+1,j  ,i) * w(k+1,j  ,i) &
               + zydx(k,j-1,i)       * dzw(k  ,j-1,i) * w(k  ,j-1,i) &
               + zydx(k,j-1,i) * two * dzw(k+1,j-1,i) * w(k+1,j-1,i) ) 
       enddo
    enddo

  end subroutine set_fluxes
 
end module mg_fluxes

