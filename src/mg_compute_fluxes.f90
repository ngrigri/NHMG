module mg_compute_fluxes

  use mg_mpi
  use mg_tictoc
  use mg_namelist
  use mg_grids
  use mg_mpi_exchange
  use mg_netcdf_out

  implicit none

contains
  !-------------------------------------------------------------------------     
  subroutine compute_fluxes(u,v,w,uf,vf)

    real(kind=rp), dimension(:,:,:), pointer, intent(in)  :: u,v,w
    real(kind=rp), dimension(:,:,:), pointer, intent(out) :: uf,vf

    integer(kind=ip):: k,j,i
    integer(kind=ip):: nx,ny,nz

    real(kind=rp), dimension(:,:)  , pointer :: dx,dy
    real(kind=rp), dimension(:,:,:), pointer :: zr,zw

    real(kind=rp), dimension(:,:,:), pointer :: dzw
    real(kind=rp), dimension(:,:,:), pointer :: Arx,Ary
    real(kind=rp), dimension(:,:,:), pointer :: zxdy,zydx
    real(kind=rp), dimension(:,:,:), pointer :: cw

    real(kind=rp), dimension(:,:,:), pointer :: uf_tmp,vf_tmp !!dirty

    real(kind=rp), parameter :: two  = 2._rp
    real(kind=rp), parameter :: one  = 1._rp
    real(kind=rp), parameter :: hlf  = 0.5_rp
    real(kind=rp), parameter :: qrt  = 0.25_rp
    real(kind=rp), parameter :: zero = 0._rp

    nx = grid(1)%nx
    ny = grid(1)%ny
    nz = grid(1)%nz

    dx => grid(1)%dx
    dy => grid(1)%dy
    zr => grid(1)%zr
    zw => grid(1)%zw

    dzw => grid(1)%dzw
    Arx => grid(1)%Arx
    Ary => grid(1)%Ary
    cw  => grid(1)%cw
    zxdy => grid(1)%zxdy
    zydx => grid(1)%zydx

    !- uf -!

    do i = 1,nx+1  
       do j = 1,ny 
!       do j = 0,ny+1

          k = 1 ! lower level
          uf(k,j,i) =  &
               qrt                                                      * & 
               ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j,i-1) - zw(k,j,i-1) ) * &
               ( dy(j,i) + dy(j,i-1) ) * u(k,j,i) &
               - qrt * ( &
               + zxdy(k,j,i  )* two * dzw(k  ,j,i  )*w(k  ,j,i  ) & 
               + zxdy(k,j,i  )*       dzw(k+1,j,i  )*w(k+1,j,i  ) &
               + zxdy(k,j,i-1)* two * dzw(k  ,j,i-1)*w(k  ,j,i-1) &
               + zxdy(k,j,i-1)*       dzw(k+1,j,i-1)*w(k+1,j,i-1) &
               )  ! umask

!!$               - qrt * ( &
!!$               + zxdy(k,j,i  )*dzw(k+1,j,i  )*w(k+1-1,j,i)   & !note the index trick for w 
!!$               + zxdy(k,j,i-1)*dzw(k+1,j,i-1)*w(k+1-1,j,i-1) ) & !note the index trick for w 
!!$               -( &
!!$               + zxdy(k,j,i  )*zxdy(k,j,i  )/(cw(k,j,i  )+cw(k+1,j,i  )) &
!!$               + zxdy(k,j,i-1)*zxdy(k,j,i-1)/(cw(k,j,i-1)+cw(k+1,j,i-1)) ) * &
!!$               (hlf * ( dx(j,i) + dx(j,i-1) )) * u(k,j,i) &
!!$               -( & 
!!$               + zxdy(k,j,i) * zydx(k,j,i)/(cw(k,j,i  )+cw(k+1,j,i  ))   &
!!$               * hlf * ( &
!!$               hlf * ( dy(j  ,i) + dy(j-1,i) )*v(k,j,i ) + &
!!$               hlf * ( dy(j+1,i) + dy(j  ,i) )*v(k,j+1,i ))   & 
!!$               + zxdy(k,j,i-1) * zydx(k,j,i-1)/(cw(k,j,i-1)+cw(k+1,j,i-1))   &
!!$               * hlf * ( &
!!$               hlf * ( dy(j  ,i-1) + dy(j-1,i-1) )*v(k,j,i-1) + &
!!$               hlf * ( dy(j+1,i-1) + dy(j  ,i-1) )*v(k,j+1,i-1)) )

          do k = 2,nz-1 !interior levels
             uf(k,j,i) =  qrt                                               * & 
                  ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j,i-1) - zw(k,j,i-1) ) * &
                  ( dy(j,i) + dy(j,i-1) ) * u(k,j,i) &
                  - qrt * ( &
                  + zxdy(k,j,i  ) * dzw(k  ,j,i  )*w(k  ,j,i  ) &
                  + zxdy(k,j,i  ) * dzw(k+1,j,i  )*w(k+1,j,i  ) &
                  + zxdy(k,j,i-1) * dzw(k  ,j,i-1)*w(k  ,j,i-1) &
                  + zxdy(k,j,i-1) * dzw(k+1,j,i-1)*w(k+1,j,i-1) &
                  )  ! umask
          enddo

          k = nz !upper level
          uf(k,j,i) = &
               qrt                                                       * & 
               ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j,i-1) - zw(k,j,i-1) ) * &
               ( dy(j,i) + dy(j,i-1) ) * u(k,j,i) &
               - qrt * ( &
               + zxdy(k,j,i  )*       dzw(k  ,j,i  )*w(k  ,j,i  ) & 
               + zxdy(k,j,i  )* two * dzw(k+1,j,i  )*w(k+1,j,i  ) & 
               + zxdy(k,j,i-1)*       dzw(k  ,j,i-1)*w(k  ,j,i-1) & 
               + zxdy(k,j,i-1)* two * dzw(k+1,j,i-1)*w(k+1,j,i-1) & 
               )  ! umask

       enddo
    enddo

    allocate(uf_tmp(1:nz,0:ny+1,0:nx+1))
    uf_tmp(:,1:ny,1:nx+1)=uf(:,1:ny,1:nx+1)
    call fill_halo(1,uf_tmp,lbc_null='u')
    uf(:,0:ny+1,1:nx+1)=uf_tmp(:,0:ny+1,1:nx+1)
    deallocate(uf_tmp)

    !- vf -!

    do i = 1,nx
!    do i = 0,nx+1
       do j = 1,ny+1

          k = 1 !lower level
          vf(k,j,i) = &
               qrt                                                       * &
               ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j-1,i) - zw(k,j-1,i) ) * &
               ( dx(j,i) + dx(j-1,i) ) * v(k,j,i) &
               - qrt * ( &
               + zydx(k,j  ,i)* two * dzw(k  ,j  ,i)*w(k  ,j  ,i) &  
               + zydx(k,j  ,i)*       dzw(k+1,j  ,i)*w(k+1,j  ,i) &
               + zydx(k,j-1,i)* two * dzw(k  ,j-1,i)*w(k  ,j-1,i) &
               + zydx(k,j-1,i)*       dzw(k+1,j-1,i)*w(k+1,j-1,i) &
               ) !* vmask(j,i)

!!$               - qrt * ( &
!!$               + zydx(k,j  ,i) * dzw(k+1,j  ,i)*w(k+1-1,j,i) & !note the index trick for w 
!!$               + zydx(k,j-1,i) * dzw(k+1,j-1,i)*w(k+1-1,j-1,i) ) & !note the index trick for w 
!!$               -( &
!!$               + zydx(k,j  ,i) * zydx(k,j  ,i)/(cw(k,j  ,i)+cw(k+1,j  ,i)) &
!!$               + zydx(k,j-1,i) * zydx(k,j-1,i)/(cw(k,j-1,i)+cw(k+1,j-1,i)) ) * &
!!$               hlf * ( dy(j,i) + dy(j-1,i) ) * v(k,j,i) &
!!$               - ( &
!!$               + zxdy(k,j  ,i) * zydx(k,j  ,i)/(cw(k,j  ,i)+cw(k+1,j  ,i))  &
!!$               * hlf * ( &
!!$               hlf * ( dx(j,i  ) + dx(j,i-1) ) * u(k,j,i) + &
!!$               hlf * ( dx(j,i+1) + dx(j,i  ) ) * u(k,j,i+1)) &
!!$               + zxdy(k,j-1,i) * zydx(k,j-1,i)/(cw(k,j-1,i)+cw(k+1,j-1,i))   &
!!$               * hlf * ( &
!!$               hlf * ( dx(j-1,i  ) + dx(j-1,i-1) ) * u(k,j-1,i) + &
!!$               hlf * ( dx(j-1,i+1) + dx(j-1,i  ) ) * u(k,j-1,i+1)) &
!!$               ) 

          do k = 2,nz-1 !interior levels
             vf(k,j,i) = &
                  qrt                                                       * &
                  ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j-1,i) - zw(k,j-1,i) ) * &
                  ( dx(j,i) + dx(j-1,i) ) * v(k,j,i) &
                  - qrt * ( &
                  + zydx(k,j  ,i) * dzw(k  ,j  ,i)*w(k  ,j  ,i) & 
                  + zydx(k,j  ,i) * dzw(k+1,j  ,i)*w(k+1,j  ,i) & 
                  + zydx(k,j-1,i) * dzw(k  ,j-1,i)*w(k  ,j-1,i) & 
                  + zydx(k,j-1,i) * dzw(k+1,j-1,i)*w(k+1,j-1,i) & 
                  )  !* vmask(j,i)

          enddo

          k = nz !upper level
          vf(k,j,i) =  &
               qrt                                                       * &
               ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j-1,i) - zw(k,j-1,i) ) * &
               ( dx(j,i) + dx(j-1,i) ) * v(k,j,i) &
               - qrt * ( &
               + zydx(k,j  ,i)*       dzw(k  ,j  ,i)*w(k  ,j  ,i) &
               + zydx(k,j  ,i)* two * dzw(k+1,j  ,i)*w(k+1,j  ,i) & 
               + zydx(k,j-1,i)*       dzw(k  ,j-1,i)*w(k  ,j-1,i) &
               + zydx(k,j-1,i)* two * dzw(k+1,j-1,i)*w(k+1,j-1,i) & 
               ) !* vmask(j,i)

       enddo
    enddo

    allocate(vf_tmp(1:nz,0:ny+1,0:nx+1))
    vf_tmp(:,1:ny+1,1:nx)=vf(:,1:ny+1,1:nx)
    call fill_halo(1,vf_tmp,lbc_null='v')
    vf(:,1:ny+1,0:nx+1)=vf_tmp(:,1:ny+1,0:nx+1)
    deallocate(vf_tmp)

  end subroutine compute_fluxes
 
end module mg_compute_fluxes

