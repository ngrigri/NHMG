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
  subroutine compute_fluxes(zr,zw,u,v,w,uf,vf)

    real(kind=rp), dimension(:,:,:), pointer, intent(in)  :: zr,zw
    real(kind=rp), dimension(:,:,:), pointer, intent(in)  :: u,v,w
    real(kind=rp), dimension(:,:,:), pointer, intent(out) :: uf,vf

    integer(kind=ip):: k,j,i
    integer(kind=ip):: nx,ny,nz

    real(kind=rp), dimension(:,:)  , pointer :: dx,dy
!    real(kind=rp), dimension(:,:,:), pointer :: zr,zw
!    real(kind=rp), dimension(:,:,:), pointer :: dzw
!    real(kind=rp), dimension(:,:,:), pointer :: zydx,zxdy
!    real(kind=rp), dimension(:,:,:), pointer :: cw

    real(kind=rp), dimension(:,:)  ,   pointer :: dxu,dyv
    real(kind=rp), dimension(:,:,:),   pointer :: Arx,Ary
    real(kind=rp), dimension(:,:)  ,   pointer :: Arz
    real(kind=rp), dimension(:,:,:),   pointer :: dz,dzw
    real(kind=rp), dimension(:,:,:),   pointer :: zy,zx
    real(kind=rp), dimension(:,:,:),   pointer :: zydx,zxdy
    real(kind=rp), dimension(:,:,:),   pointer :: zxw,zyw
    real(kind=rp), dimension(:,:,:),   pointer :: cw

    real(kind=rp), dimension(:,:,:),   pointer :: uf_tmp,vf_tmp

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

    !! Cell heights
    allocate(dz(nz,0:ny+1,0:nx+1))
    do i = 0,nx+1
       do j = 0,ny+1
          do k = 1,nz
             dz(k,j,i) = zw(k+1,j,i)-zw(k,j,i)
          enddo
       enddo
    enddo
    allocate(dzw(nz+1,0:ny+1,0:nx+1))
    do i = 0,nx+1
       do j = 0,ny+1
          dzw(1,j,i) = zr(1,j,i)-zw(1,j,i) 
          do k = 2,nz
             dzw(k,j,i) = zr(k,j,i)-zr(k-1,j,i) 
          enddo
          dzw(nz+1,j,i) = zw(nz+1,j,i)-zr(nz,j,i) 
       enddo
    enddo
    !! Cell widths
    allocate(dxu(0:ny+1,nx+1))
    do i = 1,nx+1
       do j = 0,ny+1
          dxu(j,i) = 0.5_8*(dx(j,i)+dx(j,i-1))
       enddo
    enddo
    allocate(dyv(ny+1,0:nx+1))
    do i = 0,nx+1
       do j = 1,ny+1
          dyv(j,i) = 0.5_8*(dy(j,i)+dy(j-1,i))
       enddo
    enddo
    !!  Areas
    allocate(Arx(nz,0:ny+1,nx+1))
    do i = 1,nx+1
       do j = 0,ny+1
          do k = 1,nz
             Arx(k,j,i) = 0.25_8*(dz(k,j,i)+dz(k,j,i-1))*(dy(j,i)+dy(j,i-1)) 
          enddo
       enddo
    enddo
    allocate(Ary(nz,ny+1,0:nx+1))
    do i = 0,nx+1
       do j = 1,ny+1
          do k = 1,nz
             Ary(k,j,i) = 0.25_8*(dz(k,j,i)+dz(k,j-1,i))*(dx(j,i)+dx(j-1,i)) 
          enddo
       enddo
    enddo
    allocate(Arz(0:ny+1,0:nx+1))
    do i = 0,nx+1
       do j = 0,ny+1
          Arz(j,i) = dx(j,i)*dy(j,i)
       enddo
    enddo
    !! Slopes in x- and y-direction defined at rho-points
    allocate(zx(nz,0:ny+1,0:nx+1))
    allocate(zy(nz,0:ny+1,0:nx+1))
    do i = 1,nx
       do j = 1,ny
          do k = 1,nz
             zy(k,j,i) = 0.5_8*(zr(k,j+1,i)-zr(k,j-1,i))/dy(j,i)
             zx(k,j,i) = 0.5_8*(zr(k,j,i+1)-zr(k,j,i-1))/dx(j,i)
          enddo
       enddo
    enddo
    call fill_halo(1,zy) ! copy interior value in the halo
    call fill_halo(1,zx) ! copy interior value in the halo
    allocate(zxdy(nz,0:ny+1,0:nx+1))
    allocate(zydx(nz,0:ny+1,0:nx+1))
    do k = 1,nz
       zydx(k,:,:) = zy(k,:,:)*dx(:,:)
       zxdy(k,:,:) = zx(k,:,:)*dy(:,:)
    enddo
    allocate(zyw(nz+1,0:ny+1,0:nx+1))
    allocate(zxw(nz+1,0:ny+1,0:nx+1))
    do i = 1,nx
       do j = 1,ny
          do k = 1,nz+1
             zyw(k,j,i) = 0.5_8*(zw(k,j+1,i)-zw(k,j-1,i))/dy(j,i)
             zxw(k,j,i) = 0.5_8*(zw(k,j,i+1)-zw(k,j,i-1))/dx(j,i)
          enddo
       enddo
    enddo
    call fill_halo(1,zyw) ! copy interior value in the halo
    call fill_halo(1,zxw) ! copy interior value in the halo
    allocate(cw(nz+1,0:ny+1,0:nx+1))
    do i = 0,nx+1
       do j = 0,ny+1
          do k = 1,nz+1
             cw(k,j,i) = Arz(j,i)/dzw(k,j,i) * (1._8 + zxw(k,j,i)*zxw(k,j,i)+zyw(k,j,i)*zyw(k,j,i))
          enddo
       enddo
    enddo

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

    deallocate(dz)
    deallocate(dzw)
    deallocate(dxu)
    deallocate(dyv)
    deallocate(Arx)
    deallocate(Ary)
    deallocate(Arz)
    deallocate(zx)
    deallocate(zy)
    deallocate(zxdy)
    deallocate(zydx)
    deallocate(zxw)
    deallocate(zyw)
    deallocate(cw)

  end subroutine compute_fluxes
 
end module mg_compute_fluxes

