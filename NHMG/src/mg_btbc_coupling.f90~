#include "cppdefs.h"
#ifdef NHMG
module mg_btbc_coupling

  use mg_mpi
  use mg_tictoc
  use mg_namelist
  use mg_grids
  use mg_mpi_exchange
  use mg_netcdf_out

  implicit none

contains
  !-------------------------------------------------------------------------     
  subroutine btbc_coupling(zr,zw,uf_bar,vf_bar,u,v,w,uf,vf)

    real(kind=rp), dimension(:,:),   pointer, intent(in)  :: uf_bar,vf_bar
    real(kind=rp), dimension(:,:,:), pointer, intent(in)  :: zr,zw
    real(kind=rp), dimension(:,:,:), pointer, intent(inout)  :: u,v,w
    real(kind=rp), dimension(:,:,:), pointer, optional, intent(out):: uf,vf

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

    real(kind=rp), dimension(:,:)  ,   pointer :: su_integr,sv_integr
    real(kind=rp), dimension(:,:)  ,   pointer :: uf_integr,vf_integr
!!$    real(kind=rp), dimension(:,:)  ,   pointer :: uf_integr1,vf_integr1,uf_integr2,vf_integr2
    real(kind=rp), dimension(:,:,:),   pointer :: uf_integr_tmp,vf_integr_tmp
    real(kind=rp), dimension(:,:,:),   pointer :: wc

    real(kind=rp), parameter :: two  = 2._rp
    real(kind=rp), parameter :: one  = 1._rp
    real(kind=rp), parameter :: hlf  = 0.5_rp
    real(kind=rp), parameter :: qrt  = 0.25_rp
    real(kind=rp), parameter :: zero = 0._rp

    integer(kind=ip), save :: iter_coupling_in=0
    iter_coupling_in = iter_coupling_in + 1

    nx = grid(1)%nx
    ny = grid(1)%ny
    nz = grid(1)%nz

    dx => grid(1)%dx
    dy => grid(1)%dy

!    zr => grid(1)%zr
!    zw => grid(1)%zw
!    dzw => grid(1)%dzw
!    zxdy => grid(1)%zxdy
!    zydx => grid(1)%zydx
!    cw  => grid(1)%cw

!XXX
    !! Cell heights
    allocate(dz(nz,0:ny+1,0:nx+1))
    do i = 0,nx+1
       do j = 0,ny+1
          do k = 1,nz
!             dz(k,j,i) = zw(k+1,j,i)-zw(k,j,i)
             dz(k,j,i) = zw(k,j,i)-zw(k-1,j,i) !because zw indexed as croco from 0 to nz
          enddo
       enddo
    enddo
    allocate(dzw(nz+1,0:ny+1,0:nx+1))
    do i = 0,nx+1
       do j = 0,ny+1
!          dzw(1,j,i) = zr(1,j,i)-zw(1,j,i) 
          dzw(1,j,i) = zr(1,j,i)-zw(0,j,i) !because zw indexed as croco from 0 to nz
          do k = 2,nz
             dzw(k,j,i) = zr(k,j,i)-zr(k-1,j,i) 
          enddo
!          dzw(nz+1,j,i) = zw(nz+1,j,i)-zr(nz,j,i) 
         dzw(nz+1,j,i) = zw(nz,j,i)-zr(nz,j,i) !because zw indexed as croco from 0 to nz
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
    call fill_halo(1,zy)
    call fill_halo(1,zx) 
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
!             zyw(k,j,i) = 0.5_8*(zw(k,j+1,i)-zw(k,j-1,i))/dy(j,i)
             zyw(k,j,i) = 0.5_8*(zw(k-1,j+1,i)-zw(k-1,j-1,i))/dy(j,i) !because zw indexed as croco from 0 to nz
!             zxw(k,j,i) = 0.5_8*(zw(k,j,i+1)-zw(k,j,i-1))/dx(j,i)
             zxw(k,j,i) = 0.5_8*(zw(k-1,j,i+1)-zw(k-1,j,i-1))/dx(j,i) !because zw indexed as croco from 0 to nz
          enddo
       enddo
    enddo
    call fill_halo(1,zyw)
    call fill_halo(1,zxw)
    allocate(cw(nz+1,0:ny+1,0:nx+1))
    do i = 0,nx+1
       do j = 0,ny+1
          do k = 1,nz+1
             cw(k,j,i) = Arz(j,i)/dzw(k,j,i) * (1._8 + zxw(k,j,i)*zxw(k,j,i)+zyw(k,j,i)*zyw(k,j,i))
          enddo
       enddo
    enddo
!XXX

    !------------------------------------------
    !1 - compute integrated transport anomalies

    allocate(su_integr(0:ny+1,0:nx+1))
    allocate(sv_integr(0:ny+1,0:nx+1))
!!$    allocate(uf_integr1(0:ny+1,0:nx+1))
!!$    allocate(vf_integr1(0:ny+1,0:nx+1))
!!$    allocate(uf_integr2(0:ny+1,0:nx+1))
!!$    allocate(vf_integr2(0:ny+1,0:nx+1))
    allocate(uf_integr(0:ny+1,0:nx+1))
    allocate(vf_integr(0:ny+1,0:nx+1))
    su_integr(:,:) = zero
    sv_integr(:,:) = zero
!!$    uf_integr1(:,:) = zero
!!$    vf_integr1(:,:) = zero
!!$    uf_integr2(:,:) = zero
!!$    vf_integr2(:,:) = zero
    uf_integr(:,:) = zero
    vf_integr(:,:) = zero

    do i = 1,nx+1  
       do j = 1,ny 
!       do j = 0,ny +1

          k = 1 ! lower level
          su_integr(j,i) =  &
               qrt                                                      * & 
!               ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j,i-1) - zw(k,j,i-1) ) * &
               ( zw(k,j,i) - zw(k-1,j,i) + zw(k,j,i-1) - zw(k-1,j,i-1) ) * & !because zw indexed as croco from 0 to nz
               ( dy(j,i) + dy(j,i-1) )
!!$          uf_integr1(j,i) =  &
!!$               qrt                                                      * & 
!!$!               ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j,i-1) - zw(k,j,i-1) ) * &
!!$               ( zw(k,j,i) - zw(k-1,j,i) + zw(k,j,i-1) - zw(k-1,j,i-1) ) * & !because zw indexed as croco from 0 to nz
!!$               ( dy(j,i) + dy(j,i-1) ) * u(k,j,i)
!!$          uf_integr2(j,i) =  &
!!$               - qrt * ( &
!!$               + zxdy(k,j,i  )* two * dzw(k  ,j,i  )*w(k-1,j,i) & !note the index trick for w 
!!$               + zxdy(k,j,i  )*       dzw(k+1,j,i  )*w(k+1-1,j,i) & !note the index trick for w 
!!$               + zxdy(k,j,i-1)* two * dzw(k  ,j,i-1)*w(k-1,j,i-1) & !note the index trick for w 
!!$               + zxdy(k,j,i-1)*       dzw(k+1,j,i-1)*w(k+1-1,j,i-1) & !note the index trick for w 
!!$               )
!!$          uf_integr(j,i) =  uf_integr1(j,i) + uf_integr2(j,i)
          uf_integr(j,i) =  &
               qrt                                                      * & 
!               ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j,i-1) - zw(k,j,i-1) ) * &
               ( zw(k,j,i) - zw(k-1,j,i) + zw(k,j,i-1) - zw(k-1,j,i-1) ) * & !because zw indexed as croco from 0 to nz
               ( dy(j,i) + dy(j,i-1) ) * u(k,j,i) &
               - qrt * ( &
               + zxdy(k,j,i  )* two * dzw(k  ,j,i  )*w(k-1,j,i) & !note the index trick for w 
               + zxdy(k,j,i  )*       dzw(k+1,j,i  )*w(k+1-1,j,i) & !note the index trick for w 
               + zxdy(k,j,i-1)* two * dzw(k  ,j,i-1)*w(k-1,j,i-1) & !note the index trick for w 
               + zxdy(k,j,i-1)*       dzw(k+1,j,i-1)*w(k+1-1,j,i-1) & !note the index trick for w 
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

          do k = 2,nz-1 ! interior levels
             su_integr(j,i) = su_integr(j,i) + &
                  qrt                                               * & 
!                  ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j,i-1) - zw(k,j,i-1) ) * &
                  ( zw(k,j,i) - zw(k-1,j,i) + zw(k,j,i-1) - zw(k-1,j,i-1) ) * & !because zw indexed as croco from 0 to nz
                  ( dy(j,i) + dy(j,i-1) )
!!$             uf_integr1(j,i) = uf_integr1(j,i) + &
!!$                  qrt                                               * & 
!!$!                  ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j,i-1) - zw(k,j,i-1) ) * &
!!$                  ( zw(k,j,i) - zw(k-1,j,i) + zw(k,j,i-1) - zw(k-1,j,i-1) ) * & !because zw indexed as croco from 0 to nz
!!$                  ( dy(j,i) + dy(j,i-1) ) * u(k,j,i)
!!$             uf_integr2(j,i) = uf_integr2(j,i) + &
!!$                  - qrt * ( &
!!$                  + zxdy(k,j,i  ) * dzw(k  ,j,i  ) * w(k-1,j,i) & !note the index trick for w
!!$                  + zxdy(k,j,i  ) * dzw(k+1,j,i  ) * w(k+1-1,j,i) & !note the index trick for w
!!$                  + zxdy(k,j,i-1) * dzw(k  ,j,i-1) * w(k-1,j,i-1) & !note the index trick for w
!!$                  + zxdy(k,j,i-1) * dzw(k+1,j,i-1) * w(k+1-1,j,i-1) & !note the index trick for w
!!$                  )
!!$             uf_integr(j,i) =  uf_integr1(j,i) + uf_integr2(j,i)
             uf_integr(j,i) = uf_integr(j,i) + &
                  qrt                                               * & 
!                  ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j,i-1) - zw(k,j,i-1) ) * &
                  ( zw(k,j,i) - zw(k-1,j,i) + zw(k,j,i-1) - zw(k-1,j,i-1) ) * & !because zw indexed as croco from 0 to nz
                  ( dy(j,i) + dy(j,i-1) ) * u(k,j,i) &
                  - qrt * ( &
                  + zxdy(k,j,i  ) * dzw(k  ,j,i  ) * w(k-1,j,i) & !note the index trick for w
                  + zxdy(k,j,i  ) * dzw(k+1,j,i  ) * w(k+1-1,j,i) & !note the index trick for w
                  + zxdy(k,j,i-1) * dzw(k  ,j,i-1) * w(k-1,j,i-1) & !note the index trick for w
                  + zxdy(k,j,i-1) * dzw(k+1,j,i-1) * w(k+1-1,j,i-1) & !note the index trick for w
                  )  ! umask
          enddo

          k = nz ! upper level
          su_integr(j,i) = su_integr(j,i) + &
               qrt                                                       * & 
!               ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j,i-1) - zw(k,j,i-1) ) * &
               ( zw(k,j,i) - zw(k-1,j,i) + zw(k,j,i-1) - zw(k-1,j,i-1) ) * & !because zw indexed as croco from 0 to nz
               ( dy(j,i) + dy(j,i-1) )
!!$          uf_integr1(j,i) = uf_integr1(j,i) + &
!!$               qrt                                                       * & 
!!$!               ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j,i-1) - zw(k,j,i-1) ) * &
!!$               ( zw(k,j,i) - zw(k-1,j,i) + zw(k,j,i-1) - zw(k-1,j,i-1) ) * & !because zw indexed as croco from 0 to nz
!!$               ( dy(j,i) + dy(j,i-1) ) * u(k,j,i)
!!$          uf_integr2(j,i) = uf_integr2(j,i) + &
!!$               - qrt * ( &
!!$               + zxdy(k,j,i  )*       dzw(k  ,j,i  )*w(k-1,j,i) & !note the index trick for w
!!$               + zxdy(k,j,i  )* two * dzw(k+1,j,i  )*w(k+1-1,j,i) & !note the index trick for w
!!$               + zxdy(k,j,i-1)*       dzw(k  ,j,i-1)*w(k-1,j,i-1) & !note the index trick for w
!!$               + zxdy(k,j,i-1)* two * dzw(k+1,j,i-1)*w(k+1-1,j,i-1) & !note the index trick for w
!!$               )  
!!$          uf_integr(j,i) =  uf_integr1(j,i) + uf_integr2(j,i)
          uf_integr(j,i) = uf_integr(j,i) + &
               qrt                                                       * & 
!               ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j,i-1) - zw(k,j,i-1) ) * &
               ( zw(k,j,i) - zw(k-1,j,i) + zw(k,j,i-1) - zw(k-1,j,i-1) ) * & !because zw indexed as croco from 0 to nz
               ( dy(j,i) + dy(j,i-1) ) * u(k,j,i) &
               - qrt * ( &
               + zxdy(k,j,i  )*       dzw(k  ,j,i  )*w(k-1,j,i) & !note the index trick for w
               + zxdy(k,j,i  )* two * dzw(k+1,j,i  )*w(k+1-1,j,i) & !note the index trick for w
               + zxdy(k,j,i-1)*       dzw(k  ,j,i-1)*w(k-1,j,i-1) & !note the index trick for w
               + zxdy(k,j,i-1)* two * dzw(k+1,j,i-1)*w(k+1-1,j,i-1) & !note the index trick for w
               )  ! umask

       enddo
    enddo

!!$    allocate(uf_integr_tmp(1:nz,0:ny+1,0:nx+1))
!!$    uf_integr_tmp(1,1:ny,1:nx+1)=uf_integr1(1:ny,1:nx+1)
!!$    call fill_halo(1,uf_integr_tmp,lbc_null='u')
!!$    uf_integr1(1:ny,1:nx+1)=uf_integr_tmp(1,1:ny,1:nx+1)
!!$    deallocate(uf_integr_tmp)
!!$
!!$    allocate(uf_integr_tmp(1:nz,0:ny+1,0:nx+1))
!!$    uf_integr_tmp(1,1:ny,1:nx+1)=uf_integr2(1:ny,1:nx+1)
!!$    call fill_halo(1,uf_integr_tmp,lbc_null='u')
!!$    uf_integr2(1:ny,1:nx+1)=uf_integr_tmp(1,1:ny,1:nx+1)
!!$    deallocate(uf_integr_tmp)

    allocate(uf_integr_tmp(1:nz,0:ny+1,0:nx+1))
    uf_integr_tmp(:,:,:) = zero
    uf_integr_tmp(1,1:ny,1:nx+1)=uf_integr(1:ny,1:nx+1)
    call fill_halo(1,uf_integr_tmp,lbc_null='u')
    uf_integr(0:ny+1,0:nx+1)=uf_integr_tmp(1,0:ny+1,0:nx+1)
    deallocate(uf_integr_tmp)

    if (check_output) then
!!$       call write_netcdf(uf_integr1,vname='uf1in',netcdf_file_name='coin.nc',rank=myrank,iter=iter_coupling_in)
!!$       call write_netcdf(uf_integr2,vname='uf2in',netcdf_file_name='coin.nc',rank=myrank,iter=iter_coupling_in)
       call write_netcdf(uf_integr,vname='ufin',netcdf_file_name='coin.nc',rank=myrank,iter=iter_coupling_in)
    endif

    do i = 1,nx+1  
       do j = 1,ny
!       do j = 0,ny+1
          uf_integr(j,i) = uf_integr(j,i) - uf_bar(j,i)
       enddo
    enddo

    do i = 1,nx
!    do i = 0,nx+1
       do j = 1,ny+1

          k = 1 ! lower level
          sv_integr(j,i) = &
               qrt                                                       * &
!               ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j-1,i) - zw(k,j-1,i) ) * &
               ( zw(k,j,i) - zw(k-1,j,i) + zw(k,j-1,i) - zw(k-1,j-1,i) ) * & !because zw indexed as croco from 0 to nz
               ( dx(j,i) + dx(j-1,i) )
!!$          vf_integr1(j,i) = &
!!$               qrt                                                       * &
!!$!               ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j-1,i) - zw(k,j-1,i) ) * &
!!$               ( zw(k,j,i) - zw(k-1,j,i) + zw(k,j-1,i) - zw(k-1,j-1,i) ) * & !because zw indexed as croco from 0 to nz
!!$               ( dx(j,i) + dx(j-1,i) ) * v(k,j,i)
!!$          vf_integr2(j,i) = &
!!$               - qrt * ( &
!!$               + zydx(k,j  ,i)* two * dzw(k  ,j  ,i) * w(k-1,j,i) & !note the index trick for w 
!!$               + zydx(k,j  ,i)*       dzw(k+1,j  ,i) * w(k+1-1,j,i) & !note the index trick for w 
!!$               + zydx(k,j-1,i)* two * dzw(k  ,j-1,i) * w(k-1,j-1,i) & !note the index trick for w 
!!$               + zydx(k,j-1,i)*       dzw(k+1,j-1,i) * w(k+1-1,j-1,i) & !note the index trick for w 
!!$               )
!!$          vf_integr(j,i) = vf_integr1(j,i) + vf_integr2(j,i)
          vf_integr(j,i) = &
               qrt                                                       * &
!               ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j-1,i) - zw(k,j-1,i) ) * &
               ( zw(k,j,i) - zw(k-1,j,i) + zw(k,j-1,i) - zw(k-1,j-1,i) ) * & !because zw indexed as croco from 0 to nz
               ( dx(j,i) + dx(j-1,i) ) * v(k,j,i) &
               - qrt * ( &
               + zydx(k,j  ,i)* two * dzw(k  ,j  ,i) * w(k-1,j,i) & !note the index trick for w 
               + zydx(k,j  ,i)*       dzw(k+1,j  ,i) * w(k+1-1,j,i) & !note the index trick for w 
               + zydx(k,j-1,i)* two * dzw(k  ,j-1,i) * w(k-1,j-1,i) & !note the index trick for w 
               + zydx(k,j-1,i)*       dzw(k+1,j-1,i) * w(k+1-1,j-1,i) & !note the index trick for w 
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

          do k = 2,nz-1 ! interior levels
             sv_integr(j,i) = sv_integr(j,i) + &
                  qrt                                                       * &
!                  ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j-1,i) - zw(k,j-1,i) ) * &
                  ( zw(k,j,i) - zw(k-1,j,i) + zw(k,j-1,i) - zw(k-1,j-1,i) ) * & !because zw indexed as croco from 0 to nz
                  ( dx(j,i) + dx(j-1,i) )
!!$             vf_integr1(j,i) = vf_integr1(j,i) + &
!!$                  qrt                                                       * &
!!$!                  ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j-1,i) - zw(k,j-1,i) ) * &
!!$                  ( zw(k,j,i) - zw(k-1,j,i) + zw(k,j-1,i) - zw(k-1,j-1,i) ) * & !because zw indexed as croco from 0 to nz
!!$                  ( dx(j,i) + dx(j-1,i) ) * v(k,j,i)
!!$             vf_integr2(j,i) = vf_integr2(j,i) + &
!!$                  - qrt * ( &
!!$                  + zydx(k,j  ,i) * dzw(k  ,j  ,i) * w(k-1,j,i) & !note the index trick for w
!!$                  + zydx(k,j  ,i) * dzw(k+1,j  ,i) * w(k+1-1,j,i) & !note the index trick for w
!!$                  + zydx(k,j-1,i) * dzw(k  ,j-1,i) * w(k-1,j-1,i) & !note the index trick for w
!!$                  + zydx(k,j-1,i) * dzw(k+1,j-1,i) * w(k+1-1,j-1,i) & !note the index trick for w
!!$                  )
!!$             vf_integr(j,i) = vf_integr1(j,i) + vf_integr2(j,i)
             vf_integr(j,i) = vf_integr(j,i) + &
                  qrt                                                       * &
!                  ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j-1,i) - zw(k,j-1,i) ) * &
                  ( zw(k,j,i) - zw(k-1,j,i) + zw(k,j-1,i) - zw(k-1,j-1,i) ) * & !because zw indexed as croco from 0 to nz
                  ( dx(j,i) + dx(j-1,i) ) * v(k,j,i) &
                  - qrt * ( &
                  + zydx(k,j  ,i) * dzw(k  ,j  ,i) * w(k-1,j,i) & !note the index trick for w
                  + zydx(k,j  ,i) * dzw(k+1,j  ,i) * w(k+1-1,j,i) & !note the index trick for w
                  + zydx(k,j-1,i) * dzw(k  ,j-1,i) * w(k-1,j-1,i) & !note the index trick for w
                  + zydx(k,j-1,i) * dzw(k+1,j-1,i) * w(k+1-1,j-1,i) & !note the index trick for w
                  )  !* vmask(j,i)
          enddo

          k = nz ! upper level
          sv_integr(j,i) = sv_integr(j,i) + &
               qrt                                                       * &
!               ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j-1,i) - zw(k,j-1,i) ) * &
               ( zw(k,j,i) - zw(k-1,j,i) + zw(k,j-1,i) - zw(k-1,j-1,i) ) * & !because zw indexed as croco from 0 to nz
               ( dx(j,i) + dx(j-1,i) )
!!$          vf_integr1(j,i) = vf_integr1(j,i) + &
!!$               qrt                                                       * &
!!$!               ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j-1,i) - zw(k,j-1,i) ) * &
!!$               ( zw(k,j,i) - zw(k-1,j,i) + zw(k,j-1,i) - zw(k-1,j-1,i) ) * & !because zw indexed as croco from 0 to nz
!!$               ( dx(j,i) + dx(j-1,i) ) * v(k,j,i)
!!$          vf_integr2(j,i) = vf_integr2(j,i) + &
!!$               - qrt * ( &
!!$               + zydx(k,j  ,i)*       dzw(k  ,j  ,i) * w(k-1,j,i) & !note the index trick for w
!!$               + zydx(k,j  ,i)* two * dzw(k+1,j  ,i) * w(k+1-1,j,i) & !note the index trick for w
!!$               + zydx(k,j-1,i)*       dzw(k  ,j-1,i) * w(k-1,j-1,i) & !note the index trick for w
!!$               + zydx(k,j-1,i)* two * dzw(k+1,j-1,i) * w(k+1-1,j-1,i) & !note the index trick for w
!!$               )
!!$          vf_integr(j,i) = vf_integr1(j,i) + vf_integr2(j,i)
          vf_integr(j,i) = vf_integr(j,i) + &
               qrt                                                       * &
!               ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j-1,i) - zw(k,j-1,i) ) * &
               ( zw(k,j,i) - zw(k-1,j,i) + zw(k,j-1,i) - zw(k-1,j-1,i) ) * & !because zw indexed as croco from 0 to nz
               ( dx(j,i) + dx(j-1,i) ) * v(k,j,i) &
               - qrt * ( &
               + zydx(k,j  ,i)*       dzw(k  ,j  ,i) * w(k-1,j,i) & !note the index trick for w
               + zydx(k,j  ,i)* two * dzw(k+1,j  ,i) * w(k+1-1,j,i) & !note the index trick for w
               + zydx(k,j-1,i)*       dzw(k  ,j-1,i) * w(k-1,j-1,i) & !note the index trick for w
               + zydx(k,j-1,i)* two * dzw(k+1,j-1,i) * w(k+1-1,j-1,i) & !note the index trick for w
               ) !* vmask(j,i)

       enddo
    enddo

!!$    allocate(vf_integr_tmp(1:nz,0:ny+1,0:nx+1))
!!$    vf_integr_tmp(1,1:ny+1,1:nx)=vf_integr1(1:ny+1,1:nx)
!!$    call fill_halo(1,vf_integr_tmp,lbc_null='v')
!!$    vf_integr1(1:ny+1,1:nx)=vf_integr_tmp(1,1:ny+1,1:nx)
!!$    deallocate(vf_integr_tmp)
!!$
!!$    allocate(vf_integr_tmp(1:nz,0:ny+1,0:nx+1))
!!$    vf_integr_tmp(1,1:ny+1,1:nx)=vf_integr2(1:ny+1,1:nx)
!!$    call fill_halo(1,vf_integr_tmp,lbc_null='v')
!!$    vf_integr2(1:ny+1,1:nx)=vf_integr_tmp(1,1:ny+1,1:nx)
!!$    deallocate(vf_integr_tmp)

    allocate(vf_integr_tmp(1:nz,0:ny+1,0:nx+1))
    vf_integr_tmp(:,:,:) = zero
    vf_integr_tmp(1,1:ny+1,1:nx)=vf_integr(1:ny+1,1:nx)
    call fill_halo(1,vf_integr_tmp,lbc_null='v')
    vf_integr(0:ny+1,0:nx+1)=vf_integr_tmp(1,0:ny+1,0:nx+1)
    deallocate(vf_integr_tmp)

    if (check_output) then
!!$       call write_netcdf(vf_integr1,vname='vf1in',netcdf_file_name='coin.nc',rank=myrank,iter=iter_coupling_in)
!!$       call write_netcdf(vf_integr2,vname='vf2in',netcdf_file_name='coin.nc',rank=myrank,iter=iter_coupling_in)
       call write_netcdf(vf_integr,vname='vfin',netcdf_file_name='coin.nc',rank=myrank,iter=iter_coupling_in)
    endif

    do i = 1,nx
!    do i = 0,nx+1
       do j = 1,ny+1
          vf_integr(j,i) = vf_integr(j,i) - vf_bar(j,i)
       enddo
    enddo

    !-------------------------------
    !2 - correct u,v,w at each depth

    allocate(wc(0:nz,0:ny+1,0:nx+1)) ! indexed as in croco from 0 to nz
    wc(:,:,:) = zero

    do i = 1,nx
       do j = 1,ny
          
          k = 0 ! bottom

          wc(k,j,i) = ( &

               + hlf *( &
                 (zx(k  +1,j,i)) * &!note the index trick for zx
                 ( uf_integr(j,i  )/su_integr(j,i  ) & 
                  +uf_integr(j,i+1)/su_integr(j,i+1) ) ) &
            
               + hlf *( &
                 (zy(k  +1,j,i)) * &!note the index trick for zy
                 ( vf_integr(j  ,i)/sv_integr(j  ,i) & 
                  +vf_integr(j+1,i)/sv_integr(j+1,i) ) ) &

                      )

          do k = 1,nz-1 ! interior levels

             wc(k,j,i) = ( &

               - (zw(k,j,i)-zw(0,j,i))/(zw(nz,j,i)-zw(0,j,i)) & ! zw indexed as croco from 0 to nz
                *( uf_integr(j,i+1)-uf_integr(j,i)+vf_integr(j+1,i)-vf_integr(j,i) ) / (dx(j,i)*dy(j,i)) &

               + qrt *( &
                 (zx(k-1+1,j,i)+zx(k  +1,j,i)) * &!note the index trick for zx
                 ( uf_integr(j,i  )/su_integr(j,i  ) & 
                  +uf_integr(j,i+1)/su_integr(j,i+1) ) ) &
            
               + qrt *( &
                 (zy(k-1+1,j,i)+zy(k  +1,j,i)) * &!note the index trick for zy
                 ( vf_integr(j  ,i)/sv_integr(j  ,i) & 
                  +vf_integr(j+1,i)/sv_integr(j+1,i) ) ) &
            
                          )
          enddo

          k = nz ! surface

          wc(k,j,i) = ( &

               - ( uf_integr(j,i+1)-uf_integr(j,i)+vf_integr(j+1,i)-vf_integr(j,i) ) / (dx(j,i)*dy(j,i)) &

               + hlf *( &
                 (zx(k-1+1,j,i)) * &!note the index trick for zx
                 ( uf_integr(j,i  )/su_integr(j,i  ) & 
                  +uf_integr(j,i+1)/su_integr(j,i+1) ) ) &

               + hlf *( &
                 (zy(k-1+1,j,i)) * &!note the index trick for zy
                 ( vf_integr(j  ,i)/sv_integr(j  ,i) & 
                  +vf_integr(j+1,i)/sv_integr(j+1,i) ) ) &

                                 )

       enddo
    enddo

    call fill_halo(1,wc)
    if (check_output) then
       call write_netcdf(wc,vname='wcin',netcdf_file_name='coin.nc',rank=myrank,iter=iter_coupling_in)
    endif

    do i = 1,nx+1  
       do j = 1,ny 
!       do j = 0,ny+1

          k = 1 ! lower level

             u(k,j,i) = u(k,j,i) - ( &
                 Arx(k,j,i)*uf_integr(j,i)/su_integr(j,i) &
               + qrt *( &
                 zxdy(k,j,i  )* two *dzw(k  ,j,i  )*wc(k-1,j,i  ) & !note the index trick for wc
               + zxdy(k,j,i  )      *dzw(k+1,j,i  )*wc(k+1-1,j,i) &!note the index trick for wc
               + zxdy(k,j,i-1)* two *dzw(k  ,j,i-1)*wc(k-1,j,i-1) & !note the index trick for wc
               + zxdy(k,j,i-1)      *dzw(k+1,j,i-1)*wc(k+1-1,j,i-1) ) & !note the index trick for wc        
                                   ) / Arx(k,j,i)
 
          do k = 2,nz-1 ! interior levels

             u(k,j,i) = u(k,j,i) - ( &
                 Arx(k,j,i)*uf_integr(j,i)/su_integr(j,i) &
               + qrt *( &
                 zxdy(k,j,i  ) *dzw(k  ,j,i  )*wc(k-1,j,i  ) & !note the index trick for wc
               + zxdy(k,j,i  ) *dzw(k+1,j,i  )*wc(k+1-1,j,i) &!note the index trick for wc
               + zxdy(k,j,i-1) *dzw(k  ,j,i-1)*wc(k-1,j,i-1) & !note the index trick for wc
               + zxdy(k,j,i-1) *dzw(k+1,j,i-1)*wc(k+1-1,j,i-1) ) & !note the index trick for wc
                                   ) / Arx(k,j,i)
          enddo

          k = nz !upper level

             u(k,j,i) = u(k,j,i) - ( &
                 Arx(k,j,i)*uf_integr(j,i)/su_integr(j,i) &
               + qrt *( &
                 zxdy(k,j,i  )      *dzw(k  ,j,i  )*wc(k-1,j,i  ) & !note the index trick for wc
               + zxdy(k,j,i  )* two *dzw(k+1,j,i  )*wc(k+1-1,j,i) &!note the index trick for wc
               + zxdy(k,j,i-1)      *dzw(k  ,j,i-1)*wc(k-1,j,i-1) & !note the index trick for wc
               + zxdy(k,j,i-1)* two *dzw(k+1,j,i-1)*wc(k+1-1,j,i-1) ) & !note the index trick for wc        
                                   ) / Arx(k,j,i)

       enddo
    enddo

    allocate(uf_integr_tmp(1:nz,0:ny+1,0:nx+1))
    uf_integr_tmp(:,:,:) = zero
    uf_integr_tmp(:,1:ny,1:nx+1)=u(:,1:ny,1:nx+1)
    call fill_halo(1,uf_integr_tmp,lbc_null='u')
    u(:,0:ny+1,1:nx+1)=uf_integr_tmp(:,0:ny+1,1:nx+1)
    deallocate(uf_integr_tmp)

    do i = 1,nx
!    do i = 0,nx+1
       do j = 1,ny+1

          k = 1 ! lower level

             v(k,j,i) = v(k,j,i) - ( &
                 Ary(k,j,i)*vf_integr(j,i)/sv_integr(j,i) &
               + qrt *( &
                 zydx(k,j  ,i)* two * dzw(k  ,j  ,i)*w(k-1,j,i) & !note the index trick for w 
               + zydx(k,j  ,i)*       dzw(k+1,j  ,i)*w(k+1-1,j,i) & !note the index trick for w 
               + zydx(k,j-1,i)* two * dzw(k  ,j-1,i)*w(k-1,j-1,i) & !note the index trick for w 
               + zydx(k,j-1,i)*       dzw(k+1,j-1,i)*w(k+1-1,j-1,i) ) & !note the index trick for w 
                                   ) / Ary(k,j,i)

          do k = 2,nz-1 ! interior levels

             v(k,j,i) = v(k,j,i) - ( &
                 Ary(k,j,i)*vf_integr(j,i)/sv_integr(j,i) &
               + qrt *( &
                 zydx(k,j  ,i)* dzw(k  ,j  ,i)*w(k-1,j,i) & !note the index trick for w 
               + zydx(k,j  ,i)* dzw(k+1,j  ,i)*w(k+1-1,j,i) & !note the index trick for w 
               + zydx(k,j-1,i)* dzw(k  ,j-1,i)*w(k-1,j-1,i) & !note the index trick for w 
               + zydx(k,j-1,i)* dzw(k+1,j-1,i)*w(k+1-1,j-1,i) ) & !note the index trick for w 
                                   ) / Ary(k,j,i)
          
          enddo

          k = nz !upper level

             v(k,j,i) = v(k,j,i) - ( &
                 Ary(k,j,i)*vf_integr(j,i)/sv_integr(j,i) &
               + qrt *( &
                 zydx(k,j  ,i)*       dzw(k  ,j  ,i)*w(k-1,j,i) & !note the index trick for w 
               + zydx(k,j  ,i)* two * dzw(k+1,j  ,i)*w(k+1-1,j,i) & !note the index trick for w 
               + zydx(k,j-1,i)*       dzw(k  ,j-1,i)*w(k-1,j-1,i) & !note the index trick for w 
               + zydx(k,j-1,i)* two * dzw(k+1,j-1,i)*w(k+1-1,j-1,i) ) & !note the index trick for w 
                                   ) / Ary(k,j,i)

       enddo
    enddo

    allocate(vf_integr_tmp(1:nz,0:ny+1,0:nx+1))
    vf_integr_tmp(:,:,:) = zero
    vf_integr_tmp(:,1:ny+1,1:nx)=v(:,1:ny+1,1:nx)
    call fill_halo(1,vf_integr_tmp,lbc_null='v')
    v(:,1:ny+1,0:nx+1)=vf_integr_tmp(:,1:ny+1,0:nx+1)
    deallocate(vf_integr_tmp)

    do i = 1,nx
       do j = 1,ny
          
          do k = 0,nz ! all levels, indexed as in croco from 0 to nz

             w(k,j,i) = w(k,j,i) - wc(k,j,i) 

          enddo

       enddo
    enddo

    call fill_halo(1,w)

    deallocate(su_integr)
    deallocate(sv_integr)
!!$    deallocate(uf_integr1)
!!$    deallocate(vf_integr1)
!!$    deallocate(uf_integr2)
!!$    deallocate(vf_integr2)
    deallocate(uf_integr)
    deallocate(vf_integr)
    deallocate(wc)

    if ((present(uf)).and.(present(vf))) then
    !----------------------------------------------
    !3 - compute horizontal transport at each depth

    do i = 1,nx+1  
       do j = 1,ny 
!       do j = 0,ny+1

          k = 1 ! lower level
          uf(k,j,i) =  &
               qrt                                                      * & 
!               ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j,i-1) - zw(k,j,i-1) ) * &
               ( zw(k,j,i) - zw(k-1,j,i) + zw(k,j,i-1) - zw(k-1,j,i-1) ) * & !because zw indexed as croco from 0 to nz
               ( dy(j,i) + dy(j,i-1) ) * u(k,j,i) &
               - qrt * ( &
               + zxdy(k,j,i  )* two * dzw(k  ,j,i  )*w(k-1,j,i) & !note the index trick for w 
               + zxdy(k,j,i  )*       dzw(k+1,j,i  )*w(k+1-1,j,i) & !note the index trick for w 
               + zxdy(k,j,i-1)* two * dzw(k  ,j,i-1)*w(k-1,j,i-1) & !note the index trick for w 
               + zxdy(k,j,i-1)*       dzw(k+1,j,i-1)*w(k+1-1,j,i-1) & !note the index trick for w 
               )  ! umask

!!$               - qrt * ( &
!!$               + zxdy(k,j,i  )*dzw(k+1,j,i  )*w(k+1-1,j,i)   & !note the index trick for wc
!!$               + zxdy(k,j,i-1)*dzw(k+1,j,i-1)*w(k+1-1,j,i-1) ) & !note the index trick for wc
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

          do k = 2,nz-1 ! interior levels            
             uf(k,j,i) = &
                  qrt                                               * & 
!                  ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j,i-1) - zw(k,j,i-1) ) * &
                  ( zw(k,j,i) - zw(k-1,j,i) + zw(k,j,i-1) - zw(k-1,j,i-1) ) * & !because zw indexed as croco from 0 to nz
                  ( dy(j,i) + dy(j,i-1) ) * u(k,j,i) &
                  - qrt * ( &
                  + zxdy(k,j,i  ) * dzw(k  ,j,i  )*w(k-1,j,i) & !note the index trick for wc
                  + zxdy(k,j,i  ) * dzw(k+1,j,i  )*w(k+1-1,j,i) & !note the index trick for wc
                  + zxdy(k,j,i-1) * dzw(k  ,j,i-1)*w(k-1,j,i-1) & !note the index trick for wc
                  + zxdy(k,j,i-1) * dzw(k+1,j,i-1)*w(k+1-1,j,i-1) & !note the index trick for wc
                  )  ! umask
          enddo

          k = nz ! upper level
          uf(k,j,i) = &
               qrt                                                       * & 
!               ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j,i-1) - zw(k,j,i-1) ) * &
               ( zw(k,j,i) - zw(k-1,j,i) + zw(k,j,i-1) - zw(k-1,j,i-1) ) * & !because zw indexed as croco from 0 to nz
               ( dy(j,i) + dy(j,i-1) ) * u(k,j,i) &
               - qrt * ( &
               + zxdy(k,j,i  )*       dzw(k  ,j,i  )*w(k-1,j,i) & !note the index trick for wc
               + zxdy(k,j,i  )* two * dzw(k+1,j,i  )*w(k+1-1,j,i) & !note the index trick for wc
               + zxdy(k,j,i-1)*       dzw(k  ,j,i-1)*w(k-1,j,i-1) & !note the index trick for wc
               + zxdy(k,j,i-1)* two * dzw(k+1,j,i-1)*w(k+1-1,j,i-1) & !note the index trick for wc
               )  ! umask

       enddo
    enddo

    allocate(uf_tmp(1:nz,0:ny+1,0:nx+1))
    uf_tmp(:,1:ny,1:nx+1)=uf(:,1:ny,1:nx+1)
    call fill_halo(1,uf_tmp,lbc_null='u')
    uf(:,0:ny+1,1:nx+1)=uf_tmp(:,0:ny+1,1:nx+1)
    deallocate(uf_tmp)

    do i = 1,nx
!    do i = 0,nx+1
       do j = 1,ny+1

          k = 1 ! lower level
          vf(k,j,i) = &
               qrt                                                       * &
!               ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j-1,i) - zw(k,j-1,i) ) * &
               ( zw(k,j,i) - zw(k-1,j,i) + zw(k,j-1,i) - zw(k-1,j-1,i) ) * & !because zw indexed as croco from 0 to nz
               ( dx(j,i) + dx(j-1,i) ) * v(k,j,i) &
               - qrt * ( &
               + zydx(k,j  ,i)* two * dzw(k  ,j  ,i)*w(k-1,j,i) & !note the index trick for w 
               + zydx(k,j  ,i)*       dzw(k+1,j  ,i)*w(k+1-1,j,i) & !note the index trick for w 
               + zydx(k,j-1,i)* two * dzw(k  ,j-1,i)*w(k-1,j-1,i) & !note the index trick for w 
               + zydx(k,j-1,i)*       dzw(k+1,j-1,i)*w(k+1-1,j-1,i) & !note the index trick for w 
               ) !* vmask(j,i)

!!$               - qrt * ( &
!!$               + zydx(k,j  ,i) * dzw(k+1,j  ,i)*w(k+1-1,j,i) & !note the index trick for wc
!!$               + zydx(k,j-1,i) * dzw(k+1,j-1,i)*w(k+1-1,j-1,i) ) & !note the index trick for wc        
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

          do k = 2,nz-1 ! interior levels
             vf(k,j,i) = &
                  qrt                                                       * &
!                  ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j-1,i) - zw(k,j-1,i) ) * &
                  ( zw(k,j,i) - zw(k-1,j,i) + zw(k,j-1,i) - zw(k-1,j-1,i) ) * & !because zw indexed as croco from 0 to nz
                  ( dx(j,i) + dx(j-1,i) ) * v(k,j,i) &
                  - qrt * ( &
                  + zydx(k,j  ,i) * dzw(k  ,j  ,i)*w(k-1,j,i) & !note the index trick for wc
                  + zydx(k,j  ,i) * dzw(k+1,j  ,i)*w(k+1-1,j,i) & !note the index trick for wc
                  + zydx(k,j-1,i) * dzw(k  ,j-1,i)*w(k-1,j-1,i) & !note the index trick for wc
                  + zydx(k,j-1,i) * dzw(k+1,j-1,i)*w(k+1-1,j-1,i) & !note the index trick for wc
                  )  !* vmask(j,i)
          enddo

          k = nz ! upper level
          vf(k,j,i) = &
               qrt                                                       * &
!               ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j-1,i) - zw(k,j-1,i) ) * &
               ( zw(k,j,i) - zw(k-1,j,i) + zw(k,j-1,i) - zw(k-1,j-1,i) ) * & !because zw indexed as croco from 0 to nz
               ( dx(j,i) + dx(j-1,i) ) * v(k,j,i) &
               - qrt * ( &
               + zydx(k,j  ,i)*       dzw(k  ,j  ,i)*w(k-1,j,i) & !note the index trick for wc
               + zydx(k,j  ,i)* two * dzw(k+1,j  ,i)*w(k+1-1,j,i) & !note the index trick for wc
               + zydx(k,j-1,i)*       dzw(k  ,j-1,i)*w(k-1,j-1,i) & !note the index trick for wc
               + zydx(k,j-1,i)* two * dzw(k+1,j-1,i)*w(k+1-1,j-1,i) & !note the index trick for wc
               ) !* vmask(j,i)

       enddo
    enddo

    allocate(vf_tmp(1:nz,0:ny+1,0:nx+1))
    vf_tmp(:,1:ny+1,1:nx)=vf(:,1:ny+1,1:nx)
    call fill_halo(1,vf_tmp,lbc_null='v')
    vf(:,1:ny+1,0:nx+1)=vf_tmp(:,1:ny+1,0:nx+1)
    deallocate(vf_tmp)

    !! check
    allocate(uf_integr(0:ny+1,1:nx+1))
    allocate(vf_integr(1:ny+1,0:nx+1))
    uf_integr(:,:) = zero
    vf_integr(:,:) = zero
    do i = 1,nx+1  
       do j = 1,ny 
!       do j = 0,ny+1
          do k = 1,nz
             uf_integr(j,i) = uf_integr(j,i) + uf(k,j,i)
          enddo
       enddo
    enddo
    do i = 1,nx
!    do i = 0,nx+1
       do j = 1,ny+1
          do k = 1,nz
             vf_integr(j,i) = vf_integr(j,i) + vf(k,j,i)
          enddo
       enddo
    enddo
    if (check_output) then
       call write_netcdf(uf_integr,vname='ufout',netcdf_file_name='coin.nc',rank=myrank,iter=iter_coupling_in)
       call write_netcdf(vf_integr,vname='vfout',netcdf_file_name='coin.nc',rank=myrank,iter=iter_coupling_in)
       call write_netcdf(uf_integr-uf_bar,vname='diff_uf',netcdf_file_name='coin.nc',rank=myrank,iter=iter_coupling_in)
       call write_netcdf(vf_integr-vf_bar,vname='diff_vf',netcdf_file_name='coin.nc',rank=myrank,iter=iter_coupling_in)
    endif
    deallocate(uf_integr)
    deallocate(vf_integr)
    !! check

    endif

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

  end subroutine btbc_coupling

end module mg_btbc_coupling
#else
        module mg_btbc_coupling_empty
        end module
#endif

