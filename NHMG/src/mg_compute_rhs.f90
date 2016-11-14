module mg_compute_rhs

  use mg_mpi
  use mg_tictoc
  use mg_namelist
  use mg_grids
  use mg_mpi_exchange
  use mg_netcdf_out

  implicit none

contains
  !-------------------------------------------------------------------------     
  subroutine compute_rhs(zr,zw,u,v,w)

    real(kind=rp), dimension(:,:,:), pointer, intent(in)  :: zr,zw
    real(kind=rp), dimension(:,:,:), pointer, intent(in) :: u,v,w

    integer(kind=ip):: k, j, i
    integer(kind=ip):: nx, ny, nz

    real(kind=rp), dimension(:,:)  , pointer :: dx,dy
!    real(kind=rp), dimension(:,:,:), pointer :: zr,zw
!    real(kind=rp) :: Arx, Ary
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

    real(kind=rp), dimension(:,:,:), pointer :: uf,vf,wf,rhs

    real(kind=rp), parameter :: two  = 2._rp
    real(kind=rp), parameter :: one  = 1._rp
    real(kind=rp), parameter :: hlf  = 0.5_rp
    real(kind=rp), parameter :: qrt  = 0.25_rp
    real(kind=rp), parameter :: zero = 0._rp

    integer(kind=ip), save :: iter_rhs=-2
    iter_rhs = iter_rhs + 1

    nx = grid(1)%nx
    ny = grid(1)%ny
    nz = grid(1)%nz

    dx => grid(1)%dx
    dy => grid(1)%dy

!    zr => grid(1)%zr
!    zw => grid(1)%zw

    if (myrank==0) write(*,*)'- compute rhs:'

!    !- Cell heights (computed in define matices)
!    dzw => grid(1)%dzw
!
!    !- Slopes in x- and y-direction defined at rho-points (computed in define matices)
!    zxdy => grid(1)%zxdy
!    zydx => grid(1)%zydx
!
!    !- (computed in define matices)
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

    !--------------!
    !- DIVERGENCE -!
    !--------------!
    !- Divergence is stored in rhs array pointer
    !- It is calculated progressively using first uf, 
    !- then using vf and at the e wf
    rhs => grid(1)%b
    rhs(:,:,:) = zero

    !------!
    !- UF -!
    !------!
    !- if (myrank==0) write(*,*)'- compute rhs: uf'

    uf => grid(1)%dummy3Dnz
    uf(:,:,:) = zero

    k = 1 ! lower level

    do i = 1,nx+1  
       do j = 1,ny 

          uf(k,j,i) =  &
               qrt                                                      * & 
!               ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j,i-1) - zw(k,j,i-1) ) * &
               ( zw(k,j,i) - zw(k-1,j,i) + zw(k,j,i-1) - zw(k-1,j,i-1) ) * & !because zw indexed as croco from 0 to nz
               ( dy(j,i) + dy(j,i-1) ) * u(k,j,i) &

               - qrt * ( &
               + zxdy(k,j,i  )*dzw(k+1,j,i  )*w(k+1-1,j,i)   & !note the index trick for w 
               + zxdy(k,j,i-1)*dzw(k+1,j,i-1)*w(k+1-1,j,i-1) ) & !note the index trick for w 
               -( &
               + zxdy(k,j,i  )*zxdy(k,j,i  )/(cw(k,j,i  )+cw(k+1,j,i  )) &
               + zxdy(k,j,i-1)*zxdy(k,j,i-1)/(cw(k,j,i-1)+cw(k+1,j,i-1)) ) * &
               (hlf * ( dx(j,i) + dx(j,i-1) )) * u(k,j,i) &
               -( & 
               + zxdy(k,j,i) * zydx(k,j,i)/(cw(k,j,i  )+cw(k+1,j,i  ))   &
               * hlf * ( &
               hlf * ( dy(j  ,i) + dy(j-1,i) )*v(k,j,i ) + &
               hlf * ( dy(j+1,i) + dy(j  ,i) )*v(k,j+1,i ))   & 
               + zxdy(k,j,i-1) * zydx(k,j,i-1)/(cw(k,j,i-1)+cw(k+1,j,i-1))   &
               * hlf * ( &
               hlf * ( dy(j  ,i-1) + dy(j-1,i-1) )*v(k,j,i-1) + &
               hlf * ( dy(j+1,i-1) + dy(j  ,i-1) )*v(k,j+1,i-1)) )

       enddo
    enddo

    do i = 1,nx+1  
       do j = 1,ny 

          do k = 2,nz-1 !interior levels
             
             uf(k,j,i) =  qrt                                               * & 
!                  ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j,i-1) - zw(k,j,i-1) ) * &
                  ( zw(k,j,i) - zw(k-1,j,i) + zw(k,j,i-1) - zw(k-1,j,i-1) ) * & !because zw indexed as croco from 0 to nz
                  ( dy(j,i) + dy(j,i-1) ) * u(k,j,i) &
                  - qrt * ( &
                  + zxdy(k,j,i  ) * dzw(k  ,j,i  ) * w(k  -1,j,i  ) & !note the index trick for w 
                  + zxdy(k,j,i  ) * dzw(k+1,j,i  ) * w(k+1-1,j,i  ) & !note the index trick for w 
                  + zxdy(k,j,i-1) * dzw(k  ,j,i-1) * w(k  -1,j,i-1) & !note the index trick for w 
                  + zxdy(k,j,i-1) * dzw(k+1,j,i-1) * w(k+1-1,j,i-1) & !note the index trick for w 
                  )  ! umask

          enddo

       enddo
    enddo

    k = nz !upper level

    do i = 1,nx+1  
       do j = 1,ny 

          uf(k,j,i) = &
               qrt                                                       * & 
!               ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j,i-1) - zw(k,j,i-1) ) * &
               ( zw(k,j,i) - zw(k-1,j,i) + zw(k,j,i-1) - zw(k-1,j,i-1) ) * & !because zw indexed as croco from 0 to nz
               ( dy(j,i) + dy(j,i-1) ) * u(k,j,i) &
               - qrt * ( &
               + zxdy(k,j,i  )*       dzw(k  ,j,i  )*w(k  -1,j,i  ) & !note the index trick for w 
               + zxdy(k,j,i  )* two * dzw(k+1,j,i  )*w(k+1-1,j,i  ) & !note the index trick for w 
               + zxdy(k,j,i-1)*       dzw(k  ,j,i-1)*w(k  -1,j,i-1) & !note the index trick for w 
               + zxdy(k,j,i-1)* two * dzw(k+1,j,i-1)*w(k+1-1,j,i-1) & !note the index trick for w 
               )  ! umask

       enddo
    enddo

    call fill_halo(1,uf,lbc_null='u')

!    if (check_output) then
!       if ((iter_rhs .EQ. 397) .OR. (iter_rhs .EQ. 399) .OR. &
!           (iter_rhs .EQ. 1997) .OR. (iter_rhs .EQ. 1999) .OR. &
!           (iter_rhs .EQ. 3997) .OR. (iter_rhs .EQ. 3999) .OR. &
!           (iter_rhs .GE. 6997)) then
!       call write_netcdf(uf,vname='uf',netcdf_file_name='rhs.nc',rank=myrank,iter=iter_rhs)
!       endif
!    endif

    do i = 1,nx
       do j = 1,ny 
          do k = 1,nz

             rhs(k,j,i) = uf(k,j,i+1) - uf(k,j,i) 

          enddo
       enddo
    enddo

    uf => null()

    !------!
    !- VF -!
    !------!
    !- if (myrank==0) write(*,*)'- compute rhs: vf'

    vf => grid(1)%dummy3Dnz
    vf(:,:,:) = zero

    k = 1 !lower level

    do i = 1,nx
       do j = 1,ny+1

          vf(k,j,i) = &
               qrt                                                       * &
!               ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j-1,i) - zw(k,j-1,i) ) * &
               ( zw(k,j,i) - zw(k-1,j,i) + zw(k,j-1,i) - zw(k-1,j-1,i) ) * & !because zw indexed as croco from 0 to nz
               ( dx(j,i) + dx(j-1,i) ) * v(k,j,i) &

               - qrt * ( &
               + zydx(k,j  ,i) * dzw(k+1,j  ,i)*w(k+1-1,j,i) & !note the index trick for w 
               + zydx(k,j-1,i) * dzw(k+1,j-1,i)*w(k+1-1,j-1,i) ) & !note the index trick for w 
               -( &
               + zydx(k,j  ,i) * zydx(k,j  ,i)/(cw(k,j  ,i)+cw(k+1,j  ,i)) &
               + zydx(k,j-1,i) * zydx(k,j-1,i)/(cw(k,j-1,i)+cw(k+1,j-1,i)) ) * &
               hlf * ( dy(j,i) + dy(j-1,i) ) * v(k,j,i) &
               - ( &
               + zxdy(k,j  ,i) * zydx(k,j  ,i)/(cw(k,j  ,i)+cw(k+1,j  ,i))  &
               * hlf * ( &
               hlf * ( dx(j,i  ) + dx(j,i-1) ) * u(k,j,i) + &
               hlf * ( dx(j,i+1) + dx(j,i  ) ) * u(k,j,i+1)) &
               + zxdy(k,j-1,i) * zydx(k,j-1,i)/(cw(k,j-1,i)+cw(k+1,j-1,i))   &
               * hlf * ( &
               hlf * ( dx(j-1,i  ) + dx(j-1,i-1) ) * u(k,j-1,i) + &
               hlf * ( dx(j-1,i+1) + dx(j-1,i  ) ) * u(k,j-1,i+1)) &
               ) 

       enddo
    enddo

    do i = 1,nx
       do j = 1,ny+1

          do k = 2,nz-1 !interior levels

             vf(k,j,i) = &
                  qrt                                                       * &
!                  ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j-1,i) - zw(k,j-1,i) ) * &
                  ( zw(k,j,i) - zw(k-1,j,i) + zw(k,j-1,i) - zw(k-1,j-1,i) ) * & !because zw indexed as croco from 0 to nz
                  ( dx(j,i) + dx(j-1,i) ) * v(k,j,i) &
                  - qrt * ( &
                  + zydx(k,j  ,i) * dzw(k  ,j  ,i) * w(k  -1,j,i  ) & !note the index trick for w 
                  + zydx(k,j  ,i) * dzw(k+1,j  ,i) * w(k+1-1,j,i  ) & !note the index trick for w 
                  + zydx(k,j-1,i) * dzw(k  ,j-1,i) * w(k  -1,j-1,i) & !note the index trick for w 
                  + zydx(k,j-1,i) * dzw(k+1,j-1,i) * w(k+1-1,j-1,i) & !note the index trick for w 
                  )  !* vmask(j,i)

          enddo

       enddo
    enddo

    k = nz !upper level

    do i = 1,nx
       do j = 1,ny+1

          vf(k,j,i) =  &
               qrt                                                       * &
!               ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j-1,i) - zw(k,j-1,i) ) * &
               ( zw(k,j,i) - zw(k-1,j,i) + zw(k,j-1,i) - zw(k-1,j-1,i) ) * & !because zw indexed as croco from 0 to nz
               ( dx(j,i) + dx(j-1,i) ) * v(k,j,i) &
               - qrt * ( &
               + zydx(k,j  ,i)*       dzw(k  ,j  ,i) * w(k  -1,j  ,i) & !note the index trick for w 
               + zydx(k,j  ,i)* two * dzw(k+1,j  ,i) * w(k+1-1,j  ,i) & !note the index trick for w 
               + zydx(k,j-1,i)*       dzw(k  ,j-1,i) * w(k  -1,j-1,i) & !note the index trick for w 
               + zydx(k,j-1,i)* two * dzw(k+1,j-1,i) * w(k+1-1,j-1,i) & !note the index trick for w 
               ) !* vmask(j,i)

       enddo
    enddo

    call fill_halo(1,vf,lbc_null='v')

!    if (check_output) then
!    if ((iter_rhs .EQ. 397) .OR. (iter_rhs .EQ. 399) .OR. &
!           (iter_rhs .EQ. 1997) .OR. (iter_rhs .EQ. 1999) .OR. &
!           (iter_rhs .EQ. 3997) .OR. (iter_rhs .EQ. 3999) .OR. &
!           (iter_rhs .GE. 6997)) then
!       call write_netcdf(vf,vname='vf',netcdf_file_name='rhs.nc',rank=myrank,iter=iter_rhs)
!       endif
!    endif

    do i = 1,nx
       do j = 1,ny 
          do k = 1,nz

             rhs(k,j,i) =  rhs(k,j,i) + vf(k,j+1,i) - vf(k,j,i)

          enddo
       enddo
    enddo

    vf => null()

    !------!
    !- WF -!
    !------!
    !- if (myrank==0) write(*,*)'- compute rhs: wf'

    wf => grid(1)%dummy3Dnzp
    wf(:,:,:) = zero

    k = 1 !bottom

    do i = 1,nx
       do j = 1,ny

          wf(k,j,i) = zero

       enddo
    enddo

    do i = 1,nx
       do j = 1,ny

          do k = 2,nz !interior levels

             wf(k,j,i) = cw(k,j,i) * dzw(k,j,i) * w(k-1,j,i) & !note the index trick for w
                  - qrt * hlf * ( &
                  + zxdy(k  ,j,i) * ( dx(j,i  ) + dx(j,i-1) ) * u(k  ,j,i  ) &
                  + zxdy(k  ,j,i) * ( dx(j,i+1) + dx(j,i  ) ) * u(k  ,j,i+1) &
                  + zxdy(k-1,j,i) * ( dx(j,i  ) + dx(j,i-1) ) * u(k-1,j,i  ) &
                  + zxdy(k-1,j,i) * ( dx(j,i+1) + dx(j,i  ) ) * u(k-1,j,i+1) ) &
                  - qrt * hlf * ( &
                  + zydx(k  ,j,i) * ( dy(j  ,i) + dy(j-1,i) ) * v(k  ,j  ,i) &
                  + zydx(k  ,j,i) * ( dy(j+1,i) + dy(j  ,i) ) * v(k  ,j+1,i) &
                  + zydx(k-1,j,i) * ( dy(j  ,i) + dy(j-1,i) ) * v(k-1,j  ,i) &
                  + zydx(k-1,j,i) * ( dy(j+1,i) + dy(j  ,i) ) * v(k-1,j+1,i) )
          enddo

       enddo
    enddo

    k = nz+1 !surface

    do i = 1,nx
       do j = 1,ny

          wf(k,j,i) = cw(k,j,i) * dzw(k,j,i) * w(k-1,j,i) & !note the index trick for w
               - hlf * hlf * ( &
               + zxdy(k-1,j,i) * ( dx(j,i  ) + dx(j,i-1) ) * u(k-1,j,i  ) &
               + zxdy(k-1,j,i) * ( dx(j,i+1) + dx(j,i  ) ) * u(k-1,j,i+1) ) &
               - hlf * hlf * ( &
               + zydx(k-1,j,i) * ( dy(j  ,i) + dy(j-1,i) ) * v(k-1,j  ,i) &
               + zydx(k-1,j,i) * ( dy(j+1,i) + dy(j  ,i) ) * v(k-1,j+1,i) )

       enddo
    enddo

!    if (check_output) then
!       if ((iter_rhs .EQ. 397) .OR. (iter_rhs .EQ. 399) .OR. &
!           (iter_rhs .EQ. 1997) .OR. (iter_rhs .EQ. 1999) .OR. &
!           (iter_rhs .EQ. 3997) .OR. (iter_rhs .EQ. 3999) .OR. &
!           (iter_rhs .GE. 6997)) then
!       call write_netcdf(wf,vname='wf',netcdf_file_name='rhs.nc',rank=myrank,iter=iter_rhs)
!       endif
!    endif

    do i = 1,nx
       do j = 1,ny 
          do k = 1,nz

             rhs(k,j,i) = rhs(k,j,i) + wf(k+1,j,i) - wf(k,j,i)

          enddo
       enddo
    enddo

    wf => null()

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

  end subroutine compute_rhs

end module mg_compute_rhs

