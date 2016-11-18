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
  subroutine compute_rhs(u,v,w)

    real(kind=rp), dimension(:,:,:), pointer, intent(in) :: u,v,w

    integer(kind=ip):: k, j, i
    integer(kind=ip):: nx, ny, nz

    real(kind=rp), dimension(:,:)  , pointer :: dx,dy
    real(kind=rp), dimension(:,:,:), pointer :: zr,zw

    real(kind=rp), dimension(:,:,:), pointer :: dzw
    real(kind=rp), dimension(:,:,:), pointer :: Arx,Ary
    real(kind=rp), dimension(:,:,:), pointer :: zxdy,zydx
    real(kind=rp), dimension(:,:,:), pointer :: cw
 
    real(kind=rp), dimension(:,:,:), pointer :: uf,vf,wf
    real(kind=rp), dimension(:,:,:), pointer :: rhs

    real(kind=rp), parameter :: two  = 2._rp
    real(kind=rp), parameter :: one  = 1._rp
    real(kind=rp), parameter :: hlf  = 0.5_rp
    real(kind=rp), parameter :: qrt  = 0.25_rp
    real(kind=rp), parameter :: zero = 0._rp

    integer(kind=ip), save :: iter_rhs=-2
    iter_rhs = iter_rhs + 1

    if (myrank==0) write(*,*)'- compute rhs:',iter_rhs

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
               ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j,i-1) - zw(k,j,i-1) ) * &
               ( dy(j,i) + dy(j,i-1) ) * u(k,j,i) &

               - qrt * ( &
               + zxdy(k,j,i  )*dzw(k+1,j,i  )*w(k+1,j,i)   &
               + zxdy(k,j,i-1)*dzw(k+1,j,i-1)*w(k+1,j,i-1) ) &
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
                  ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j,i-1) - zw(k,j,i-1) ) * &
                  ( dy(j,i) + dy(j,i-1) ) * u(k,j,i) &
                  - qrt * ( &
                  + zxdy(k,j,i  ) * dzw(k  ,j,i  ) * w(k  ,j,i  ) & 
                  + zxdy(k,j,i  ) * dzw(k+1,j,i  ) * w(k+1,j,i  ) & 
                  + zxdy(k,j,i-1) * dzw(k  ,j,i-1) * w(k  ,j,i-1) & 
                  + zxdy(k,j,i-1) * dzw(k+1,j,i-1) * w(k+1,j,i-1) & 
                  )  ! umask

          enddo

       enddo
    enddo

    k = nz !upper level

    do i = 1,nx+1  
       do j = 1,ny 

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
               ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j-1,i) - zw(k,j-1,i) ) * &
               ( dx(j,i) + dx(j-1,i) ) * v(k,j,i) &

               - qrt * ( &
               + zydx(k,j  ,i) * dzw(k+1,j  ,i)*w(k+1,j,i) & 
               + zydx(k,j-1,i) * dzw(k+1,j-1,i)*w(k+1,j-1,i) ) & 
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
                  ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j-1,i) - zw(k,j-1,i) ) * &
                  ( dx(j,i) + dx(j-1,i) ) * v(k,j,i) &
                  - qrt * ( &
                  + zydx(k,j  ,i) * dzw(k  ,j  ,i) * w(k  ,j,i  ) & 
                  + zydx(k,j  ,i) * dzw(k+1,j  ,i) * w(k+1,j,i  ) &
                  + zydx(k,j-1,i) * dzw(k  ,j-1,i) * w(k  ,j-1,i) & 
                  + zydx(k,j-1,i) * dzw(k+1,j-1,i) * w(k+1,j-1,i) &
                  )  !* vmask(j,i)

          enddo

       enddo
    enddo

    k = nz !upper level

    do i = 1,nx
       do j = 1,ny+1

          vf(k,j,i) =  &
               qrt                                                       * &
               ( zw(k+1,j,i) - zw(k,j,i) + zw(k+1,j-1,i) - zw(k,j-1,i) ) * &
               ( dx(j,i) + dx(j-1,i) ) * v(k,j,i) &
               - qrt * ( &
               + zydx(k,j  ,i)*       dzw(k  ,j  ,i) * w(k  ,j  ,i) & 
               + zydx(k,j  ,i)* two * dzw(k+1,j  ,i) * w(k+1,j  ,i) &
               + zydx(k,j-1,i)*       dzw(k  ,j-1,i) * w(k  ,j-1,i) & 
               + zydx(k,j-1,i)* two * dzw(k+1,j-1,i) * w(k+1,j-1,i) & 
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

             wf(k,j,i) = cw(k,j,i) * dzw(k,j,i) * w(k,j,i) & 
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

          wf(k,j,i) = cw(k,j,i) * dzw(k,j,i) * w(k,j,i) & 
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

  end subroutine compute_rhs

end module mg_compute_rhs

