program mg_testseamount

  use mg_mpi
  use mg_tictoc
  use mg_setup_tests
  use mg_mpi_exchange_ijk
  use nhydro

  implicit none

  integer(kind=ip):: nxg    ! global x dimension
  integer(kind=ip):: nyg    ! global y dimension
  integer(kind=ip):: nzg    ! z dimension

  integer(kind=ip):: npxg   ! number of processes in x
  integer(kind=ip):: npyg   ! number of processes in y

  integer(kind=ip) :: nx, ny, nz  ! local dimensions

  real(kind=rp), dimension(:,:,:), allocatable, target :: u,v,w
  real(kind=rp), dimension(:,:,:), allocatable, target :: tmp_rnd
  real(kind=rp), dimension(:,:,:), allocatable, target :: tmp_rnd2
  real(kind=rp), dimension(:,:,:), pointer :: up,vp,wp
  real(kind=rp) :: Lx, Ly, Htot

  integer(kind=ip) :: np, ierr, rank

  integer(kind=ip) :: it, nit

  real(kind=rp), dimension(:,:), pointer :: dx, dy, zeta, h
  real(kind=rp), dimension(:,:), pointer :: rmask
  real(kind=rp) :: hc, theta_b, theta_s

  integer(kind=ip) :: pi, pj
  integer(kind=ip) :: ib, ie, jb, je, kb, ke

  call tic(1,'mg_bench_seamount')

  nit = 1

  !---------------------!
  !- Global/local dim  -!
  !---------------------!
  nxg   = 64
  nyg   = 64
  nzg   = 64

  npxg  = 2
  npyg  = 2

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world, rank, ierr)
  call mpi_comm_size(mpi_comm_world, np, ierr)

  if (np /= (npxg*npyg)) then
     write(*,*) "Error: in number of processes !", np, npxg, npyg
     stop -1
  endif

  nx = nxg / npxg
  ny = nyg / npyg
  nz = nzg

  !---------------------!
  !- Initialise nhydro -!
  !---------------------!
  if (rank == 0) write(*,*)'Initialise nhydro grids'

  call nhydro_init(nx,ny,nz,npxg,npyg)

  !---------------------!
  !- Setup seamount    -!
  !---------------------!
  if (rank == 0) write(*,*)'Initialise seamount bench'

  Lx   =  1.e4_rp
  Ly   =  1.e4_rp
  Htot =  4.e3_rp 

  hc      = 4.e3_rp
  theta_b =  0._rp
  theta_s =  0._rp

  if (myrank==0) then
     write(*,*)''
     write(*,*)'Lx, Ly, Htot:',Lx, Ly, Htot
     write(*,*)'hc, theta_b, theta_s:',hc, theta_b, theta_s
  endif

  allocate(   dx(0:ny+1,0:nx+1))
  allocate(   dy(0:ny+1,0:nx+1))
  allocate( zeta(0:ny+1,0:nx+1))
  allocate(    h(0:ny+1,0:nx+1))

  call setup_seamount(nx,ny,nz,npxg,npyg,Lx,Ly,Htot,dx,dy,zeta,h)

  allocate(rmask(0:ny+1,0:nx+1))
  rmask(:,:) = 1._rp
  if (bmask) then
     !- Mask the boundaries
     call fill_halo_2D_bmask(1,rmask)   !- Generic n procs
  endif

  call nhydro_matrices(dx,dy,zeta,h,rmask,hc,theta_b,theta_s)

  !-------------------------------------!
  !- U,V,W initialisation (model vars) -!
  !-------------------------------------!
  if (rank == 0) write(*,*)'Initialise u, v, w'

  allocate(u(1:nx+1,0:ny+1,1:nz))
  allocate(v(0:nx+1,1:ny+1,1:nz))
  allocate(w(0:nx+1,0:ny+1,0:nz))

  if (myrank==0) then
     write(*,*)''
     write(*,*)'U=0, V=0 and W=-1 except at bottom'
  endif
  u(:,:,:)    =  0._8
  v(:,:,:)    =  0._8

  w(:,:,0)    =  0._8
  w(:,:,1:nz) = -1._8

  up => u
  vp => v
  wp => w

!!$  if (myrank==0) then
!!$     write(*,*)''
!!$     write(*,*)'U, V and W are initalized with random numbers /= on each process'
!!$  endif
!!$  pj = myrank/npxg
!!$  pi = mod(myrank,npxg)
!!$
!!$  ib = 1 + pi * nx
!!$  jb = 1 + pj * ny
!!$
!!$  ie = ib + nx - 1
!!$  je = jb + ny - 1
!!$
!!$  !  allocate(tmp_rnd (1:nx,1:ny,1:nz))
!!$  !  allocate(tmp_rnd2(1:ny,1:ny,0:nz))
!!$
!!$  allocate(tmp_rnd (1:nxg,1:nyg,1:nzg))
!!$  allocate(tmp_rnd2(1:nxg,1:nyg,0:nzg))
!!$
     do it = 1, nit
!!$
!!$     kb = 1
!!$     ke = nz
!!$
!!$     call random_number(tmp_rnd)
!!$     tmp_rnd = 2._8 * tmp_rnd - 1._8
!!$     u(1:nx,1:ny,1:nz) = tmp_rnd(ib:ie,jb:je,kb:ke)
!!$     up => u
!!$     call fill_halo_ijk(nx,ny,up,'u') ! depend of mg_grids for MPI neighbours !
!!$
!!$     call random_number(tmp_rnd)
!!$     tmp_rnd = 2._8 * tmp_rnd - 1._8
!!$     v(1:nx,1:ny,1:nz) = tmp_rnd(ib:ie,jb:je,kb:ke)
!!$     vp => v
!!$     call fill_halo_ijk(nx,ny,vp,'v') ! depend of mg_grids for MPI neighbours !
!!$
!!$     kb = 0
!!$     ke = nz
!!$
!!$     call random_number(tmp_rnd2)
!!$     tmp_rnd2 = 2._8 * tmp_rnd2 - 1._8
!!$     w(1:nx,1:ny,0:nz) = tmp_rnd2(ib:ie,jb:je,kb:ke)
!!$     wp => w
!!$     call fill_halo_ijk(nx,ny,wp,'w') ! depend of mg_grids for MPI neighbours !

     if (netcdf_output) then
        call write_netcdf(u,vname='u',netcdf_file_name='u.nc',rank=myrank,iter=it)
        call write_netcdf(v,vname='v',netcdf_file_name='v.nc',rank=myrank,iter=it)
        call write_netcdf(w,vname='w',netcdf_file_name='w.nc',rank=myrank,iter=it)
     endif

     !----------------------!
     !- Call nhydro solver -!
     !----------------------!
     if (rank == 0) write(*,*)'Call nhydro solver'

     call nhydro_solve(nx,ny,nz,rmask,u,v,w)

     if (netcdf_output) then
        call write_netcdf(u,vname='uc',netcdf_file_name='uc.nc',rank=myrank,iter=it)
        call write_netcdf(v,vname='vc',netcdf_file_name='vc.nc',rank=myrank,iter=it)
        call write_netcdf(w,vname='wc',netcdf_file_name='wc.nc',rank=myrank,iter=it)
     endif

     !------------------------------------------------------------!
     !- Check if nh correction is correct                        -!
     !------------------------------------------------------------!
     if (rank == 0) write(*,*)'Check nondivergence'

     call nhydro_check_nondivergence(nx,ny,nz,rmask,u,v,w)

     if (netcdf_output) then
        call write_netcdf(grid(1)%b,vname='bc',netcdf_file_name='bc.nc',rank=myrank,iter=it)
     endif

  enddo

!!$  deallocate(tmp_rnd)
!!$  deallocate(tmp_rnd2)

  !---------------------!
  !- Deallocate memory -!
  !---------------------!
  if (rank == 0) write(*,*)'Clean memory before to finish the program.'
  call nhydro_clean()

  !----------------------!
  !- End Bench-seamount -!
  !----------------------!
  call mpi_finalize(ierr)

  call toc(1,'mg_bench_seamount')
  if(myrank == 0) call print_tictoc(myrank)

end program mg_testseamount

