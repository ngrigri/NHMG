module nhmg

  use mg_cst
  use mg_mpi
  use mg_grids
  use mg_namelist
  use mg_tictoc
  use mg_mpi_exchange
  use mg_horiz_grids
  use mg_vert_grids
  use mg_bbc
  use mg_fluxes
  use mg_coupling
  use mg_projection
  use mg_solvers
  use mg_netcdf_out

  implicit none

contains

  !--------------------------------------------------------------
  subroutine nhmg_write(nx,ny,nz,ua,va,wa)

    integer(kind=ip), intent(in) :: nx, ny, nz

    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: ua
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: va
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz+1), target, intent(in) :: wa

    real(kind=rp), dimension(:,:,:), pointer :: u,v,w

    integer(kind=ip), save :: iter_write=0
    iter_write = iter_write + 1

    u => ua
    v => va
    w => wa

!    if ((iter_write .EQ. 1)) then
       call write_netcdf(u,vname='u',netcdf_file_name='wrt.nc',rank=myrank,iter=iter_write)
       call write_netcdf(v,vname='v',netcdf_file_name='wrt.nc',rank=myrank,iter=iter_write)
       call write_netcdf(w,vname='w',netcdf_file_name='wrt.nc',rank=myrank,iter=iter_write)
!    endif

  end subroutine nhmg_write

  !--------------------------------------------------------------
  subroutine nhmg_init(nx,ny,nz,npxg,npyg)
      
    integer(kind=ip), intent(in) :: nx, ny, nz
    integer(kind=ip), intent(in) :: npxg, npyg

    call mg_mpi_init()

    if (myrank==0) write(*,*)' nhmg_init:'

    call read_nhnamelist(vbrank=myrank)

    call define_grids(npxg,npyg,nx,ny,nz)

    call define_neighbours()

    call print_grids()

  end subroutine nhmg_init

  !--------------------------------------------------------------
  subroutine nhmg_fluxes(nx,ny,nz,ua,va,wa,ufa,vfa)

    integer(kind=ip), intent(in) :: nx, ny, nz

    !real(kind=rp), dimension(1:nx+1,0:ny+1,1:nz  ), target, intent(in) :: ua
    !real(kind=rp), dimension(0:nx+1,1:ny+1,1:nz  ), target, intent(in) :: va
    !real(kind=rp), dimension(0:nx+1,0:ny+1,1:nz+1), target, intent(in) :: wa
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz  ), target, intent(in) :: ua !ND 11/01
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz  ), target, intent(in) :: va !ND 11/01
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz+1), target, intent(in) :: wa !ND 11/01

    real(kind=rp), dimension(1:nx+1,0:ny+1,1:nz  ),   target, intent(out):: ufa
    real(kind=rp), dimension(0:nx+1,1:ny+1,1:nz  ),   target, intent(out):: vfa

    real(kind=rp), dimension(:,:,:), pointer :: u,v,w
    real(kind=rp), dimension(:,:,:), pointer :: uf,vf

!!! dirty reshape arrays indexing ijk -> kji !!!
    real(kind=rp), dimension(:,:,:), allocatable, target :: ub,vb,wb
    real(kind=rp), dimension(:,:,:), allocatable, target :: ufb,vfb
!!! dirty reshape arrays indexing ijk -> kji !!!

    integer(kind=ip) :: i,j,k

    integer(kind=ip), save :: iter_fluxes=-1
    iter_fluxes = iter_fluxes + 1

    if (myrank==0) write(*,*)' nhmg_fluxes:',iter_fluxes

    call tic(1,'nhmg_fluxes')

!!! dirty reshape arrays indexing ijk -> kji !!!
    !allocate(ub(1:nz  ,0:ny+1,0:nx+1))
    !allocate(vb(1:nz  ,0:ny+1,0:nx+1))
    !allocate(wb(1:nz+1,0:ny+1,0:nx+1))
    allocate(ub(1:nz  ,-1:ny+2,-1:nx+2)) ! ND 11/01
    allocate(vb(1:nz  ,-1:ny+2,-1:nx+2)) ! ND 11/01
    allocate(wb(1:nz+1,-1:ny+2,-1:nx+2)) ! ND 11/01
    !do i = 1,nx+1
    !  do j = 0,ny+1
    do i = -1,nx+2
      do j = -1,ny+2
        do k = 1,nz
          ub(k,j,i) = ua(i,j,k)
        enddo
      enddo
    enddo
    !do j = 0,ny+1
    !   do k = 1,nz
    !      ub(k,j,0) = zero
    !   enddo
    !enddo
    !do i = 0,nx+1
    !  do j = 1,ny+1
    do i = -1,nx+2
      do j = -1,ny+2
        do k = 1,nz
          vb(k,j,i) = va(i,j,k)
        enddo
      enddo
    enddo
    !do i = 0,nx+1
    !   do k = 1,nz
    !      vb(k,0,i) = zero
    !   enddo
    !enddo
    !do i = 0,nx+1
    !  do j = 0,ny+1
    do i = -1,nx+2
      do j = -1,ny+2
        do k = 1,nz+1
          wb(k,j,i) = wa(i,j,k)
        enddo
      enddo
    enddo
    u => ub
    v => vb
    w => wb
    !allocate(ufb(1:nz,0:ny+1,0:nx+1))
    !allocate(vfb(1:nz,0:ny+1,0:nx+1))
    allocate(ufb(1:nz,0:ny+1,1:nx+1))
    allocate(vfb(1:nz,1:ny+1,0:nx+1))
    uf => ufb
    vf => vfb
!!! dirty reshape arrays indexing ijk -> kji !!!

!    u => ua
!    v => va
!    w => wa
!    uf => ufa
!    vf => vfa

    if (check_output) then
       !if ((iter_fluxes .EQ. 1)) then
          call write_netcdf(u,vname='uin',netcdf_file_name='fl.nc',rank=myrank,iter=iter_fluxes)
          call write_netcdf(v,vname='vin',netcdf_file_name='fl.nc',rank=myrank,iter=iter_fluxes)
          call write_netcdf(w,vname='win',netcdf_file_name='fl.nc',rank=myrank,iter=iter_fluxes)
       !endif
    endif

    call set_fluxes(u,v,w,uf,vf)

    if (check_output) then
       !if ((iter_fluxes .EQ. 1)) then
          call write_netcdf(uf,vname='uf',netcdf_file_name='fl.nc',rank=myrank,iter=iter_fluxes)
          call write_netcdf(vf,vname='vf',netcdf_file_name='fl.nc',rank=myrank,iter=iter_fluxes)
       !endif
    endif

!!! dirty reshape arrays indexing kji -> ijk !!!
   do i = 1,nx+1
      do j = 0,ny+1
        do k = 1,nz
          ufa(i,j,k) = uf(k,j,i)
        enddo
      enddo
    enddo
    do i = 0,nx+1
      do j = 1,ny+1
        do k = 1,nz
          vfa(i,j,k) = vf(k,j,i)
        enddo
      enddo
    enddo
!!! dirty reshape arrays indexing kji -> ijk !!!

    if (associated(u)) u => null()
    if (associated(v)) v => null()
    if (associated(w)) w => null()
    if (associated(uf)) uf => null()
    if (associated(vf)) vf => null()

!!! dirty reshape arrays indexing kji -> ijk !!!
    deallocate(ub)
    deallocate(vb)
    deallocate(wb)
    deallocate(ufb)
    deallocate(vfb)
!!! dirty reshape arrays indexing kji -> ijk !!!

    call toc(1,'nhmg_fluxes')	
 
  end subroutine nhmg_fluxes

  !--------------------------------------------------------------
  subroutine nhmg_matrices(nx,ny,nz,zra,Hza,dxa,dya)

    integer(kind=ip), intent(in) :: nx, ny, nz

    real(kind=rp), dimension(0:nx+1,0:ny+1,1:nz)          , intent(in) :: zra
    real(kind=rp), dimension(0:nx+1,0:ny+1,1:nz)          , intent(in) :: Hza
    real(kind=rp), dimension(0:nx+1,0:ny+1)     , optional, intent(in) :: dxa
    real(kind=rp), dimension(0:nx+1,0:ny+1)     , optional, intent(in) :: dya

    real(kind=rp), dimension(:,:,:), pointer :: zr,Hz
    real(kind=rp), dimension(:,:)  , pointer :: dx,dy

!!! dirty reshape arrays indexing ijk -> kji !!!
    integer(kind=ip) :: i,j,k
    real(kind=rp), dimension(1:nz,0:ny+1,0:nx+1), target :: zrb
    real(kind=rp), dimension(1:nz,0:ny+1,0:nx+1), target :: Hzb 
    real(kind=rp), dimension(0:ny+1,0:nx+1),      target :: dxb,dyb
!!!

    integer(kind=ip), save :: iter_matrices=0
    iter_matrices = iter_matrices + 1

    if (myrank==0) write(*,*)' nhmg_matrices: ',iter_matrices

    !--------------------!
    !- Horizontal grids -!
    !--------------------!
    if (present(dxa) .and. present(dya)) then

       dxb = transpose(dxa)
       dyb = transpose(dya)
       dx => dxb
       dy => dyb

       call set_horiz_grids(dx,dy)

       if (check_output) then
          call write_netcdf(grid(1)%dx,vname='dx',netcdf_file_name='mat.nc',rank=myrank,iter=iter_matrices)
          call write_netcdf(grid(1)%dy,vname='dy',netcdf_file_name='mat.nc',rank=myrank,iter=iter_matrices)
       endif

       if (associated(dx)) dx => null()
       if (associated(dy)) dy => null()

    end if

    !------------------!
    !- Vertical grids -!
    !------------------!
!!! dirty reshape arrays indexing ijk -> kji !!!
    do i = 0,nx+1
       do j = 0,ny+1
          do k = 1,nz
             zrb(k,j,i) = zra(i,j,k)
             Hzb(k,j,i) = Hza(i,j,k)
          enddo
       enddo
    enddo
    zr => zrb
    Hz => Hzb
!!!

    call set_vert_grids(zr,Hz)

    if (check_output) then
       !if ((iter_matrices .EQ. 1) .OR. (iter_matrices .EQ. 2)) then
          call write_netcdf(grid(1)%zr,vname='zr',netcdf_file_name='mat.nc',rank=myrank,iter=iter_matrices)
          call write_netcdf(grid(1)%dz,vname='dz',netcdf_file_name='mat.nc',rank=myrank,iter=iter_matrices)
       !endif
    endif

    if (associated(zr)) zr => null()
    if (associated(Hz)) Hz => null()

    !------------!
    !- matrices -!
    !------------!
    call set_matrices()

    if (check_output) then
       !if ((iter_matrices .EQ. 1) .OR. (iter_matrices .EQ. 2)) then
          call write_netcdf(grid(1)%cA,vname='cA',netcdf_file_name='mat.nc',rank=myrank,iter=iter_matrices)
       !endif
    endif

  end subroutine nhmg_matrices

  !--------------------------------------------------------------
  subroutine nhmg_bbc(nx,ny,nz,ua,va,wa)

    integer(kind=ip), intent(in) :: nx, ny, nz

    !real(kind=rp), dimension(1:nx+1,0:ny+1,1:nz  ), target, intent(in) :: ua
    !real(kind=rp), dimension(0:nx+1,1:ny+1,1:nz  ), target, intent(in) :: va
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz  ), target, intent(in) :: ua !ND 13/01
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz  ), target, intent(in) :: va !ND 13/01
    real(kind=rp), dimension( 0:nx+1,0:ny+1,1:nz+1)  , target, intent(inout) :: wa

    real(kind=rp), dimension(:,:,:), pointer :: u,v,w

!!! dirty reshape arrays indexing ijk -> kji !!!
    real(kind=rp), dimension(:,:,:), allocatable, target :: ub,vb,wb
    integer(kind=ip) :: i,j,k
!!!

    integer(kind=ip), save :: iter_bbc=0
    iter_bbc = iter_bbc + 1

    if (myrank==0) write(*,*)' nhmg_bbc:',iter_bbc

!!! dirty reshape arrays indexing ijk -> kji !!!
    !allocate( ub(1:nz,0:ny+1,0:nx+1))
    !allocate( vb(1:nz,0:ny+1,0:nx+1))
    allocate( ub(1:nz  ,-1:ny+2,-1:nx+2)) !!ND 12/01
    allocate( vb(1:nz  ,-1:ny+2,-1:nx+2)) !ND 12/01
    allocate( wb(1:nz+1, 0:ny+1, 0:nx+1))
    do i = -1,nx+2
      do j = -1,ny+2
        do k = 1,nz
          ub(k,j,i) = ua(i,j,k)
        enddo
      enddo
    enddo
    !i = 0
    !do j = 0,ny+1
    !   do k = 1,nz
    !      ub(k,j,i) = zero
    !   enddo
    !enddo
    do i = -1,nx+2
      do j = -1,ny+2
        do k = 1,nz
          vb(k,j,i) = va(i,j,k)
        enddo
      enddo
    enddo
    !j = 0
    !do i = 0,nx+1
    !   do k = 1,nz
    !      vb(k,j,i) = zero
    !   enddo
    !enddo
    do i = 0,nx+1
      do j = 0,ny+1
        do k = 1,nz+1
          wb(k,j,i) = wa(i,j,k)
        enddo
      enddo
    enddo
    u => ub
    v => vb
    w => wb
!!!

    !call fill_halo(1,u,lbc_null='u') !ND 12/01
    !call fill_halo(1,v,lbc_null='v') !ND 12/01
    !call fill_halo(1,w)              !ND 12/01

    call set_bbc(u,v,w)

!!! dirty reshape arrays indexing kji -> ijk !!!
    do i = 0,nx+1
      do j = 0,ny+1
        do k = 1,nz+1
          wa(i,j,k) = w(k,j,i)
        enddo
      enddo
    enddo
!!!

    if (associated(u)) u => null()
    if (associated(v)) v => null()
    if (associated(w)) w => null()

!!! dirty reshape arrays indexing kji -> ijk !!!
    deallocate(ub)
    deallocate(vb)
    deallocate(wb)
!!!

  end subroutine nhmg_bbc

  !--------------------------------------------------------------
  subroutine nhmg_coupling(nx,ny,nz,uf_bara,vf_bara,ua,va,wa,ufa,vfa)

    integer(kind=ip), intent(in) :: nx, ny, nz

    real(kind=rp), dimension(1:nx+1,0:ny+1),        target, intent(in) :: uf_bara
    real(kind=rp), dimension(0:nx+1,1:ny+1),        target, intent(in) :: vf_bara

    real(kind=rp), dimension(1:nx+1,0:ny+1,1:nz  ), target, intent(inout) :: ua
    real(kind=rp), dimension(0:nx+1,1:ny+1,1:nz  ), target, intent(inout) :: va
    real(kind=rp), dimension(0:nx+1,0:ny+1,1:nz+1), target, intent(inout) :: wa

    real(kind=rp), dimension(1:nx+1,0:ny+1,1:nz  ), target, optional, intent(out):: ufa
    real(kind=rp), dimension(0:nx+1,1:ny+1,1:nz  ), target, optional, intent(out):: vfa

    real(kind=rp), dimension(:,:),   pointer :: uf_bar,vf_bar
    real(kind=rp), dimension(:,:,:), pointer :: u,v,w
    real(kind=rp), dimension(:,:,:), pointer :: uf,vf

!!! dirty reshape arrays indexing ijk -> kji !!!
    real(kind=rp), dimension(:,:),   allocatable, target :: uf_barb,vf_barb 
    real(kind=rp), dimension(:,:,:), allocatable, target :: ub,vb,wb
    real(kind=rp), dimension(:,:,:), allocatable, target :: ufb,vfb
!!! dirty reshape arrays indexing ijk -> kji !!!

    integer(kind=ip) :: i,j,k

    integer(kind=ip), save :: iter_coupling=0
    iter_coupling = iter_coupling + 1

    if (myrank==0) write(*,*)' nhmg_coupling:',iter_coupling

    call tic(1,'nhmg_coupling')

!!! dirty reshape arrays indexing ijk -> kji !!!
    !allocate(uf_barb(0:ny+1,0:nx+1))   
    !allocate(vf_barb(0:ny+1,0:nx+1))  
    !allocate(ub(1:nz  ,0:ny+1,0:nx+1))
    !allocate(vb(1:nz  ,0:ny+1,0:nx+1))
    !allocate(wb(1:nz+1,0:ny+1,0:nx+1))
    allocate(uf_barb(0:ny+1,1:nx+1))  !ND 11/01
    allocate(vf_barb(1:ny+1,0:nx+1))  !ND 11/01
    allocate(ub(1:nz  ,0:ny+1,1:nx+1))!ND 11/01
    allocate(vb(1:nz  ,1:ny+1,0:nx+1))!ND 11/01
    allocate(wb(1:nz+1,0:ny+1,0:nx+1))!ND 11/01
    do i = 1,nx+1
      do j = 0,ny+1
          uf_barb(j,i) = uf_bara(i,j)
      enddo
    enddo
    do i = 0,nx+1
      do j = 1,ny+1
          vf_barb(j,i) = vf_bara(i,j)
      enddo
    enddo
    do i = 1,nx+1
      do j = 0,ny+1
        do k = 1,nz
          ub(k,j,i) = ua(i,j,k)
        enddo
      enddo
    enddo
    do i = 0,nx+1
      do j = 1,ny+1
        do k = 1,nz
          vb(k,j,i) = va(i,j,k)
        enddo
      enddo
    enddo
    do i = 0,nx+1
      do j = 0,ny+1
        do k = 1,nz+1
          wb(k,j,i) = wa(i,j,k)
        enddo
      enddo
    enddo
    uf_bar => uf_barb
    vf_bar => vf_barb
    u => ub
    v => vb
    w => wb
!!! dirty reshape arrays indexing ijk -> kji !!!

!    uf_bar => uf_bara
!    vf_bar => vf_bara
!    u => ua
!    v => va
!    w => wa

!    call fill_halo(1,uf_bar,lbc_null='u') !ND 11/01
!    call fill_halo(1,vf_bar,lbc_null='v') !ND 11/01
!    call fill_halo(1,u,lbc_null='u')      !ND 11/01
!    call fill_halo(1,v,lbc_null='v')      !ND 11/01
!    call fill_halo(1,w)                   !ND 11/01

    if (check_output) then
       !if ((iter_coupling .EQ. 1) .OR. (iter_coupling .EQ. 2)) then
          call write_netcdf(uf_bar,vname='uf_bar',netcdf_file_name='co.nc',rank=myrank,iter=iter_coupling)
          call write_netcdf(vf_bar,vname='vf_bar',netcdf_file_name='co.nc',rank=myrank,iter=iter_coupling)
          call write_netcdf(u,vname='uin',netcdf_file_name='co.nc',rank=myrank,iter=iter_coupling)
          call write_netcdf(v,vname='vin',netcdf_file_name='co.nc',rank=myrank,iter=iter_coupling)
          call write_netcdf(w,vname='win',netcdf_file_name='co.nc',rank=myrank,iter=iter_coupling)
       !endif
    endif

    if ((present(ufa)).and.(present(vfa))) then
 
       !allocate(ufb(1:nz,0:ny+1,0:nx+1))
       !allocate(vfb(1:nz,0:ny+1,0:nx+1))
       allocate(ufb(1:nz,0:ny+1,1:nx+1)) 
       allocate(vfb(1:nz,1:ny+1,0:nx+1)) 
       uf => ufb
       vf => vfb

!       uf => ufa
!       vf => vfa

       call bt2bc_coupling(uf_bar,vf_bar,u,v,w,uf,vf)

    else

       call bt2bc_coupling(uf_bar,vf_bar,u,v,w)

    endif

    if (check_output) then
       !if ((iter_coupling .EQ. 1) .OR. (iter_coupling .EQ. 2)) then
          call write_netcdf(u,vname='uout',netcdf_file_name='co.nc',rank=myrank,iter=iter_coupling)
          call write_netcdf(v,vname='vout',netcdf_file_name='co.nc',rank=myrank,iter=iter_coupling)
          call write_netcdf(w,vname='wout',netcdf_file_name='co.nc',rank=myrank,iter=iter_coupling)
          if ((present(ufa)).and.(present(vfa))) then
             call write_netcdf(uf,vname='uf',netcdf_file_name='co.nc',rank=myrank,iter=iter_coupling)
             call write_netcdf(vf,vname='vf',netcdf_file_name='co.nc',rank=myrank,iter=iter_coupling)
          endif
       !endif
    endif
    
    !!- check non-divergence of the corrected u,v,w
    !call set_rhs(u,v,w)
    !
    !if (check_output) then
    !   !if ((iter_coupling .EQ. 0) .OR. (iter_coupling .EQ. 1)) then
    !      call write_netcdf(grid(1)%b,vname='bout',netcdf_file_name='co.nc',rank=myrank,iter=iter_coupling)
    !   !endif
    !endif

!!! dirty reshape arrays indexing kji -> ijk !!!
    do i = 1,nx+1
      do j = 0,ny+1
        do k = 1,nz
          ua(i,j,k) = u(k,j,i)
        enddo
      enddo
    enddo
    do i = 0,nx+1
      do j = 1,ny+1
        do k = 1,nz
          va(i,j,k) = v(k,j,i)
        enddo
      enddo
    enddo
    do i = 0,nx+1
      do j = 0,ny+1
        do k = 1,nz+1
          wa(i,j,k) = w(k,j,i)
        enddo
      enddo
    enddo
    if ((present(ufa)).and.(present(vfa))) then
    do i = 1,nx+1
      do j = 0,ny+1
        do k = 1,nz
          ufa(i,j,k) = uf(k,j,i)
        enddo
      enddo
    enddo
    do i = 0,nx+1
      do j = 1,ny+1
        do k = 1,nz
          vfa(i,j,k) = vf(k,j,i)
        enddo
      enddo
    enddo
    deallocate(ufb)
    deallocate(vfb)
    endif
    deallocate(uf_barb)
    deallocate(vf_barb)
    deallocate(ub)
    deallocate(vb)
    deallocate(wb)
!!! dirty reshape arrays indexing kji -> ijk !!!

    call toc(1,'nhmg_coupling')

  end subroutine nhmg_coupling

  !--------------------------------------------------------------
  subroutine nhmg_solve(nx,ny,nz,ua,va,wa,rua,rva,rwa,dt,rufrca,rvfrca)

    integer(kind=ip), intent(in) :: nx, ny, nz

    real(kind=rp), dimension(1:nx+1,0:ny+1,1:nz  ), target, intent(inout) :: ua
    real(kind=rp), dimension(0:nx+1,1:ny+1,1:nz  ), target, intent(inout) :: va
    real(kind=rp), dimension(0:nx+1,0:ny+1,1:nz+1), target, intent(inout) :: wa

    real(kind=rp), dimension(1:nx+1,0:ny+1,1:nz  ), target, optional, intent(in) :: rua
    real(kind=rp), dimension(0:nx+1,1:ny+1,1:nz  ), target, optional, intent(in) :: rva
    real(kind=rp), dimension(0:nx+1,0:ny+1,1:nz+1), target, optional, intent(in) :: rwa
    real(kind=rp),                                          optional, intent(in) :: dt
    real(kind=rp), dimension(1:nx+1,0:ny+1),        target, optional, intent(out):: rufrca
    real(kind=rp), dimension(0:nx+1,1:ny+1),        target, optional, intent(out):: rvfrca


    real(kind=rp), dimension(:,:,:), pointer :: u,v,w
    real(kind=rp), dimension(:,:,:), pointer :: ru,rv,rw
    real(kind=rp), dimension(:,:)  , pointer :: rufrc,rvfrc

!!! dirty reshape arrays indexing
    real(kind=rp), dimension(:,:,:), allocatable, target :: ub,vb,wb
    real(kind=rp), dimension(:,:,:), allocatable, target :: rub,rvb,rwb
    real(kind=rp), dimension(:,:,:), allocatable, target :: ruc,rvc
    real(kind=rp), dimension(:,:),   allocatable, target :: rufrcb,rvfrcb
    integer(kind=ip) :: i,j,k
!!! 

    integer(kind=ip), save :: iter_solve=0
    iter_solve = iter_solve + 1

    if (myrank==0) write(*,*)' nhmg_solve:',iter_solve

    call tic(1,'nhmg_solve')

!!! dirty reshape arrays indexing ijk -> kji !!!
    allocate(ub(1:nz  ,0:ny+1,1:nx+1))
    allocate(vb(1:nz  ,1:ny+1,0:nx+1))
    allocate(wb(1:nz+1,0:ny+1,0:nx+1))
    do i = 1,nx+1
      do j = 0,ny+1
        do k = 1,nz
          ub(k,j,i) = ua(i,j,k)
        enddo
      enddo
    enddo
    do i = 0,nx+1
      do j = 1,ny+1
        do k = 1,nz
          vb(k,j,i) = va(i,j,k)
        enddo
      enddo
    enddo
    do i = 0,nx+1
       do j = 0,ny+1
          do k = 1,nz+1
             wb(k,j,i) = wa(i,j,k)
          enddo
       enddo
    enddo
    u => ub
    v => vb
    w => wb
!!! dirty reshape arrays indexing ijk -> kji !!!

    !u => ua
    !v => va
    !w => wa
   
    if (check_output) then
       !if ((iter_solve .EQ. 1) .OR. (iter_solve .EQ. 2)) then
          call write_netcdf(u,vname='uin',netcdf_file_name='so.nc',rank=myrank,iter=iter_solve)
          call write_netcdf(v,vname='vin',netcdf_file_name='so.nc',rank=myrank,iter=iter_solve)
          call write_netcdf(w,vname='win',netcdf_file_name='so.nc',rank=myrank,iter=iter_solve)
       !endif
    endif

    !- step 1 - 
    call set_rhs(u,v,w)

    if (check_output) then
       !if ((iter_solve .EQ. 1) .OR. (iter_solve .EQ. 2)) then
          call write_netcdf(grid(1)%b,vname='b',netcdf_file_name='so.nc',rank=myrank,iter=iter_solve)
       !endif
    endif

    !- step 2 -
    call solve_p()

    if (check_output) then
       !if ((iter_solve .EQ. 1) .OR. (iter_solve .EQ. 2)) then
          call write_netcdf(grid(1)%p,vname='p',netcdf_file_name='so.nc',rank=myrank,iter=iter_solve)
          call write_netcdf(grid(1)%r,vname='r',netcdf_file_name='so.nc',rank=myrank,iter=iter_solve)
       !endif
    endif

    !- step 3 -
    call correct_uvw(u,v,w)

    if (check_output) then
       !if ((iter_solve .EQ. 1) .OR. (iter_solve .EQ. 2)) then
          call write_netcdf(u,vname='uout',netcdf_file_name='so.nc',rank=myrank,iter=iter_solve)
          call write_netcdf(v,vname='vout',netcdf_file_name='so.nc',rank=myrank,iter=iter_solve)
          call write_netcdf(w,vname='wout',netcdf_file_name='so.nc',rank=myrank,iter=iter_solve)
       !endif
    endif

    !- check - non-divergence of the projected u,v,w
    !call set_rhs(u,v,w)
    !
    !if (check_output) then
    !   !if ((iter_solve .EQ. 1) .OR. (iter_solve .EQ. 2)) then
    !      call write_netcdf(grid(1)%b,vname='bout',netcdf_file_name='so.nc',rank=myrank,iter=iter_solve)
    !   !endif
    !endif

    !- step 4 - bc2bt coupling
    if ((present(rua)).and.(present(rva)).and.(present(rwa))) then

       !!! dirty reshape arrays indexing ijk -> kji !!!
       !allocate(rub(1:nz  ,0:ny+1,1:nx+1))
       !allocate(rvb(1:nz  ,1:ny+1,0:nx+1))
       allocate(rub(1:nz  ,0:ny+1,0:nx+1)) !ND 17/01
       allocate(rvb(1:nz  ,0:ny+1,0:nx+1)) !ND 17/01
       allocate(rwb(1:nz+1,0:ny+1,0:nx+1))
       do i = 1,nx+1
          do j = 0,ny+1
             do k = 1,nz
                rub(k,j,i) = rua(i,j,k)
             enddo
          enddo
       enddo
       do i = 0,nx+1
          do j = 1,ny+1
             do k = 1,nz
                rvb(k,j,i) = rva(i,j,k)
             enddo
          enddo
       enddo
       do i = 0,nx+1
          do j = 0,ny+1
             do k = 1,nz+1
                rwb(k,j,i) = rwa(i,j,k)
             enddo
          enddo
       enddo
       ru => rub
       rv => rvb
       rw => rwb
       call fill_halo(1,ru) !ND 17/01
       call fill_halo(1,rv) !ND 17/01
       call fill_halo(1,rw) !ND 17/01
       !
       allocate(ruc(1:nz  ,0:ny+1,1:nx+1)) !ND 17/01
       allocate(rvc(1:nz  ,1:ny+1,0:nx+1)) !ND 17/01
       do i = 1,nx+1
          do j = 0,ny+1
             do k = 1,nz
                ruc(k,j,i) = rub(k,j,i)
             enddo
          enddo
       enddo
       do i = 0,nx+1
          do j = 1,ny+1
             do k = 1,nz
                rvc(k,j,i) = rvb(k,j,i)
             enddo
          enddo
       enddo
       if (associated(ru)) ru => null()
       if (associated(rv)) rv => null()
       ru => ruc !ND 17/01
       rv => rvc !ND 17/01
       !
       allocate(rufrcb(0:ny+1,1:nx+1))
       allocate(rvfrcb(1:ny+1,0:nx+1))
       rufrc => rufrcb
       rvfrc => rvfrcb
       !!! 
   
       if (check_output) then
          !if ((iter_solve .EQ. 1) .OR. (iter_solve .EQ. 2)) then
             call write_netcdf(ru,vname='ruin',netcdf_file_name='so.nc',rank=myrank,iter=iter_solve)
             call write_netcdf(rv,vname='rvin',netcdf_file_name='so.nc',rank=myrank,iter=iter_solve)
             call write_netcdf(rw,vname='rwin',netcdf_file_name='so.nc',rank=myrank,iter=iter_solve)
          !endif
       endif

       call bc2bt_coupling(ru,rv,rw,dt,rufrc,rvfrc)

       if (check_output) then
          !if ((iter_solve .EQ. 1) .OR. (iter_solve .EQ. 2)) then
             call write_netcdf(rufrc,vname='rufrc',netcdf_file_name='so.nc',rank=myrank,iter=iter_solve)
             call write_netcdf(rvfrc,vname='rvfrc',netcdf_file_name='so.nc',rank=myrank,iter=iter_solve)
          !endif
       endif

    endif

    if (associated(u)) u => null()
    if (associated(v)) v => null()
    if (associated(w)) w => null()

!!! dirty reshape arrays indexing kji -> ijk !!!
    do i = 1,nx+1
       do j = 0,ny+1
          do k = 1,nz
             ua(i,j,k) = ub(k,j,i)
          enddo
       enddo
    enddo
    do i = 0,nx+1
       do j = 1,ny+1
          do k = 1,nz
             va(i,j,k) = vb(k,j,i)
          enddo
       enddo
    enddo
    do i = 0,nx+1
       do j = 0,ny+1
          do k = 1,nz+1
             wa(i,j,k) = wb(k,j,i)
          enddo
       enddo
    enddo
    deallocate(ub)
    deallocate(vb)
    deallocate(wb)
    if ((present(rua)).and.(present(rva))) then
    deallocate(rub)
    deallocate(rvb)
    deallocate(ruc)
    deallocate(rvc)
    deallocate(rwb)
    do i = 1,nx+1
       do j = 0,ny+1      
          rufrca(i,j) = rufrcb(j,i)
       enddo
    enddo
    do i = 0,nx+1
       do j = 1,ny+1
          rvfrca(i,j) = rvfrcb(j,i) 
       enddo
    enddo
    deallocate(rufrcb)
    deallocate(rvfrcb)
    endif
!!! dirty reshape arrays indexing kji -> ijk !!!

    call toc(1,'nhmg_solve')

  end subroutine nhmg_solve

  !--------------------------------------------------------------
  subroutine nhmg_clean()

    real(kind=rp) :: tstart,tend,perf

    call grids_dealloc()

    call print_tictoc()

  end subroutine nhmg_clean

end module nhmg
