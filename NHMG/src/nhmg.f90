!!
!! nhmg module is here just to ceate a nhmg.mod to call the libnhmg.la library
!!
module nhmg

  implicit none

  private

  ! Overloads

!!$  interface define_matrices
!!$     module procedure           &
!!$          define_matrices_topo
!!$  end interface define_matrices

!!$  interface gather
!!$     module procedure &
!!$          gather_2D,  &
!!$          gather_3D
!!$  end interface gather

!!$  interface fill_halo
!!$     module procedure   &
!!$          fill_halo_2D, &
!!$          fill_halo_3D, &
!!$          fill_halo_4D
!!$  end interface fill_halo
!!$
!!$  interface write_netcdf
!!$     module procedure                &
!!$          sub_netcdf_write_fast_r2D, &
!!$          sub_netcdf_write_fast_r3D, &
!!$          sub_netcdf_write_fast_r3D_p, &
!!$          sub_netcdf_write_fast_r4D
!!$  end interface write_netcdf

!!$  interface setup_zr_zw
!!$     module procedure           &
!!$          setup_zr_zw_seamount, &
!!$          setup_zr_zw_croco
!!$  end interface setup_zr_zw


  ! Visibility

!!$  public :: read_nhnamelist
!!$
!!$  public :: btbc_coupling
!!$  public :: compute_barofrc
!!$  public :: compute_fluxes
!!$  public :: compute_rhs
!!$  public :: correct_uvw
!!$  public :: define_grids
!!$  public :: mg_mpi_init
!!$  public :: set_bbc
!!$  public :: solve_p

end module nhmg
