!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !===================================================================!
  !- All constants of the program have to be declared in this module -!
  !===================================================================!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module mg_cst

  implicit none

  !======================================!
  !- DECLARATIONS of GLOBAL VARIABLE(S) -!
  !======================================!
  integer(kind = 4), parameter :: ip = 4
  integer(kind = 4), parameter :: rp = 8

  real(kind = rp)   , parameter :: rZero   =  0.0_rp
  real(kind = rp)   , parameter :: rHalf   =  0.5_rp
  real(kind = rp)   , parameter :: rQuarter=  0.25_rp
  real(kind = rp)   , parameter :: rOne    =  1.0_rp
  real(kind = rp)   , parameter :: rTwo    =  2.0_rp
  real(kind = rp)   , parameter :: rThree  =  3.0_rp
  real(kind = rp)   , parameter :: rFour   =  4.0_rp
  real(kind = rp)   , parameter :: rFive   =  5.0_rp
  real(kind = rp)   , parameter :: rTen    = 10.0_rp
  real(kind = rp)   , parameter :: rOneHundred =  100.0_rp
  real(kind = rp)   , parameter :: rcircle = 360_rp
  real(kind = rp)   , parameter :: rOneThousan = 1000.0_rp
  real(kind = rp)   , parameter :: rMax    = 1.e20_rp
  

  integer(kind = ip), parameter :: iZero   =  0_ip
  integer(kind = ip), parameter :: iOne    =  1_ip
  integer(kind = ip), parameter :: iTwo    =  2_ip
  integer(kind = ip), parameter :: iThree  =  3_ip
  integer(kind = ip), parameter :: iFour   =  4_ip
  integer(kind = ip), parameter :: iFive   =  5_ip
  integer(kind = ip), parameter :: iSix    =  6_ip
  integer(kind = ip), parameter :: iSeven  =  7_ip
  integer(kind = ip), parameter :: iEight  =  8_ip
  integer(kind = ip), parameter :: iNine   =  9_ip
  integer(kind = ip), parameter :: iTen    = 10_ip
  integer(kind = ip), parameter :: iTwenty = 20_ip
  integer(kind = ip), parameter :: iHigh   = 9999999_ip
  integer(kind = ip), parameter :: iMega   = 1000000_ip
  integer(kind = ip), parameter :: iMaxSubDomEddies = iMega

  integer(kind = ip), parameter :: nb_max_tracers = 50_ip
  integer(kind = ip), parameter :: pcat = iFour ! simon's cat www.simonscat.com

  real(kind = rp), parameter :: mask_value = 1.e20_rp
  real(kind = rp), parameter :: missing_value = 1.e20_rp

  integer(kind=ip), parameter :: i_defp =  iOne
  integer(kind=ip), parameter :: i_defn = -iOne
  integer(kind=ip), parameter :: i_high = 9999999_ip
  real(kind=rp)   , parameter :: r_def  = rZero

  character(len = 4) , parameter :: c_def  = 'NONE'
  logical            , parameter :: l_def  = .false.

  real(kind=rp)   , parameter :: pi = acos(-1._rp)


end module mg_cst
