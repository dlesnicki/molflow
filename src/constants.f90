
module constants

  use kinds

  implicit none

  real(kind=dp), parameter :: aufs=2.41888432649994E-02_dp
  real(kind=dp), parameter :: bohr=0.529177249_dp
  real(kind=dp), parameter :: lightspeed=2.99792458000000E+08_dp
  real(kind=dp), parameter :: pi=3.14159265358979323846264338_dp
  real(kind=dp), parameter :: ry=13.60569193_dp
  real(kind=dp), parameter :: evtokel=11604.505_dp
  real(kind=dp), parameter :: scmass=1822.888485_dp
  real(kind=dp), parameter :: au_kb=294210.1080_dp!315774.6821_dp
  real(kind=dp), parameter :: kb_au=1.0_dp/au_kb
  real(kind=dp), parameter :: au_kjm=2.62549962505E3_dp
  real(kind=dp), parameter :: au_kcm=6.275094706142E2_dp
  real(kind=dp), parameter :: kbolz=kb_au*bohr*bohr/aufs/aufs/scmass
  real(kind=dp), parameter :: elchg=1.602176462E-19_dp
  real(kind=dp), parameter :: esu=3.33564E-10_dp
  real(kind=dp), parameter :: eAtoDebye=elchg/esu*1E-8_dp*1E18_dp


end module constants
