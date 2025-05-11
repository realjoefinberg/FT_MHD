! parameters_mhd1d_1.f90
module parameters
  implicit none
  integer, parameter :: dp    = kind(1.0d0)
  integer, parameter :: nx    = 400        ! grid cells
  real(dp), parameter :: x0   = 0.0_dp
  real(dp), parameter :: x1   = 1.0_dp
  real(dp), parameter :: dx   = (x1 - x0) / nx
  real(dp), parameter :: CFL  = 0.4_dp
  real(dp), parameter :: gamma = 5.0_dp/3.0_dp
  integer, parameter :: nvar  = 8          ! # of conserved vars
end module parameters
