module parameters
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: nx = 400        ! number of grid cells
  real(dp), parameter :: x0 = 0.0_dp, x1 = 1.0_dp
  real(dp), parameter :: dx = (x1-x0)/nx
  real(dp), parameter :: CFL = 0.4_dp
  real(dp), parameter :: gamma = 5.0_dp/3.0_dp
  integer             :: nvar = 8       ! ρ,ρvₓ,ρv_y,ρv_z,Bₓ,B_y,B_z,E
end module parameters

module mhd_types
  use parameters, only: dp, nvar
  implicit none
  type :: state
    real(dp), dimension(nvar) :: U
  end type state
end module mhd_types

module mhd_flux
  use parameters, only: dp, gamma, nvar
  implicit none
contains
  subroutine compute_flux(Ui, Fi)
    ! Compute flux Fi(:) from state Ui(:)
    real(dp), intent(in)  :: Ui(nvar)
    real(dp), intent(out) :: Fi(nvar)
    !--- YOUR TASK: implement ideal‐MHD flux here ---
  end subroutine compute_flux
end module mhd_flux

program mhd1d
  use parameters
  use mhd_types
  use mhd_flux
  implicit none

  type(state), dimension(nx) :: grid
  real(dp)                   :: dt, t, t_end
  real(dp), dimension(nx, nvar) :: U_old, U_new, F
  integer :: i

  ! Initialize time
  t = 0.0_dp
  t_end = 0.2_dp

  ! 1) Initialize Brio–Wu shock tube: left state vs. right state
  do i=1,nx
    if (i .le. nx/2) then
      ! e.g. ρ=1, p=1, Bx=0.75, By=1, vz=0
      ! pack into U_old(:,i) ...
    else
      ! right state: ρ=0.125, p=0.1, Bx=0.75, By=-1, vz=0
    end if
  end do

  ! Main time‑stepping loop
  do while (t < t_end)
    ! 2) Compute dt from CFL condition using max wave speed
    ! 3) Compute fluxes F(:,i) at each cell
    ! 4) Update conserved variables:
    !    U_new(:,i) = U_old(:,i) - dt/dx*(F(:,i+1)-F(:,i))
    ! 5) Apply periodic (or reflecting) BCs on U_new
    ! 6) Swap U_old <‑> U_new, t = t + dt
  end do

  ! 7) Output final grid to file
  open(unit=10, file="mhd1d_final.dat", status="replace")
  do i=1,nx
    write(10,'(1X,8(1PE12.5))') U_old(:,i)
  end do
  close(10)

end program mhd1d
