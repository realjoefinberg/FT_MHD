module parameters
  implicit none
  integer, parameter :: dp    = kind(1.0d0)
  integer, parameter :: nx    = 400
  integer, parameter :: nvar  = 8        ! ρ, ρvₓ, ρv_y, ρv_z, Bₓ, B_y, B_z, E
  real(dp),   parameter :: x0  = 0.0_dp
  real(dp),   parameter :: x1  = 1.0_dp
  real(dp),   parameter :: dx  = (x1 - x0) / nx
  real(dp),   parameter :: CFL = 0.4_dp
  real(dp),   parameter :: gamma = 5.0_dp/3.0_dp
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
    !-------------------------------------------------------------------
    ! Compute 1D ideal‑MHD flux Fi(:) from conserved state Ui(:)
    !-------------------------------------------------------------------
    real(dp), intent(in)  :: Ui(nvar)
    real(dp), intent(out) :: Fi(nvar)

    !— unpack Ui into (rho, vx, …, vdotB) —!
    include 'flux_unpack_mhd1d_1.inc'

    !— assemble the flux components —!
    include 'flux_compute_mhd1d_1.inc'

  end subroutine compute_flux

end module mhd_flux

program mhd1d
  use parameters
  use mhd_types
  use mhd_flux
  implicit none

  real(dp), dimension(nvar,nx) :: U_old, U_new, F
  real(dp) :: dt, t, t_end
  integer   :: i

  ! Initialize time
  t     = 0.0_dp
  t_end = 0.2_dp

  ! 1) Initialize Brio–Wu shock tube
  do i = 1, nx
    if (i <= nx/2) then
      ! left state: ρ=1, v=(0,0,0), B=(0.75,1,0), p=1
      call set_brio_wu(U_old(:,i), 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                       0.75_dp, 1.0_dp, 0.0_dp, 1.0_dp, gamma)
    else
      ! right state: ρ=0.125, v=(0,0,0), B=(0.75,-1,0), p=0.1
      call set_brio_wu(U_old(:,i), 0.125_dp, 0.0_dp, 0.0_dp, 0.0_dp,&
                       0.75_dp, -1.0_dp, 0.0_dp, 0.1_dp, gamma)
    end if
  end do

  ! Main time‑stepping loop
  do while (t < t_end)
    ! 2) Compute dt from CFL condition (you'll need a wave‐speed routine)
    !    dt = CFL*dx / max_speed(U_old)
    ! 3) Compute fluxes in each cell
    do i = 1, nx
      call compute_flux(U_old(:,i), F(:,i))
    end do
    ! 4) Update U_new (finite‐volume update using F(:,i))
    ! 5) Apply boundary conditions on U_new
    ! 6) Swap U_old <-> U_new; t = t + dt
    exit  ! remove once steps 2–6 are implemented
  end do

  ! 7) Output final grid
  open(unit=10, file="mhd1d_final.dat", status="replace")
  do i = 1, nx
    write(10,'(1X,8(1PE12.5))') U_old(:,i)
  end do
  close(10)

end program mhd1d
