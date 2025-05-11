! mhd1d.f90
program mhd1d
  use parameters       ! provides dp, nx, nvar, dx, CFL, gamma
  use mhd_types        ! provides type(state)
  use mhd_init         ! provides set_brio_wu
  use mhd_flux         ! provides compute_flux
  implicit none

  !-------------------------------------------------------------------
  ! Arrays: U_old/U_new(:,i) hold conserved vars at cell i
  !         F(:,i)      holds flux at cell i
  !-------------------------------------------------------------------
  real(dp), dimension(nvar,nx) :: U_old, U_new, F
  real(dp)                     :: dt, t, t_end, max_speed
  integer                       :: i

  ! Initialize time
  t     = 0.0_dp
  t_end = 0.2_dp

  !-------------------------------------------------------------------
  ! 1) Initialize Brio–Wu shock tube
  !    Left state for i<=nx/2, right state for i>nx/2
  !-------------------------------------------------------------------
  do i = 1, nx
    if (i <= nx/2) then
      call set_brio_wu(U_old(:,i), &
           1.0_dp,  0.0_dp, 0.0_dp, 0.0_dp, &  ! rho, vx, vy, vz
           0.75_dp, 1.0_dp, 0.0_dp,           &  ! Bx, By, Bz
           1.0_dp)                              ! gas pressure
    else
      call set_brio_wu(U_old(:,i), &
           0.125_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
           0.75_dp, -1.0_dp, 0.0_dp,
           0.1_dp)
    end if
  end do

  !-------------------------------------------------------------------
  ! 2–6) Main time‑stepping loop (skeleton)
  !    - compute dt via CFL
  !    - compute fluxes
  !    - update U_new
  !    - apply BCs
  !    - swap
  !-------------------------------------------------------------------
  do while (t < t_end)

    ! 2) Compute dt from CFL: need max wave speed routine (placeholder)
    !    For now, set a fixed small dt:
    dt = 0.5_dp * dx / 1.0_dp  

    ! 3) Compute flux F(:,i) at each cell
    do i = 1, nx
      call compute_flux(U_old(:,i), F(:,i))
    end do

    ! 4) Update U_new by simple upwind (1st order) / finite‑volume:
    !    U_new_i = U_old_i - (dt/dx)*(F_i - F_{i-1})
    !    here we do a backward difference in x
    do i = 2, nx
      U_new(:,i) = U_old(:,i) - (dt/dx) * (F(:,i) - F(:,i-1))
    end do

    ! 5) Boundary conditions (e.g., periodic)
    U_new(:,1)  = U_new(:,nx)
    U_new(:,nx) = U_new(:,1)

    ! 6) Advance
    U_old = U_new
    t     = t + dt

  end do

  !-------------------------------------------------------------------
  ! 7) Write final state to file
  !-------------------------------------------------------------------
  open(unit=10, file="mhd1d_final.dat", status="replace")
  do i = 1, nx
    write(10,'(1X,8(1PE12.5))') U_old(:,i)
  end do
  close(10)

end program mhd1d
