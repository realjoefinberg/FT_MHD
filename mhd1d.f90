! mhd1d.f90
program mhd1d
  use parameters         ! provides dp, nx, nvar, dx, CFL, gamma
  use mhd_types          ! provides type(state)
  use mhd_init           ! provides set_brio_wu
  use mhd_flux           ! provides compute_flux
  use mhd_wave      ! provides max_wave_speed(U_old, nx)

  implicit none

  real(dp), dimension(nvar,nx) :: U_old, U_new, F
  real(dp)                     :: dt, t, t_end
  integer                       :: i

  ! Initialize time
  t     = 0.0_dp
  t_end = 0.2_dp

  !-------------------------------------------------------------------
  ! 1) Initialize Brio–Wu shock tube
  !-------------------------------------------------------------------
  do i = 1, nx
    if (i <= nx/2) then
      call set_brio_wu(U_old(:,i), &
           1.0_dp,   0.0_dp, 0.0_dp, 0.0_dp, &  ! rho, vx, vy, vz
           0.75_dp,  1.0_dp, 0.0_dp, &          ! Bx, By, Bz
           1.0_dp)                             ! p
    else
      call set_brio_wu(U_old(:,i), &
           0.125_dp, 0.0_dp, 0.0_dp, 0.0_dp, &  ! rho, vx, vy, vz
           0.75_dp, -1.0_dp, 0.0_dp, &          ! Bx, By, Bz
           0.1_dp)                             ! p
    end if
  end do

    !-------------------------------------------------------------------
  ! DEBUG: print the very first cell and its flux, then stop
  !-------------------------------------------------------------------
  print *, 'U_old(:,1) =', U_old(:,1)
  call compute_flux( U_old(:,1), F(:,1) )
  print *, 'F(:,1)     =', F(:,1)
  stop


  !-------------------------------------------------------------------
  ! 2–6) Main time‑stepping loop (skeleton)
  !-------------------------------------------------------------------
  do while (t < t_end)

    ! 2) CFL time step
    dt = CFL * dx / max_wave_speed(U_old, nx) 

    ! 3) Compute flux at each cell
    do i = 1, nx
      call compute_flux(U_old(:,i), F(:,i))
    end do

    ! 4) Finite-volume update (simple backward difference)
    do i = 2, nx
      U_new(:,i) = U_old(:,i) - (dt/dx) * (F(:,i) - F(:,i-1))
    end do

    ! 5) Periodic boundary conditions
    U_new(:,1)  = U_new(:,nx)
    U_new(:,nx) = U_new(:,1)

    ! 6) Advance
    U_old = U_new
    t     = t + dt

  end do

  !-------------------------------------------------------------------
  ! 7) Write final state
  !-------------------------------------------------------------------
  open(unit=10, file="mhd1d_final.dat", status="replace")
  do i = 1, nx
    write(10,'(1X,8(1PE12.5))') U_old(:,i)
  end do
  close(10)

end program mhd1d
