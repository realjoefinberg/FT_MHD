module mhd_wave
  use parameters, only: dp, gamma, nvar
  implicit none
contains

  ! Compute the max characteristic speed in one cell
  real(dp) function max_wave_speed_cell(U) result(s_max)
    real(dp), intent(in) :: U(nvar)
    real(dp) :: rho, vx, p, Bx, cs, vAx

    rho = U(1)
    vx  = U(2) / rho
    Bx  = U(5)
    ! recover gas pressure
    p = (gamma - 1._dp) * (U(8) - 0.5_dp*rho*vx*vx - 0.5_dp*Bx*Bx)

    cs   = sqrt(gamma * p / rho)
    vAx  = Bx / sqrt(rho)
    ! fast‐magnetosonic in 1D: sqrt(cs^2 + vAx^2)
    s_max = abs(vx) + sqrt(cs*cs + vAx*vAx)
  end function max_wave_speed_cell

  ! Loop over all cells to find the global max speed
  real(dp) function max_wave_speed(U, nx) result(max_s)
    real(dp), intent(in) :: U(nvar, nx)
    integer,  intent(in) :: nx
    integer :: i
    real(dp) :: s

    max_s = 0._dp
    do i = 1, nx
      s = max_wave_speed_cell(U(:,i))
      if (s > max_s) max_s = s
    end do
  end function max_wave_speed

end module mhd_wave
