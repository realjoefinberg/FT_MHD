module mhd_init
  use parameters, only: dp, gamma, nvar
  implicit none
contains

  subroutine set_brio_wu(U, rho, vx, vy, vz, Bx, By, Bz, p)
    real(dp), intent(out) :: U(nvar)
    ! primitive inputs:
    real(dp), intent(in ) :: rho, vx, vy, vz, Bx, By, Bz, p
    real(dp) :: kinetic, magnetic

    ! pack conserved:
    U(1) = rho
    U(2) = rho * vx
    U(3) = rho * vy
    U(4) = rho * vz

    U(5) = Bx
    U(6) = By
    U(7) = Bz

    kinetic  = 0.5_dp * rho * (vx*vx + vy*vy + vz*vz)
    magnetic = 0.5_dp * (Bx*Bx + By*By + Bz*Bz)
    U(8) = p/(gamma-1._dp) + kinetic + magnetic
  end subroutine set_brio_wu

end module mhd_init
