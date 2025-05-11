module mhd_init
  use parameters, only: dp, gamma, nvar
  implicit none
contains

  subroutine set_brio_wu(U, rho, vx, vy, vz, Bx, By, Bz, p)
    !-------------------------------------------------------------------
    ! Pack a Brio–Wu state into the conserved vector U(:)
    ! Inputs:
    !   rho, vx,vy,vz : mass density and velocity components
    !   Bx,By,Bz      : magnetic field components
    !   p             : gas pressure
    ! Outputs:
    !   U(1:nvar)     : conserved variables [ρ,ρv, B, E]
    !-------------------------------------------------------------------
    real(dp), intent(out) :: U(nvar)
    real(dp), intent(in)  :: rho, vx, vy, vz, Bx, By, Bz, p
    real(dp)             :: Ekin, Emag, Eint

    Ekin = 0.5_dp * rho * (vx*vx + vy*vy + vz*vz)
    Emag = 0.5_dp * (Bx*Bx + By*By + Bz*Bz)
    Eint = p / (gamma - 1._dp)

    U(1) = rho
    U(2) = rho * vx
    U(3) = rho * vy
    U(4) = rho * vz
    U(5) = Bx
    U(6) = By
    U(7) = Bz
    U(8) = Ekin + Emag + Eint
  end subroutine set_brio_wu

end module mhd_init
