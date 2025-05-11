module mhd_types
  use parameters, only: dp, nvar
  implicit none

  !------------------------------------------------------------------
  ! Derived type to hold the conserved variables for one grid cell:
  !   U(1) = ρ        (mass density)
  !   U(2) = ρ v_x    (momentum x)
  !   U(3) = ρ v_y    (momentum y)
  !   U(4) = ρ v_z    (momentum z)
  !   U(5) = B_x      (magnetic field x)
  !   U(6) = B_y      (magnetic field y)
  !   U(7) = B_z      (magnetic field z)
  !   U(8) = E        (total energy)
  !------------------------------------------------------------------
  type :: state
    real(dp), dimension(nvar) :: U
  end type state

end module mhd_types
