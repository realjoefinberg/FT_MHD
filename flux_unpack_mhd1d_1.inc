! flux_unpack_mhd1d_1.inc
! Unpack conserved MHD state Ui(1:nvar) into physical variables

  real(dp) :: rho, vx, vy, vz, Bx, By, Bz, E
  real(dp) :: B2, p, vdotB

  rho   = Ui(1)
  vx    = Ui(2) / rho
  vy    = Ui(3) / rho
  vz    = Ui(4) / rho
  Bx    = Ui(5)
  By    = Ui(6)
  Bz    = Ui(7)
  E     = Ui(8)

  B2    = Bx**2 + By**2 + Bz**2
  p     = (gamma - 1._dp) * (E - 0.5_dp*rho*(vx**2 + vy**2 + vz**2) \
                             - 0.5_dp*B2)
  vdotB = vx*Bx + vy*By + vz*Bz
