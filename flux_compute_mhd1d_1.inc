! flux_compute_mhd1d_1.inc
! Assemble the 1D ideal‑MHD flux components into Fi(1:8)
  Fi(1) = rho*vx
  Fi(2) = rho*vx*vx + p + 0.5_dp*B2 - Bx*Bx
  Fi(3) = rho*vx*vy    - Bx*By
  Fi(4) = rho*vx*vz    - Bx*Bz
  Fi(5) = 0._dp
  Fi(6) = vx*By - vy*Bx
  Fi(7) = vx*Bz - vz*Bx
  Fi(8) = (E + p + 0.5_dp*B2)*vx - vdotB*Bx
