module mhd_flux
  use parameters, only: dp, gamma, nvar
  implicit none
contains

  subroutine compute_flux(Ui, Fi)
    !-------------------------------------------------------------------
    ! Compute the 1D ideal‑MHD flux vector Fi(:) from state Ui(:)
    !-------------------------------------------------------------------
    real(dp), intent(in)  :: Ui(nvar)
    real(dp), intent(out) :: Fi(nvar)

    !— pull in the unpacking boilerplate —!
    include 'flux_unpack_mhd1d_1.inc'

    !— pull in the flux‑assembly boilerplate —!
    include 'flux_compute_mhd1d_1.inc'

  end subroutine compute_flux

end module mhd_flux
