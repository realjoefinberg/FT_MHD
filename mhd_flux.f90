module mhd_flux
  use parameters, only: dp, gamma, nvar
  implicit none
contains

  subroutine compute_flux(Ui, Fi)
    real(dp), intent(in)  :: Ui(nvar)
    real(dp), intent(out) :: Fi(nvar)

    include 'flux_unpack_mhd1d_1.inc'
    include 'flux_compute_mhd1d_1.inc'
  end subroutine compute_flux

end module mhd_flux
