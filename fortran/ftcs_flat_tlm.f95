subroutine ftcs_flat_tlm(m0, m0d, n0, n0d, phi0, phi0d, md, nd, &
& phid, f, dx, dt, k, length)
  ! This subroutine is the tangent linear model of the discrete
  ! shallow water model from ftcs_flat.f95.
  !
  ! Inputs:
  !     m0, n0, phi0, m0d, n0d, phi0d, f, k, dx, dt
  !
  !     The tangent linear input variables are m0d, n0d, phi0d
  !
  ! Output:
  !     Tangent linear output variables md, nd, phid
  !
  ! Author:             Jeremy Shaw
  ! Institution:        Portland State University
  ! Date Created:       5 November 2016
  ! Last Modified Date: 8 November 2016  
  implicit none

  integer, intent(in) :: length

  double precision, dimension(length), intent(in) :: m0
  double precision, dimension(length), intent(in) :: m0d
  double precision, dimension(length), intent(in) :: n0
  double precision, dimension(length), intent(in) :: n0d
  double precision, dimension(length), intent(in) :: phi0
  double precision, dimension(length), intent(in) :: phi0d
  double precision, dimension(length), intent(out) :: md
  double precision, dimension(length), intent(out) :: nd
  double precision, dimension(length), intent(out) :: phid
  double precision, intent(in) :: f
  double precision, intent(in) :: K
  double precision, intent(in) :: dx
  double precision, intent(in) :: dt

  integer :: j
  double precision :: coeff1
  double precision :: coeff2
  double precision :: coeff3
  double precision, dimension(length) :: u0
  double precision, dimension(length) :: u0d
  double precision, dimension(length) :: v0
  double precision, dimension(length) :: v0d

  u0 = m0 / phi0
  v0 = n0 / phi0
  u0d = (m0d - phi0d * u0) / phi0
  v0d = (n0d - phi0d * v0) / phi0
  
  coeff1 = dt / (4.d0 * dx)
  coeff2 = dt / (2.d0 * dx)
  coeff3 = K * dt / (dx ** 2)

  md(1) = m0d(1) - coeff1 * ((u0d(2) + u0d(1)) * (m0(2) + m0(1)) + &
&   (u0(2) + u0(1)) * (m0d(2) + m0d(1)) - (u0d(1) + u0d(length)) * &
&   (m0(1) + m0(length)) - (u0(1) + u0(length)) * (m0d(1) + &
&   m0d(length)) + 2.d0 * (phi0(2) * phi0d(2) - phi0(length) * &
&   phi0d(length))) + dt * f * n0d(1) + coeff3 * (m0d(length) - 2.d0 &
&   * m0d(1) + m0d(2))

  nd(1) = n0d(1) - coeff1 * ((v0d(2) + v0d(1)) * (m0(2) + m0(1)) + &
&   (v0(2) + v0(1)) * (m0d(2) + m0d(1)) - (v0d(1) + v0d(length)) * &
&   (m0(1) + m0(length)) - (v0(1) + v0(length)) * (m0d(1) + &
&   m0d(length))) - dt * f * m0d(1) + coeff3 * (n0d(length) - 2.d0 * &
&   n0d(1) + n0d(2))

  phid(1) = phi0d(1)  - coeff2 * (m0d(2) - m0d(length)) + coeff3 * &
&   (phi0d(length) - 2.d0 * phi0d(1) + phi0d(2))

  do j = 2, length - 1
     md(j) = m0d(j) - coeff1 * ((u0d(j + 1) + u0d(j)) * (m0(j + 1) + &
&       m0(j)) + (u0(j + 1) + u0(j)) * (m0d(j + 1) + m0d(j)) - &
&       (u0d(j) + u0d(j - 1)) * (m0(j) + m0(j - 1)) - (u0(j) + &
&       u0(j - 1)) * (m0d(j) + m0d(j - 1)) + 2.d0 * (phi0(j + 1) * &
&       phi0d(j + 1) - phi0(j - 1) * phi0d(j - 1))) + dt * f * &
&       n0d(j) + coeff3 * (m0d(j - 1) - 2.d0 * m0d(j) + m0d(j + 1))

     nd(j) = n0d(j) - coeff1 * ((v0d(j + 1) + v0d(j)) * (m0(j + 1) + &
&       m0(j)) + (v0(j + 1) + v0(j)) * (m0d(j + 1) + m0d(j)) - &
&       (v0d(j) + v0d(j - 1)) * (m0(j) + m0(j - 1)) - (v0(j) + &
&       v0(j - 1)) * (m0d(j) + m0d(j - 1))) - dt * f * m0d(j) + &
&       coeff3 * (n0d(j - 1) - 2.d0 * n0d(j) + n0d(j + 1))

     phid(j) = phi0d(j) - coeff2 * (m0d(j + 1) - m0d(j - 1)) + &
&       coeff3 * (phi0d(j - 1) - 2.d0 * phi0d(j) + phi0d(j + 1))
  end do

  md(length) = m0d(length) - coeff1 * ((u0d(1) + u0d(length)) * &
&   (m0(1) + m0(length)) + (u0(1) + u0(length)) * (m0d(1) + &
&   m0d(length)) - (u0d(length) + u0d(length - 1)) * (m0(length) + &
&   m0(length - 1)) - (u0(length) + u0(length - 1)) * (m0d(length) + &
&   m0d(length - 1)) + 2.d0 * (phi0(1) * phi0d(1) - phi0(length - 1) &
&   * phi0d(length - 1))) + dt * f * n0d(length) + coeff3 * &
&   (m0d(length - 1) - 2.d0 * m0d(length) + m0d(1))

  nd(length) = n0d(length) - coeff1 * ((v0d(1) + v0d(length)) * &
&   (m0(1) + m0(length)) + (v0(1) + v0(length)) * (m0d(1) + &
&   m0d(length)) - (v0d(length) + v0d(length - 1)) * (m0(length) + &
&   m0(length - 1)) - (v0(length) + v0(length - 1)) * (m0d(length) + &
&   m0d(length - 1))) - dt * f * m0d(length) + coeff3 * &
&   (n0d(length - 1) - 2.d0 * n0d(length) + n0d(1))

  phid(length) = phi0d(length)  - coeff2 * (m0d(1) - m0d(length - 1))&
&   + coeff3 * (phi0d(length - 1) - 2.d0 * phi0d(length) + phi0d(1))
end subroutine ftcs_flat_tlm
