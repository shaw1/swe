subroutine ftcs_flat(m0, n0, phi0, m, n, phi, f, dx, dt, K, length)
  implicit none

  integer, intent(in) :: length

  double precision, dimension(length), intent(in) :: m0
  double precision, dimension(length), intent(in) :: n0
  double precision, dimension(length), intent(in) :: phi0
  double precision, dimension(length), intent(out) :: m
  double precision, dimension(length), intent(out) :: n
  double precision, dimension(length), intent(out) :: phi
  double precision, intent(in) :: f
  double precision, intent(in) :: K
  double precision, intent(in) :: dx
  double precision, intent(in) :: dt

  integer :: j
  double precision :: coeff1
  double precision :: coeff2
  double precision :: coeff3
  double precision, dimension(length) :: u0
  double precision, dimension(length) :: v0

  u0 = m0 / phi0
  v0 = n0 / phi0
  
  coeff1 = dt / (4.d0 * dx)
  coeff2 = dt / (2.d0 * dx)
  coeff3 = K * dt / (dx ** 2)

  m(1) = m0(1) - coeff1 * ((u0(2) + u0(1)) * (m0(2) + m0(1)) - &
       (u0(1) + u0(length)) * (m0(1) + m0(length)) + phi0(2) ** 2 - &
       phi0(length) ** 2) + dt * f * n0(1) + coeff3 * (m0(length) - &
       2.d0 * m0(1) + m0(2))

  n(1) = n0(1) - coeff1 * ((v0(2) + v0(1)) * (m0(2) + m0(1)) - &
       (v0(1) + v0(length)) * (m0(1) + m0(length))) - dt * f * m0(1) &
       + coeff3 * (n0(length) - 2.d0 * n0(1) + n0(2))

  phi(1) = phi0(1) - coeff2 * (m0(2) - m0(length)) + coeff3 * (phi0(length) - &
       2.d0 * phi0(1) + phi0(2))

  do j = 2, length - 1
     m(j) = m0(j) - coeff1 * ((u0(j + 1) + u0(j)) * (m0(j + 1) + m0(j)) - &
          (u0(j) + u0(j - 1)) * (m0(j) + m0(j - 1)) + phi0(j + 1) ** 2 - &
          phi0(j - 1) ** 2) + dt * f * n0(j) + coeff3 * (m0(j - 1) - &
          2.d0 * m0(j) + m0(j + 1))

     n(j) = n0(j) - coeff1 * ((v0(j + 1) + v0(j)) * (m0(j + 1) + m0(j)) - &
          (v0(j) + v0(j - 1)) * (m0(j) + m0(j - 1))) - dt * f * m0(j) &
          + coeff3 * (n0(j - 1) - 2.d0 * n0(j) + n0(j + 1))

     phi(j) = phi0(j) - coeff2 * (m0(j + 1) - m0(j - 1)) + coeff3 * &
            (phi0(j - 1) - 2.d0 * phi0(j) + phi0(j + 1))
  end do

  m(length) = m0(length) - coeff1 * ((u0(1) + u0(length)) * &
       (m0(1) + m0(length)) - (u0(length) + u0(length - 1)) * (m0(length) + &
       m0(length - 1)) + phi0(1) ** 2 - phi0(length - 1) ** 2) + &
       dt * f * n0(length) + coeff3 * (m0(length - 1) - 2.d0 * m0(length) &
       + m0(1))

  n(length) = n0(length) - coeff1 * ((v0(1) + v0(length)) * (m0(1) + &
       m0(length)) - (v0(length) + v0(length - 1)) * (m0(length) + &
       m0(length - 1))) - dt * f * m0(length) + coeff3 * (n0(length - 1) &
       - 2.d0 * n0(length) + n0(1))

  phi(length) = phi0(length) - coeff2 * (m0(1) - m0(length - 1)) + &
       coeff3 * (phi0(length - 1) - 2.d0 * phi0(length) + phi0(1))
end subroutine ftcs_flat
