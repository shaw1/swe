subroutine ftcs_flat(m0, n0, phi0, m, n, phi, f, dx, dt, K, length)
  ! Consider the shallow water equations with rotation and flat bottom
  ! topography
  !
  !               u_t + u * u_x + phi_x = f * v                   (1a)
  !                       v_t + u * v_x = -f * u                  (1b)
  !       phi_t + u * phi_x + phi * u_x = 0                       (1c)
  !
  ! where u and v are the eastward and northward wind of a fluid, and
  ! phi is the geopotential, i.e. phi(x, t) = g * h(x, t), where g is
  ! the acceleration due to gravity. These equations are the
  ! two-dimensional shallow water equations derived in [2] with the
  ! additional assumption that all functions are independent of y.
  ! Writing equations (1) in flux form, the equations become
  !
  ! m_t + [m * u + (phi ** 2) / 2]_x = f * n                    (2a)
  !                  n_t + (m * v)_x = -f * m                   (2b)
  !                      phi_t + m_x = 0                        (2c)
  !
  ! by letting m = phi * u and n = phi * v. Both [1] and [3] add
  ! diffusion terms to (2). The discrete model uses forward time and
  ! centered space finite differences.
  !
  ! Inputs:
  !    m0, n0, phi0, f, dx, dt, K, length
  !
  ! Outputs:
  !    m, n, phi
  !
  ! References:
  !     [1] Griffith, A. and Nichols, N. (2000). Adjoint methods in
  !         data assimilation for estimating model error. Flow,
  !         Turbulence, and Combustion, 65: 469--488.
  !
  !     [2] Pedlosky, J. (1987). Geophysical Fluid Dynamics. Second
  !         edition.
  !
  !     [3] Parrett, C. A., & Cullen, M. J. P. (1984). Simulation of
  !         hydraulic jumps in the presence of rotation and mountains.
  !         Quarterly Journal of the Royal Meteorological Society,
  !         147--165.
  !
  ! Author:             Jeremy Shaw
  ! Institution:        Portland State University
  ! Date Created:       Spring 2014
  ! Last Modified Date: 10 November 2016
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

  phi(1) = phi0(1) - coeff2 * (m0(2) - m0(length)) + coeff3 * &
       (phi0(length) - 2.d0 * phi0(1) + phi0(2))

  do j = 2, length - 1
     m(j) = m0(j) - coeff1 * ((u0(j + 1) + u0(j)) * (m0(j + 1) + &
          m0(j)) - (u0(j) + u0(j - 1)) * (m0(j) + m0(j - 1)) + &
          phi0(j + 1) ** 2 - phi0(j - 1) ** 2) + dt * f * n0(j) + &
          coeff3 * (m0(j - 1) - 2.d0 * m0(j) + m0(j + 1))

     n(j) = n0(j) - coeff1 * ((v0(j + 1) + v0(j)) * (m0(j + 1) + &
          m0(j)) - (v0(j) + v0(j - 1)) * (m0(j) + m0(j - 1))) - dt &
          * f * m0(j) + coeff3 * (n0(j - 1) - 2.d0 * n0(j) + &
          n0(j + 1))

     phi(j) = phi0(j) - coeff2 * (m0(j + 1) - m0(j - 1)) + coeff3 * &
            (phi0(j - 1) - 2.d0 * phi0(j) + phi0(j + 1))
  end do

  m(length) = m0(length) - coeff1 * ((u0(1) + u0(length)) * &
       (m0(1) + m0(length)) - (u0(length) + u0(length - 1)) * &
       (m0(length) + m0(length - 1)) + phi0(1) ** 2 - &
       phi0(length - 1) ** 2) + dt * f * n0(length) + coeff3 * &
       (m0(length - 1) - 2.d0 * m0(length) + m0(1))

  n(length) = n0(length) - coeff1 * ((v0(1) + v0(length)) * (m0(1) + &
       m0(length)) - (v0(length) + v0(length - 1)) * (m0(length) + &
       m0(length - 1))) - dt * f * m0(length) + coeff3 * &
       (n0(length - 1) - 2.d0 * n0(length) + n0(1))

  phi(length) = phi0(length) - coeff2 * (m0(1) - m0(length - 1)) + &
       coeff3 * (phi0(length - 1) - 2.d0 * phi0(length) + phi0(1))
end subroutine ftcs_flat
