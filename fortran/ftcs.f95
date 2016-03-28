subroutine ftcs(m0, n0, phi0, m, n, phi, H, f, g, dx, dt, K, length)
  ! The two-dimensional shallow water equations with bottom topography and
  ! rotation as derived in [2] are implemented with the additional assumption
  ! that the functions u, v, and phi are independent of y. These equations
  ! appear in [2], and are solved with a forward time and centered space
  ! finite difference scheme. Diffusion terms are included.
  !
  ! References:
  !   [1] Griffith, A. and Nichols, N. (2000). Adjoint methods in data
  !       assimilation for estimating model error. Flow, Turbulence, and
  !       Combustion, 65: 469--488.
  !
  !   [2] Pedlosky, J. (1987). Geophysical Fluid Dynamics. Second edition.
  !
  !   [3] Parrett, C. A., & Cullen, M. J. P. (1984). Simulation of hydraulic
  !       jumps in the presence of rotation and mountains. Quarterly Journal
  !       of the Royal Meteorological Society, 147--165.
  !
  ! Author:             Jeremy Shaw
  ! Institution:        Portland State University
  ! Date Created:       17 May 2015
  ! Last Modified Date: 12 October 2015
  implicit none

  integer, intent(in) :: length     ! Size of state vectors

  double precision, dimension(length), intent(in) :: m0     ! Initial m-state
  double precision, dimension(length), intent(in) :: n0     ! Initial n-state
  double precision, dimension(length), intent(in) :: phi0   ! Initial phi-state
  double precision, dimension(length), intent(out) :: m     ! Updated m-state
  double precision, dimension(length), intent(out) :: n     ! Updated n-state
  double precision, dimension(length), intent(out) :: phi   ! Updated phi-state
  double precision, dimension(length), intent(in) :: H      ! Bottom topography
  double precision, intent(in) :: g    ! Acceleration due to gravity
  double precision, intent(in) :: f    ! Coriolis coefficient
  double precision, intent(in) :: K    ! Diffusion coefficient
  double precision, intent(in) :: dx   ! Spatial step-size
  double precision, intent(in) :: dt   ! Temporal step-size

  ! Local variables
  integer :: j          ! Loop control variable
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
       phi0(length) ** 2 + g * ((phi0(2) + phi0(1)) * (H(2) - H(1)) &
       + (phi0(1) + phi0(length)) * (H(1) - H(length)))) + dt * f * n0(1) &
       + coeff3 * (m0(length) - 2.d0 * m0(1) + m0(2))

  n(1) = n0(1) - coeff1 * ((v0(2) + v0(1)) * (m0(2) + m0(1)) - &
       (v0(1) + v0(length)) * (m0(1) + m0(length))) - dt * f * m0(1) &
       + coeff3 * (n0(length) - 2.d0 * n0(1) + n0(2))

  phi(1) = phi0(1) - coeff2 * (m0(2) - m0(length)) + coeff3 * (phi0(length) - &
       2.d0 * phi0(1) + phi0(2))

  do j = 2, length - 1
     m(j) = m0(j) - coeff1 * ((u0(j + 1) + u0(j)) * (m0(j + 1) + m0(j)) - &
          (u0(j) + u0(j - 1)) * (m0(j) + m0(j - 1)) + phi0(j + 1) ** 2 - &
          phi0(j - 1) ** 2 + g * ((phi0(j + 1) + phi0(j)) * (H(j + 1) - H(j)) &
          + (phi0(j) + phi0(j - 1)) * (H(j) - H(j - 1)))) + dt * f * n0(j) &
          + coeff3 * (m0(j - 1) - 2.d0 * m0(j) + m0(j + 1))

     n(j) = n0(j) - coeff1 * ((v0(j + 1) + v0(j)) * (m0(j + 1) + m0(j)) - &
          (v0(j) + v0(j - 1)) * (m0(j) + m0(j - 1))) - dt * f * m0(j) &
          + coeff3 * (n0(j - 1) - 2.d0 * n0(j) + n0(j + 1))

     phi(j) = phi0(j) - coeff2 * (m0(j + 1) - m0(j - 1)) + coeff3 * &
            (phi0(j - 1) - 2.d0 * phi0(j) + phi0(j + 1))
  end do

  m(length) = m0(length) - coeff1 * ((u0(1) + u0(length)) * &
       (m0(1) + m0(length)) - (u0(length) + u0(length - 1)) * (m0(length) + &
       m0(length - 1)) + phi0(1) ** 2 - phi0(length - 1) ** 2 + g * &
       ((phi0(1) + phi0(length)) * (H(1) - H(length)) + &
       (phi0(length) + phi0(length - 1)) * (H(length) - H(length - 1)))) + &
       dt * f * n0(length) + coeff3 * (m0(length - 1) - 2.d0 * m0(length) &
       + m0(1))

  n(length) = n0(length) - coeff1 * ((v0(1) + v0(length)) * (m0(1) + &
       m0(length)) - (v0(length) + v0(length - 1)) * (m0(length) + &
       m0(length - 1))) - dt * f * m0(length) + coeff3 * (n0(length - 1) &
       - 2.d0 * n0(length) + n0(1))

  phi(length) = phi0(length) - coeff2 * (m0(1) - m0(length - 1)) + &
       coeff3 * (phi0(length - 1) - 2.d0 * phi0(length) + phi0(1))
end subroutine ftcs
