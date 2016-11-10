subroutine ftcs_flat_adj(m0, m0b, n0, n0b, phi0, phi0b, mb, nb, &
& phib, f, dx, dt, k, length)
  ! This subroutine is the adjoint of the discrete shallow water model
  ! from ftcs_flat.f95.
  !
  ! Inputs:
  !     m0, n0, phi0, mb, nb, phib, f, k, dx, dt
  !
  !     The input adjoint variables are mb, nb, phib
  !
  ! Output:
  !     Adjoint variables m0b, n0b, phi0b
  !
  ! Author:             Jeremy Shaw
  ! Institution:        Portland State University
  ! Date Created:       6 November 2016
  ! Last Modified Date: 8 November 2016
  implicit none
  
  integer, intent(in) :: length
  double precision, dimension(length), intent(in) :: m0
  double precision, dimension(length), intent(out) :: m0b
  double precision, dimension(length), intent(in) :: n0
  double precision, dimension(length), intent(out) :: n0b
  double precision, dimension(length), intent(in) :: phi0
  double precision, dimension(length), intent(out) :: phi0b
  double precision, dimension(length), intent(in) :: mb
  double precision, dimension(length), intent(in) :: nb
  double precision, dimension(length), intent(in) :: phib
  double precision, intent(in) :: f
  double precision, intent(in) :: k
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
  coeff3 = k * dt / dx ** 2

  m0b(1) = mb(1) - coeff1 * (((m0(2) + m0(1)) / phi0(1) + &
&    (u0(2) + u0(1))) * (mb(1) - mb(2)) + ((m0(1) + &
&    m0(length)) / phi0(1) + (u0(1) + u0(length))) * (mb(length) - &
&    mb(1))) + coeff3 * (mb(length) - 2.d0 * mb(1) + mb(2)) + &
&    coeff1 * ((v0(2) + v0(1)) * (nb(2) - nb(1)) + &
&    (v0(1) + v0(length)) * (nb(1) - nb(length))) - dt * f * nb(1) &
&    + coeff2 * (phib(2) - phib(length))

  n0b(1) = nb(1) - coeff1 * ((m0(2) + m0(1)) * (nb(1) - &
&    nb(2)) + (m0(1) + m0(length)) * (nb(length) - nb(1))) &
&    / phi0(1) + coeff3 * (nb(length) - 2.d0 * nb(1) + nb(2)) &
&    + dt * f * mb(1)

  phi0b(1) = coeff1 * ((m0(2) + m0(1)) * (mb(2) - mb(1)) &
&    + (m0(1) + m0(length)) * (mb(1) - mb(length))) * (-u0(1) / &
&    phi0(1)) + coeff2 * phi0(1) * (mb(2) - mb(length)) - &
&    coeff1 * ((m0(2) + m0(1)) * (nb(2) - nb(1)) + (m0(1) &
&    + m0(length)) * (nb(1) - nb(length))) * (v0(1) / phi0(1)) + &
&    phib(1) + coeff3 * (phib(length) - 2.d0 * phib(1) + phib(2))

  do j = 2, length - 1
     m0b(j) = mb(j) - coeff1 * (((m0(j + 1) + m0(j)) / phi0(j) + &
&       (u0(j + 1) + u0(j))) * (mb(j) - mb(j + 1)) + ((m0(j) + &
&       m0(j - 1)) / phi0(j) + (u0(j) + u0(j - 1))) * (mb(j - 1) - &
&       mb(j))) + coeff3 * (mb(j - 1) - 2.d0 * mb(j) + mb(j + 1)) + &
&       coeff1 * ((v0(j + 1) + v0(j)) * (nb(j + 1) - nb(j)) + &
&       (v0(j) + v0(j - 1)) * (nb(j) - nb(j - 1))) - dt * f * nb(j) &
&       + coeff2 * (phib(j + 1) - phib(j - 1))

     n0b(j) = nb(j) - coeff1 * ((m0(j + 1) + m0(j)) * (nb(j) - &
&       nb(j + 1)) + (m0(j) + m0(j - 1)) * (nb(j - 1) - nb(j))) &
&       / phi0(j) + coeff3 * (nb(j - 1) - 2.d0 * nb(j) + nb(j + 1)) &
&       + dt * f * mb(j)

     phi0b(j) = phib(j) - coeff1 * ((m0(j + 1) + m0(j)) * (mb(j + 1) &
&       - mb(j)) + (m0(j) + m0(j - 1)) * (mb(j) - mb(j - 1))) * &
&       (u0(j) / phi0(j)) + coeff2 * phi0(j) * (mb(j + 1) - &
&       mb(j - 1)) - coeff1 * ((m0(j + 1) + m0(j)) * (nb(j + 1) - &
&       nb(j)) + (m0(j) + m0(j - 1)) * (nb(j) - nb(j - 1))) * (v0(j) &
&       / phi0(j)) + coeff3 * (phib(j - 1) - 2.d0 * phib(j) + &
&       phib(j + 1))
  end do

  m0b(length) = mb(length) - coeff1 * (((m0(1) + m0(length)) / &
&    phi0(length) + (u0(1) + u0(length))) * (mb(length) - mb(1)) + &
&    ((m0(length) + m0(length - 1)) / phi0(length) + (u0(length) + &
&    u0(length - 1))) * (mb(length - 1) - mb(length))) + coeff3 * &
&    (mb(length - 1) - 2.d0 * mb(length) + mb(1)) + coeff1 * ((v0(1) &
&    + v0(length)) * (nb(1) - nb(length)) + (v0(length) + &
&    v0(length - 1)) * (nb(length) - nb(length - 1))) - dt * f * &
&    nb(length) + coeff2 * (phib(1) - phib(length - 1))

  n0b(length) = nb(length) - coeff1 * ((m0(1) + m0(length)) * &
&    (nb(length) - nb(1)) + (m0(length) + m0(length - 1)) * &
&    (nb(length - 1) - nb(length))) / phi0(length) + coeff3 * &
&    (nb(length - 1) - 2.d0 * nb(length) + nb(1)) + dt * f * &
&    mb(length)

  phi0b(length) = coeff1 * ((m0(1) + m0(length)) * (mb(1) - &
&    mb(length)) + (m0(length) + m0(length - 1)) * (mb(length) - &
&    mb(length - 1))) * (-u0(length) / phi0(length)) + coeff2 * &
&    phi0(length) * (mb(1) - mb(length - 1)) - coeff1 * ((m0(1) + &
&    m0(length)) * (nb(1) - nb(length)) + (m0(length) + &
&    m0(length - 1)) * (nb(length) - nb(length - 1))) * (v0(length) &
&    / phi0(length)) + phib(length) + coeff3 * (phib(length - 1) - &
&    2.d0 * phib(length) + phib(1))
end subroutine ftcs_flat_adj
