!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.10 (r5717) - 30 Jul 2015 16:03
!
!  Differentiation of ftcs in forward (tangent) mode:
!   variations   of useful results: m n phi
!   with respect to varying inputs: m0 phi0 n0
!   RW status of diff variables: m:out n:out m0:in phi0:in n0:in
!                phi:out
SUBROUTINE FTCS_TOPOG_TLM(m0, m0d, n0, n0d, phi0, phi0d, md, nd, phid&
& , h, f, g, dx, dt, k, length)
  ! This subroutine is the tangent linear model of the discrete
  ! shallow water model from ftcs_topog.f95, as generated by Tapenade
  ! 3.10
  !
  ! Inputs:
  !     m0, n0, phi0, m0d, n0d, phi0d, f, k, dx, dt
  !
  !     The tangent linear input variables are m0d, n0d, phi0d
  !
  ! Output:
  !     Tangent linear output variables md, nd, phid
  !
  ! Author:             Tapenade 3.10
  ! Edited By:          Jeremy Shaw
  ! Institution:        Portland State University
  ! Date Created:       Spring 2016
  ! Last Modified Date: 10 November 2016
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: length
  DOUBLE PRECISION, DIMENSION(length), INTENT(IN) :: m0
  DOUBLE PRECISION, DIMENSION(length), INTENT(IN) :: m0d
  DOUBLE PRECISION, DIMENSION(length), INTENT(IN) :: n0
  DOUBLE PRECISION, DIMENSION(length), INTENT(IN) :: n0d
  DOUBLE PRECISION, DIMENSION(length), INTENT(IN) :: phi0
  DOUBLE PRECISION, DIMENSION(length), INTENT(IN) :: phi0d
  DOUBLE PRECISION, DIMENSION(length), INTENT(OUT) :: md
  DOUBLE PRECISION, DIMENSION(length), INTENT(OUT) :: nd
  DOUBLE PRECISION, DIMENSION(length), INTENT(OUT) :: phid
  DOUBLE PRECISION, DIMENSION(length), INTENT(IN) :: h
  DOUBLE PRECISION, INTENT(IN) :: g
  DOUBLE PRECISION, INTENT(IN) :: f
  DOUBLE PRECISION, INTENT(IN) :: k
  DOUBLE PRECISION, INTENT(IN) :: dx
  DOUBLE PRECISION, INTENT(IN) :: dt
  INTEGER :: j
  DOUBLE PRECISION :: coeff1
  DOUBLE PRECISION :: coeff2
  DOUBLE PRECISION :: coeff3
  DOUBLE PRECISION, DIMENSION(length) :: u0
  DOUBLE PRECISION, DIMENSION(length) :: u0d
  DOUBLE PRECISION, DIMENSION(length) :: v0
  DOUBLE PRECISION, DIMENSION(length) :: v0d
  u0d = (m0d*phi0-m0*phi0d)/phi0**2
  u0 = m0/phi0
  v0d = (n0d*phi0-n0*phi0d)/phi0**2
  v0 = n0/phi0
  coeff1 = dt/(4.d0*dx)
  coeff2 = dt/(2.d0*dx)
  coeff3 = k*dt/dx**2
  md = 0.D0
  md(1) = m0d(1) - coeff1*((u0d(2)+u0d(1))*(m0(2)+m0(1))+(u0(2)+u0(1))*(&
&   m0d(2)+m0d(1))-(u0d(1)+u0d(length))*(m0(1)+m0(length))-(u0(1)+u0(&
&   length))*(m0d(1)+m0d(length))+2*phi0(2)*phi0d(2)-2*phi0(length)*&
&   phi0d(length)+g*((h(2)-h(1))*(phi0d(2)+phi0d(1))+(h(1)-h(length))*(&
&   phi0d(1)+phi0d(length)))) + dt*f*n0d(1) + coeff3*(m0d(length)-2.d0*&
&   m0d(1)+m0d(2))

  nd = 0.D0
  nd(1) = n0d(1) - coeff1*((v0d(2)+v0d(1))*(m0(2)+m0(1))+(v0(2)+v0(1))*(&
&   m0d(2)+m0d(1))-(v0d(1)+v0d(length))*(m0(1)+m0(length))-(v0(1)+v0(&
&   length))*(m0d(1)+m0d(length))) - dt*f*m0d(1) + coeff3*(n0d(length)-&
&   2.d0*n0d(1)+n0d(2))

  phid = 0.D0
  phid(1) = phi0d(1) - coeff2*(m0d(2)-m0d(length)) + coeff3*(phi0d(&
&   length)-2.d0*phi0d(1)+phi0d(2))

  DO j=2,length-1
    md(j) = m0d(j) - coeff1*((u0d(j+1)+u0d(j))*(m0(j+1)+m0(j))+(u0(j+1)+&
&     u0(j))*(m0d(j+1)+m0d(j))-(u0d(j)+u0d(j-1))*(m0(j)+m0(j-1))-(u0(j)+&
&     u0(j-1))*(m0d(j)+m0d(j-1))+2*phi0(j+1)*phi0d(j+1)-2*phi0(j-1)*&
&     phi0d(j-1)+g*((h(j+1)-h(j))*(phi0d(j+1)+phi0d(j))+(h(j)-h(j-1))*(&
&     phi0d(j)+phi0d(j-1)))) + dt*f*n0d(j) + coeff3*(m0d(j-1)-2.d0*m0d(j&
&     )+m0d(j+1))

    nd(j) = n0d(j) - coeff1*((v0d(j+1)+v0d(j))*(m0(j+1)+m0(j))+(v0(j+1)+&
&     v0(j))*(m0d(j+1)+m0d(j))-(v0d(j)+v0d(j-1))*(m0(j)+m0(j-1))-(v0(j)+&
&     v0(j-1))*(m0d(j)+m0d(j-1))) - dt*f*m0d(j) + coeff3*(n0d(j-1)-2.d0*&
&     n0d(j)+n0d(j+1))

    phid(j) = phi0d(j) - coeff2*(m0d(j+1)-m0d(j-1)) + coeff3*(phi0d(j-1)&
&     -2.d0*phi0d(j)+phi0d(j+1))

  END DO
  md(length) = m0d(length) - coeff1*((u0d(1)+u0d(length))*(m0(1)+m0(&
&   length))+(u0(1)+u0(length))*(m0d(1)+m0d(length))-(u0d(length)+u0d(&
&   length-1))*(m0(length)+m0(length-1))-(u0(length)+u0(length-1))*(m0d(&
&   length)+m0d(length-1))+2*phi0(1)*phi0d(1)-2*phi0(length-1)*phi0d(&
&   length-1)+g*((h(1)-h(length))*(phi0d(1)+phi0d(length))+(h(length)-h(&
&   length-1))*(phi0d(length)+phi0d(length-1)))) + dt*f*n0d(length) + &
&   coeff3*(m0d(length-1)-2.d0*m0d(length)+m0d(1))

  nd(length) = n0d(length) - coeff1*((v0d(1)+v0d(length))*(m0(1)+m0(&
&   length))+(v0(1)+v0(length))*(m0d(1)+m0d(length))-(v0d(length)+v0d(&
&   length-1))*(m0(length)+m0(length-1))-(v0(length)+v0(length-1))*(m0d(&
&   length)+m0d(length-1))) - dt*f*m0d(length) + coeff3*(n0d(length-1)-&
&   2.d0*n0d(length)+n0d(1))

  phid(length) = phi0d(length) - coeff2*(m0d(1)-m0d(length-1)) + coeff3*&
&   (phi0d(length-1)-2.d0*phi0d(length)+phi0d(1))
END SUBROUTINE FTCS_TOPOG_TLM
