PROGRAM QMC
!***********************************************************************
! This program was made by Professor:  Michel CAFFAREL and also share
! by him to his MSc Students. Some changes was done me. But still be
! the same idea. If you want to know more about the calculations is
! attache a pdf explainig the theory.
!
!		VCastor, 2021
!***********************************************************************

!-----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE

PARAMETER(m_max=1000)
INTEGER                         :: m_max
INTEGER                         :: i, j, l, n_paths, i_path
REAL (KIND=8)                   :: twopi, tau, prod, elocal, alpha, gam, gauss
REAL (KIND=8), DIMENSION(m_max) :: h, s, w, eloc
REAL (KIND=8), DIMENSION(3,2)   :: xold, xnew, b
!-----------------------------------------------------------------------------------------------------------------------------------

twopi = 2.d0*DACOS(1.d0)

PRINT*,'N?'
READ(*,*) n_paths

tau   = 0.0075d0
alpha = 0.35d0
gam   = 2.d0

h(:) = 0.d0; s(:) = 0.d0

DO i_path = 1, n_paths

    DO i = 1, 2
      DO l = 1, 3
        xold(l,i) = 1. + RAND()
      ENDDO
    ENDDO

    DO j = 1, m_max
      CALL drift(xold, b, gam, alpha)
        DO i = 1, 2
          DO l = 1, 3
            xnew(l,i) = xold(l,i) + b(l,i)*tau + DSQRT(tau)*gauss()
          ENDDO
        ENDDO
      eloc(j) = elocal(xnew,gam,alpha)
      w(j)    = EXP(-tau*eloc(j))
      xold    = xnew
    ENDDO

    prod = 1.d0
    DO j = 1, m_max
        prod = prod*w(j)
        h(j) = h(j) + prod*eloc(j)
        s(j) = s(j) + prod
    ENDDO

ENDDO

OPEN(11,FILE='He_Montecarlo.out')
DO j = 1, m_max
    WRITE(11,*) 'time', j*tau, 'Energy = ', h(j)/s(j)
ENDDO
WRITE(*,*) '***************************************************************'
WRITE(*,*) 'THE JOB WAS DONE AND ITS WRITTEN IN THE FILE: He_Montecarlo.out'
WRITE(*,*) '***************************************************************'
WRITE(*,*) 'Sale, bye /(o.o)/'

ENDPROGRAM QMC

!***********************************************************************
!                       SUBROUTINES AND FUNCTIONS
!***********************************************************************

SUBROUTINE DRIFT(x, b, gam, alpha)

IMPLICIT NONE

INTEGER                       :: i, l
REAL (KIND=8)                 :: gam, alpha, r12, u
REAL (KIND=8), DIMENSION(2)   :: r
REAL (KIND=8), DIMENSION(3,2) :: x, b

DO i = 1, 2
    r(i) = DSQRT(x(1,i)*x(1,i) + x(2,i)*x(2,i) + x(3,i)*x(3,i))
    DO l = 1, 3
        b(l,i) = -gam*x(l,i)/r(i)
    ENDDO
ENDDO

r12 = DSQRT((x(1,1) - x(1,2))**2 + (x(2,1) - x(2,2))**2 + (x(3,1) - x(3,2))**2)
u   = alpha/(1.d0 + alpha*r12)

DO l = 1, 3
    b(l,1) = b(l,1) + u*(x(l,1) - x(l,2))/r12
    b(l,2) = b(l,2) + u*(x(l,2) - x(l,1))/r12
ENDDO

ENDSUBROUTINE DRIFT

!***********************************************************************

FUNCTION elocal(x, gam, alpha)

IMPLICIT NONE
INTEGER                       :: l
REAL (KIND=8)                 :: elocal, alpha, gam, r12, u, prod
REAL (KIND=8), DIMENSION(2)   :: r
REAL (KIND=8), DIMENSION(3,2) :: x

r(1)   = SQRT(x(1,1)**2 + x(2,1)**2 + x(3,1)**2)
r(2)   = SQRT(x(1,2)**2 + x(2,2)**2 + x(3,2)**2)
r12    = SQRT((x(1,1) - x(1,2) )**2 + (x(2,1)-x(2,2))**2 + (x(3,1)-x(3,2))**2)
elocal = (gam-2.d0)*(1.d0/r(1) + 1.d0/r(2)) + 1.d0/r12 - gam*gam
u      = alpha/(1.d0 + alpha*r12)

prod = 0.d0

DO l = 1, 3
    prod = prod + (x(l,1) - x(l,2))/r12*(x(l,1)/r(1) - x(l,2)/r(2))
ENDDO

elocal = elocal - 2.d0*u/r12 + gam*u*prod

ENDFUNCTION

!***********************************************************************
!    BOX MULLER GAUSSIAN GENERATOR

FUNCTION gauss()

u = ALOG(RAND())
v = RAND()

gauss = DSQRT(-2.d0*u)*DCOS(twopi*v)

ENDFUNCTION

!***********************************************************************
!    This program does not have easter eggs on the source code, sorry :c
