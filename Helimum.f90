PROGRAM QMC
!================================================================================
!
!This program was made by Professor: Michel CAFFAREL and also share by him to
!his MSc Students.  Some changes was made by me. But is still the main idea. If
!you want to know more about the calculations read the file attached in the main
!branch.
!
!VCastor
!
!================================================================================

IMPLICIT NONE

INTEGER m_max
PARAMETER(m_max=1000)
INTEGER :: i_path, i, j, l, n_paths
REAL :: twopi, tau, prod, elocal, alpha, gam, xold(3,2), xnew(3,2), b(3,2), gauss
REAL, DIMENSION(m_max) :: h, s, w, eloc

twopi = 2.*ACOS(-1.)

PRINT*,'N?'
READ(*,*) n_paths
tau = 0.0075
alpha = 0.35
gam = 2.

h(:) = 0.
s(:) = 0.

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
                    xnew(l,i) = xold(l,i) + b(l,i)*tau + SQRT(tau)*gauss()
                ENDDO
            ENDDO
        eloc(j) = elocal(xnew,gam,alpha)
        w(j) = EXP(-tau*eloc(j))
        xold = xnew
    ENDDO

    prod = 1.
    DO j = 1, m_max
        prod = prod*w(j)
        h(j) = h(j) + prod*eloc(j)
        s(j) = s(j) + prod
    ENDDO
ENDDO

DO j = 1, m_max
    WRITE(*,*) 'time', j*tau, 'Energy = ', h(j)/s(j)
!    WRITE(*,*) h(1), s(1)
!    WRITE(*,*) prod
ENDDO

END

!====================================================================
!     SUBROUTINES AND FUNCTIONS
!====================================================================

SUBROUTINE DRIFT(x, b, gam, alpha)

IMPLICIT NONE

INTEGER i, l
REAL x(3,2), b(3,2), gam, alpha, r(2), r12, u

DO i = 1, 2
    r(i) = SQRT(x(1,i)**2 + x(2,i)**2 + x(3,i)**2)
    DO l = 1, 3
        b(l,i) = -gam*x(l,i)/r(i)
    ENDDO
ENDDO

r12 = SQRT((x(1,1) - x(1,2))**2 + (x(2,1) - x(2,2))**2 + (x(3,1) - x(3,2))**2)
u = alpha/(1. + alpha*r12)

DO l = 1, 3
    b(l,1) = b(l,1) + u*(x(l,1) - x(l,2))/r12
    b(l,2) = b(l,2) + u*(x(l,2) - x(l,1))/r12
ENDDO

END

!================================================================================

FUNCTION elocal(x, gam, alpha)

IMPLICIT NONE

REAL elocal, x(3,2), alpha, gam, r(2), r12, u, prod
INTEGER l

r(1)   = SQRT(x(1,1)**2 + x(2,1)**2 + x(3,1)**2)
r(2)   = SQRT(x(1,2)**2 + x(2,2)**2 + x(3,2)**2)
r12    = SQRT((x(1,1) - x(1,2) )**2 + (x(2,1)-x(2,2))**2 + (x(3,1)-x(3,2))**2)
elocal = (gam-2.)*(1./r(1) + 1./r(2)) + 1./r12 - gam**2
u      = alpha/(1. + alpha*r12)

prod = 0.

DO l = 1, 3
    prod = prod + (x(l,1) - x(l,2))/r12*(x(l,1)/r(1) - x(l,2)/r(2))
ENDDO

elocal = elocal - 2.*u/r12 + gam*u*prod

END

!================================================================================
!    BOX MULLER GAUSSIAN GENERATOR

FUNCTION gauss()

u = ALOG(RAND())
v = RAND()

gauss = SQRT(-2.*u)*COS(twopi*v)

END
