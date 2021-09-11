      program pson_mtg

      INTEGER N, MAX_IT, ITER_TOL
      REAL, ALLOCATABLE :: X(:),Y(:),Z(:)
      REAL, ALLOCATABLE :: S(:,:,:),  S2(:,:,:)  ,S3(:,:,:)
      REAL, ALLOCATABLE :: P(:,:,:),  P2(:,:,:)  ,P3(:,:,:)
      REAL, ALLOCATABLE :: RES(:,:,:),RES2(:,:,:),RES3(:,:,:)
      COMMON DL
      REAL, EXTERNAL :: COMPUTE_RES

      PI = ATAN(1.)*4.
      DL = PI*2.
      OPEN(1, FILE = 'poisson.in')
      READ(1,*) N        ! finest grid
      READ(1,*) TOL      ! residual at convergence
      READ(1,*) MAX_IT   ! the maxmun iteration number
      READ(1,*) 
      READ(1,*) ITNUM    ! the number of iteration on each grid
      CLOSE(1)

      WRITE(6,*) "N       =", N
      WRITE(6,*) "TOL     =", TOL
      WRITE(6,*) "MAX_IT  =", MAX_IT
      WRITE(6,*) "ITNUM   =", ITNUM

      ALLOCATE (X(N+1),Y(N+1),Z(N+1))
      ALLOCATE (S(N+1,N+1,N+1),P(N+1,N+1,N+1))
      ALLOCATE (S2(N/2+1,N/2+1,N/2+1),S3(N/4+1,N/4+1,N/4+1))
      ALLOCATE (P2(N/2+1,N/2+1,N/2+1),P3(N/4+1,N/4+1,N/4+1))
      ALLOCATE (RES(N+1,N+1,N+1),RES2(N/2+1,N/2+1,N/2+1))
      ALLOCATE (RES3(N/4+1,N/4+1,N/4+1))

      DS = DL/FLOAT(N)
      DO I=1,N+1
         X(I) = (I-1)*DS
         Y(I) = X(I)
         Z(I) = X(I)
      END DO
      P = 0.
      RES = 0.
      RES2 =0.
      RES3 =0.

      RES_MAX = 1.              !residual
      ITER_TOL = 0              !iteration on the whole circle

      DO  I = 1,N+1
         DO J = 1,N+1
            DO L = 1,N+1
               S(I,J,L) = 0.
               DO K = 1,4
                  S(I,J,L)=S(I,J,L)+sin(X(I)*K)*sin(Y(J)*K)*sin(Z(L)*K)
               END DO
            END DO
         END DO
      END DO

      DO WHILE ((RES_MAX.GT.TOL) .AND. (ITER_TOL.LT.MAX_IT))
c- iterate on N grid
         CALL ADI(S,ITNUM,N,P,RES)
         RES_MAX = COMPUTE_RES(RES,N)
c-go down N/2 grid 
         CALL INJECT(RES,N,S2)
         P2 = 0.
         CALL ADI(S2,ITNUM,N/2,P2,RES2)
c-go down N/4 grid 
         CALL INJECT(RES2,N/2,S3)
         P3 = 0.
         CALL ADI(S3,ITNUM,N/4,P3,RES3)

c-go up to N/2 grid 
         CALL INTERP(P3,N/2,P2)
         CALL ADI(S2,ITNUM,N/2,P2,RES2)
c-go up to N grid 
         CALL INTERP(P2,N,P)
         call check_err(S,N,P,RES)

      ITER_TOL = ITER_TOL +1
      WRITE(6,20) ITER_TOL,RES_MAX
      END DO
 20   FORMAT(1x,4X,'loop iteration =',I4,5X,'residual = ',F10.6)

c output the results
      OPEN (1,FILE = 'psonout.m',STATUS = 'OLD')
      DO K = 1,7
      WRITE(1,49) K
         DO J = 1,N+1
            WRITE(1,50) (P(K*N/8+1,J,L), L = 1,N+1) 
         END DO
      WRITE(1,*) '];'
      END DO
      CLOSE(1)
 49   FORMAT(1x, 'p',I1,' =[')
 50   FORMAT(1x,100(F10.5,1x))

      STOP
      END
c------||--------||---------||---------||----------||----------
        SUBROUTINE INJECT(P1,N,P2)
c this routine injects P1 of (N+1)^3 onto P2 of (N/2+1)^3
      INTEGER N
      REAL P1(N+1,N+1,N+1),P2(N/2+1,N/2+1,N/2+1)
         DO  I = 1,N/2+1
            DO J = 1,N/2+1
               DO L = 1,N/2+1
                  P2(I,J,L)=P1(2*I-1,2*J-1,2*L-1)
               END DO
            END DO
         END DO
      END SUBROUTINE INJECT 
c------||--------||---------||---------||----------||----------
        SUBROUTINE INTERP(P1,N,P2)
c this routine interpolates P1 of (N/2+1)^3 onto P2 of (N+1)^3
c note that P2 has elements on each grid point, so the result is
c P2 plus the interpolation
      INTEGER N
      REAL P1(N/2+1,N/2+1,N/2+1),P2(N+1,N+1,N+1),TMP(N+1,N+1,N+1)

      DO  I = 1,N/2+1
         DO J = 1,N/2+1
            DO L = 1,N/2+1
               TMP(2*I-1,2*J-1,2*L-1) = P1(I,J,L)
            END DO
         END DO
      END DO
      DO  I = 1,N+1,2
         DO J = 1,N+1,2
            DO L = 2,N,2
               TMP(I,J,L)=(TMP(I,J,L-1) + TMP(I,J,L+1))/2.                            
            END DO
         END DO
      END DO
      DO  I = 1,N+1,2
         DO L = 1,N+1
            DO J = 2,N,2
               TMP(I,J,L)=(TMP(I,J-1,L) + TMP(I,J+1,L))/2. 
            END DO
         END DO
      END DO
      DO J = 1,N+1
         DO L = 1,N+1
            DO I = 2,N,2
               TMP(I,J,L)=(TMP(I-1,J,L) + TMP(I+1,J,L))/2. 
            END DO
         END DO
      END DO
      P2 = P2 + TMP
      END SUBROUTINE INTERP

c------||--------||---------||---------||----------||----------
      FUNCTION COMPUTE_RES(RES,N)
c--compute the 1 norm residual
      INTEGER N
      REAL RES(N+1,N+1,N+1),COMPUTE_RES
      COMPUTE_RES= 0.
      DO  I = 1,N+1
         DO L = 1,N+1
            DO J = 1,N+1
               TMP = ABS(RES(I,J,L))
               IF (TMP.GT.COMPUTE_RES) COMPUTE_RES = TMP
            END DO
         END DO
      END DO
      END FUNCTION COMPUTE_RES
c------||--------||---------||---------||----------||----------
        SUBROUTINE ADI(S,ITNUM,N,DP,RES)
C this subroutine iterate the linear system (L)*dp = S, where (L) is
c the operator--poisson 2nd order discretization with N grid, using the ADI method.
c the routine iterate ITNUM times and return the residual RES at each point
c dp is the initial guess and also the return approx. solution
      INTEGER N,ITNUM
      REAL S(N+1,N+1,N+1),DP(N+1,N+1,N+1),RES(N+1,N+1,N+1)
      REAL d1(N-1),d2(N-1),d3(N-1),b(N-1)
      COMMON DL

      DS  = DL/float(N)
      DS2 = DS**2

      ITER = 0
      DO I = 1,N-1
         d1(I) = 1./DS2
         d2(I) = -2./DS2*3.
         d3(I) = 1./DS2         
      END DO

      WRITE(6,35) N,N,N
 35   FORMAT(5X,'now I am iterating on ',I2,'X',I2,'X',I2,' grid....')
      DO WHILE (ITER.LT.ITNUM)
c--sweep in z direction
      DO  I = 2,N
         DO J = 2,N
            DO L = 2,N
               b(L-1)=S(I,J,L)-(DP(I+1,J,L)+DP(I-1,J,L))/DS2
               b(L-1)=b(L-1) - (DP(I,J+1,L)+DP(I,J-1,L))/DS2
            END DO
            CALL THOMAS(d1,d2,d3,b,N-1)
            DP(I,J,2:N) = b(1:N-1)
         END DO
      END DO

c--sweep in x direction
      DO  J = 2,N
         DO L = 2,N
            DO I = 2,N
               b(I-1) = S(I,J,L)-(DP(I,J,L+1)+DP(I,J,L-1))/DS2
               b(I-1) = b(I-1) - (DP(I,J+1,L)+DP(I,J-1,L))/DS2
            END DO
            CALL THOMAS(d1,d2,d3,b,N-1)
            DP(2:N,J,L) = b(1:N-1)
         END DO
      END DO

c--sweep in y direction
      DO  I = 2,N
         DO L = 2,N
            DO J = 2,N
               b(J-1) = S(I,J,L)-(DP(I+1,J,L)+DP(I-1,J,L))/DS2
               b(J-1) = b(J-1) - (DP(I,J,L+1)+DP(I,J,L-1))/DS2
            END DO
            CALL THOMAS(d1,d2,d3,b,N-1)
            DP(I,2:N,L) = b(1:N-1)
         END DO
      END DO

c--compute the residual
      RES_MAX = 0.
      DO  I = 2,N
         DO L = 2,N
            DO J = 2,N
               TMP = S(I,J,L)-(DP(I+1,J,L)-2.*DP(I,J,L)+DP(I-1,J,L))/DS2
               TMP = TMP     -(DP(I,J+1,L)-2.*DP(I,J,L)+DP(I,J-1,L))/DS2
               TMP = TMP     -(DP(I,J,L+1)-2.*DP(I,J,L)+DP(I,J,L-1))/DS2
               RES(I,J,L) = TMP
               IF(ABS(TMP).GT.RES_MAX) RES_MAX = ABS(TMP)
            END DO
         END DO
      END DO
      ITER = ITER +1
      WRITE(6,40) ITER,RES_MAX
      END DO
 40   FORMAT(1x,'iter = ',I4,5X,'residual = ',F10.6)
      END SUBROUTINE ADI

c------||--------||---------||---------||----------||----------
        SUBROUTINE check_err(S,N,DP,RES)
C this subroutine iterate the linear system (L)*dp = S, where (L) is
c the operator--poisson 2nd order discretization with N grid, using the ADI method.
c the routine iterate ITNUM times and return the residual RES at each point
c dp is the initial guess and also the return approx. solution
      INTEGER N,ITNUM
      REAL S(N+1,N+1,N+1),DP(N+1,N+1,N+1),RES(N+1,N+1,N+1)
      REAL d1(N-1),d2(N-1),d3(N-1),b(N-1)
      COMMON DL

      DS  = DL/float(N)
      DS2 = DS**2

      RES_MAX = 0.
      DO  I = 2,N
         DO L = 2,N
            DO J = 2,N
               TMP = S(I,J,L)-(DP(I+1,J,L)-2.*DP(I,J,L)+DP(I-1,J,L))/DS2
               TMP = TMP     -(DP(I,J+1,L)-2.*DP(I,J,L)+DP(I,J-1,L))/DS2
               TMP = TMP     -(DP(I,J,L+1)-2.*DP(I,J,L)+DP(I,J,L-1))/DS2
               RES(I,J,L) = TMP
               IF(ABS(TMP).GT.RES_MAX) RES_MAX = ABS(TMP)
            END DO
         END DO
      END DO
      WRITE(6,50) RES_MAX
 50   FORMAT(1x,'the residual = ',F10.6)
      END SUBROUTINE check_err
c------||--------||---------||---------||----------||----------
        SUBROUTINE THOMAS(a0,b0,c0,g,Nx)

c Solves the system Ax=g for x using the Thomas algorithm,
c assuming A is tridiagonal and diagonally dominant.  It is
c assumed that (a,b,c,g) are previously-defined vectors of
c length n, where a is the subdiagonal, b is the main diagonal,
c and c is the superdiagonal of the matrix A.  The vectors
c (a,b,c) are replaced by the m_i and U on exit, and the vector
c g is replaced by the solution x of the original system.  

        INTEGER Nx
        REAL a0(Nx),b0(Nx),c0(Nx)
        REAL a(Nx),b(Nx),c(Nx),g(Nx)

c -----keep a0, b0, c0 unchanged
	a = a0
	b = b0
	c = c0

c  -------------- FORWARD SWEEP --------------

        DO I = 1,Nx-1

c Compute m_(j+1).  Note that we can put m_(j+1) in the location
c (below the diagonal!) that a_(j+1) used to sit without disrupting
c the rest of the algorithm, as a_(j+1) is set to zero by construction
c during this iteration.

           a(I+1) = - a(I+1) / b(I);

c Add m_(j+1) times the upper triangular part of the j'th row of
c the augmented matrix to the (j+1)'th row of the augmented
c matrix.

           b(I+1) = b(I+1) + a(I+1) * c(I);
           g(I+1) = g(I+1) + a(I+1) * g(I);
        END DO

c ------------ BACK SUBSTITUTION ------------

        g(Nx) = g(Nx) / b(Nx);
        DO I = Nx-1,1,-1
           g(I) = ( g(I) - c(I) * g(I+1) ) / b(I);
        END DO
        END SUBROUTINE THOMAS




