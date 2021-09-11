      program pson_cg

      INTEGER N, MAX_IT, ITER, PRECONDTION
      REAL, ALLOCATABLE :: X(:),Y(:),Z(:)
      REAL, ALLOCATABLE :: S(:,:,:),P(:,:,:),r(:,:,:),g(:,:,:)
      REAL, EXTERNAL :: VEC_VEC
      COMMON DS,DS2
      COMMON PRECONDTION

      PI = ATAN(1.)*4.
      OPEN(1, FILE = 'poisson.in')
      READ(1,*) N
      READ(1,*) TOL
      READ(1,*) MAX_IT
      READ(1,*) PRECONDTION
      CLOSE(1)

      WRITE(6,*) "N          =", N
      WRITE(6,*) "TOL        =", TOL
      WRITE(6,*) "MAX_IT     =", MAX_IT
      WRITE(6,*) "PRECONDTION=", PRECONDTION

      ALLOCATE (X(N+1),Y(N+1),Z(N+1))
      ALLOCATE (S(N+1,N+1,N+1),P(N+1,N+1,N+1))
      ALLOCATE (r(N+1,N+1,N+1),g(N+1,N+1,N+1))

      DS = 2.*PI/FLOAT(N)
      DO I=1,N+1
         X(I) = (I-1)*DS
         Y(I) = X(I)
         Z(I) = X(I)
      END DO
      P = 0.
      r = 0.
      g = 0.
      DS2 = DS**2
      RES = 1.                  !residual, 1-norm
      RES2 = 1.                 !residual, 2-norm
      ITER = 0

      DO  I =1,N+1
         DO J = 1,N+1
            DO L = 1,N+1
               S(I,J,L) = 0.
               DO K = 1,4
                  S(I,J,L)=S(I,J,L)+sin(X(I)*K)*sin(Y(J)*K)*sin(Z(L)*K)
               END DO
            END DO
         END DO
      END DO
c - preconditioning the matrix S
      IF(PRECONDTION.EQ.1) CALL PRECOND(S,N)

      DO WHILE ((RES.GT.TOL) .AND. (ITER.LT.MAX_IT))
         CALL COMPUTE_AX(P,N,r)
         r = S - r
         RES2_SAVE = RES2
         RES2 = VEC_VEC(r,r,N)   !here the RES2 is the residual r'*r
         IF (ITER .EQ. 1) THEN
            g = r
         ELSE
            BETA = RES2/RES2_SAVE
            g = r + BETA*g
         END IF
         CALL COMPUTE_AX(g,N,r)  !here we don't need r anymore, so it is overwritten

         ALFA = RES2/VEC_VEC(g,r,N)
         P = P + ALFA*g

c--compute the residual
         RES_MAX = 0.
         DO  I = 2,N
            DO L = 2,N
               DO J = 2,N
                  RES = ABS(r(I,J,L))
                  IF (RES.GT.RES_MAX) RES_MAX = RES
               END DO
            END DO
         END DO
         RES = RES_MAX

         ITER = ITER +1
         WRITE(6,40) ITER,RES
      END DO
 40   FORMAT(1x,'iter = ',I4,5X,'residual = ',F10.6)

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
 50   FORMAT(1x,50(F10.5,1x))

      STOP
      END

c------||--------||---------||---------||----------||----------

      FUNCTION VEC_VEC(X,Y,M)
c compute the inner product of two "3-D" vectors X and Y
      INTEGER M
      REAL X(M+1,M+1,M+1), Y(M+1,M+1,M+1), VEC_VEC

      VEC_VEC = 0.
      DO  I = 2,M
         DO L = 2,M
            DO J = 2,M
               VEC_VEC = VEC_VEC + X(I,J,L)*Y(I,J,L) 
            END DO
         END DO
      END DO
      END FUNCTION VEC_VEC

c------||--------||---------||---------||----------||----------
      SUBROUTINE  COMPUTE_AX(X,N,Y)
c compute the product of the matrix A and the vector x
c the result is saved in the vetor y
c note that here the vetors x and y are "3-D" vetors
      INTEGER N, PRECONDTION
      REAL X(N+1,N+1,N+1), Y(N+1,N+1,N+1)
      COMMON DS,DS2
      COMMON PRECONDTION

      DO  I = 2,N
         DO J = 2,N
            DO L = 2,N
               TMP =       (X(I+1,J,L)-2.*X(I,J,L)+X(I-1,J,L))/DS2
               TMP = TMP + (X(I,J+1,L)-2.*X(I,J,L)+X(I,J-1,L))/DS2
               TMP = TMP + (X(I,J,L+1)-2.*X(I,J,L)+X(I,J,L-1))/DS2
               Y(I,J,L) = TMP
            END DO
         END DO
      END DO
c -precondtioning Y
      IF(PRECONDTION.EQ.1) CALL PRECOND(Y,N)
      END SUBROUTINE  COMPUTE_AX

c------||--------||---------||---------||----------||----------
      SUBROUTINE  PRECOND(X,N)
c compute the product of the matrix inv(M) and the vector x
c that is, solve M*y = x, the result is saved in the vetor x
c note that here the vetor x is "3-D" vetors
      INTEGER N
      REAL X(N+1,N+1,N+1)  
      COMMON DS,DS2
      REAL XT((N-1)**3),YT((N-1)**3),d1((N-1)**3)
      REAL d2((N-1)**3),d3((N-1)**3)

c -shift 3-D vector X to 1-D vector XT
      DO  I = 2,N
         DO J = 2,N
            DO L = 2,N
               IT=(I-2)*(N-1)**2+(J-2)*(N-1)+L-1
               XT(IT) = X(I,J,L)
            END DO
         END DO
      END DO

c -compute the preconditioner matrix T
      DO I = 1,(N-1)**3
         d2(I) = 1./DS2      !here the matrix T is the tridiagonal of matrix A
      END DO
      d1 = d2
      d3 = d2
      DO I = 1,(N-1)**2
         d1((I-1)*(N-1)+1) = 0.  
         d3(I*(N-1)) = 0. 
      END DO

      CALL THOMAS(d1,-6.*d2,d3,XT,(N-1)**3)
c -shift  1-D vector XT to 3-D vector X
      DO  I = 2,N
         DO J = 2,N
            DO L = 2,N
               IT=(I-2)*(N-1)**2+(J-2)*(N-1)+L-1
               X(I,J,L) = XT(IT)
            END DO
         END DO
      END DO
      END SUBROUTINE PRECOND

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

