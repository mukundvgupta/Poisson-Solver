      program poisson

      INTEGER N, MAX_IT, ITER
      double precision, ALLOCATABLE :: d1(:),d2(:),d3(:),X(:),Y(:),b(:)
      double precision, ALLOCATABLE :: S(:,:),P(:,:)
      PI = ATAN(1.)*4.
      OPEN(1, FILE = 'pson2d.in')
      READ(1,*) N
      READ(1,*) TOL
      READ(1,*) MAX_IT
      READ(1,*) beta
      CLOSE(1)

      WRITE(6,*) "N       =", N
      WRITE(6,*) "TOL     =", TOL
      WRITE(6,*) "MAX_IT  =", MAX_IT
      WRITE(6,*) "beta    =", beta

      ALLOCATE (d1(N-1),d2(N-1),d3(N-1),X(N+1),Y(N+1),b(N-1))
      ALLOCATE (S(N+1,N+1),P(N+1,N+1))

      DS = 2.*PI/FLOAT(N)
      DO I=1,N+1
         X(I) = (I-1)*DS
         Y(I) = X(I)
      END DO
      P  = 0.
      DS2 = DS**2
      RES = 1.                  !residual
      ITER = 1

      DO I = 1,N-1
         d1(I) = 1./DS2
         d2(I) = -(2+beta)/DS2  !-2./DS2 two small and -2./DS2*2. two large
         d3(I) = 1./DS2         
      END DO

c	WRITE(*,10) (d1(I),I=1,N)
c	WRITE(*,10) (d2(I),I=1,N)
c	WRITE(*,10) (d3(I),I=1,N)
c10 	format(1x,10(F10.5,1x))

      DO  I =1,N+1
         DO J = 1,N+1
            S(I,J) = 0.
            DO K = 1,4
               S(I,J)=S(I,J)+sin(X(I)*K)*sin(Y(J)*K)
            END DO
         END DO
      END DO

      DO WHILE ((RES.GT.TOL) .AND. (ITER.LT.MAX_IT))
c--sweep in x direction
      DO  J = 2,N
         DO I = 2,N
            b(I-1) = S(I,J)  - (P(I,J+1)+P(I,J-1) - (2-beta)*P(I,J))/DS2
         END DO
         CALL THOMAS(d1,d2,d3,b,N-1)
         P(2:N,J) = b(1:N-1)
      END DO

c--sweep in y direction
      DO  I = 2,N
         DO J = 2,N
            b(J-1) = S(I,J)-(P(I+1,J)+P(I-1,J) - (2-beta)*P(I,J))/DS2
         END DO
         CALL THOMAS(d1,d2,d3,b,N-1)
         P(I,2:N) = b(1:N-1)
      END DO

c--compute the residual
      RES_MAX = 0.
      DO  I = 2,N
         DO J = 2,N
            RES = S(I,J)-(P(I+1,J)-2.*P(I,J)+P(I-1,J))/DS2
            RES = RES     -(P(I,J+1)-2.*P(I,J)+P(I,J-1))/DS2
            RES = ABS(RES)
            IF (RES.GT.RES_MAX) RES_MAX = RES
         END DO
      END DO
      RES = RES_MAX
      ITER = ITER +1
      WRITE(6,40) ITER,RES
      END DO
 40   FORMAT(1x,'iter = ',I4,5X,'residual = ',es10.4)
      
c--output the results
      OPEN (1,FILE = 'psonout.m',STATUS = 'OLD')
      WRITE(1,*) 'p=['
      DO K = 1,N+1
         WRITE(1,50) (S(K,J), J = 1,N+1) 
      END DO
      WRITE(1,*) '];'
      CLOSE(1)
 50   FORMAT(1x,50(F10.5,1x))

      STOP
      END

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
        double precision a0(Nx),b0(Nx),c0(Nx)
        double precision a(Nx),b(Nx),c(Nx),g(Nx)

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




