! SUBROUTINE to calculate the inverse of a matrix
!------------------------------------------------
!Modified/Created (2017) M.Drolia
!--------------------------------
! Given matrix A, finds the inverse of ANew(NA,NA) which is the NA*NA 
! matrix taken as the first NA*NA elements from top left row/column.
!
! Uses LAPACK library.
!Input Arguments:~~~~~~~~~~~~~~~~~~~~~~~~~~
!A	Input matrix, supposed to be square
!IA	Dimension of the matrix A = (IA,IA)
!NA	Dimension of the submatrix that needs to be inverted. NA <= IA
!Output Arguments:~~~~~~~~~~~~~~~~~~~~~~~~~
!Ainv	Inverted matrix of dimension (NA,NA)
!_________________________________________

SUBROUTINE inv(A,IA,Ainv,NA)
double complex A,ASaved,Ainv
integer IA,IPIV, INFO, index_i,NA,index_j

double complex, allocatable :: ANew(:,:),B(:,:)

Dimension A(IA,IA),ASaved(IA,IA),IPIV(NA),Ainv(NA,NA)



! External procedures defined in LAPACK
  external ZGESV

!save a copy of A
ASaved = A

allocate (ANew(NA,NA),B(NA,NA))
!Create identity matrix of size NA
B = 0
Do index_i = 1,NA
	B(index_i,index_i) = 1
enddo

!Create the sub-matrix ANew
Do index_i = 1,NA
	Do index_j = 1,NA
		ANew(index_i,index_j) = A(index_i,index_j)
	enddo
enddo
!write(*,*)'--------------------------------------------------'
!write(*,*)ANew
!write(*,*)'--------------------------------------------------'

!Invert the matrix ANew
CALL ZGESV( NA, NA, ANew, NA, IPIV, B, NA, INFO )
!Error check
      IF( INFO.GT.0 ) THEN
         WRITE(*,*)'The diagonal element of the triangular factor of A,'
         WRITE(*,*)'U(',INFO,',',INFO,') is zero, so that'
         WRITE(*,*)'A is singular; the solution could not be computed.'
         STOP
      END IF
!Set the corresponding arrays
      Ainv = B
      A = ASaved
Return
end !function inv


SUBROUTINE MATMULCPLX(A, IA, JA, B, IB, JB, C, IC, JC, L, M, N)
!
!***MATRIX MULTIPLICATION
!***(c) Peter and Jacqueline A. Bettess, 1986
!**Modified (2017) M.Drolia
!-----------------------------------------------------------------------
! PURPOSE
!      Post multiplies matrix A by matrix B to give matrix C.  All are
!      as 2 dimensional arrays
! HISTORY
!      Written June 1986
!
! ARGUMENTS IN
!      A      Matrix A
!      IA     First dimension of matrix A
!      JA     Second dimension of matrix A
!      B      Matix B
!      IB     First dimension of matrix B
!      JB     Second dimension of matrix B
!      IC     Number of rows in C
!      JC     Number of columns in C
!      L      Number of rows used in A and C
!      M      Number of columns used in B and C
!      N      Number of columns used in A and rows used in B
!
! ARGUMENTS OUT
!      C      Product matrix, C = A * B
!
!***********************************************************************
      DOUBLE complex A, B, C
      INTEGER IA, IB, IC, IL, IM, I_N, JA, JB, JC, L, M, N
      DIMENSION A(IA,JA), B(IB,JB), C(IC,JC)
!
!***process rows in A and C
!
      DO 30 IL = 1, L
!
!***process columns in B and C
!
        DO 20 IM = 1, M
          C(IL,IM) = DCMPLX(0,0)
!
!***form inner product
!
          DO 10 I_N = 1, N
            C(IL,IM) = C(IL,IM) + A(IL,I_N) * B(I_N,IM)
   10     CONTINUE
   20   CONTINUE
   30 CONTINUE
      RETURN
      END
!************************************************************************
SUBROUTINE modquo(a,b,q,r)
integer a,b,q,r

!Compute the quotient
q = floor(real(a)/real(b))
!Compute the remainder
r = mod(a,b)

Return
end !function modquo



! SUBROUTINE to calculate the product of a block matrix 
! with another regular matrix (or vector)
!------------------------------------------------
!Modified/Created (2017) M.Drolia
!--------------------------------
!_________________________________________

SUBROUTINE MATMULCPLX_block(A,IA,NANGL,B,IB,JB,C,TOTNOD)
double complex A,B,C
integer IA,NANGL,IB,JB,TOTNOD,iNODE,iROW,IL,IM,I_N,i,j
double complex, allocatable :: ANew(:,:),BNew(:,:),CTemp(:,:)

Dimension A(IA,NANGL),B(IB,JB),C(IA,JB)

allocate (ANew(NANGL,NANGL),BNew(NANGL,JB),CTemp(NANGL,JB))

Do iNODE = 1,TOTNOD
	iROW = (iNODE - 1)*NANGL !This is the node for which the corresponding block is being multiplied
	!Create the sub-matrices ANew,BNew
	ANew = A(iROW+1:iROW+NANGL,1:NANGL)
	BNew = B(iROW+1:iROW+NANGL,1:JB)
	!Multiply the sub-matrices
	DO IL = 1,NANGL
		DO IM = 1,JB
			CTemp(IL,IM) = DCMPLX(0,0)
			Do I_N = 1,NANGL
				CTemp(IL,IM) = CTemp(IL,IM) + ANew(IL,I_N)*BNew(I_N,IM)
			enddo
		enddo
	enddo		
	!Update the product in the subsequent section of C (output) matrix
	C(iROW+1:iROW+NANGL,1:JB) = CTemp(1:NANGL,1:JB)
enddo !iNODE
Return
end !function MATMULCPLX_block

! SUBROUTINE to calculate the inverse of a (Block) matrix
!------------------------------------------------
!Modified/Created (2017) M.Drolia
!--------------------------------
! Given matrix A, finds the inverse in blocks of ANew(NA,NA) 
! which is the NANGL*NANGL matrix taken as the 
! first NANGL*NANGL elements from top left row/column, for each 
! of the TOTNOD nodes.
!
! Uses LAPACK library.
!Input Arguments:~~~~~~~~~~~~~~~~~~~~~~~~~~
!A	Input matrix, supposed to be square
!IA	Dimension of the matrix A = (IA,NANGL), IA = TOTDOF = NANGL*TOTNOD
!NANGL	Enrichments
!Output Arguments:~~~~~~~~~~~~~~~~~~~~~~~~~
!Ainv	Inverted matrix of dimension (IA,NANGL)
!_________________________________________

SUBROUTINE inv_Block(A,IA,NANGL,Ainv,TOTNOD)
double complex A,ASaved,Ainv
integer IA,IPIV, INFO, index_i,NANGL
integer iNODE,iROW,TOTNOD,i,j

double complex, allocatable :: ANew(:,:),B(:,:)

Dimension A(IA,NANGL),ASaved(IA,NANGL),IPIV(NANGL),Ainv(IA,NANGL)


! External procedures defined in LAPACK
  external ZGESV

!save a copy of A
ASaved = A

allocate (ANew(NANGL,NANGL),B(NANGL,NANGL))

Do iNODE = 1,TOTNOD
	iROW = (iNODE - 1)*NANGL	!This is the node for which the corresponding block is being inverted
	!Create identity matrix of size NANGL
	B = 0
	Do index_i = 1,NANGL
		B(index_i,index_i) = 1
	enddo
	!Create the sub-matrix ANew
	ANew = A(iROW+1:iROW+NANGL,1:NANGL)
	!Invert the matrix ANew
	CALL ZGESV( NANGL, NANGL, ANew, NANGL, IPIV, B, NANGL, INFO )
	!Error check
	      IF( INFO.GT.0 ) THEN
		 WRITE(*,*)'The diagonal element of the triangular factor of A,'
		 WRITE(*,*)'U(',INFO,',',INFO,') is zero, so that'
		 WRITE(*,*)'A is singular; the solution could not be computed.'
		 STOP
	      END IF
	!Set the corresponding arrays
	Ainv(iROW+1:iROW+NANGL,1:NANGL) = B
enddo !iNODE
A = ASaved
Return
end !function inv_Block

! SUBROUTINE to calculate the Gauss-Legendre abscissae and weights
!------------------------------------------------
!Modified/Created (2020) M.Drolia
!Description:
! gauleg.f90     P145 Numerical Recipes in Fortran
! compute x(i) and w(i)  i=1,n  Legendre ordinates and weights
! on interval -1.0 to 1.0 (length is 2.0)
! use ordinates and weights for Gauss Legendre integration
!	In: (x1,x2) == range of integration
!	  : n       == dimension of abscissa and weight vectors
!	OUT: (x,w)  == abscissae and weights
!	NB:PI=4.0D0*DATAN(1.0D0) gives 
!	erroneous x,w for the middle point when n=odd.
!--------------------------------
subroutine gaulegf(x1, x2, x, w, n)
  implicit none
  integer n
  double precision x1, x2
  double precision x(n), w(n)
  integer i, j, m
  double precision p1, p2, p3, pp, xl, xm, z, z1
  double precision, parameter :: eps=3.d-14
  !definitions extra
  double precision PI
  PI = 3.141592654d0!4.0D0*DATAN(1.0D0) 
      
  m = (n+1)/2
  xm = 0.5d0*(x2+x1)
  xl = 0.5d0*(x2-x1)
  do i=1,m
    z = cos(PI*(i-0.25d0)/(n+0.5d0))
    z1 = 0.0
    do while(abs(z-z1) .gt. eps)
      p1 = 1.0d0
      p2 = 0.0d0
      do j=1,n
        p3 = p2
        p2 = p1
        p1 = ((2.0d0*j-1.0d0)*z*p2-(j-1.0d0)*p3)/j
      end do
      pp = n*(z*p1-p2)/(z*z-1.0d0)
      z1 = z
      z = z1 - p1/pp
    end do
    x(i) = xm - xl*z
    x(n+1-i) = xm + xl*z
    w(i) = (2.0d0*xl)/((1.0d0-z*z)*pp*pp)
    w(n+1-i) = w(i)
  end do
end subroutine gaulegf
