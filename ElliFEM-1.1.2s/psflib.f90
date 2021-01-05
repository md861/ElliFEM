!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getSF(XI, ETA, SF, NODEL, SFDL, NDIME,NORDER)

           !Description: This subroutine creates the array to store 
           		!shape functions and their derivatives
           !Arguments out
           		!SF: 	shape functions
           		!SFDL:	Local derivatives

  IMPLICIT NONE
  double precision XI, ETA, SF, SFDL
  double precision XI_i,ETA_i,LPx,LPy,XI_j,ETA_j
  double precision XI_k,ETA_k
  double precision Sx,Sy,LP_Sx,LP_Sy
  double precision RefNODES(NDIME,NORDER+1)
  integer NODEL, NDIME,NORDER
  integer ELNODMAP(NDIME,NODEL)
  integer I,J,K
  DIMENSION SF(NODEL),SFDL(NDIME,NODEL)
  
  
  !get (default equi-spaced) local co-ordinates
  call getRefNODES(NORDER,NDIME,RefNODES)
  !get MAP that gives index for local co-ordinates in RefNODES
  call getELNODMAP(NORDER,NODEL,NDIME,ELNODMAP)
  
  Do I = 1,NODEL       
       !LPx = Lagrange Polynomials in xi (order = NORDER)
       XI_i = RefNODES(1,ELNODMAP(1,I))
       LPx = 1.0D0      
       !LPy = Lagrange Polynomials in eta (order = NORDER)
       ETA_i = RefNODES(2,ELNODMAP(2,I))
       LPy = 1.0D0
       
       Sx = 0.0D0	!Sum to store LP'x
       Sy = 0.0D0	!Sum to store LP'y
       
       DO J = 1,NORDER+1    
         !for XI     
         if(J .NE. ELNODMAP(1,I))THEN !{j .NE. i}   
           !Build Lagrange polynomials
           XI_j = RefNODES(1,J)
           LPx  = LPx * ((XI - XI_j)/(XI_i - XI_j))
           !Build derivatives for Lagrange polynomials
           LP_Sx = 1.0D0
           Do K = 1,NORDER+1
             if((K .NE. ELNODMAP(1,I)).AND.(K .NE. J))THEN !{K .NE. i}and{K .NE. j}  
               XI_k  = RefNODES(1,K) 
               LP_Sx = LP_Sx * ((XI - XI_k)/(XI_i - XI_k))
             endif
           endDO
           LP_Sx = LP_Sx * (1.0D0/(XI_i - XI_j))
           Sx = Sx + LP_Sx
         endif   
         
         !for ETA
         if(J .NE. ELNODMAP(2,I))THEN !{j .NE. i}
           !Build Lagrange Polynomials LPy
           ETA_j = RefNODES(2,J)
           LPy = LPy * ((ETA - ETA_j)/(ETA_i - ETA_j))
           !Build derivatives for LPy
           LP_Sy = 1.0D0
           Do K = 1,NORDER+1
             if((K .NE. ELNODMAP(2,I)).AND.(K .NE. J))THEN !{K .NE. i}and{K .NE. j} 
               ETA_k = RefNODES(2,K)
               LP_Sy = LP_Sy * ((ETA - ETA_k)/(ETA_i - ETA_k)) 
             endif
           endDO
           LP_Sy = LP_Sy * (1.0D0/(ETA_i - ETA_j))
           Sy = Sy + LP_Sy
         endif         
       endDO !{j}
       
       !Calculate the shape function          
       SF(I)     = LPx*LPy
       !Calculate local derivatives
       SFDL(1,I) = Sx*LPy
       SFDL(2,I) = LPx*Sy      
  endDO  
end subroutine getSF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getSFN(X,Y,SF,NODEL,K_W,NANGL,BETA,SFN,ENR_Type)
           !Description: This subroutine creates the array to store 
           		!enriched shape functions
           !Arguments out
           		!SFN: Enriched shape functions
           		
IMPLICIT NONE
double precision X,Y,SF(NODEL),K_W
integer NANGL,NODEL
  !definitions for ENR_Type
  integer ENR_type		! 1=ENR2, 0=ENR1
  !changes for ENR_Type
  double precision BETA(NANGL-ENR_type)
double complex SFN(NODEL*NANGL)
double complex CmplxTerm
INTEGER temp,I,q
double precision Z(NANGL)

	CmplxTerm = DCMPLX(1,0)
	temp = NANGL - ENR_type
	DO I = 1,temp
	  Z(I) = X*DCOS(BETA(I))+Y*DSIN(BETA(I))
	endDO
	if(ENR_Type .EQ. 1)Then
	  Z(NANGL)=0
	endif

	DO I = 1,NODEL
	  DO q = 1,temp
	    SFN(NANGL*(I-1)+q) = SF(I)*CDEXP(DCMPLX(0.0D0,K_W*Z(q)))
	  endDO
	  if(ENR_Type .EQ. 1)Then
	    SFN(NANGL*(I-1)+NANGL) = CmplxTerm*SF(I)
	  endif
	endDO
end subroutine getSFN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getSFDGN(X, Y,SF,NODEL,SFDG,NDIME,K_W,& 
                    NANGL,BETA,SFDGN,ENR_type)
                    
IMPLICIT NONE
double precision X,Y,SF(NODEL),&
                 SFDG(NDIME,NODEL),&
                 K_W
  !definitions for ENR_Type
  integer ENR_type		! 1=ENR2, 0=ENR1
  double precision BETA(NANGL-ENR_type)
integer NANGL,NODEL,NDIME
double complex SFDGN(NDIME,NODEL*NANGL)
double complex CmplxTerm
integer I,temp,q
double precision Z(NANGL)
                    
      CmplxTerm = DCMPLX(1,0)
      temp = NANGL - ENR_Type
      DO I = 1,temp
         Z(I) = X*DCOS(BETA(I))+Y*DSIN(BETA(I))
      endDO
      if(ENR_Type .EQ. 1)Then
        Z(NANGL)=0
      endif
 
 
      DO I = 1,NODEL
      DO  q = 1,NANGL-ENR_Type
         SFDGN(1,NANGL*(I-1)+q)=(SFDG(1,I)+SF(I)*DCOS(BETA(q))*&
           DCMPLX(0.0D0,K_W)) * CDEXP(DCMPLX(0.0D0,K_W*Z(q)))
         SFDGN(2,NANGL*(I-1)+q)=(SFDG(2,I)+SF(I)*DSIN(BETA(q))*&
           DCMPLX(0.0D0,K_W)) * CDEXP(DCMPLX(0.0D0,K_W*Z(q)))
      endDO
         if(ENR_Type .EQ. 1)Then
           SFDGN(1,NANGL*(I-1)+NANGL)=CmplxTerm*SFDG(1,I)
           SFDGN(2,NANGL*(I-1)+NANGL)=CmplxTerm*SFDG(2,I)
         endif
     endDO
                    
end subroutine getSFDGN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getELNODMAP(NORDER,NODEL,NDIME,ELNODMAP)

           !Description: This subroutine creates the array to store 
           		!local co-ordinates with some nodal 
           		!spacing. (default is equispaced)
           !Arguments out
           		!ELNODMAP: Table that stores the 
           			!local co-ordinates of elemental
           			!nodes.
           			
  IMPLICIT NONE         			
  integer  NORDER ,NODEL ,NDIME
  integer ELNODMAP
  DIMENSION ELNODMAP(NDIME,NODEL)
  integer I,J,counter
  
    if((NORDER) .EQ. 1) THEN!for 1st-order elements {Historic indexing}
      !4----3
      !|    |
      !|    |
      !1----2
      ELNODMAP(1,1) = 1
      ELNODMAP(1,2) = 2
      ELNODMAP(1,3) = 2
      ELNODMAP(1,4) = 1
      
      ELNODMAP(2,1) = 1
      ELNODMAP(2,2) = 1
      ELNODMAP(2,3) = 2
      ELNODMAP(2,4) = 2
    else!for higher-order elements
      !7--8--9		!13--14--15--16
      !|     |		!|            |
      !4  5  6	   ,	!9   10  11  12    , ...
      !|     |		!|            |
      !1--2--3		!5   6   7    8 
      			!|            |
      			!1---2---3----4
      counter = 1			
      Do I = 1,NORDER+1
        DO J = 1,NORDER+1
          ELNODMAP(1,counter) = J
          ELNODMAP(2,counter) = I
          counter = counter + 1
        endDO
      endDO
      !write(*,*)'counter = ',counter
    endif
end subroutine getELNODMAP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getRefNODES(NORDER,NDIME,RefNODES)
           !Description: This subroutine creates the array to store 
           		!local co-ordinates in the (XI,ETA) plane.
           !Arguments out
           		!RefNODES: Local co-ordinates of element nodes.
  IMPLICIT NONE
  integer NORDER,NDIME
  double precision RefNODES(NDIME,NORDER+1),h
  integer I
  
  h = (1.0D0 - (-1.0D0))/(NORDER)  
  !Create equi-spaced local nodes in [-1,1]
  DO I = 1,NORDER+1
    RefNODES(1,I) = -1.0D0 + (h)*(I-1)
    RefNODES(2,I) = -1.0D0 + (h)*(I-1)
  endDO
end subroutine getRefNODES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine GetELNODE(ELNODDAT,NODEL,NORDER)

           !Description: This subroutine creates the array to store 
           		!referencing indices for local nodes, from 
           		!the global array (ELDAT).
           !Arguments out
           		!ELNODDAT: Table that stores the order of 
           			!global nodes in ELDAT that 
           			!correspond to the local index
           			!in the reference element.
           			!This order is (for high order)
           			!Bottom->Top, Left->Right
           			
  IMPLICIT NONE         			
  integer ELNODDAT, NODEL, NORDER
  DIMENSION ELNODDAT(NODEL)

    if((NORDER) .EQ. 1) THEN!for 1st-order elements
      !4----3
      !|    |
      !|    |
      !1----2
      ELNODDAT = (/2,3,4,1/)
    else if((NORDER) .EQ. 2) THEN!for 2nd-order elements
      !7--8--9
      !|     |
      !4  5  6
      !|     |
      !1--2--3
      ELNODDAT = (/2,6,3,5,9,7,1,8,4/)
    else if((NORDER) .EQ. 3) THEN!for 3rd-order elements
      !13--14--15--16
      !|            |
      !9   10  11  12
      !|            |
      !5   6   7    8 
      !|            |
      !1---2---3----4
      ELNODDAT = (/2,7,8,3,6,14,15,9,5,13,16,10,1,12,11,4/)
    else if((NORDER) .EQ. 4) THEN!for 4th-order elements
      !21--22--23--24----25
      !|                 |
      !16  17  18  19    20 
      !|                 |
      !11  12  13  14    15    
      !|                 |
      !6   7   8    9    10
      !|                 |
      !1---2---3----4----5
      ELNODDAT = (/2,8,9,10,3,7,18,22,19,11,6,21,25,23,12,5,17,24,20,13,1,16,15,14,4/)
    endif
end subroutine GetELNODE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
