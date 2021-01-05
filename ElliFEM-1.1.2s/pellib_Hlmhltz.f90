!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE getHlmhltz_M_S_RHS(IELEM,ELCOR,NODEL,NDIME,&
		   NGAUSS,NORDER,K_W,Omega,&
		   ELM,ELK,PTLD_Phi,&!PTLD_Vlcty,&
		   NANGL,NODEF,BETA,&
		   ENR_Type,&
		   EFACT,NMARK,ELTYPE,&
		   Kp,Theta_p)
           !Description: This subroutine builds elementary 
           		!matrices for Mass, Stiffness, and 
           		!vector for load. 
           		!This was modified over DIC subroutine.
           !Arguments out
           		!ELM: Elementary Mass 
           		!ELK: Elementary stiffness 
           		!PTLD: Elementary load vector

IMPLICIT NONE
double precision ELCOR
integer NODEL,NDIME,NORDER
DIMENSION ELCOR(NODEL,NDIME)

INTEGER NGAUSS, IGAUSS, JGAUSS,I,J, IELEM
DOUBLE PRECISION WG(NGAUSS), XG(NGAUSS),&
                 XI, ETA, WTX, WTY, WTXY,&
                 X,Y
                 
double precision SF(NODEL), SFDL(NDIME,NODEL)
double precision JCBT(NDIME,NDIME),&
                 JCBTI(NDIME,NDIME), DETJ


double precision n(NDIME)
integer EDG_TYPE
integer NBC_pos

!definitions for DIC,M,S
double precision SFDG(NDIME,NODEL),&
                 r_mod,K_W,Omega
double complex	 ELM(NODEF,NODEF),&
		 ELK(NODEF,NODEF),&
		 FZ_Phi,&!FZ_Vlcty,&
		 PTLD_Phi(NODEF)!,&
		 !PTLD_Vlcty(NODEF)
		 
!definitions for PUFEM
integer NANGL,NODEF
double complex SFN(NODEF),SFDGN(NDIME,NODEF)

    !definitions for ENR_Type
    integer ENR_type		! 1=ENR2, 0=ENR1
    !changes for ENR_Type
    double precision BETA(NANGL-ENR_type)
    
!definitions for Hlmhltz
integer NMARK,ELTYPE,EFACT(NMARK,ELTYPE)
double complex GZ_Phi,U_p,dUdx,dUdy
double precision Kp,Theta_p

  
  !Initialize elementary matrices and vectors
  ELM = dcmplx(0.0D0,0.0D0)
  ELK = dcmplx(0.0D0,0.0D0)
  PTLD_Phi = dcmplx(0.0D0,0.0D0)
  !PTLD_Vlcty = dcmplx(0.0D0,0.0D0)
  
  !Get integration points.
  !call GAULEG(NGAUSS, XG, NGAUSS, WG, NGAUSS)
  call gaulegf(-1.0d0,1.0d0, XG, WG, NGAUSS)
  
  !Integrate inside domain.
    DO IGAUSS = 1,NGAUSS
      XI = XG(IGAUSS)
      WTX = WG(IGAUSS)
      DO JGAUSS = 1,NGAUSS
        ETA = XG(JGAUSS)
        WTY = WG(JGAUSS)
        WTXY = WTX*WTY
        !Get shape functions
        call getSF(XI, ETA, SF, NODEL, SFDL, NDIME,NORDER)
        !Get Jacobians.
          !Get JCBT = transpose(Jacobian)
          CALL MATMUL(SFDL, NDIME, NODEL, ELCOR, NODEL, NDIME, JCBT,&
          NDIME, NDIME, NDIME, NDIME, NODEL)
          !Invert JCBT
          CALL invJCBT(JCBT, NDIME, JCBTI, DETJ)


        !Update weights
        WTXY = WTXY*DETJ
        
        !Get global partial derivatives, i.e. (d/dx,d/dy) 
        call MATMUL(JCBTI,NDIME,NDIME,SFDL,NDIME,NODEL,SFDG,&
        	    NDIME,NODEL,NDIME,NODEL,NDIME)
        
        !Get (X,Y) from (XI,ETA)
        X = 0.0D0
        Y = 0.0D0
        DO I=1,NODEL
          X = X + SF(I)*ELCOR(I,1)
          Y = Y + SF(I)*ELCOR(I,2)
        ENDDO        
        
        !get Enriched shape functions and derivatives
          !get SFN
          call getSFN(X,Y,SF,NODEL,K_W,NANGL,BETA,SFN,ENR_type)
          !get SFDGN
          call getSFDGN(X, Y,SF,NODEL,SFDG,NDIME,K_W,& 
                    NANGL,BETA,SFDGN,ENR_type)
          !write(*,*)'X = ',x,'Y = ',y
          !Do I = 1,NODEF
          !   write(*,*)'SFDGN (1,',I,') =',SFDGN(1,I)
          !  write(*,*)'SFDGN (2,',I,') =',SFDGN(2,I)
          !endDO

        !Build Mass and Stiffness
        DO I = 1,NODEF
          DO J = 1,NODEF
          ELM(J,I) = ELM(J,I) + SFN(J)*SFN(I)*WTXY
          ELK(J,I) = ELK(J,I) + (SFDGN(1,J)*SFDGN(1,I) + &
                                 SFDGN(2,J)*SFDGN(2,I))*WTXY
          endDO
        endDO

        !Build RHS for DIC
          !Calculate load          
          r_mod = sqrt(X*X + Y*Y)
          FZ_Phi = dcmplx(0.0D0,0.0D0)!CDEXP(DCMPLX(0.0D0,1.0D0)*(K_W*r_mod&
				      !	 - Omega*(0.0D0)))
          !FZ_Vlcty = dcmplx(0.0D0,0.0D0)!CDEXP(DCMPLX(0.0D0,1.0D0)*(K_W*r_mod&
	  !				!- Omega*(0.0D0)))&
	  !		!*(DCMPLX(0.0D0,-1.0D0)*Omega)
					 
          !Calculate inner product of load with weight functions
          DO J = 1,NODEF
            PTLD_Phi(J) = PTLD_Phi(J) + FZ_Phi*SFN(J)*WTXY
            !PTLD_Vlcty(J) = PTLD_Vlcty(J) + FZ_Vlcty*SFN(J)*WTXY
          endDO
        
      endDO
    endDO
 
  !Integrate over boundary. {For non-zero BCs}
  NBC_pos = 1
  !edge 1
   if(EFACT(NBC_pos,1).EQ. 1)THEN   
   EDG_TYPE = 1
   !write(*,*)'-----------------'
   !write(*,*) 'ELEMENT = ',IELEM,' EDGE = 1'
     Do IGAUSS = 1,NGAUSS
       XI = XG(IGAUSS)
       WTX = WG(IGAUSS)
       ETA = -1.0
       WTY = 1
       WTXY = WTX*WTY
       !Get shape functions
       call getSF(XI, ETA, SF, NODEL, SFDL, NDIME,NORDER)
       !Get Jacobians.
         !Get JCBT = transpose(Jacobian)
         CALL MATMUL(SFDL, NDIME, NODEL, ELCOR, NODEL, NDIME, JCBT,&
         NDIME, NDIME, NDIME, NDIME, NODEL)      
       
       !Update weights
         !get Determinant for line integral = dsqrt((dx/dXI)^2 + (dy/dXI)^2)
         DETJ = dsqrt((JCBT(1,1)**2) + (JCBT(1,2)**2))
         !Update weight
         WTXY = WTXY*DETJ 
         !write(*,*) 'DETJ = ',DETJ
       
       !Get (X,Y) from (XI,ETA)
       X = 0.0D0
       Y = 0.0D0
       DO I=1,NODEL
         X = X + SF(I)*ELCOR(I,1)
         Y = Y + SF(I)*ELCOR(I,2)
       ENDDO 
       
       !get Enriched shape functions
         !get SFN
         call getSFN(X,Y,SF,NODEL,K_W,NANGL,BETA,SFN,ENR_type)
       
       !Get normal to edge
       call getNormVec(EDG_TYPE,JCBT,NDIME,NODEL,X,Y,ELCOR,n,NORDER)
       !write(*,*)'n = ',n  
       
       !Compute load     
       U_p = CDEXP(DCMPLX(0.0D0,Kp)*(X*DCOS(Theta_p)+Y*DSIN(Theta_p)))
       dUdx = U_p*DCMPLX(0.0D0,Kp*DCOS(Theta_p))
       dUdy = U_p*DCMPLX(0.0D0,Kp*DSIN(Theta_p))
       GZ_Phi = dcmplx(0.0D0,0.0D0)
       GZ_Phi = (-1.0D0)*(dUdx*n(1) + dUdy*n(2))
       
       !Compute inner product of load with weight functions
       DO J = 1,NODEF
           PTLD_Phi(J) = PTLD_Phi(J) + GZ_Phi*SFN(J)*WTXY            
       endDO
     endDO
     !write(*,*)'-----------------'
   endif
  
   !edge 2
   if(EFACT(NBC_pos,2).EQ. 1)THEN  
   EDG_TYPE = 0 
   !write(*,*)'-----------------'
   !write(*,*) 'ELEMENT = ',IELEM,' EDGE = 2'
     Do IGAUSS = 1,NGAUSS
       ETA = XG(IGAUSS)
       WTY = WG(IGAUSS)
       XI = 1.0
       WTX = 1
       WTXY = WTX*WTY
       !Get shape functions
       call getSF(XI, ETA, SF, NODEL, SFDL, NDIME,NORDER)
       !Get Jacobians.
         !Get JCBT = transpose(Jacobian)
         CALL MATMUL(SFDL, NDIME, NODEL, ELCOR, NODEL, NDIME, JCBT,&
         NDIME, NDIME, NDIME, NDIME, NODEL)          
       
       
       !Update weights
         !get Determinant for line integral = dsqrt((dx/dETA)^2 + (dy/dETA)^2)
         DETJ = dsqrt((JCBT(2,1)**2) + (JCBT(2,2)**2))
         !Update weight
         WTXY = WTXY*DETJ  
         !write(*,*) 'DETJ = ',DETJ
!
       
       !Get (X,Y) from (XI,ETA)
       X = 0.0D0
       Y = 0.0D0
       DO I=1,NODEL
         X = X + SF(I)*ELCOR(I,1)
         Y = Y + SF(I)*ELCOR(I,2)
       ENDDO 
       
       !get Enriched shape functions
         !get SFN
         call getSFN(X,Y,SF,NODEL,K_W,NANGL,BETA,SFN,ENR_type)
       
       !Get normal to edge
       call getNormVec(EDG_TYPE,JCBT,NDIME,NODEL,X,Y,ELCOR,n,NORDER)
       !write(*,*)'n = ',n  
       
       !Compute load
       U_p = CDEXP(DCMPLX(0.0D0,Kp)*(X*DCOS(Theta_p)+Y*DSIN(Theta_p)))
       dUdx = U_p*DCMPLX(0.0D0,Kp*DCOS(Theta_p))
       dUdy = U_p*DCMPLX(0.0D0,Kp*DSIN(Theta_p))
       GZ_Phi = dcmplx(0.0D0,0.0D0)
       GZ_Phi = (-1.0D0)*(dUdx*n(1) + dUdy*n(2))
       
       !Compute inner product of load with weight functions
       DO J = 1,NODEF
           PTLD_Phi(J) = PTLD_Phi(J) + GZ_Phi*SFN(J)*WTXY            
       endDO
     endDO
   !write(*,*)'-----------------'
   endif

   !edge 3
   if(EFACT(NBC_pos,3).EQ. 1)THEN   
   EDG_TYPE = 1
   !write(*,*)'-----------------'
   !write(*,*) 'ELEMENT = ',IELEM,' EDGE = 3'
     Do IGAUSS = 1,NGAUSS
       XI = XG(IGAUSS)
       WTX = WG(IGAUSS)
       ETA = 1.0
       WTY = 1
       WTXY = WTX*WTY
       !Get shape functions
       call getSF(XI, ETA, SF, NODEL, SFDL, NDIME,NORDER)
       
       !Get Jacobians.
         !Get JCBT = transpose(Jacobian)
         CALL MATMUL(SFDL, NDIME, NODEL, ELCOR, NODEL, NDIME, JCBT,&
         NDIME, NDIME, NDIME, NDIME, NODEL)
       
       !Update weights
         !get Determinant for line integral = dsqrt((dx/dXI)^2 + (dy/dXI)^2)
         DETJ = dsqrt((JCBT(1,1)**2) + (JCBT(1,2)**2))
         !Update weight
         WTXY = WTXY*DETJ        
         !write(*,*) 'DETJ = ',DETJ
       
       !Get (X,Y) from (XI,ETA)
       X = 0.0D0
       Y = 0.0D0
       DO I=1,NODEL
         X = X + SF(I)*ELCOR(I,1)
         Y = Y + SF(I)*ELCOR(I,2)
       ENDDO 
       
       !get Enriched shape functions
         !get SFN
         call getSFN(X,Y,SF,NODEL,K_W,NANGL,BETA,SFN,ENR_type)
       
       !Get normal to edge  
       call getNormVec(EDG_TYPE,JCBT,NDIME,NODEL,X,Y,ELCOR,n,NORDER)
       !write(*,*)'n = ',n  
         
       !Compute load
       U_p = CDEXP(DCMPLX(0.0D0,Kp)*(X*DCOS(Theta_p)+Y*DSIN(Theta_p)))
       dUdx = U_p*DCMPLX(0.0D0,Kp*DCOS(Theta_p))
       dUdy = U_p*DCMPLX(0.0D0,Kp*DSIN(Theta_p))
       GZ_Phi = dcmplx(0.0D0,0.0D0)
       GZ_Phi = (-1.0D0)*(dUdx*n(1) + dUdy*n(2))
       
       !Compute inner product of load with weight functions
       DO J = 1,NODEF
           PTLD_Phi(J) = PTLD_Phi(J) + GZ_Phi*SFN(J)*WTXY            
       endDO
     endDO
   !write(*,*)'-----------------'
   endif
   
   !edge 4
   if(EFACT(NBC_pos,4).EQ. 1)THEN   
   EDG_TYPE = 0
   !write(*,*)'-----------------'
   !write(*,*) 'ELEMENT = ',IELEM,' EDGE = 4'
     Do IGAUSS = 1,NGAUSS
       ETA = XG(IGAUSS)
       WTY = WG(IGAUSS)
       XI = -1.0
       WTX = 1
       WTXY = WTX*WTY
       !Get shape functions
       call getSF(XI, ETA, SF, NODEL, SFDL, NDIME,NORDER)
       !Get Jacobians.
         !Get JCBT = transpose(Jacobian)
         CALL MATMUL(SFDL, NDIME, NODEL, ELCOR, NODEL, NDIME, JCBT,&
         NDIME, NDIME, NDIME, NDIME, NODEL)
         
       
       !Update weights
         !get Determinant for line integral = dsqrt((dx/dETA)^2 + (dy/dETA)^2)
         DETJ = dsqrt((JCBT(2,1)**2) + (JCBT(2,2)**2))
         !Update weight
         WTXY = WTXY*DETJ  
         !write(*,*) 'DETJ = ',DETJ
         
       
       !Get (X,Y) from (XI,ETA)
       X = 0.0D0
       Y = 0.0D0
       DO I=1,NODEL
         X = X + SF(I)*ELCOR(I,1)
         Y = Y + SF(I)*ELCOR(I,2)
       ENDDO 
       
       !get Enriched shape functions
         !get SFN
         call getSFN(X,Y,SF,NODEL,K_W,NANGL,BETA,SFN,ENR_type)
!
       !Get normal to edge
       call getNormVec(EDG_TYPE,JCBT,NDIME,NODEL,X,Y,ELCOR,n,NORDER)
       !write(*,*)'n = ',n  
       
       !Compute load
       U_p = CDEXP(DCMPLX(0.0D0,Kp)*(X*DCOS(Theta_p)+Y*DSIN(Theta_p)))
       dUdx = U_p*DCMPLX(0.0D0,Kp*DCOS(Theta_p))
       dUdy = U_p*DCMPLX(0.0D0,Kp*DSIN(Theta_p))
       GZ_Phi = dcmplx(0.0D0,0.0D0)
       GZ_Phi = (-1.0D0)*(dUdx*n(1) + dUdy*n(2))
       !write(*,*)'GZ_phi = ',GZ_phi
       
       !Compute inner product of load with weight functions
       DO J = 1,NODEF
           PTLD_Phi(J) = PTLD_Phi(J) + GZ_Phi*SFN(J)*WTXY       
           !write(*,*)'PTLD_Phi(J) = ',PTLD_Phi(J)     
       endDO
     endDO
   !write(*,*)'-----------------'
   endif

end SUBROUTINE getHlmhltz_M_S_RHS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
