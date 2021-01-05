!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine plotPara(IStep,NPOINTPLOT,TOTELS,ELDAT,COORD,ofst_ELDAT,NODEL,&
                   TOTCOORD,NDIME,NORDER,NANGL,TOTDOF,Phi,GLDF,&!Vlcty,&
                   K_W,Omega,dt,&
                   NODEF,BETA,&
                   ENR_Type,&
		   Kp,Theta_p)
  
  IMPLICIT NONE
  integer IStep, NPOINTPLOT, TOTELS, ELDAT,ofst_ELDAT,NODEL,TOTCOORD,NDIME,NORDER
  double precision COORD
  integer IELEM, P, Q, I
  character (len=20) temp_char
  
  double precision XI,ETA,SF, SFDL, ELCOR, X, Y
  integer ELNODDAT, ELNDS
  
  DIMENSION ELDAT(TOTELS, ofst_ELDAT + NODEL), COORD(TOTCOORD,NDIME)
  
  DIMENSION SF(NODEL), SFDL(NDIME,NODEL), ELNODDAT(NODEL), ELNDS(NODEL),&
            ELCOR(NODEL,NDIME)
  !definitions for DIC,M,S
  integer TOTDOF,NANGL,NODEF,&
          GLDF(TOTELS,NODEL),ELDF(NODEL*NANGL)
  double complex Phi(TOTDOF),&!Vlcty(TOTDOF),&
                 PhiNum,PhiExa!,VlctyNum
  !definitions for Implicit
  double precision r_mod,K_W,Omega,T,dt
  
  !definitions for PUFEM
  double complex SFN(NODEF)
  
    !definitions for ENR_Type
    integer ENR_type		! 1=ENR2, 0=ENR1
    !changes for ENR_Type
    double precision BETA(NANGL-ENR_type)
  
  !definitions for Hlmhltz
  double precision Kp,Theta_p
  
  write(temp_char,*)IStep
  open(17,file='PLOTS/pot-Re-Im'//trim(adjustl(temp_char))//'.dat',form='formatted')
  write(17,*)'VARIABLES = "X" , "Y" , "ReNum" , "ImNum" , "ReAna" ,&
      "ImAna" '
  

  T = IStep*dt
  
  DO IELEM = 1,TOTELS
    call GetELNODE(ELNODDAT,NODEL,NORDER)
    ELNDS = ELDAT(IELEM, ofst_ELDAT + ELNODDAT)
    !write(*,*)ELNDS
    DO P = 1,NODEL
      DO Q = 1,NDIME
        ELCOR(P,Q) = COORD(ELNDS(P),Q)
      endDO
    endDO
    !write(*,*)ELCOR(:,1)
    !write(*,*)ELCOR(:,2)
    !write(*,*)COORD(ELNDS,1)
    !write(*,*)COORD(ELNDS,2)
    
    call getELDF(IELEM,GLDF,TOTELS,NODEL,ELNODDAT,NORDER,ELDF,NANGL)
    !write(*,*)IELEM,'ELDF1 = ',ELDF
     
    write(17,*)'ZONE T="ZONE',IELEM,'" I=',NPOINTPLOT+1,' J=',NPOINTPLOT+1,&
    ' F=POINT'
    
    DO P = 1,NPOINTPLOT+1
      XI = -1.0d0+(DBLE(P-1)/DBLE(NPOINTPLOT))*2.0d0
      DO Q = 1,NPOINTPLOT+1
        ETA = -1.0d0+(DBLE(Q-1)/DBLE(NPOINTPLOT))*2.0d0
        
        call getSF(XI, ETA, SF, NODEL, SFDL, NDIME,NORDER)
        X = 0.0D0
        Y = 0.0D0
        DO I=1,NODEL
          X = X + SF(I)*ELCOR(I,1)
          Y = Y + SF(I)*ELCOR(I,2)
        endDO
        
        !get Enriched shape functions
          !get SFN
          call getSFN(X,Y,SF,NODEL,K_W,NANGL,BETA,SFN,ENR_Type)
        
        PhiNum = dcmplx(0.0D0,0.0D0)
        !VlctyNum = dcmplx(0.0D0,0.0D0)
        !Interpolate numerical
        Do I=1,NODEF
          PhiNum = PhiNum + Phi(ELDF(I))*SFN(I)
          !VlctyNum = VlctyNum + Vlcty(ELDF(I))*SFN(I)
        endDO
        !Calculate exact
        PhiExa = dcmplx(0.0D0,0.0D0)
        r_mod = sqrt(X*X + Y*Y)
        PhiExa = CDEXP(DCMPLX(0.0D0,Kp)*&
        	(X*DCOS(Theta_p)+Y*DSIN(Theta_p)))
        !PhiExa = CDEXP(DCMPLX(0.0D0,1)*(K_W*r_mod&
	!		-Omega*T))
        
        write(17,201)x,y,dreal(PhiNum),dimag(PhiNum),&
                     dreal(PhiExa),dimag(PhiExa)
        
      endDO
    endDO
  endDO
  
  201 FORMAT (7(E14.5,2X))
  close(17)
  
end subroutine plotPara
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
