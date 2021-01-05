use module_data_read

implicit none
!Definitions
  integer NDIME,TOTELS,NODEL,TOTCOORD,NODEDGE,NORDER,NMARK,NGAUSS
  integer, allocatable :: ELDAT(:,:)
  integer ofst_ELDAT,ELTYPE
  double precision, allocatable :: COORD(:,:),ELCOR(:,:)
  integer, allocatable :: ELEDGDAT(:,:), ELNODDAT(:), ELNDS(:)
  integer, allocatable :: FACT(:,:,:)!,EFACT(:,:)  
  integer IStep, NPOINTPLOT, IELEM, I, J,counter
  
  integer TOTNOD
  integer, allocatable ::  GLDF(:,:),MapCoord2Coef(:),ELDF(:)
  integer NANGL
  !definitions for book-keeping 
  character (len=52) char_name
  integer n_record_IStep
  !definitions for DIC,M,S
  integer TOTDOF
  double precision K_W,Omega,PI
  double complex, allocatable :: ELM(:,:),ELK(:,:),&!Elementary Mass & Stiffness
  				 GM(:,:),GS(:,:),&  !Global Mass & Stiffness
  				 InvGM(:,:),&
  				 PTLD_Phi(:),&!PTLD_Vlcty(:),&
  				 RHS_Phi(:),&!RHS_Vlcty(:),&
  				 Phi(:)!,Vlcty(:)
  !definitions for getRHS
  double precision Cp	!Wave speed
  double precision dt
  double complex, allocatable :: RHS(:)
  !definitions for Implicit formulation
  double complex, allocatable :: PrevPhi(:),&
                                 temp_vec(:),InvLHS(:,:),&
                                 LHS(:,:)
  integer n_istep
  
  !definitions for other postprocessing
  double precision error1,error2
  
  !definitions for PUFEM
  integer NODEF 		!NODEF : NODEL*NANGL
  double precision, allocatable :: BETA(:)
  
  !definitions for marked BC nodes  
  integer, allocatable :: MrkdNds_BC(:,:)   !This matrix stores the indices 
					    !of BC nodes for each boundary type (NMARK)
  integer, allocatable :: TOTNODS_BC(:)     !This vector stores TOTal number 
                        		    !of NODeS marked for a given boundary type
  integer,allocatable :: NonDBCNodes(:),&   !Vector to store indices of Non-DBC coefficients	
  			 DBCNodes(:) 	    !Vector to store indices of DBC marked coefficients.
  integer DBC_pos !This is the position of DBC in NMARK.	
  double complex,allocatable :: LHS_rdcd(:,:),&
    		  		InvLHS_rdcd(:,:),&
    		  		RHS_rdcd(:),&
    		  		Phi_rdcd(:)
  integer TotNonDBCNODS,TotDBCNODS
  integer,allocatable :: temp_NonDBCNodes(:),&   !dummy vector
  			 temp_DBCNodes(:) 	 !dummy vector
    !definitions for DBC2
    integer DBC2_pos
    integer TotNonDBC2NODS,TotDBC2NODS
    integer,allocatable :: temp_NonDBC2Nodes(:),&   
  			   temp_DBC2Nodes(:) 	 
    integer,allocatable :: NonDBC2Nodes(:),&   
  			   DBC2Nodes(:) 	
    double complex, allocatable :: InvGM_DBC2(:,:) !Inverted Mass for the DBC2 type coefficients
  			   
    integer TotUnMrkdNOD    
    integer,allocatable :: MergeFlag(:),&
                           dummy_unMarkedNodes(:),&
                           temp_unMarkedNodes(:),&
                           unMarkedNodes(:),&
                           dummy_MarkedNodes(:),&
                           temp_MarkedNodes(:),&
                           MarkedNodes(:)
                           
    double complex, allocatable :: Phi_DBC2(:),&
                                   RHS_DBC2(:),&
                                   temp_vec_DBC2(:)
                              
  !definitions for ENR_Type
  integer ENR_type		! 1=ENR2, 0=ENR1    
  
  !definitions for Hlmhltz
  integer, allocatable :: EFACT(:,:)
  double precision Kp, Theta_p
  
  !definitions for condition number
  integer info
  double complex, allocatable :: work(:),sys_mat(:,:)
  double precision, allocatable :: rwork(:)
  integer, allocatable :: I_PIV(:)
  double precision row_sum, max_sum,rcond
  logical flag_cndtn
  
  !definitions for time-stamping
  double precision Tt_strt, Ta_strt, Ta_stp, Tp_strt, Tp_stp, Tt_stp,&
  		   Tinv_strt,Tinv_stp, Time_total,Time_assmb,Time_Procs
  
  !definitions for Eigenvalues and determinants
  logical Slct, flag_dtrmnnt
  integer SDIM,lwork_e,info_e
  double complex, allocatable :: EigVals(:),UntryMat(:,:),work_e(:)
  double precision, allocatable :: rwork_e(:)
  character SORT,JOBVS
  logical, allocatable :: bwork_e(:)
  double complex detVal

 
  !Time stamping - Tt_strt marker  
  call cpu_time(Tt_strt)
  !Get input data  
    call getData(NDIME,TOTELS,NODEL,TOTCOORD,NODEDGE,NORDER,NMARK,&
               ELDAT,ofst_ELDAT,COORD,ELEDGDAT,FACT,ELTYPE,GLDF,&
               TOTNOD,MapCoord2Coef,&
               MrkdNds_BC,TOTNODS_BC)  

               
  !book-keeping
  char_name = "case_default"
  call system("mkdir "//char_name)	!Used for storing all results
  call system("mkdir PLOTS")		!Used for plots
  open(1002, file='logfile.txt', status="old", position="append", action="write")
  
               
  !Set-up problem related parameters---------------------------------Starts
  PI = 4.0D0*DATAN(1.0D0)
  K_W = 30.0D0*PI
    !Changes/Addition in parameters for PUFEM
    NGAUSS = 30!10!NGAUSS + CEILING(dsqrt(1.0D0/TOTELS)*K_W*10.0D0/(2.0D0*PI))
    NANGL = 15!9	!Must be replaced to 1 + n_Beta if ENR_type = 1 (ENR2)
    NODEF = NODEL*NANGL
    !Changes for ENR_Type
    ENR_type = 1! 1=ENR2, 0=ENR1
    DO I = 1,NANGL - ENR_type
      if(I.EQ.1)THEN
        allocate (BETA(NANGL-ENR_type))
      endif
      BETA(I) = 2.0d0*PI*dble(i)/(nangl-ENR_type)
      write(*,*)'BETA(',I,') = ',BETA(I)
    endDO
  TOTDOF = TOTNOD*NANGL 
  write(*,*)'NGAUSS = ',NGAUSS,' This must be (h/Lambda)*10'
  write(1002,*)'NGAUSS = ',NGAUSS  
  write(1002,*)'NANGL = ',NANGL,' is 1 for FEM'   
  write(1002,*)'TOTDOF = ',TOTDOF
  Omega = 2.0D0*PI
  Cp = Omega/K_W
  dt = 0.001  !this is redundant variable
  ! n_istep = 2500
  NPOINTPLOT = 30	!# of integration points for plots
  !For Hlmhltz problem
  Kp = K_W!1.0D0*PI
  Theta_p = PI/(3.0D0)
  !For condition number
  flag_cndtn = .FALSE.!.TRUE. 
  !For Eigenvalues and determinants
  flag_dtrmnnt = .FALSE.!.TRUE.
  !Set-up problem related parameters---------------------------------Stops
  
  !Build marked indices for DBC--------------------------------------Starts
    DBC_pos = 2				!NB: DBC at pos = 2 is always zero Dirichlet 
    					!condition. (although this is changeable as 
    					!long as we are consistent in GMSH and here 
    					!for DBC_pos
    TotDBCNODS = TOTNODS_BC(DBC_pos)	!# of DBC nodes p-FEM
    TotNonDBCNODS = TOTNOD - TotDBCNODS	!# of nonDBC nodes p-FEM    
    allocate(temp_NonDBCNodes(TotNonDBCNODS))
    allocate(temp_DBCNodes(TotDBCNODS))
    allocate(NonDBCNodes(TotNonDBCNODS*NANGL))
    allocate(DBCNodes(TotDBCNODS*NANGL))
    call getDBCMarker(MrkdNds_BC,TOTNODS_BC,TOTCOORD,NMARK,&
         TOTNOD,temp_NonDBCNodes,temp_DBCNodes,DBC_pos)
      !write(*,*)'temp_DBCNodes = ',temp_DBCNodes
      !write(*,*)'temp_NonDBCNodes = ',temp_NonDBCNodes
    !update the marked indices for PUFEM coefficients.
    call getPUFEM_DBCMarker(temp_NonDBCNodes,temp_DBCNodes,&
                            TotDBCNODS,TOTNOD,NANGL,&
                            NonDBCNodes,DBCNodes)
      !write(*,*)'DBCNodes = ',DBCNodes
      !write(*,*)'NonDBCNodes = ',NonDBCNodes
    !update # of DBC/NonDBC coefficients
    TotDBCNODS = TotDBCNODS*NANGL
    TotNonDBCNODS = TotNonDBCNODS*NANGL
  !Build marked indices for DBC--------------------------------------Stops
    !get unmarked nodes
      allocate(MergeFlag(NMARK))
      allocate(dummy_unMarkedNodes(TOTNOD))
      allocate(dummy_MarkedNodes(TOTNOD))
      MergeFlag = 0
      call MergeMarkers(MrkdNds_BC,TOTNODS_BC,TOTCOORD,NMARK,&
        	         TOTNOD,MergeFlag,&
      		         dummy_unMarkedNodes,TotUnMrkdNOD,&
      		         dummy_MarkedNodes)
      allocate(temp_unMarkedNodes(TotUnMrkdNOD))
      allocate(temp_MarkedNodes(TOTNOD-TotUnMrkdNOD))
      temp_unMarkedNodes = dummy_unMarkedNodes(1:TotUnMrkdNOD)
      temp_MarkedNodes = dummy_MarkedNodes(1:(TOTNOD-TotUnMrkdNOD))
      write(*,*)'TotUnMrkdNOD = ',TotUnMrkdNOD
      !write(*,*)'temp_unMarkedNodes = ',temp_unMarkedNodes
      !write(*,*)'temp_MarkedNodes = ',temp_MarkedNodes
      !update unmarked indices for PUFEM coefficients
      allocate(unMarkedNodes(TotUnMrkdNOD*NANGL))
      allocate(MarkedNodes((TOTNOD-TotUnMrkdNOD)*NANGL))
      call getPUFEM_DBCMarker(temp_MarkedNodes,temp_unMarkedNodes,&
                            TotUnMrkdNOD,TOTNOD,NANGL,&
                            MarkedNodes,unMarkedNodes)
      !write(*,*)'unMarkedNodes = ',unMarkedNodes
      !update # of unMarked coefficients
      TotUnMrkdNOD = TotUnMrkdNOD*NANGL
      
    !!Allocate dummy vectors to store LHS(DBC)*Phi(DBC)  
    !allocate(temp_vec_DBC2(TotUnMrkdNOD))
  !Build marked indices for DBC2--------------------------------------Stops
  
  !allocate variables--------------------------------------------Starts
    allocate(ELNODDAT(NODEL),ELNDS(NODEL),ELCOR(NODEL,NDIME))
    allocate(EFACT(NMARK,ELTYPE))
    allocate(ELDF(NODEF))
    !allocate for DIC,M,S
    allocate (ELM(NODEF,NODEF),ELK(NODEF,NODEF))
    allocate (PTLD_Phi(NODEF))!,PTLD_Vlcty(NODEF))
    allocate (GM(TOTDOF,TOTDOF),InvGM(TOTDOF,TOTDOF))
    allocate (GS(TOTDOF,TOTDOF))
    allocate (RHS_Phi(TOTDOF))!,RHS_Vlcty(TOTDOF))
    allocate (Phi(TOTDOF))!,Vlcty(TOTDOF))
    !allocate for getRHS
    allocate (RHS(TOTDOF))
    !allocate for Implicit formulation
    allocate (PrevPhi(TOTDOF))
    allocate (temp_vec(TOTDOF))
    allocate (LHS(TOTDOF,TOTDOF))
    allocate (InvLHS(TOTDOF,TOTDOF))
    !allocate for DBC
      !none
  !allocate variables--------------------------------------------Stops
  
  !~~~~~~~~~~~~~~~~~~~~Helmholtz_M,S,RHS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!Starts
  write(*,*)'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  write(*,*)'Processing Helmholtz_M,S,RHS'
  write(*,*)'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'  
  !Initialize matrices and vectors
  GM = dcmplx(0.0D0,0.0D0)
  GS = dcmplx(0.0D0,0.0D0)
  RHS_Phi = dcmplx(0.0D0,0.0D0)
  !RHS_Vlcty = dcmplx(0.0D0,0.0D0)
  Phi = dcmplx(0.0D0,0.0D0)
  !Vlcty = dcmplx(0.0D0,0.0D0)
  !Assemble
  call cpu_time(Ta_strt)
  Do IELEM = 1,TOTELS
    !get element co-ordinates
    call GetELNODE(ELNODDAT,NODEL,NORDER)
    ELNDS = ELDAT(IELEM, ofst_ELDAT + ELNODDAT)
    DO I = 1,NODEL
      DO J = 1,NDIME
        ELCOR(I,J) = COORD(ELNDS(I),J)
      endDO
    endDO
    
    !get element edge-factors
    !write(*,*) 'IELEM = ',IELEM
    !write(*,*) 'Elcord = ',ELCOR(:,:)
    DO I = 1,NMARK
      EFACT(I,1:ELTYPE)=FACT(IELEM,I,1:ELTYPE)
      !write(*,*) ' EFACT (I,:) = ',EFACT(I,1:ELTYPE)
    endDO
    
        
    !Build local system    
    call getHlmhltz_M_S_RHS(IELEM,ELCOR,NODEL,NDIME,&
		   NGAUSS,NORDER,K_W,Omega,&
		   ELM,ELK,PTLD_Phi,&!PTLD_Vlcty,&
		   NANGL,NODEF,BETA,&
		   ENR_Type,&
		   EFACT,NMARK,ELTYPE,&
		   Kp,Theta_p)
		       
    !Build global system
      !get indices that map from ELememnt Degree of Freedom {ELDF} to GLobal Degree of Freedom {GLDF table}
      call getELDF(IELEM,GLDF,TOTELS,NODEL,ELNODDAT,NORDER,ELDF,NANGL)
      !transfer ELementary matrices/vectors to Global matrices/vectors
      DO I = 1,NODEF
        DO J = 1,NODEF
          GM(ELDF(I),ELDF(J)) = GM(ELDF(I),ELDF(J)) + ELM(I,J)
          GS(ELDF(I),ELDF(J)) = GS(ELDF(I),ELDF(J)) + ELK(I,J)
        endDO
        RHS_Phi(ELDF(I)) = RHS_Phi(ELDF(I)) + PTLD_Phi(I)
        !RHS_Vlcty(ELDF(I)) = RHS_Vlcty(ELDF(I)) + PTLD_Vlcty(I)
      endDO    
  endDO  
  call cpu_time(Ta_stp)  

  write(*,*)'~~~~~~~~~~~~~~~~~~~'
  write(*,*)'LHS for Helmholtz'
  write(*,*)'~~~~~~~~~~~~~~~~~~~'  
  
  !!Update Stiffness for wave problem
  !GS = (Cp*Cp)*GS
  
  !!Compute LHS for implicit (Backward-Euler)
  !  !LHS = [M] + dt*dt*[S]
  !  LHS = GM + (dt*dt*GS)
  
  !Compute LHS for Hlmhltz
    !LHS = k^2[M] - [S]
    LHS = (Kp*Kp*GM) - GS
    
    !Changes for DBC LHS
      !allocate DBC matrices/vectors            
      allocate(LHS_rdcd(TotUnMrkdNOD,TotUnMrkdNOD))
      allocate(InvLHS_rdcd(TotUnMrkdNOD,TotUnMrkdNOD))
      allocate(RHS_rdcd(TotUnMrkdNOD))
      allocate(Phi_rdcd(TotUnMrkdNOD))        
	
      !call inv(LHS,TOTDOF,InvLHS,TOTDOF)	!The bigger system is no longer valid.
      !Reduce LHS
      LHS_rdcd = LHS(unMarkedNodes,unMarkedNodes)
      !invert reduced system
      call cpu_time(Tinv_strt)
      call inv(LHS_rdcd,TotUnMrkdNOD,InvLHS_rdcd,TotUnMrkdNOD)
      call cpu_time(Tinv_stp)
      	 !Compute condition number
      	 if(flag_cndtn .eqv. .TRUE.)THEN
      	   write(*,*)'~~~~~~~~~~~~'
      	   write(*,*)'Condition # '
      	   !allocate for condition #
      	     allocate(sys_mat(TotUnMrkdNOD,TotUnMrkdNOD))
             allocate(I_PIV(TotUnMrkdNOD))
	           allocate(work(2*TotUnMrkdNOD),rwork(2*TotUnMrkdNOD))
      	   !get LU decomposition !NB: the system matrix is copied, so LHS is unchanged.
      	   sys_mat = LHS_rdcd
      	   call zgetrf(TotUnMrkdNOD,TotUnMrkdNOD,sys_mat,TotUnMrkdNOD,&
      	               I_PIV, info)
      	   !get 1-norm
      	   Do J = 1,TotUnMrkdNOD
            row_sum = 0.0D0
            Do I = 1,TotUnMrkdNOD
              row_sum = row_sum + cdabs(sys_mat(I,J))
            endDO
            if(row_sum .GT. max_sum)THEN
              max_sum = row_sum
            endif
           endDO
      	   !get reciprocal condition number
      	   call zgecon( '1', TotUnMrkdNOD, sys_mat, TotUnMrkdNOD,& 
      	   	       max_sum, rcond, work, rwork, info )
      	   write(*,*)'rcond = ',rcond  
      	   open(1120, file = 'rcond.txt')
      	   write(1120,*)rcond
      	   close(1120)
      	   write(*,*)'~~~~~~~~~~~~'
      	 endif !if to calculate conditioning
      	 
      	 !Compute Eigenvalues and determinant
      	 if(flag_dtrmnnt .eqv. .TRUE.)THEN
      	   write(*,*)'~~~~~~~~~~~~'
      	   write(*,*)'det(sys) calculations '
      	   !allocate for Eigenvalues and determinant
	          allocate(EigVals(TotUnMrkdNOD))
	         !allocate(UntryMat(TotUnMrkdNOD,TotUnMrkdNOD))
            allocate(work_e(2*TotUnMrkdNOD),rwork_e(TotUnMrkdNOD))
            allocate(bwork_e(TotUnMrkdNOD))
      	   JOBVS = 'N' !V = compute Schur form. N = don't
      	   	       !NB: UntryMat is not referenced for Jobvs=N
	         SORT = 'N'  !N = don't, S= do sort EigenVals
           SDIM = 0    !0 if sort=N, otherwise =# eigenvals for which slct returns TRUE if sort=S
           lwork_e = 2*TotUnMrkdNOD !Has to be >=max(1,2*n_dim)
	   	                             !NB: Slct = not referenced for sort=N 
      	   call zgees( JOBVS, SORT, Slct, TotUnMrkdNOD, LHS_rdcd,& 
      	   	          TotUnMrkdNOD, SDIM, EigVals, UntryMat,&
                      TotUnMrkdNOD, WORK_e, LWORK_e, RWORK_e,&
                      BWORK_e, INFO_e )
           detVal = Product(EigVals) 
           write(*,*) '{det(sys)} = ',detVal
           write(*,*) 'abs{det(sys)} = ',abs(detVal)
           open(1123, file = 'det_sys.txt')
      	   write(1123,*)'{det(sys)} = ',detVal
           write(1123,*)'abs{det(sys)} = ',abs(detVal)
      	   close(1123)
           write(*,*)'~~~~~~~~~~~~'
         endif!<- if to calculate eigen-stuff  
  write(*,*)'~~~~~~~~~~~~~~~~~~~'
  write(*,*)'LHS (Helmholtz) Reduced to implement DBC'
  write(*,*)'~~~~~~~~~~~~~~~~~~~' 
  write(*,*)'~~~~~~~~~~~~~~~~~~~'
  write(*,*)'Processing solution vector'
  write(*,*)'~~~~~~~~~~~~~~~~~~~' 
  
  !Compute DBC co-efficients
    !It is assumed Phi_DBC1 = 0
    !Compute Phi_DBC2 
  !Solver for Phi (and Vlcty)    
    !compute Phi(t)
      !Changes for DBC
        !Reduce RHS
        RHS_rdcd = RHS_Phi(unMarkedNodes)
        !Update RHS(reduced) = RHS(reduced) - LHS(DBC2)Phi(DBC2)
                               
        !Compute reduced Solution
        call cpu_time(Tp_strt)
        call MATMULCPLX(InvLHS_rdcd,TotUnMrkdNOD,TotUnMrkdNOD,&
                        RHS_rdcd,TotUnMrkdNOD,1,&
                        Phi_rdcd,TotUnMrkdNOD,1,&
                        TotUnMrkdNOD,1,TotUnMrkdNOD)

        !Update total solution
        Phi = dcmplx(0.0D0,0.0D0)	!Zero DBC1
        !Phi(DBC2Nodes) = Phi_DBC2       !Non-Zero DBC2
        Phi(unMarkedNodes) = Phi_rdcd
        call cpu_time(Tp_stp)
        
  !Post Processing     
      ISTEP = 1!This is redundant for Hlmhltz
      ! !Errors
      ! call getLn_norm(IStep,TOTELS,ELDAT,COORD,ofst_ELDAT,NODEL,&
      !              TOTCOORD,NDIME,NORDER,NANGL,TOTDOF,Phi,GLDF,&
      !              ELCOR,NGAUSS,&
      !              K_W,Omega,dt,Cp,error1,error2,&
      !              NODEF,BETA,ENR_type,&
		  !  Kp,Theta_p)
      ! write(*,*)'L2 error % =',error2*100
      ! write(*,*)'L1 error % =',error1*100
      !Plots            
      !n_record_IStep = CEILING(dble(n_istep)/200.0D0)
      !if(((IStep .EQ. 1).OR.(mod(IStep,n_record_IStep).EQ. 0)).OR.(IStep .EQ. n_istep))Then
        call plotPara(IStep,NPOINTPLOT,TOTELS,ELDAT,COORD,ofst_ELDAT,NODEL,&
                   TOTCOORD,NDIME,NORDER,NANGL,TOTDOF,Phi,GLDF,&!Vlcty,&
                   K_W,Omega,dt,&
                   NODEF,BETA,ENR_type,&
      		   Kp,Theta_p)
      !endif 
        
        
  !~~~~~~~~~~~~~~~~~~~~Helmholtz_M,S,RHS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!Stops  

  call cpu_time(Tt_stp)
  !Total time
  Time_total = Tt_stp - Tt_strt
  !Assembly time
  Time_assmb = Ta_stp - Ta_strt
  !Processing time
  Time_Procs = (Tinv_stp - Tinv_strt)+(Tp_stp - Tp_strt)
  write(*,*)'Total time = ',Time_total
  write(*,*)'Assmb time = ',Time_assmb
  write(*,*)'Procs time = ',Time_Procs
  open(1220, file = 'Time.txt')
  write(1220,*)'Total time = ',Time_total
  write(1220,*)'Assmb time = ',Time_assmb
  write(1220,*)'Procs time = ',Time_Procs
  close(1220)
  
  !book-keeping	
  call system("cp dat "//char_name)
  call system("cp -r PLOTS "//char_name)
  call system("cp error_1_data "//char_name)
  call system("cp error_2_data "//char_name)
  call system("cp logfile.txt "//char_name)
  call system("cp rcond.txt "//char_name)
  call system("cp Time.txt "//char_name)
  call system("cp det_sys.txt "//char_name)
  !call system("rm dat")
  call system("rm -r PLOTS")
  call system("rm error_1_data")
  call system("rm error_2_data")
  call system("rm logfile.txt")
  call system("rm rcond.txt")
  call system("rm Time.txt")
  call system("rm det_sys.txt")
                   
END
