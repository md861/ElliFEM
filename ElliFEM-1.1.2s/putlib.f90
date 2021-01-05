!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getDBCMarker(MrkdNds_BC,TOTNODS_BC,TOTCOORD,NMARK,&
         TOTNOD,NonDBCNodes,DBCNodes,DBC_pos)
         
           !Description: This subroutine computes the sorted sequence of
           		!indices that correspond to DBC and non-DBC nodes
           !Arguments out
           		!NonDBCNodes: Vector that stores the 
           		            !indices of non-DBC nodes in a sorted
           		            !fashion.
           		!DBCNodes: Vector that stores the 
           			    !indices of DBC marked nodes in a sorted 
           			    !fashion.

IMPLICIT NONE
integer TOTCOORD,NMARK,TOTNOD,DBC_pos,&
        MrkdNds_BC(NMARK,TOTCOORD),TOTNODS_BC(NMARK),&
        NonDBCNodes(TOTNOD -TOTNODS_BC(DBC_pos)),&
        DBCNodes(TOTNODS_BC(DBC_pos))
integer MapOfIndex(TOTNOD)
        
integer I,I_nDBC,I_DBC

  MapOfIndex = 0
  I_nDBC = 0
  I_DBC = 0
  MapOfIndex(MrkdNds_BC(DBC_pos,1:TOTNODS_BC(DBC_pos))) = -1
  Do I = 1,TOTNOD
    if(MapOfIndex(I).EQ. 0)THEN !NonDBC node encountered
      I_nDBC = I_nDBC + 1
      NonDBCNodes(I_nDBC) = I
    else 		        !DBC node encountered
      I_DBC = I_DBC + 1
      DBCNodes(I_DBC) = I
    endif    
  endDO
  !write(*,*)I_DBC,TOTNODS_BC(DBC_pos)
  !write(*,*)'DBCNodes = ',DBCNodes
  !write(*,*)I_nDBC,TOTNOD-TOTNODS_BC(DBC_pos)
  !write(*,*)'NonDBCNodes = ',NonDBCNodes
end subroutine getDBCMarker
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getPUFEM_DBCMarker(temp_NonDBCNodes,temp_DBCNodes,&
                            TotDBCNODS,TOTNOD,NANGL,&
                            NonDBCNodes,DBCNodes)
                            
           !Description: This subroutine computes the sorted sequence of
           		!indices that correspond to DBC and non-DBC coefficients
           !Arguments out
           		!NonDBCNodes: Vector that stores the 
           		            !indices of non-DBC coefficients.
           		!DBCNodes: Vector that stores the 
           			    !indices of DBC marked coefficients.
           	        !The (sorted) order of input arguments
           	        !*temp_vectors* is preserved.
IMPLICIT NONE                          
integer temp_NonDBCNodes(TOTNOD-TotDBCNODS),&
        temp_DBCNodes(TotDBCNODS),&
        TotDBCNODS,TOTNOD,NANGL,&
        NonDBCNodes((TOTNOD-TotDBCNODS)*NANGL),&
        DBCNodes(TotDBCNODS*NANGL)  
        
integer I,J,counter
  NonDBCNodes = -1
  DBCNodes = -1  
  
  !Populate DBC coefficients
  counter = 0
  Do I = 1,TotDBCNODS
    Do J = 1,NANGL
      counter = counter + 1
      DBCNodes(counter) = (temp_DBCNodes(I) - 1)*NANGL + J
    endDO
  endDO
  
  !Populate NonDBC coefficients
  counter = 0
  Do I = 1,TOTNOD - TotDBCNODS
    Do J = 1,NANGL
      counter = counter + 1
      NonDBCNodes(counter) = (temp_NonDBCNodes(I) - 1)*NANGL + J
    endDO
  endDO
  
        
end subroutine getPUFEM_DBCMarker
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine MergeMarkers(MrkdNds_BC,TOTNODS_BC,TOTCOORD,NMARK,&
		       TOTNOD,MergeFlag,&
		       unMarkedNodes,TotUnMrkdNOD,&
		       MarkedNodes)
         
           !Description: This subroutine computes the sorted sequence of
           		!indices that not correspond to any marked nodes
           		!as referred to by the MergeFlag.
           !Arguments out
           		!unMarkedNodes: Vector that stores the 
           		              !indices of unmarked nodes in a sorted
           		              !fashion.
           		              

IMPLICIT NONE
integer TOTCOORD,NMARK,TOTNOD,&
        MrkdNds_BC(NMARK,TOTCOORD),TOTNODS_BC(NMARK),&
        MergeFlag(NMARK),&
        TotUnMrkdNOD,&
        unMarkedNodes(TOTNOD),&
        MarkedNodes(TOTNOD)
integer MapOfIndex(TOTNOD)
        
integer I,I_unMrkd,I_Mrkd

  MapOfIndex = 0
  I_unMrkd = 0
  I_Mrkd = 0
  Do I = 1,NMARK
    if(MergeFlag(I) .EQ. 1)THEN !Flag is 1 for current BC type
      MapOfIndex(MrkdNds_BC(I,1:TOTNODS_BC(I))) = -1
    endif
  endDO
  Do I = 1,TOTNOD
    if(MapOfIndex(I).EQ. 0)THEN !unMarked node encountered
      I_unMrkd = I_unMrkd + 1
      unMarkedNodes(I_unMrkd) = I
    else 		        !Marked node encountered
      I_Mrkd = I_Mrkd + 1
      MarkedNodes(I_Mrkd) = I
    endif    
  endDO
  TotUnMrkdNOD = I_unMrkd
  !write(*,*)I_unMrkd,TotUnMrkdNOD
  !write(*,*)'unMarkedNodes = ',unMarkedNodes
end subroutine MergeMarkers	       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getELDF(IELEM,GLDF,TOTELS,NODEL,ELNODDAT,NORDER,ELDF,NANGL)
IMPLICIT NONE
integer TOTELS,NODEL,NORDER,NANGL,IELEM
integer GLDF(TOTELS,NODEL),ELNODDAT(NODEL),ELDF(NODEL*NANGL)
integer temp_vec(NODEL)
integer I,J,counter
           !Description: This subroutine computes the sequence of
           		!indices that map from local (sequence)
           		!of co-efficients to global (sequence) of 
           		!co-efficients. This sequence used to assemble
           		!from local to global system, and also for 
           		!post processing.
           !Arguments out
           		!ELDF: Vector that stores the 
           		    !sequence of global indices that 
           		    !correspond to the local sequence
           		    !of indices pertaining to co-efficients.
counter = 0
temp_vec = GLDF(IELEM,ELNODDAT)
DO I = 1,NODEL
  DO J = 1,NANGL
    counter = counter + 1
    ELDF(counter) = (temp_vec(I) - 1)*NANGL + J
  endDO
endDO
!write(*,*)'ELDF = ',ELDF
end subroutine getELDF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE invJCBT(JCBT, NDIME, JCBTI, DETJ)
IMPLICIT NONE
integer NDIME
double precision JCBT(NDIME,NDIME),&
                 JCBTI(NDIME,NDIME),DETJ

  DETJ = JCBT(1,1) * JCBT(2,2) - JCBT(1,2) * JCBT(2,1)
  JCBTI(1,1) = JCBT(2,2) / DETJ
  JCBTI(1,2) = -JCBT(1,2) / DETJ
  JCBTI(2,1) = -JCBT(2,1) / DETJ
  JCBTI(2,2) = JCBT(1,1) / DETJ
end subroutine invJCBT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getNormVec(EDG_TYPE,JCBT,NDIME,NODEL,X,Y,ELCOR,n,NORDER)

           !Description: This subroutine computes the outward normal 
           		!vector at the edge at global point X,Y.
           		!It assumes that 
           		!1)the mid-point of element,i.e.
           		!  the point in the element for XI=ETA=0, does 
           		!  not fall inside the inner-cusp formed by any 
           		!  two adjacent nodes on the edge. See reference 
           		!  document for more clarification.
           		!2)Non-zero area is spanned by the element.
           		!3)Quadrangle elements are present only.
           		!4)NDIME = 2
           		!Method: First calculates 
           		!        m = dy/dx = (dy/dXI) / (dx/dXI), EDG_TYPE = 1
           		!          = (dy/dETA) / (dx/dETA), otherwise
           		!        Then predicts
           		!	 n = (m,-1) {This is normalised}
           		!	 Then calculates
           		! 	 m0= (X0-X,Y0-Y) {This is normalised}
           		!            as direction cosine to (X0,Y0),
           		!	     where (X0,Y0) is the mid-point in the element.
           		!	 Then corrects n as (n) or (-n), based on
           		! 	 <n,m0> which should be negative (i.e. outward normal).
           !Arguments out
           		!n: Vector that stores the 
           		    !directional cosines of the outward
           		    !normal, at the given global point X,Y.
IMPLICIT NONE
integer NDIME,NODEL,EDG_TYPE,NORDER
double precision JCBT(NDIME,NDIME),X,Y,ELCOR(NODEL,NDIME)
double precision m,temp_dbl,n(NDIME),m0(NDIME)
double precision SF_temp(NODEL), SFDL_temp(NDIME,NODEL),&
                 X_temp,Y_temp
integer I

  !Compute n for edges 1 and 3.
  if(EDG_TYPE .EQ. 1)THEN
        !predict n based on m        
          if(JCBT(1,1) .EQ. 0.0D0)THEN !m = infinity
            !predict n
            n = (/1.0D0,0.0D0/)
            !write(*,*)'<debugger>m = infinity'
            !write(*,*)'<debugger>n = ',n
          else
            !get m =  dy/dXI / dx/dXI
            m = JCBT(1,2)/JCBT(1,1)
            temp_dbl = dsqrt(1 + (m**2))
            !predict n
            n = (/ m/temp_dbl, -1.0D0/temp_dbl/)
            !write(*,*)'<debugger>m = ',m
            !write(*,*)'<debugger>n = ',n
          endif
        !correct n based on middle node
            !get middle node
            call getSF(0.0D0, 0.0D0, SF_temp, NODEL, SFDL_temp, NDIME,NORDER)
            X_temp = 0.0D0
            Y_temp = 0.0D0
            DO I=1,NODEL
              X_temp = X_temp + SF_temp(I)*ELCOR(I,1)
              Y_temp = Y_temp + SF_temp(I)*ELCOR(I,2)
            ENDDO 
            !get m0
            if((X_temp - X) .EQ. 0.0D0)THEN 			!m0 = infinity
              if((Y_temp - Y) .GT. 0.0D0)THEN			!m0 = +infinity
                m0 = (/0.0D0,1.0D0/)
              else						!m0 = -infinity, NB:Both dy and dx can't be 0
                m0 = (/0.0D0,-1.0D0/)
              endif
            else
              temp_dbl = (Y_temp - Y)/(X_temp - X)
              m0 = (/ (X_temp - X), (Y_temp - Y) /)
              temp_dbl = dsqrt(((X_temp - X)**2)+((Y_temp - Y)**2))
              m0 = m0/temp_dbl
            endif
            !write(*,*)'<debugger>m0 = ',m0
            !write(*,*)'<debugger>dx = ',(X_temp - X),'dy = ',(Y_temp - Y)
            !update n
            temp_dbl = dot_product(n,m0)
            if(temp_dbl .GT. 0.0D0)THEN	
              n = -n            
            endif	!NB: <n,m0> .NE. 0. Also if <n,m0> .LT. 0 then no need to change n.
            !write(*,*)'<debugger>n = ',n
            
  !Compute n for edges 2 and 4.
  else
        !predict n based on m        
          if(JCBT(2,1) .EQ. 0.0D0)THEN !m = infinity
            !predict n
            n = (/1.0D0,0.0D0/)
            !write(*,*)'<debugger>m = infinity'
            !write(*,*)'<debugger>n = ',n
          else
            !get m =  dy/dETA / dx/dETA
            m = JCBT(2,2)/JCBT(2,1)
            temp_dbl = dsqrt(1 + (m**2))
            !predict n
            n = (/ m/temp_dbl, -1.0D0/temp_dbl/)
            !write(*,*)'<debugger>m = ',m
            !write(*,*)'<debugger>n = ',n
          endif
        !correct n based on middle node
            !get middle node
            call getSF(0.0D0, 0.0D0, SF_temp, NODEL, SFDL_temp, NDIME,NORDER)
            X_temp = 0.0D0
            Y_temp = 0.0D0
            DO I=1,NODEL
              X_temp = X_temp + SF_temp(I)*ELCOR(I,1)
              Y_temp = Y_temp + SF_temp(I)*ELCOR(I,2)
            ENDDO 
            !get m0
            if((X_temp - X) .EQ. 0.0D0)THEN 			!m0 = infinity
              if((Y_temp - Y) .GT. 0.0D0)THEN			!m0 = +infinity
                m0 = (/0.0D0,1.0D0/)
              else						!m0 = -infinity, NB:Both dy and dx can't be 0
                m0 = (/0.0D0,-1.0D0/)
              endif
            else
              temp_dbl = (Y_temp - Y)/(X_temp - X)
              m0 = (/ (X_temp - X), (Y_temp - Y) /)
              temp_dbl = dsqrt(((X_temp - X)**2)+((Y_temp - Y)**2))
              m0 = m0/temp_dbl
            endif
            !write(*,*)'<debugger>m0 = ',m0
            !write(*,*)'<debugger>dx = ',(X_temp - X),'dy = ',(Y_temp - Y)
            !update n
            temp_dbl = dot_product(n,m0)
            if(temp_dbl .GT. 0.0D0)THEN	
              n = -n            
            endif	!NB: <n,m0> .NE. 0. Also if <n,m0> .LT. 0 then no need to change n.
            !write(*,*)'<debugger>n = ',n
  endif
            
end subroutine getNormVec
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MATMUL(A, IA, JA, B, IB, JB, C, IC, JC, L, M, N)
!
!***MATRIX MULTIPLICATION
!***(c) Peter and Jacqueline A. Bettess, 1986
!
!-----------------------------------------------------------------------
! PURPOSE
!      Post multiplies matrix A by matrix B to give matrix C.  All are
!      as 2 dimensional arrays
!
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
      IMPLICIT NONE
      DOUBLE PRECISION A, B, C
      INTEGER IA, IB, IC, IL, IM, IN, JA, JB, JC, L, M, N
      DIMENSION A(IA,JA), B(IB,JB), C(IC,JC)
!
!***process rows in A and C
!
      DO 30 IL = 1, L
!
!***process columns in B and C
!
        DO 20 IM = 1, M
          C(IL,IM) = 0.0
!
!***form inner product
!
          DO 10 IN = 1, N
            C(IL,IM) = C(IL,IM) + A(IL,IN) * B(IN,IM)
   10     CONTINUE
   20   CONTINUE
   30 CONTINUE
      RETURN
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GAULEG(N, X, IX, W, IW)
!                                                                      *
!  SETS UP ARRAY OF N GAUSS POINTS, WITH CORRECT SIGNS FOR ABSCISSA    *
!                                                                      *
!                                                                      *
!  N           NUMBER OF GAUSS POINTS                                  *
!  X(IX)       ARRAY OF ABSCISSAS                                      *
!  W(IW)       ARRAY OF WEIGHTS                                        *
!                                                                      *
  IMPLICIT NONE
  INTEGER IW, IX, N
  DOUBLE PRECISION W(IW), X(IX)
  
  double precision x_dummy(IX),w_dummy(IW)
  integer b,q,r, I
  PARAMETER (b=2)
  
  
  call modquo(N,b,q,r)
  call getINTGRNPNTS(N,x_dummy,IX,w_dummy,IW)
  Do I = 1,r	!This is to update at X = 0 (i.e. the mid-point for odd N)
    X(q+I) = x_dummy(I)
    W(q+I) = w_dummy(I)
  endDO
  Do I = 1,q
    X(q+r+I)= x_dummy(r+I)
    W(q+r+I)= w_dummy(r+I)
    
    X(q-I+1)= -x_dummy(r+I)
    W(q-I+1)= w_dummy(r+I)
  endDo
  
  !write(*,*)'<debugger>X = ',X
  !write(*,*)'<debugger>W = ',W
  
end subroutine GAULEG
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getINTGRNPNTS(N,X,IX,W,IW)
  IMPLICIT NONE
  INTEGER IW, IX, N
  DOUBLE PRECISION W(IW), X(IX)
  
  if(N .EQ. 1)THEN
       X( 1)=   0.00000000000000000000
       W( 1)=   2.00000000000000000000
  elseif(N .EQ. 2)THEN
       X( 1)=   0.57735026918962576451
       W( 1)=   1.00000000000000000000
  elseif(N .EQ. 3)THEN
       X( 1)=   0.00000000000000000000
       X( 2)=   0.77459666924148337704
       W( 1)=   0.88888888888888888889
       W( 2)=   0.55555555555555555556
  elseif(N .EQ. 4)THEN
       X( 1)=   0.33998104358485626480
       X( 2)=   0.86113631159405257522
       W( 1)=   0.65214515486254614263
       W( 2)=   0.34785484513745385737
  elseif(N .EQ. 5)THEN
       X( 1)=   0.00000000000000000000
       X( 2)=   0.53846931010568309104
       X( 3)=   0.90617984593866399280
       W( 1)=   0.56888888888888888889
       W( 2)=   0.47862867049936646804
       W( 3)=   0.23692688505618908751
  elseif(N .EQ. 6)THEN
       X( 1)=   0.23861918608319690863
       X( 2)=   0.66120938646626451366
       X( 3)=   0.93246951420315202781
       W( 1)=   0.46791393457269104739
       W( 2)=   0.36076157304813860757
       W( 3)=   0.17132449237917034504
  elseif(N .EQ. 7)THEN
       X( 1)=   0.00000000000000000000
       X( 2)=   0.40584515137739716691
       X( 3)=   0.74153118559939443986
       X( 4)=   0.94910791234275852453
       W( 1)=   0.41795918367346938776
       W( 2)=   0.38183005050511894495
       W( 3)=   0.27970539148927666790
       W( 4)=   0.12948496616886969327
  elseif(N .EQ. 8)THEN
       X( 1)=   0.18343464249564980494
       X( 2)=   0.52553240991632898582
       X( 3)=   0.79666647741362673959
       X( 4)=   0.96028985649753623168
       W( 1)=   0.36268378337836198297
       W( 2)=   0.31370664587788728734
       W( 3)=   0.22238103445337447054
       W( 4)=   0.10122853629037625915
  elseif(N .EQ. 9)THEN
       X( 1)=   0.00000000000000000000
       X( 2)=   0.32425342340380892904
       X( 3)=   0.61337143270059039731
       X( 4)=   0.83603110732663579430
       X( 5)=   0.96816023950762608984
       W( 1)=   0.33023935500125976316
       W( 2)=   0.31234707704000284007
       W( 3)=   0.26061069640293546232
       W( 4)=   0.18064816069485740406
       W( 5)=   0.08127438836157441197
  elseif(N .EQ. 10)THEN
       X( 1)=   0.14887433898163121088
       X( 2)=   0.43339539412924719080
       X( 3)=   0.67940956829902440623
       X( 4)=   0.86506336668898451073
       X( 5)=   0.97390652851717172008
       W( 1)=   0.29552422471475287017
       W( 2)=   0.26926671930999635509
       W( 3)=   0.21908636251598204400
       W( 4)=   0.14945134915058059315
       W( 5)=   0.06667134430868813759  
  elseif(N .EQ. 11)THEN
       X( 1)=   0.00000000000000000000
       X( 2)=   0.26954315595234497233
       X( 3)=   0.51909612920681181593
       X( 4)=   0.73015200557404932409
       X( 5)=   0.88706259976809529908
       X( 6)=   0.97822865814605699280
       W( 1)=   0.27292508677790063071
       W( 2)=   0.26280454451024666218
       W( 3)=   0.23319376459199047992
       W( 4)=   0.18629021092773425143
       W( 5)=   0.12558036946490462463
       W( 6)=   0.05566856711617366648
  elseif(N .EQ. 12)THEN
       X( 1)=   0.12523340851146891547
       X( 2)=   0.36783149899818019375
       X( 3)=   0.58731795428661744730
       X( 4)=   0.76990267419430468704
       X( 5)=   0.90411725637047485668
       X( 6)=   0.98156063424671925069
       W( 1)=   0.24914704581340278500
       W( 2)=   0.23349253653835480876
       W( 3)=   0.20316742672306592175
       W( 4)=   0.16007832854334622633
       W( 5)=   0.10693932599531843096
       W( 6)=   0.04717533638651182719
  elseif(N .EQ. 13)THEN
       X( 1)=   0.00000000000000000000
       X( 2)=   0.23045831595513479407
       X( 3)=   0.44849275103644685288
       X( 4)=   0.64234933944034022064
       X( 5)=   0.80157809073330991279
       X( 6)=   0.91759839922297796521
       X( 7)=   0.98418305471858814947
       W( 1)=   0.23255155323087391019
       W( 2)=   0.22628318026289723841
       W( 3)=   0.20781604753688850231
       W( 4)=   0.17814598076194573828
       W( 5)=   0.13887351021978723846
       W( 6)=   0.09212149983772844791
       W( 7)=   0.04048400476531587952
  elseif(N .EQ. 14)THEN
       X( 1)=   0.10805494870734366207
       X( 2)=   0.31911236892788976044
       X( 3)=   0.51524863635815409197
       X( 4)=   0.68729290481168547015
       X( 5)=   0.82720131506976499319
       X( 6)=   0.92843488366357351734
       X( 7)=   0.98628380869681233884
       W( 1)=   0.21526385346315779020
       W( 2)=   0.20519846372129560397
       W( 3)=   0.18553839747793781374
       W( 4)=   0.15720316715819353457
       W( 5)=   0.12151857068790318469
       W( 6)=   0.08015808715976020981
       W( 7)=   0.03511946033175186303
  elseif(N .EQ. 15)THEN
       X( 1)=   0.00000000000000000000
       X( 2)=   0.20119409399743452230
       X( 3)=   0.39415134707756336990
       X( 4)=   0.57097217260853884754
       X( 5)=   0.72441773136017004742
       X( 6)=   0.84820658341042721620
       X( 7)=   0.93727339240070590431
       X( 8)=   0.98799251802048542849
       W( 1)=   0.20257824192556127288
       W( 2)=   0.19843148532711157646
       W( 3)=   0.18616100001556221103
       W( 4)=   0.16626920581699393355
       W( 5)=   0.13957067792615431445
       W( 6)=   0.10715922046717193501
       W( 7)=   0.07036604748810812471
       W( 8)=   0.03075324199611726835
  elseif(N .EQ. 16)THEN
       X( 1)=   0.09501250983763744019
       X( 2)=   0.28160355077925891323
       X( 3)=   0.45801677765722738634
       X( 4)=   0.61787624440264374845
       X( 5)=   0.75540440835500303390
       X( 6)=   0.86563120238783174388
       X( 7)=   0.94457502307323257608
       X( 8)=   0.98940093499164993260
       W( 1)=   0.18945061045506849629
       W( 2)=   0.18260341504492358887
       W( 3)=   0.16915651939500253819
       W( 4)=   0.14959598881657673208
       W( 5)=   0.12462897125553387205
       W( 6)=   0.09515851168249278481
       W( 7)=   0.06225352393864789286
       W( 8)=   0.02715245941175409485
  elseif(N .EQ. 17)THEN
       X( 1)=   0.00000000000000000000
       X( 2)=   0.17848418149584785585
       X( 3)=   0.35123176345387631530
       X( 4)=   0.51269053708647696789
       X( 5)=   0.65767115921669076585
       X( 6)=   0.78151400389680140693
       X( 7)=   0.88023915372698590212
       X( 8)=   0.95067552176876776122
       X( 9)=   0.99057547531441733568
       W( 1)=   0.17944647035620652546
       W( 2)=   0.17656270536699264633
       W( 3)=   0.16800410215645004451
       W( 4)=   0.15404576107681028808
       W( 5)=   0.13513636846852547329
       W( 6)=   0.11188384719340397109
       W( 7)=   0.08503614831717918088
       W( 8)=   0.05545952937398720113
       W( 9)=   0.02414830286854793196
  elseif(N .EQ. 18)THEN
       X( 1)=   0.08477501304173530124
       X( 2)=   0.25188622569150550959
       X( 3)=   0.41175116146284264604
       X( 4)=   0.55977083107394753461
       X( 5)=   0.69168704306035320787
       X( 6)=   0.80370495897252311568
       X( 7)=   0.89260246649755573921
       X( 8)=   0.95582394957139775518
       X( 9)=   0.99156516842093094673
       W( 1)=   0.16914238296314359184
       W( 2)=   0.16427648374583272299
       W( 3)=   0.15468467512626524493
       W( 4)=   0.14064291467065065120
       W( 5)=   0.12255520671147846018
       W( 6)=   0.10094204410628716556
       W( 7)=   0.07642573025488905653
       W( 8)=   0.04971454889496979645
       W( 9)=   0.02161601352648331031
  elseif(N .EQ. 19)THEN
       X( 1)=   0.00000000000000000000
       X( 2)=   0.16035864564022537587
       X( 3)=   0.31656409996362983199
       X( 4)=   0.46457074137596094572
       X( 5)=   0.60054530466168102347
       X( 6)=   0.72096617733522937862
       X( 7)=   0.82271465653714282498
       X( 8)=   0.90315590361481790164
       X( 9)=   0.96020815213483003085
       X(10)=   0.99240684384358440319
       W( 1)=   0.16105444984878369598
       W( 2)=   0.15896884339395434765
       W( 3)=   0.15276604206585966678
       W( 4)=   0.14260670217360661178
       W( 5)=   0.12875396253933622768
       W( 6)=   0.11156664554733399472
       W( 7)=   0.09149002162244999946
       W( 8)=   0.06904454273764122658
       W( 9)=   0.04481422676569960033
       W(10)=   0.01946178822972647704
  elseif(N .EQ. 20)THEN
       X( 1)=   0.07652652113349733375
       X( 2)=   0.22778585114164507808
       X( 3)=   0.37370608871541956067
       X( 4)=   0.51086700195082709800
       X( 5)=   0.63605368072651502545
       X( 6)=   0.74633190646015079261
       X( 7)=   0.83911697182221882339
       X( 8)=   0.91223442825132590587
       X( 9)=   0.96397192727791379127
       X(10)=   0.99312859918509492479
       W( 1)=   0.15275338713072585070
       W( 2)=   0.14917298647260374679
       W( 3)=   0.14209610931838205133
       W( 4)=   0.13168863844917662690
       W( 5)=   0.11819453196151841731
       W( 6)=   0.10193011981724043504
       W( 7)=   0.08327674157670474872
       W( 8)=   0.06267204833410906357
       W( 9)=   0.04060142980038694133
       W(10)=   0.01761400713915211831
  elseif(N .EQ. 21)THEN
       X( 1)=   0.00000000000000000000
       X( 2)=   0.14556185416089509094
       X( 3)=   0.28802131680240109660
       X( 4)=   0.42434212020743878357
       X( 5)=   0.55161883588721980706
       X( 6)=   0.66713880419741231931
       X( 7)=   0.76843996347567790862
       X( 8)=   0.85336336458331728365
       X( 9)=   0.92009933415040082879
       X(10)=   0.96722683856630629432
       X(11)=   0.99375217062038950026
       W( 1)=   0.14608113364969042719
       W( 2)=   0.14452440398997005906
       W( 3)=   0.13988739479107315472
       W( 4)=   0.13226893863333746178
       W( 5)=   0.12183141605372853420
       W( 6)=   0.10879729916714837766
       W( 7)=   0.09344442345603386155
       W( 8)=   0.07610011362837930202
       W( 9)=   0.05713442542685720828
       W(10)=   0.03695378977085249380
       W(11)=   0.01601722825777433332
  elseif(N .EQ. 22)THEN
       X( 1)=   0.06973927331972222121
       X( 2)=   0.20786042668822128548
       X( 3)=   0.34193582089208422516
       X( 4)=   0.46935583798675702641
       X( 5)=   0.58764040350691159296
       X( 6)=   0.69448726318668278005
       X( 7)=   0.78781680597920816200
       X( 8)=   0.86581257772030013654
       X( 9)=   0.92695677218717400052
       X(10)=   0.97006049783542872712
       X(11)=   0.99429458548239929207
       W( 1)=   0.13925187285563199338
       W( 2)=   0.13654149834601517135
       W( 3)=   0.13117350478706237073
       W( 4)=   0.12325237681051242429
       W( 5)=   0.11293229608053921839
       W( 6)=   0.10041414444288096493
       W( 7)=   0.08594160621706772741
       W( 8)=   0.06979646842452048809
       W( 9)=   0.05229333515268328594
       W(10)=   0.03377490158481415479
       W(11)=   0.01462799529827220068
  elseif(N .EQ. 23)THEN
       X( 1)=   0.00000000000000000000
       X( 2)=   0.13325682429846611093
       X( 3)=   0.26413568097034493053
       X( 4)=   0.39030103803029083142
       X( 5)=   0.50950147784600754969
       X( 6)=   0.61960987576364615639
       X( 7)=   0.71866136313195019446
       X( 8)=   0.80488840161883989215
       X( 9)=   0.87675235827044166738
       X(10)=   0.93297108682601610235
       X(11)=   0.97254247121811523196
       X(12)=   0.99476933499755212352
       W( 1)=   0.13365457218610617535
       W( 2)=   0.13246203940469661737
       W( 3)=   0.12890572218808214998
       W( 4)=   0.12304908430672953047
       W( 5)=   0.11499664022241136494
       W( 6)=   0.10489209146454141007
       W( 7)=   0.09291576606003514748
       W( 8)=   0.07928141177671895492
       W( 9)=   0.06423242140852585213
       W(10)=   0.04803767173108466857
       W(11)=   0.03098800585697944431
       W(12)=   0.01341185948714177208
  elseif(N .EQ. 24)THEN
       X( 1)=   0.06405689286260562609
       X( 2)=   0.19111886747361630916
       X( 3)=   0.31504267969616337439
       X( 4)=   0.43379350762604513849
       X( 5)=   0.54542147138883953566
       X( 6)=   0.64809365193697556925
       X( 7)=   0.74012419157855436424
       X( 8)=   0.82000198597390292195
       X( 9)=   0.88641552700440103421
       X(10)=   0.93827455200273275852
       X(11)=   0.97472855597130949820
       X(12)=   0.99518721999702136018
       W( 1)=   0.12793819534675215697
       W( 2)=   0.12583745634682829612
       W( 3)=   0.12167047292780339120
       W( 4)=   0.11550566805372560135
       W( 5)=   0.10744427011596563478
       W( 6)=   0.09761865210411388827
       W( 7)=   0.08619016153195327592
       W( 8)=   0.07334648141108030573
       W( 9)=   0.05929858491543678075
       W(10)=   0.04427743881741980617
       W(11)=   0.02853138862893366318
       W(12)=   0.01234122979998719955
  elseif(N .EQ. 25)THEN
       X( 1)=   0.00000000000000000000
       X( 2)=   0.12286469261071039639
       X( 3)=   0.24386688372098843205
       X( 4)=   0.36117230580938783774
       X( 5)=   0.47300273144571496052
       X( 6)=   0.57766293024122296772
       X( 7)=   0.67356636847346836449
       X( 8)=   0.75925926303735763058
       X( 9)=   0.83344262876083400142
       X(10)=   0.89499199787827536885
       X(11)=   0.94297457122897433941
       X(12)=   0.97666392145951751150
       X(13)=   0.99555696979049809791
       W( 1)=   0.12317605372671545120
       W( 2)=   0.12224244299031004169
       W( 3)=   0.11945576353578477223
       W( 4)=   0.11485825914571164834
       W( 5)=   0.10851962447426365312
       W( 6)=   0.10053594906705064420
       W( 7)=   0.09102826198296364981
       W( 8)=   0.08014070033500101801
       W( 9)=   0.06803833381235691721
       W(10)=   0.05490469597583519193
       W(11)=   0.04093915670130631266
       W(12)=   0.02635498661503213726
       W(13)=   0.01139379850102628795
  elseif(N .EQ. 26)THEN
       X( 1)=   0.05923009342931320709
       X( 2)=   0.17685882035689018397
       X( 3)=   0.29200483948595689514
       X( 4)=   0.40305175512348630648
       X( 5)=   0.50844071482450571770
       X( 6)=   0.60669229301761806323
       X( 7)=   0.69642726041995726486
       X( 8)=   0.77638594882067885619
       X( 9)=   0.84544594278849801880
       X(10)=   0.90263786198430707422
       X(11)=   0.94715906666171425014
       X(12)=   0.97838544595647099110
       X(13)=   0.99588570114561692900
       W( 1)=   0.11832141527926227652
       W( 2)=   0.11666044348529658204
       W( 3)=   0.11336181654631966655
       W( 4)=   0.10847184052857659066
       W( 5)=   0.10205916109442542324
       W( 6)=   0.09421380035591414846
       W( 7)=   0.08504589431348523921
       W( 8)=   0.07468414976565974589
       W( 9)=   0.06327404632957483554
       W(10)=   0.05097582529714781200
       W(11)=   0.03796238329436276395
       W(12)=   0.02441785109263190879
       W(13)=   0.01055137261734300716
  elseif(N .EQ. 27)THEN
       X( 1)=   0.00000000000000000000
       X( 2)=   0.11397258560952996693
       X( 3)=   0.22645936543953685886
       X( 4)=   0.33599390363850889973
       X( 5)=   0.44114825175002688059
       X( 6)=   0.54055156457945689490
       X( 7)=   0.63290797194649514093
       X( 8)=   0.71701347373942369929
       X( 9)=   0.79177163907050822714
       X(10)=   0.85620790801829449030
       X(11)=   0.90948232067749110430
       X(12)=   0.95090055781470500685
       X(13)=   0.97992347596150122286
       X(14)=   0.99617926288898856694
       W( 1)=   0.11422086737895698905
       W( 2)=   0.11347634610896514862
       W( 3)=   0.11125248835684519267
       W( 4)=   0.10757828578853318721
       W( 5)=   0.10250163781774579867
       W( 6)=   0.09608872737002850757
       W( 7)=   0.08842315854375695019
       W( 8)=   0.07960486777305777126
       W( 9)=   0.06974882376624559298
       W(10)=   0.05898353685983359911
       W(11)=   0.04744941252061506270
       W(12)=   0.03529705375741971102
       W(13)=   0.02268623159618062320
       W(14)=   0.00979899605129436026
  elseif(N .EQ. 28)THEN
       X( 1)=   0.05507928988403427043
       X( 2)=   0.16456928213338077128
       X( 3)=   0.27206162763517807768
       X( 4)=   0.37625151608907871022
       X( 5)=   0.47587422495511826103
       X( 6)=   0.56972047181140171931
       X( 7)=   0.65665109403886496122
       X( 8)=   0.73561087801363177203
       X( 9)=   0.80564137091717917145
       X(10)=   0.86589252257439504894
       X(11)=   0.91563302639213207387
       X(12)=   0.95425928062893819725
       X(13)=   0.98130316537087275369
       X(14)=   0.99644249757395444995
       W( 1)=   0.11004701301647519628
       W( 2)=   0.10871119225829413525
       W( 3)=   0.10605576592284641791
       W( 4)=   0.10211296757806076981
       W( 5)=   0.09693065799792991585
       W( 6)=   0.09057174439303284094
       W( 7)=   0.08311341722890121839
       W( 8)=   0.07464621423456877902
       W( 9)=   0.06527292396699959579
       W(10)=   0.05510734567571674543
       W(11)=   0.04427293475900422784
       W(12)=   0.03290142778230437998
       W(13)=   0.02113211259277125975
       W(14)=   0.00912428259309451774
  elseif(N .EQ. 29)THEN
       X( 1)=   0.00000000000000000000
       X( 2)=   0.10627823013267923017
       X( 3)=   0.21135228616600107451
       X( 4)=   0.31403163786763993495
       X( 5)=   0.41315288817400866389
       X( 6)=   0.50759295512422764210
       X( 7)=   0.59628179713822782038
       X( 8)=   0.67821453760268651516
       X( 9)=   0.75246285173447713391
       X(10)=   0.81818548761525244499
       X(11)=   0.87463780492010279042
       X(12)=   0.92118023295305878509
       X(13)=   0.95728559577808772580
       X(14)=   0.98254550526141317487
       X(15)=   0.99667944226059658616
       W( 1)=   0.10647938171831424425
       W( 2)=   0.10587615509732094141
       W( 3)=   0.10407331007772937391
       W( 4)=   0.10109127375991496612
       W( 5)=   0.09696383409440860630
       W( 6)=   0.09173775713925876335
       W( 7)=   0.08547225736617252755
       W( 8)=   0.07823832713576378383
       W( 9)=   0.07011793325505127857
       W(10)=   0.06120309065707913854
       W(11)=   0.05159482690249792391
       W(12)=   0.04140206251868283610
       W(13)=   0.03074049220209362264
       W(14)=   0.01973208505612270598
       W(15)=   0.00851690387874640965
  elseif(N .EQ. 30)THEN
       X( 1)=   0.05147184255531769583
       X( 2)=   0.15386991360858354696
       X( 3)=   0.25463692616788984644
       X( 4)=   0.35270472553087811347
       X( 5)=   0.44703376953808917678
       X( 6)=   0.53662414814201989926
       X( 7)=   0.62052618298924286114
       X( 8)=   0.69785049479331579693
       X( 9)=   0.76777743210482619492
       X(10)=   0.82956576238276839744
       X(11)=   0.88256053579205268154
       X(12)=   0.92620004742927432588
       X(13)=   0.96002186496830751222
       X(14)=   0.98366812327974720997
       X(15)=   0.99689348407464954027
       W( 1)=   0.10285265289355884034
       W( 2)=   0.10176238974840550460
       W( 3)=   0.09959342058679526706
       W( 4)=   0.09636873717464425964
       W( 5)=   0.09212252223778612872
       W( 6)=   0.08689978720108297980
       W( 7)=   0.08075589522942021535
       W( 8)=   0.07375597473770520627
       W( 9)=   0.06597422988218049513
       W(10)=   0.05749315621761906648
       W(11)=   0.04840267283059405290
       W(12)=   0.03879919256962704960
       W(13)=   0.02878470788332336935
       W(14)=   0.01846646831109095914
       W(15)=   0.00796819249616660562
  endif

       
end subroutine getINTGRNPNTS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         
  
  
