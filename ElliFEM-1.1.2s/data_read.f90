module module_data_read
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getData(NDIME,TOTELS,NODEL,TOTCOORD,NODEDGE,NORDER,NMARK,&
               ELDAT,ofst_ELDAT,COORD,ELEDGDAT,FACT,ELTYPE,GLDF,&
               TOTNOD,MapCoord2Coef,&
               MrkdNds_BC,TOTNODS_BC)
IMPLICIT NONE
!Definitions for Read data starts
integer lngth_txt	!Buffer size to read the texts from input data.
integer offset_ELDAT	!To leave some room in memory for extra element 
			!data (e.g. material properties). By default
			!each element is considered to have same number
			!of nodes and, asumed to be of same material.
integer ETYPE		!Defines the number of edges in an element. E.g
			!ETYPE = 4 for Quadangles, 3 for Triangles, etc.
			
integer ofst_ELDAT	!dummy variable
integer ELTYPE  	!dummy variable
			
PARAMETER (lngth_txt = 500,offset_ELDAT=0,ETYPE=4) 
integer temp_int1,temp_int2 !Buffers for redundant data from input file.

character(len=lngth_txt) temp_txt
integer counter,Iflag,IELEM,INODE,IBC
integer NDIME, NODEL, TOTELS, TOTCOORD, NMARK,&
        NMARK_ELEMS, NODEDGE, NORDER
integer, allocatable :: ELDAT(:,:)
double precision, allocatable :: COORD(:,:)
integer, allocatable :: FACT(:,:,:)   	!FACT(Element#,BoundaryTag,Flag_WhichEdges)
integer, allocatable :: ELEDGDAT(:,:) 	!ELEDGDAT(ETYPE,NODEDGE) stores the ordering 
					!of nodes on each edge of a reference element.		
					
!Definitions for Read data stops
!--------------------------------
!Definitions for build GLDF starts
integer I
integer TOTNOD
integer, allocatable :: GLDF(:,:)	!GLDF is the table that contains element nodes
					!and their connectivities. This would have been
					!un-necessary if the TOTCOORD = TOTNOD.
					
integer, allocatable :: MapCoord2Coef(:)!This is a table (for future/debugging use) that MAPs
					!the CO-ORDinates index to CO-EFficients index.
integer, allocatable :: MapCoef2Coord(:),Map_temp(:)
!Definitions for build GLDF stops

!deifnitions for marking BC nodes
integer, allocatable :: MrkdNds_BC(:,:)   !This matrix stores the indices 
					  !of BC nodes for each boundary type (NMARK)
integer, allocatable :: TOTNODS_BC(:)     !This vector stores TOTal number 
                        		  !of NODeS marked for a given boundary type


ofst_ELDAT = offset_ELDAT
ELTYPE = ETYPE

open(1002,file='logfile.txt')

!Read data_starts
open(1001, file ='dat')
  !Read NDIME
  read(1001,*)temp_txt, NDIME
  write(1002,*)'NDIME	= ',NDIME
  !Read Total number of elements {TOTELS}
  read(1001,*)temp_txt, TOTELS
  write(1002,*)'TOTELS	= ',TOTELS
  
  !Start reading Element Data {ELDAT}
    do IELEM = 1,TOTELS
      !Read the first line to estimate number of NODes in the ELement {NODEL}
      if(IELEM .EQ. 1) THEN
        read(1001,'(A)')temp_txt          
        counter = 1
        Iflag = 1
        NODEL = 0
        do while (counter .LE. LEN_TRIM(temp_txt))
          if (temp_txt(counter:counter) .EQ. '') THEN
  	    if(Iflag .NE. 0) THEN
	      Iflag = 0
              NODEL = NODEL + 1
            endif
	  else
	    Iflag = 1
	  endif
	  counter = counter + 1
        enddo
        NODEL = NODEL - 1
        write(1002,*)'NODEL	= ',NODEL
        NODEDGE = int(sqrt(real(NODEL)))	!NODes per EDGE {NODEDGE}
        write(1002,*)'NODEDGE	= ',NODEDGE
        NORDER = NODEDGE - 1			!order of polynomials
        write(1002,*)'NORDER	= ',NORDER
        !Allocate the correct size of ELDAT array, and populate its first entry.
        allocate(ELDAT(TOTELS, offset_ELDAT + NODEL))        
        read(temp_txt,*)temp_int1,(ELDAT(IELEM,offset_ELDAT + counter),counter=1,NODEL),temp_int2
        ELDAT(IELEM,offset_ELDAT + 1:NODEL) = ELDAT(IELEM, offset_ELDAT + 1:NODEL) + 1
        write(1002,*)IELEM,(ELDAT(IELEM,offset_ELDAT + counter),counter=1,NODEL)
      else
        !Populate the line as an entry in the ELDAT      
        read(1001,*)temp_int1,(ELDAT(IELEM,offset_ELDAT + counter),counter=1,NODEL),temp_int2
        ELDAT(IELEM,offset_ELDAT + 1:NODEL) = ELDAT(IELEM,offset_ELDAT + 1:NODEL) + 1
        write(1002,*)IELEM,(ELDAT(IELEM,offset_ELDAT + counter),counter=1,NODEL)
      endif
    enddo
  !Stop  reading Element Data {ELDAT}
  
  !Read TOTal number of COORDinates {TOTCOORD}
  read(1001,*)temp_txt, TOTCOORD
  write(1002,*)'TOTCOORD	= ',TOTCOORD
  allocate(COORD(TOTCOORD,NDIME))
  !Start reading COORDinate data {COORD}
    do INODE = 1,TOTCOORD
      read(1001,*)(COORD(INODE,counter),counter=1,NDIME),temp_int1
      write(1002,*)INODE,COORD(INODE,1:2)
    enddo
  !Stop reading COORDinate data {COORD}
  
  !Construct ordering of nodes on element edges. Starts
    allocate(ELEDGDAT(ETYPE,NODEDGE))
    call GetELEDGE(ELEDGDAT,ETYPE,NODEDGE,NORDER)
  !Construct ordering of nodes on element edges. Stops
  
  !Read Total Number of MARKed boundaries {NMARK}. 
  read(1001,*)temp_txt, NMARK
  write(1002,*)'NMARK	= ',NMARK
  
  
  !Update the global table for edge boundary conditions. Starts
  allocate(FACT(TOTELS,NMARK,ETYPE))	!Allocate the array for global boundary data.
  FACT = 0				!Reset global edgefactor table.
  allocate(MrkdNds_BC(NMARK,TOTCOORD),TOTNODS_BC(NMARK)) !These matrices/vectors are for marking BC nodes.
  MrkdNds_BC = -1 !-1 represents unchanged 
  TOTNODS_BC = -1 !
    !Read Nodes that belong to the given boundaries. Starts
    do IBC = 1,NMARK
      !Read boundary tag
      read(1001,*)temp_txt, temp_txt
      !Read number of marked edges for the current boundary type (IBC)
      read(1001,*)temp_txt, NMARK_ELEMS
      write(1002,*)'NMARK_ELEMS	= ',NMARK_ELEMS
      !Populate Boundary Markers
      call BndryMark_sub(1001,NMARK_ELEMS,ELEDGDAT,NODEDGE,ETYPE,&
           TOTELS,FACT,NMARK,IBC,ELDAT,offset_ELDAT,NODEL,&
           MrkdNds_BC,TOTNODS_BC,TOTCOORD)
    
    enddo
    !Read Nodes that belong to the given boundaries. Stops
  !Update the global table for edge boundary conditions. Starts
  write(1002,*)'<debugger>-------------FACT	starts'
  do IELEM = 1,TOTELS
    do INODE = 1,NMARK
      write(1002,*)'FACT(',IELEM,',',INODE,',1:ETYPE) = ',FACT(IELEM,INODE,1:ETYPE)
    enddo
  enddo
  write(1002,*)'<debugger>-------------FACT	stops'
  
!Read data_stops
 close(1001) !close dat file

!build GLobal Degree of Freedom table {GLDF} starts
   write(1002,*)'<debugger>-------------GLDF	starts'
   allocate(GLDF(TOTELS,NODEL))
   allocate(MapCoord2Coef(TOTCOORD))
   allocate(Map_temp(TOTCOORD))
   
   MapCoord2Coef = -1 !Initialize the CO-ORDinates index to CO-EFficients index map.
   Map_temp = -1
   counter = 0
   DO IELEM = 1,TOTELS   
     DO I = 1,NODEL
       INODE = MapCoord2Coef(ELDAT(IELEM, ofst_ELDAT + I))
       if(INODE .EQ. -1)THEN !New (uncounted) node encountered
         counter = counter + 1
         MapCoord2Coef(ELDAT(IELEM, ofst_ELDAT + I)) = counter
         GLDF(IELEM,I) = counter
         Map_temp(counter) = ELDAT(IELEM, ofst_ELDAT + I)
       else !Old (Counted) node encountered
         GLDF(IELEM,I) = INODE
       endif
     endDO   
     write(1002,*)IELEM,GLDF(IELEM,1:NODEL)  
   endDO   
   TOTNOD = counter
   write(1002,*)'TOTNOD = ',TOTNOD
   write(1002,*)'<debugger>-------------GLDF	stops'
   allocate(MapCoef2Coord(TOTNOD))
   MapCoef2Coord = Map_temp(1:TOTNOD)
!build GLobal Degree of Freedom table {GLDF} stops
   !Update Marked BC nodes to reflect Coeff indices
   Do IBC = 1,NMARK
     MrkdNds_BC(IBC,1:TOTNODS_BC(IBC)) = &
       MapCoord2Coef(MrkdNds_BC(IBC,1:TOTNODS_BC(IBC)))
   endDO
  close(1002) !close log file
end subroutine getData
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine BndryMark_sub(fid,NMARK_ELEMS,ELEDGDAT,NODEDGE,&
           ETYPE,TOTELS,FACT,NMARK,iNMARK,ELDAT,offset_ELDAT,NODEL,&
           MrkdNds_BC,TOTNODS_BC,TOTCOORD)
           
           !Description: This subroutine reads the marked nodes for 
           		!a given tagged boundary. The marked nodes 
           		!are stored in a non-repeating array, named 
           		!as Marked_Nodes. This array is then used 
           		!to update the FACT table that stores the 
           		!global data for which of the edges in an 
           		!element have what type of tagged boundary
           		!present (as a flag = 1 means present).
           !Arguments out
           		!FACT   :Global table to store the flag bit 
           			!of which boundary condition must 
           			!be applied to a given edge of an
           			!element.
  IMPLICIT NONE
  integer NMARK_ELEMS, NODEDGE, ETYPE, ELEDGDAT,&
          TOTELS,FACT,NMARK, ELDAT, offset_ELDAT, NODEL 
  integer fid,counter,index_loop,INODE,iflag, IBC, IELEM, iNMARK
  integer ,allocatable :: Marked_Nodes(:), temp_int(:),flagNODE(:)
  integer temp_nodes
  DIMENSION ELEDGDAT(ETYPE,NODEDGE), temp_nodes(NODEDGE),&
            FACT(TOTELS,NMARK,ETYPE),&
            ELDAT(TOTELS, offset_ELDAT + NODEL)
            
  !definitions for marking BC nodes
  integer TOTCOORD,MrkdNds_BC(NMARK,TOTCOORD),TOTNODS_BC(NMARK)
  
  
  allocate(Marked_Nodes(NMARK_ELEMS * NODEDGE))
  allocate(temp_int(1+NODEDGE))
  allocate(flagNODE(NODEDGE))  
  counter = 0  
  !Start reading marked nodes for a given boundary tag
  do INODE = 1,NMARK_ELEMS
    flagNODE(1:NODEDGE) = 1	!Reset flags to count them both in
    !Read first line
    read(fid,*)temp_int
    !Set flags to check which marked node to keep (based on if already saved or not)
    do index_loop = 1,counter
      do iflag = 1,NODEDGE
        if(temp_int(iflag+1)==Marked_Nodes(index_loop))THEN
          flagNODE(iflag) = 0
        endif
      enddo
    enddo    
    !Add the Marked nodes if they were not previously saved
    do iflag = 1,NODEDGE
      if(flagNODE(iflag).EQ. 1)THEN
        counter = counter + 1
        Marked_Nodes(counter) = temp_int(iflag+1)
      endif
    enddo
  enddo
  write(1002,*)'<debugger>Total counted marked nodes =',counter
  Marked_Nodes(1:counter) = Marked_Nodes(1:counter) + 1
  write(1002,*)'<debugger>Marked_Nodes =',Marked_Nodes(1:counter)
  !Save these marked indices
  TOTNODS_BC(iNMARK) = counter !These many nodes are marked as type iNMARK
  MrkdNds_BC(iNMARK,1:counter) = Marked_Nodes(1:counter)
  !Stop reading marked nodes for a given boundary tag
  
  !Create edgeFactor data starts
    do IELEM = 1,TOTELS
      !Get edge nodes for each edge
      do IBC = 1,ETYPE
        !Get reference nodes for the current element edge
        temp_nodes = ELEDGDAT(IBC,1:NODEDGE)        
        !Get global nodes for the current element edge
        temp_nodes = ELDAT(IELEM,offset_ELDAT + temp_nodes)
                
        flagNODE = 0	!Reset flag
        !Check if marked for all nodes of the given edge
        do iflag = 1,NODEDGE
          !loop through all marked nodes for the given BC tag.
          do index_loop = 1,counter
            if(temp_nodes(iflag) == Marked_Nodes(index_loop))THEN
              flagNODE(iflag) = 1
            endif
          enddo	!index_loop ends
        enddo	!iflag loop ends
        !Update global edgefactor table if all nodes of given edge were marked
        if(ALL(flagNODE .EQ. 1))THEN
          FACT(IELEM,iNMARK,IBC) = 1
        endif
      enddo 	!IBC loop ends
    enddo	!IELEM loop ends
  !Create edgeFactor data stops
  
  DEALLOCATE (Marked_Nodes,temp_int,flagNODE)  
end subroutine BndryMark_sub
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine GetELEDGE(ELEDGDAT,ETYPE,NODEDGE,NORDER)

           !Description: This subroutine creates the array to store 
           		!local reference nodes for each edge of the 
           		!reference element. The indices of nodes for 
           		!each edge, when dereferenced from the global 
           		!array ELDAT, provide the global co-ordinates 
           		!for the nodes on the respective elements for 
           		!the elements. NB: This ordering is dependent 
           		!on the GMSH. It has been confirmed to be 
           		!constant for a given order of p-FEM 
           		!Quadrilateral based meshes.
           !Arguments out
           		!ELEDGDAT: Table that stores the order of 
           			!indices for all the nodes 
           			!corresponding to a given edge 
           			!of an element. 
           			!The indexing scheme for the 
           			!implemented code follows:
           			!Bottom->Top, Left->right.
           			
  IMPLICIT NONE         			
  integer ELEDGDAT, ETYPE, NODEDGE, NORDER
  DIMENSION ELEDGDAT(ETYPE,NODEDGE)

    if((NORDER) .EQ. 1) THEN!for 1st-order elements
      ELEDGDAT(1,1) = 2
      ELEDGDAT(1,2) = 3
      
      ELEDGDAT(2,1) = 3
      ELEDGDAT(2,2) = 4
      
      ELEDGDAT(3,1) = 4
      ELEDGDAT(3,2) = 1
      
      ELEDGDAT(4,1) = 1
      ELEDGDAT(4,2) = 2
    else if((NORDER) .EQ. 2) THEN!for 2nd-order elements
      ELEDGDAT(1,1) = 2
      ELEDGDAT(1,2) = 6
      ELEDGDAT(1,3) = 3
      
      ELEDGDAT(2,1) = 3
      ELEDGDAT(2,2) = 7
      ELEDGDAT(2,3) = 4
      
      ELEDGDAT(3,1) = 4
      ELEDGDAT(3,2) = 8
      ELEDGDAT(3,3) = 1
      
      ELEDGDAT(4,1) = 1
      ELEDGDAT(4,2) = 5
      ELEDGDAT(4,3) = 2
    else if((NORDER) .EQ. 3) THEN!for 3rd-order elements
      ELEDGDAT(1,1) = 2
      ELEDGDAT(1,2) = 7
      ELEDGDAT(1,3) = 8
      ELEDGDAT(1,4) = 3
      
      ELEDGDAT(2,1) = 3
      ELEDGDAT(2,2) = 9
      ELEDGDAT(2,3) = 10
      ELEDGDAT(2,4) = 4
      
      ELEDGDAT(3,1) = 4 
      ELEDGDAT(3,2) = 11
      ELEDGDAT(3,3) = 12
      ELEDGDAT(3,4) = 1
      
      ELEDGDAT(4,1) = 1
      ELEDGDAT(4,2) = 5
      ELEDGDAT(4,3) = 6
      ELEDGDAT(4,4) = 2
    else if((NORDER) .EQ. 4) THEN!for 4th-order elements
      ELEDGDAT(1,1) = 2
      ELEDGDAT(1,2) = 8
      ELEDGDAT(1,3) = 9
      ELEDGDAT(1,4) = 10
      ELEDGDAT(1,5) = 3
      
      ELEDGDAT(2,1) = 3
      ELEDGDAT(2,2) = 11
      ELEDGDAT(2,3) = 12
      ELEDGDAT(2,4) = 13
      ELEDGDAT(2,5) = 4
      
      ELEDGDAT(3,1) = 4 
      ELEDGDAT(3,2) = 14
      ELEDGDAT(3,3) = 15
      ELEDGDAT(3,4) = 16
      ELEDGDAT(3,5) = 1
      
      ELEDGDAT(4,1) = 1
      ELEDGDAT(4,2) = 5
      ELEDGDAT(4,3) = 6
      ELEDGDAT(4,4) = 7
      ELEDGDAT(4,5) = 2
    endif
end subroutine GetELEDGE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module module_data_read
