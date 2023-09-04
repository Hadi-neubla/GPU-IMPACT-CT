!**************************************************************************************************                                           
!* by Hadi Zolfaghari, ARTORG Center, University of Bern, (zolfaghari.haadi@gmail.com)            *      
!* October 2015 - March 2020                                                                      *            
!**************************************************************************************************

  
!===========================================================================================================
  
!=== GPU host allocataions (zolfagha) ======================================================================
 
!===========================================================================================================

allocate(west_ghost(1:N2-1,1:N3-1))
allocate(east_ghost(1:N2-1,1:N3-1))
allocate(front_ghost(1:N1-1,1:N2-1))
allocate(rear_ghost(1:N1-1,1:N2-1))
allocate(lower_ghost(1:N1-1,1:N3-1))
allocate(upper_ghost(1:N1-1,1:N3-1))


allocate(west_chunk(1:N2-1,1:N3-1))
allocate(east_chunk(1:N2-1,1:N3-1))
allocate(front_chunk(1:N1-1,1:N2-1))
allocate(rear_chunk(1:N1-1,1:N2-1))
allocate(lower_chunk(1:N1-1,1:N3-1))
allocate(upper_chunk(1:N1-1,1:N3-1))
