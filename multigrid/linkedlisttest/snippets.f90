! simple 1 level GS
subroutine mg_solve                                                                                   
    type(grid), pointer :: new_grid
    type(grid), pointer :: current
                        
    real(num) :: L2, L2_old
                        
    L2_old = 1e6_num 
    current => head  
                        
    do               
      call relax(current) 
      call residual(current)
                        
      L2 = sqrt(sum(abs(current%residue)**2)/real(current%nx*current%ny,num))
      if (abs(L2-L2_old) < 1e-12_num) exit ! should replace with user chosen tol eventually
      L2_old = L2       
    enddo               
                     
  end subroutine mg_solve

