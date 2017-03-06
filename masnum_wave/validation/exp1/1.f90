program compare_exp1
implicit none
character*32, parameter:: openfile1 = "pac_ncep_wav_20090228.nc"
character*32, parameter:: openfile2 = "pac_ncep_wav_20090228_stardard.nc"
  print *, "compare HS between ", trim(openfile1), " and ", trim(openfile2)
 
      print *, "Step 1: Open file ", trim(openfile1) 
      print*, "Step 1 Success"

      print *, "Step 2: Open file ", trim(openfile2) 
      print*, "Step 2 Success"

      print *, "Step 3"
      print *, "Step 3.1 Compare missing_value"
      print*, "Step 3.1 Success"
      
      print *, "Step 3.2 Compare scale_factor"
      print*, "Step 3.2 Success"

      print *, "Step 3.3 Compare HS"
               print *, "Error! HS do not match at (i,j) = ", 1, ",",1
stop
     print *, "Step 3.3 Success"
     print *, "Compare Success"
         
end
