program compare_exp1
implicit none

  include 'netcdf.inc'

  integer, parameter :: IM = 307
  integer, parameter :: JM = 301
  integer, parameter :: year_i = 2009 
  integer, parameter :: month_i = 2 
  integer, parameter :: day_i = 28
  real, parameter :: eps = 1e-36 
  real, parameter :: eps_r = 1e-2 
  character * 32,parameter :: prefix="pac_ncep_wav_"
  integer * 4 :: rcode
  integer :: status
  integer :: NCID
  integer :: VLONID, VLATID
  integer :: VHSID
  integer :: vland_1, vland_2
  real :: scale_f_1, scale_f_2

  integer :: hs_i_1(IM, JM), hs_i_2(IM, JM)
  real :: hs_r_1(IM, JM), hs_r_2(IM, JM)
  real :: lon_1(IM), lat_1(JM)
  real :: lon_2(IM), lat_2(JM)

  integer :: i, j
  integer :: yy,mon
  character * 4 :: year_s
  character * 2 :: month_s
  character * 2 :: day_s 
  character * 256 :: openfile1, openfile2

  write(year_s,'(i4.4)') year_i
  write(month_s,'(i2.2)') month_i
  write(day_s,'(i2.2)') day_i
  openfile1 = "../../exp/exp1/"//trim(prefix)//year_s//month_s//day_s//".nc"
  openfile2 = "./"//trim(prefix)//year_s//month_s//day_s//"_stardard.nc"
  print *, "compare HS between ", trim(openfile1), " and ", trim(openfile2)
 
      print *, "Step 1: Open file ", trim(openfile1) 
      status = NF_OPEN (trim(openfile1), NF_NOWRITE, NCID)
      if (status .NE. NF_NOERR) CALL Handle_err(status)

      status = NF_INQ_VARID (NCID, 'lon', VLONID)
      if (status .NE. NF_NOERR) CALL Handle_err(status)

      status = NF_INQ_VARID (NCID, 'lat', VLATID)
      if (status .NE. NF_NOERR) CALL Handle_err(status)

      status = NF_INQ_VARID (NCID, 'hs', VHSID)
      if (status .NE. NF_NOERR) CALL Handle_err(status)

      status = NF_GET_ATT_REAL(NCID, VHSID, 'missing_value',vland_1)
      if (status .NE. NF_NOERR) CALL Handle_err(status)

      status = NF_GET_ATT_REAL(NCID, VHSID, 'scale_factor',scale_f_1)
      if (status .NE. NF_NOERR) CALL Handle_err(status)

      status = NF_GET_VAR_REAL (NCID, VLONID, lon_1)
      if (status .NE. NF_NOERR) CALL Handle_err(status)

      status = NF_GET_VAR_REAL (NCID, VLATID, lat_1)
      if (status .NE. NF_NOERR) CALL Handle_err(status)

      status = NF_GET_VAR_INT (NCID, VHSID, hs_i_1)
      if (status .NE. NF_NOERR) CALL Handle_err(status)

      status = NF_CLOSE(NCID)
      if (status .NE. NF_NOERR) CALL Handle_err(status)
      print*, "Step 1 Success"

      print *, "Step 2: Open file ", trim(openfile2) 
      status = NF_OPEN (trim(openfile2), NF_NOWRITE, NCID)
      if (status .NE. NF_NOERR) CALL Handle_err(status)

      status = NF_INQ_VARID (NCID, 'lon', VLONID)
      if (status .NE. NF_NOERR) CALL Handle_err(status)

      status = NF_INQ_VARID (NCID, 'lat', VLATID)
      if (status .NE. NF_NOERR) CALL Handle_err(status)

      status = NF_INQ_VARID (NCID, 'hs', VHSID)
      if (status .NE. NF_NOERR) CALL Handle_err(status)

      status = NF_GET_ATT_REAL(NCID, VHSID, 'missing_value',vland_2)
      if (status .NE. NF_NOERR) CALL Handle_err(status)

      status = NF_GET_ATT_REAL(NCID, VHSID, 'scale_factor',scale_f_2)
      if (status .NE. NF_NOERR) CALL Handle_err(status)

      status = NF_GET_VAR_REAL (NCID, VLONID, lon_2)
      if (status .NE. NF_NOERR) CALL Handle_err(status)

      status = NF_GET_VAR_REAL (NCID, VLATID, lat_2)
      if (status .NE. NF_NOERR) CALL Handle_err(status)

      status = NF_GET_VAR_INT (NCID, VHSID, hs_i_2)
      if (status .NE. NF_NOERR) CALL Handle_err(status)

      status = NF_CLOSE(NCID)
      if (status .NE. NF_NOERR) CALL Handle_err(status)
      print*, "Step 2 Success"

      print *, "Step 3"
      print *, "Step 3.1 Compare missing_value"
      if (vland_1 /= vland_2) then
         print *, "Error! missing_vales do not match"
         stop
      endif
      print*, "Step 3.1 Success"
      
      print *, "Step 3.2 Compare scale_factor"
      if ( abs(scale_f_1 - scale_f_2) > eps ) then
          print *, "Error! scale_factors do not match"
          stop
      endif
      print*, "Step 3.2 Success"

      print *, "Step 3.3 Compare HS"
      do j = 1, JM
         do i = 1, IM
          if (hs_i_1(i,j) /= vland_1) then
            if (hs_i_2(i,j) /= vland_2) then
              hs_r_1(i,j) = hs_i_1(i,j) * scale_f_1
              hs_r_2(i,j) = hs_i_2(i,j) * scale_f_2
              if (abs(hs_r_1(i,j)-hs_r_2(i,j))/hs_r_2(i,j) > eps_r) then
                print *, "Error! HS do not match at (i,j) = ", i, ",",j
                stop
              endif 
            else
              print *, "Error! HS do not match at (i,j) = ", i, ",",j
             stop 
          endif  
         endif
       enddo
     enddo
     print *, "Step 3.3 Success"
     print *, "Compare Success"
         
end
