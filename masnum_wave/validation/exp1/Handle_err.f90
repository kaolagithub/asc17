subroutine Handle_err(status)

  include 'netcdf.inc'

  integer :: status

  if (status .NE. NF_NOERR) then
    print *, NF_STRERROR(status)
	stop 'Stopped'
  endif


end