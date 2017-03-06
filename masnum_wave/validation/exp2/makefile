#Makefile for building compare_exp2 
# Written by Zhenya.Song, 01/06/2016
NETCDF_PATH = /usr/sw-mpp/apps/lib/netcdf
FC = sw5f90 -host
FFLAGS = 
LFOPTS = -L${NETCDF_PATH}/lib -lnetcdf 
INC_PATH = -I$(NETCDF_PATH)/include

SRC = compare_exp2.f90 Handle_err.f90

EXEC=compare_exp2

OBJS = ${SRC:.f90=.o}
INC = \
# netcdf.inc

.SUFFIXES: .f90 .o
.f90.o:
	$(FC) ${INC_PATH} $(FFLAGS) -c $< 

$(EXEC) :$(INC) $(OBJS)  
	$(FC) $(FFLAGS) -o  $@ $(OBJS) $(LFOPTS)

$(OBJS) :$(INC)


clean:
#	rm -f  $(EXEC) $(OBJS) *.mod 
	rm -f  $(OBJS) *.mod 



