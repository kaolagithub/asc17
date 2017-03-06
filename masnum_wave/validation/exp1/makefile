#Makefile for building compare_exp1 
# Written by Zhenya.Song, 01/06/2016
#NETCDF_PATH = /WORK/home/fio_climate1/tools/netcdf/netcdf-3.6.2
NETCDF_PATH = /usr/local/netcdf
FC = ifort -host
FFLAGS = 
LFOPTS = -L${NETCDF_PATH}/lib -lnetcdff
INC_PATH = -I$(NETCDF_PATH)/include

SRC = compare_exp1.f90 Handle_err.f90

EXEC=compare_exp1

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



