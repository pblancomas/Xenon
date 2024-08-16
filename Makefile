prefix = /usr/local
exec_prefix = ${prefix}
libdir = ${exec_prefix}/lib

local_CFLAGS = -g -O4
GSL_CFLAGS = $(shell gsl-config --cflags)
GSL_LIBS = $(shell gsl-config --libs)

BIN = main 
OBJ = main.o 
all: $(BIN)
	
main: main.o core.o 8B_flux.o interpolators.o FFM.o
	gcc $(CFLAGS) $(GSL_CFLAGS) -o main main.o core.o 8B_flux.o interpolators.o FFM.o $(GSL_LIBS)


%.o: %.c
	gcc $(CFLAGS) $(GSL_CFLAGS) -c $<

.PHONY: clean
clean:
	rm -f $(BIN) $(OBJ)
