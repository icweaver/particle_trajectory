CFLAGS = -fdefault-real-8

sho:
	gfortran-5 ${CFLAGS} sho.f90 main.f90 -o main
	./main
	
