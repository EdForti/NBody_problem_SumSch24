# Copyright (c) 2020 NVIDIA Corporation.  All rights reserved. 

FC := nvfortran
FLAGS := -O3 -Minfo=all

Serial_NBody: Serial_NBody.F90
	${FC} ${FLAGS} -o executable.exe Serial_NBody.F90 
clean:
	rm -f *.o *.exe
