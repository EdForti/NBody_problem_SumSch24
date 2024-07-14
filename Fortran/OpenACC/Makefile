# Copyright (c) 2020 NVIDIA Corporation.  All rights reserved. 

FC := nvfortran -acc -O3 -ta=tesla -Minfo=all -acc=noautopar
FLAGS := -mp 
# ACCFLAGS := -Minfo=all

OpenAcc_NBody: OpenACC_NBody.F90
	${FC} ${FLAGS} ${ACCFLAGS} -o executable.exe OpenACC_NBody.F90 
clean:
	rm -f *.o executable.exe
