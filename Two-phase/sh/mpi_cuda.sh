#!/bin/bash
if [ -z "$2" ]; then
    echo Use: $0 processors_count time_limit_in_minutes
    exit
fi
#
for i in $( ls | grep .img )
do
 rm $i
done;
#
  for i in $( ls |grep -v .sh|grep -v .1 | grep -v .px | grep -v machines. | grep -v result. ); 
  do 
    dos2unix $i;
  done;
  ARCH='20'
#
nvcc -D TWO_PHASE -c -arch sm_$ARCH gpu.o ../../gpu.cu
mpiCC  -D TWO_PHASE -L/common/cuda/lib64 -lcudart ../../main.cpp ../../mpi.cpp gpu.o -o ../Debug/mpi_cuda.px
mpirun -ppn 3 -np $1 -maxtime $2 ../Debug/mpi_cuda.px
