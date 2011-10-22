#!/bin/bash
if [ -z "$1" ]; then
    echo Use: $0 time_limit_in_minutes
    exit
fi
#
for i in $( ls | grep machines. )
do
 rm $i
done;
#
for i in $( ls | grep result. )
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
nvcc -D B_L -c -arch sm_$ARCH gpu.o ../../gpu.cu
mpiCC  -D B_L -L/common/cuda/lib64 -lcudart ../../main.cpp ../../no_communication.cpp gpu.o -o ../Debug/cuda.px
mpirun -np 1 -maxtime $1 ../Debug/cuda.px
