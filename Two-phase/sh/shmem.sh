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
shmemcc -D TWO_PHASE ../../main.cpp ../../shmem.cpp ../../cpu.cpp ../two-phase.cpp ../../shared_test.cpp -o ../Debug/shmem.px
shmemrun -np $1 -maxtime $2 ../Debug/shmem.px
