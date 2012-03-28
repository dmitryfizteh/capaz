#!/bin/bash
if [ -z "$1" ]; then
    echo Use: $0 time_limit_in_minutes
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
#
mpicc -D THREE_PHASE ../../shared_test.cpp ../../no_communication.cpp ../../cpu.cpp ../three-phase.cpp ../../main.cpp -o ../Debug/cpu.px
mpirun -np 1 -maxtime $1 ../Debug/cpu.px
