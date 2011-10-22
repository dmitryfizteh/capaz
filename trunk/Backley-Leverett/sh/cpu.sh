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
#
mpiCC -D B_L ../../main.cpp ../../no_communication.cpp ../../cpu.cpp ../two-phase.cpp -o ../Debug/cpu.px
mpirun -np 1 -maxtime $1 ../Debug/cpu.px
