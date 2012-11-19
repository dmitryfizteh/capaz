#!/bin/bash
# Two cpu and one gpu are needed
# 

debug_name="measuring"

hostname=`hostname`
#echo $hostname

if [ "$hostname" = "mvse" ]
then
    ARCH=13
    PPN="-ppn 1"
	lib_path="-L/common/cuda/lib -lcudart"
    maxtime="-maxtime $5"
else 	if [ "$hostname" = "k100" ]
	then
	    ARCH=20
	    PPN="-ppn 1"
	    lib_path="-L/common/cuda/lib64 -lcudart"
	    maxtime="-maxtime $5"
	else
	    ARCH=20
	    PPN="-ppn 1"
	    lib_path="-L/usr/local/cuda/lib64 -lcudart"
	fi
fi

echo "nvcc -c -arch sm_$ARCH ../Debug/measuring.o ./measuring.cu"
	  nvcc -c -arch sm_$ARCH ../Debug/measuring.o ./measuring.cu
arch_file="../Debug/measuring.o"

if [ "$hostname" = "k100" ]
then
      compilator="mpicxx"
else 	if [ "$hostname" = "mvse" ]
	then		
			compilator="mpicc"
	else
		compilator="mpicxx"
	fi
fi

mkdir ../Debug

echo "$compilator $lib_path ./measuring.cpp $arch_file -o ../Debug/$debug_name.px"
      $compilator $lib_path ./measuring.cpp $arch_file -o ../Debug/$debug_name.px

cd ../Debug
echo "mpirun $PPN -np 2 10 ./$debug_name.px"
mpirun $PPN -np 2 10 ./$debug_name.px

exit

