#!/bin/bash
# Two cpu and one gpu are needed
# 

debug_name = "measuring"

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

echo "nvcc -c -arch sm_$ARCH gpu.o ./measuring.cu"
	  nvcc -c -arch sm_$ARCH gpu.o ./measuring.cu
arch_file="gpu.o"

if [ "$hostname" = "k100" ]
then
      compilator="mpicxx"
   fi
else 	if [ "$hostname" = "mvse" ]
	then		
			compilator="mpicc"
	    fi
	else
		compilator="mpicxx"
	fi
fi

mkdir ../Debug

echo "$compilator $lib_path ./measuring.cpp $comm_file ../shared_test.cpp $arch_file -o ../Debug/$2_$3$debug_name.px"
      $compilator $lib_path ./measuring.cpp $arch_file -o ../Debug/$debug_name.px

cd ../Debug
echo "mpirun $PPN -np 2 10 ./$debug_name.px"
mpirun $PPN -np 2 10 ./$debug_name.px

exit
