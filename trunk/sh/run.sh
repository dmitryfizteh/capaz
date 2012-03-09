#!/bin/bash
# ADD: remove temp files
# ADD: test script parameters

# Описание параметров запуска скрипта
if [ -z "$5" ]; then
    echo Use: ./$0 task_name architecture communication processors_count time_limit_in_minutes [debug]
    echo task_name: 2ph,3ph,bl
    echo architecture: cpu,gpu
    echo communication: no,mpi,shmem
    exit
fi

# Режим отладки
if [ "$6" = "debug" ]
then
    debug="-D MY_TEST"
    debug_name="_$6"
fi

hostname=`hostname`
#echo $hostname

if [ "$hostname" = "mvse" ]
then
    ARCH=13
else
    ARCH=20
fi

if [ "$hostname" = "k100" ]
then
    PPN=3
	lib_path="-L/common/cuda/lib64 -lcudart"
else
    PPN=1
fi

if [ "$1" = "2ph" ]
then
    task_name="-D TWO_PHASE"
    project_file="two-phase.cpp"
    project_folder="Two-phase"
fi

if [ "$1" = "3ph" ]
then
    task_name="-D THREE_PHASE"
    project_file="three-phase.cpp"
    project_folder="Three-phase"
fi

if [ "$1" = "bl" ]
then
    task_name="-D B_L"
    project_file="b-l.cpp"
    project_folder="Backley-Leverett"
fi

if [ "$2" = "gpu" ]
then
    echo "nvcc $task_name $debug -c -arch sm_$ARCH gpu.o ../$project_folder/gpu.cu"
	nvcc $task_name $debug -c -arch sm_$ARCH gpu.o ../gpu.cu
    arch_file="gpu.o"
else
    arch_file="../cpu.cpp ../$project_folder/$project_file"
	lib_path=""
fi

if [ "$3" = "no" ]
then
    comm_file="../no_communication.cpp"
fi

if [ "$3" = "mpi" ]
then
    comm_file="../mpi.cpp"
fi

if [ "$3" = "shmem" ]
then
    comm_file="../shmem.cpp"
    compilator="shmemcc"
else
    compilator="mpicxx"
fi

echo "$compilator $task_name $debug $lib_path ../main.cpp $comm_file ../shared_test.cpp $arch_file -o ../$project_folder/Debug/$2_$3$debug_name.px"
$compilator $task_name $debug $lib_path ../main.cpp $comm_file ../shared_test.cpp $arch_file -o ../$project_folder/Debug/$2_$3$debug_name.px

cd ../$project_folder/Debug
echo "mpirun -ppn $PPN -np $4 -maxtime $5 ../$project_folder/Debug/$2_$3$debug_name.px"
mpirun -ppn $PPN -np $4 -maxtime $5 ./$2_$3$debug_name.px

exit