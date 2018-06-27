#!/bin/bash

# run intel compiler config
source /data/yiyang.zhang/bin/intel/compilers_and_libraries_2018.3.222/linux/bin/compilervars.sh intel64

# intel mpi include and bin
export PATH="/data/yiyang.zhang/bin/intel/compilers_and_libraries_2018.3.222/linux/mpi/intel64/include:$PATH"
export PATH="/data/yiyang.zhang/bin/intel/compilers_and_libraries_2018.3.222/linux/mpi/intel64/bin:$PATH"

# user installed packages
# export PATH="~/.local/bin:$PATH"

cmake -DCMAKE_C_COMPILER=/data/yiyang.zhang/bin/intel/compilers_and_libraries_2018.3.222/linux/bin/intel64/icc \
-DCMAKE_CXX_COMPILER=/data/yiyang.zhang/bin/intel/compilers_and_libraries_2018.3.222/linux/bin/intel64/icpc \
-H. -Bbuild