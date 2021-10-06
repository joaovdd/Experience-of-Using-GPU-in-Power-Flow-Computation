# Set preprocessor variable
sed -i "s/^#define DOUBLE_MODE .*/#define DOUBLE_MODE false/" opcoesDeCompilacao.h
# Buiild parallel executable
nvcc -gencode arch=compute_61,code=sm_61 -ccbin g++ -L"/usr/include/mkl/intel64" -lmkl_intel_lp64 -lmkl_gnu_thread -lgomp -lmkl_core -lpthread -w -lm -ldl -DMKL_LP64 -m64 -I"/usr/include/mkl" main.cu solverCPU.cpp -o flFloat -lcublas -lcusolver -lcusparse -O0