# Experience of Using GPU in Power Flow Computation
## Build Binaries
### Dependencies
This project uses Nvidia Cuda Tooklit and Intel MKL. Note that (when applicable) the libraries' include paths in the build scripts should match the ones where the respective libraries' are present on the system.

### CPU Versions

Run buildParallelDPFP.sh, buildParallelSPFP.sh, buildSequentialDPFP.sh and buildSequentialSPFP.sh on the CPU folder.

### GPU Versions

Run buildDPFP.sh and buildSPFP.sh on the GPU folder.

## Run Tests

The programs' execution options are set in the CPU/data.txt and GPU/data.txt files.

### CPU

Main options:

* metodo (sets the method to be used):
    * esparso (sparse routines),
    * denso (dense routines)

* openmp:
    * ON - enable openmp parallelism (option used for parallell tests),
    * OFF - disable openmp parallelism (option used for sequential tests)

* arq_entrada (sets the system data file to be read):
    * ieee14.cdf
    * ieee30.cdf
    * ieee57.cdf
    * ieee118.cdf
    * ieee118x2.cdf
    * ieee118x4.cdf
    * ieee118x8.cdf
    * ieee118x16.cdf
    * ieee118x32.cdf
    * ieee118x64.cdf
    * ieee118x128.cdfx

### GPU

Main options:

* metodo (sets the method to be used):
    * hibridoC (hybrid CPU-GPU routines),
    * esparso (sparse routines),
    * denso (dense routines)

* arq_entrada (sets the system data file to be read):
    * ieee14.cdf
    * ieee30.cdf
    * ieee57.cdf
    * ieee118.cdf
    * ieee118x2.cdf
    * ieee118x4.cdf
    * ieee118x8.cdf
    * ieee118x16.cdf
    * ieee118x32.cdf
    * ieee118x64.cdf
    * ieee118x128.cdfx