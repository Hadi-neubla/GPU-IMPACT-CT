#!/bin/bash
. /etc/profile.d/modules.sh
module purge
module load rhel8/default-amp
module load hdf5
export HDF5_DIR=/usr/local/software/spack/spack-rhel8-20210927/opt/spack/linux-centos8-zen2/gcc-9.4.0/hdf5-1.10.7-krjrm53dyfdtyaqwab47n22s6xdrne2d/

export HDF5_DIR=/usr/local/software/spack/spack-views/rocky8-icelake-20220710/hdf5-1.10.8/intel-2021.6.0/intel-oneapi-mpi-2021.6.0/h75adcalc32r6k5hkaajovxvzvqa6rec/

export HDF5ROOT=/usr/local/software/spack/spack-rhel8-20210927/opt/spack/linux-centos8-zen2/gcc-9.4.0/hdf5-1.10.7-krjrm53dyfdtyaqwab47n22s6xdrne2d/

# unload all the intel modules
#module unload intel/compilers/2017.4
#module unload intel/mkl/2017.4
#module unload intel/impi/2017.4/intel
#module unload intel/libs/idb/2017.4
#module unload intel/libs/tbb/2017.4
#module unload intel/libs/ipp/2017.4
#module unload intel/libs/daal/2017.4 
#module unload intel/bundles/complib/2017.4
#
# load openmpi for skylake
#$module load openmpi-1.10.7-gcc-5.4.0-jdc7f4f
