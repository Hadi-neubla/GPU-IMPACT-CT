#*************************************************************************************************************
# GPU targets by Hadi Zolfaghari, ARTORG Center, University of Bern, (zolfaghari.haadi@gmail.com)            *
#* October 2015 - March 2020                                                                                 *
#*************************************************************************************************************

ifeq ($(machine),artorg)
  ALL = local
else 
ifeq ($(machine),daint)
  ALL = xc30
else
ifeq ($(machine),login)
  ALL = wilkes3
else
  ALL = portable
endif
endif
endif


execute:	wilkes3

# target for local system (PC)
local:	COMP = mpifort
local:	INCL = -I/usr/include -I/opt/hdf5_gcc_native_shared/include -I/opt/openmpi_gcc_native_shared/include
local:	LIBS = -L/usr/lib/x86_64-linux-gnu -lm -ldl -llapack -L/opt/openmpi_gcc_native_shared/lib -lmpi_mpifh -lmpi_usempif08 -L/opt/hdf5_gcc_native_shared/lib -lhdf5_fortran -lhdf5hl_fortran -lhdf5_hl -lhdf5 -lz -lm -ldl 
local:	OPT1 = -ffree-form -ffree-line-length-none -cpp -fdefault-real-8 -fdefault-double-8 \
-O3 -ffast-math -funroll-all-loops -fpeel-loops -ftracer -funswitch-loops \
#-ftree-vectorize -march=native -funsafe-math-optimizations -fomit-frame-pointer -fno-math-errno -fno-stack-limit 
	#-fopenmp # -J$(DST) 
local:	OPT2 = $(OPT1) -xf95-cpp-input
local:	$(BUILD_TARGET)
#local:  ftopy

# target for laptop
portable: COMP = mpifort
portable: INCL = -I/usr/include -I/opt/hdf5_gcc_moose/include -I/opt/moose/openmpi/openmpi-1.8.4/gcc-opt/include
portable: LIBS = -L/usr/include -lm -ldl -llapack -L/opt/moose/openmpi/openmpi-1.8.4/gcc-opt/lib -lmpi_mpifh -lmpi_usempi -L/opt/hdf5_gcc_moose/lib -lhdf5_fortran -lhdf5_hl -lhdf5 -lz -lm -ldl
portable:	OPT1 = -ffree-form -ffree-line-length-none -cpp -fdefault-real-8 \
-g
#-O3 -ffast-math -funroll-all-loops -fpeel-loops -ftracer -funswitch-loops \
#-ftree-vectorize -funsafe-math-optimizations -fomit-frame-pointer -fno-math-errno
portable:	OPT2 = $(OPT1) -xf95-cpp-input
portable:	$(BUILD_TARGET)
#portable: f2py

# target for Cray XC30 (Piz Daint, 5272 8-core SandyBridge 64-bit CPU compute nodes, 5272 NVidia Tesla K20X w/ 6GB memory)
xc30:	COMP = ftn
xc30:   nvcc = nvcc
xc30:   NVCC = nvcc -Wno-deprecated-gpu-targets -lineinfo
xc30:	INCL = -I$(HDF5ROOT)/include
#xc30:   LIBCUDA = -L$(shell dirname $(shell which nvcc))/  -lcudart -lcuda
#xc30:   LIBCUDA = -L/opt/nvidia/cudatoolkit10/10.0.130_3.22-7.0.1.0_5.2__gdfb4ce5/extras/CUPTI/lib64 -lcudart -L/opt/cray/nvidia/default/lib64 -lcuda
#xc30:   LIBCUDA = -L/opt/nvidia/cudatoolkit10/10.0.130_3.22-7.0.1.0_5.2__gdfb4ce5/lib64  
xc30:   LIBCUDA = -lcuda -lcudart
xc30:	LIBS = -lstdc++ -L$(HDF5ROOT)/lib -lhdf5_fortran -lhdf5_hl -lhdf5 -lz -lm -ldl
xc30:	OPT1 = -ffree-form -ffree-line-length-none -cpp -fdefault-real-8 \
	-g -fopenmp
#	-O3 -ffast-math -funroll-all-loops -fpeel-loops -ftracer -funswitch-loops \
#	-ftree-vectorize -funsafe-math-optimizations -fomit-frame-pointer -fno-math-errno
xc30:	OPT2 = $(OPT1) -xf95-cpp-input
xc30:	$(BUILD_TARGET)
#xc30:	f2py


# target for Wilkes-2 (Cambridge Computing Services, Cambridge Service for Data Driven Discovery (CSD3))
wilkes:	COMP = mpifort
wilkes:   nvcc = nvcc
wilkes:   NVCC = nvcc   -gencode=arch=compute_87,code=sm_87 \
           -Wno-deprecated-gpu-targets -lineinfo
wilkes:	INCL = -I$(HDF5ROOT)/include 
#xc30:   LIBCUDA = -L$(shell dirname $(shell which nvcc))/  -lcudart -lcuda
#xc30:   LIBCUDA = -L/opt/nvidia/cudatoolkit10/10.0.130_3.22-7.0.1.0_5.2__gdfb4ce5/extras/CUPTI/lib64 -lcudart -L/opt/cray/nvidia/default/lib64 -lcuda
#xc30:   LIBCUDA = -L/opt/nvidia/cudatoolkit10/10.0.130_3.22-7.0.1.0_5.2__gdfb4ce5/lib64  
wilkes:   LIBCUDA = -lcuda -lcudart
wilkes:	LIBS = -lstdc++ -L$(HDF5ROOT)/lib -lhdf5_fortran -lhdf5_hl -lhdf5 -lz -lm -ldl 
wilkes:	OPT1 = -ffree-form -ffree-line-length-none -cpp -fdefault-real-8 \
	-g -fopenmp
#	-O3 -ffast-math -funroll-all-loops -fpeel-loops -ftracer -funswitch-loops \
#	-ftree-vectorize -funsafe-math-optimizations -fomit-frame-pointer -fno-math-errno
wilkes:	OPT2 = $(OPT1) -xf95-cpp-input
wilkes:	$(BUILD_TARGET)
#xc30:	f2py


# target for Wilkes3 (Cambridge Computing Services, Cambridge Service for Data Driven Discovery (CSD3))
wilkes3:   COMP = mpifort
wilkes3:   nvcc = nvcc
wilkes3:   NVCC = nvcc -Wno-deprecated-gpu-targets -lineinfo
wilkes3:   INCL = -I$(HDF5ROOT)/include 
#xc30:   LIBCUDA = -L$(shell dirname $(shell which nvcc))/  -lcudart -lcuda
#xc30:   LIBCUDA = -L/opt/nvidia/cudatoolkit10/10.0.130_3.22-7.0.1.0_5.2__gdfb4ce5/extras/CUPTI/lib64 -lcudart -L/opt/cray/nvidia/default/lib64 -lcuda
#xc30:   LIBCUDA = -L/opt/nvidia/cudatoolkit10/10.0.130_3.22-7.0.1.0_5.2__gdfb4ce5/lib64  
wilkes3:   LIBCUDA = -lcuda -lcudart
wilkes3:   LIBS = -lstdc++ -L$(HDF5ROOT)/lib -lhdf5_fortran -lhdf5_hl -lhdf5 -lz -lm -ldl 
wilkes3:   OPT1 = -ffree-form -ffree-line-length-none -cpp \
	-g -fopenmp
wilkes3:   OPT3 = -ffree-form -ffree-line-length-none -cpp -g -fopenmp 
#	-O3 -ffast-math -funroll-all-loops -fpeel-loops -ftracer -funswitch-loops \
#	-ftree-vectorize -funsafe-math-optimizations -fomit-frame-pointer -fno-math-errno
wilkes3:   OPT2 = $(OPT1) -xf95-cpp-input
wilkes3:   $(BUILD_TARGET)
#xc30:	f2py
