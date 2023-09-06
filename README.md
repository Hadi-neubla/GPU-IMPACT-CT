# GPU-IMPACT-CT
A data and task parallel framework for image-based simulations of blood flow in patient-specific geometries

Compile:
- Edit src/startup.sh according to your HDF5 setup and load GPU-related environmental modules.
- cd src && source startup.sh
- make

Run:
- cd ../prog && open config.txt
- change the simulation inputs (refer to src/usr_config.f90 file for info).
- copy relevant geometry input (vtk) and waveform input (csv) files from preprocessing directory.
- enable CUDA multiprocess service
- mpirun -n $NUM_MPI_THREADS ./impact.exe
  
