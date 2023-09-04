!**************************************************************************************************                                           
!* by Hadi Zolfaghari, ARTORG Center, University of Bern, (zolfaghari.haadi@gmail.com)            *      
!* October 2015 - March 2020                                                                      *            
!**************************************************************************************************
module mod_cudafor
! Use Fortran standard module iso_c_binding,
! which contains types and functions for
! C compatibility

  use iso_c_binding
! Define dim3 type, which is needed for compute grids
! in CUDA C
type, bind(C) :: dim3
  integer(c_int) :: x, y, z
end type dim3



type, bind(C) :: cudaStream_t
   integer(c_int) :: stream
end type cudaStream_t 

! Define interface with CUDA runtime functions, that
! are going to be used in our GPU-enabled Fortran application.
interface

  ! Declare cudaMalloc function interface
  function cudaMalloc(ptr, length) bind(C, name = "cudaMalloc")
  use iso_c_binding
  implicit none

  ! "ptr" argument has a type of "C pointer" ~ void* in C,
  ! passed by address itself => void**
  type(c_ptr) :: ptr
  ! "length" argument has a type of "size_t", attribute "value"
  ! means that it is passed by value, not by address
  integer(c_size_t), value :: length
  ! Type of function return value
  integer(c_int) :: cudaMalloc
  end function

  function cudaMemcpy(dst, src, length, dir) bind(C, name = "cudaMemcpy")
  use iso_c_binding
  implicit none

  ! "dst" and "src" arguments are pointers -- "void*",
  ! passed by "value" => still void*, contrary to "ptr" above
  type(c_ptr), value :: dst, src
  integer(c_size_t), value :: length
  integer(c_int), value :: dir
  integer(c_int) :: cudaMemcpy
  end function

  function cudaHostRegister(arg, length, flag) bind(C, name = "cudaHostRegister")
  use iso_c_binding
  implicit none

  ! "dst" and "src" arguments are pointers -- "void*",
  ! passed by "value" => still void*, contrary to "ptr" above
  type(c_ptr)  :: arg
  integer(c_int) :: flag 
  integer(c_size_t), value :: length
  integer(c_int) :: cudaHostRegister
  end function





  function cudaMemcpyAsync(dst, src, length, dir, stream) bind(C, name = "cudaMemcpyAsync")
  use iso_c_binding
  import :: cudaStream_t
  implicit none

  ! "dst" and "src" arguments are pointers -- "void*",
  ! passed by "value" => still void*, contrary to "ptr" above
  type(c_ptr), value :: dst, src
  integer(c_size_t), value :: length
  integer(c_int), value :: dir
  type(cudaStream_t), value  ::  stream
 

  integer(c_int) ::  cudaMemcpyAsync
  end function

  function cudaStreamCreate(stream) bind(C, name = "cudaStreamCreate")
  use iso_c_binding
  import cudaStream_t
  implicit none

!  type(cudaStream_t), value :: stream

  type(cudaStream_t), value :: stream
  integer(c_int) :: cudaStreamCreate
  end function

  function cudaStreamSynchronize(stream) bind(C, name = "cudaStreamsynchronize")
  use iso_c_binding

  import :: cudaStream_t
  implicit none

!  type(cudaStream_t), value :: stream
  type(cudaStream_t), value :: stream

  integer(c_int) :: cudaStreamSynchronize
  end function  


 function cudaFree(ptr) bind(C, name = "cudaFree")
  use iso_c_binding
  implicit none

  type(c_ptr), value :: ptr
  integer(c_int) :: cudaFree
  end function

  function cudaGetErrorString(istat) bind(C, name = "cudaGetErrorString")
  use iso_c_binding
  implicit none

  ! cudaGetErrorString returns "const char*" in C,
  ! but here we make it to return "void*"
  type(c_ptr) :: cudaGetErrorString
  integer(c_int), value :: istat
  end function

  function cudaGetLastError() bind(C, name = "cudaGetLastError")
  use iso_c_binding
  implicit none

  integer(c_int) :: cudaGetLastError
  end function

  function cudaThreadSynchronize() bind(C, name = "cudaThreadSynchronize")
  use iso_c_binding
  implicit none

  integer(c_int) :: cudaThreadSynchronize
  end function

  function cudaConfigureCall(gridDim, blockDim, sharedMem, stream) &
    bind(C, name = "cudaConfigureCall")
  use iso_c_binding
  ! Import dim3 type here
  import :: dim3
  import :: cudaStream_t
  implicit none

  type(dim3), value :: gridDim, blockDim
  integer(c_size_t), value :: sharedMem
  type(c_ptr), value :: stream
  integer(c_int) :: cudaConfigureCall
  end function

  function cudaLaunch(handle) bind(C, name = "cudaLaunch_")
  use iso_c_binding
  implicit none

  character(c_char) :: handle
  integer(c_int) :: cudaLaunch
  end function

end interface

contains

!  function cudaSetupScalarArgument(arg, length, offset)
!
!  use iso_c_binding
!  implicit none
!
!  interface
!    function cudaSetupArgument(arg, length, offset) &
!      bind(C, name = "cudaSetupArgument")
!    use iso_c_binding
!    implicit none
!
!    type(c_ptr), value :: arg
!    integer(c_size_t), value :: length
!    integer(c_size_t), value :: offset
!    integer(c_int) :: cudaSetupArgument
!    end function
!
!  end interface
!
!  type(c_ptr), value :: arg
!  integer(c_size_t), value :: length
!  integer(c_size_t), value :: offset
!  integer(c_int) :: cudaSetupScalarArgument
!
!  cudaSetupScalarArgument = cudaSetupArgument(arg, length, offset)
!
!  end function cudaSetupScalarArgument



!  function cudaSetupArrayArgument(arg, length, offset)
!
!  use iso_c_binding
!  implicit none
!
!  interface
!
!    function cudaSetupArgument(arg, length, offset) &
!      bind(C, name = "cudaSetupArgument")
!    use iso_c_binding
!    implicit none
!
!    type(c_ptr) :: arg
!    integer(c_size_t), value :: length
!    integer(c_size_t), value :: offset
!    integer(c_int) :: cudaSetupArgument
!    end function
!
!  end interface
!
!  type(c_ptr) :: arg
!  integer(c_size_t), value :: length
!  integer(c_size_t), value :: offset
!  integer(c_int) :: cudaSetupArrayArgument
!
!  cudaSetupArrayArgument = cudaSetupArgument(arg, length, offset)
!
!  end function cudaSetupArrayArgument
end module mod_cudafor
