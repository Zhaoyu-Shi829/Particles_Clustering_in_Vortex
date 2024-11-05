#!/bin/bash
#
# Build options for MGLET
#
# USAGE:
#   1. Place this file somewhere it is easily reached, for example in your home
#      folder. It might very well be a hidden file in your home directory, i.e.
#        ~/.mglet-buildconf
#   2. Update the variables in this file according to your computer and build
#      environment.
#   3. Source the file from your ~/.bashrc file, i.e. insert the following:
#        source ~/.mglet-buildconf
#      This will always make your build configuration available to the makefile.
#   4. If you want, you might create multiple copies of the file, with different
#      settings. This can be for building with different compilers, debug- and
#      release builds etc.
#


################################################################################
# (NORMAL) USER EDITABLE PART:

# Path to MPI include files (where to find mpif.h), including '-I'
MGLET_MPI_INC="${MPI_INC}"

# Compilation flags for Debug and Release builds
# -check bounds,pointers,uninit -gen-interfaces nosource -warn interfaces
if [ "$MGLET_BUILD_TYPE" == "Debug" ]; then
    export MGLET_CFLAGS="-c -O0 -g -Wall"
    export MGLET_FFLAGS="-c -O0 -g -check all -fpe0 -init=snan,arrays -traceback -fp-stack-check -heap-arrays -Iqpack4 -I.."
else
    export MGLET_CFLAGS="-c -g -O3 -xHost -no-prec-div -qopenmp-simd"
    export MGLET_FFLAGS="-c -g -O3 -xHost -no-prec-div -qopenmp-simd -heap-arrays -diag-disable 8291,8293"
fi

if [ "$MGLET_PRECISION" == "Double" ]; then
    export MGLET_FFLAGS="$MGLET_FFLAGS -real-size 64"
fi

# Extra flags for Fortran 90 code
export MGLET_F90FLAGS="-free"

# Extra flags for linking MGLET executable
export MGLET_LFLAGS="-shared-intel"

# Extra libraries to link MGLET executable against (QPack is automatically included)
export MGLET_LLIBS=""

# C compiler
export MGLET_CC="mpiicc"

# Fortran compiler
export MGLET_FC="mpiifort"

# END OF (NORMAL) USER EDITABLE PART
################################################################################

# Detect if HDF5 is present and set correct options
function mglet_set_hdf5 {
    H5PATH=$1
    export MGLET_HDF5="Yes"
    
    # Find HDF5 include path where the compiler can find "hdf5.mod"
    if [ -e "$H5PATH/include/hdf5.mod" ] || [ -e "$H5PATH/include/HDF5.mod" ] ; then
        H5_INCPATH="${H5PATH}/include"
    elif [ -e "$H5PATH/include/shared/hdf5.mod" ] || [ -e "$H5PATH/include/shared/HDF5.mod" ] ; then
        H5_INCPATH="${H5PATH}/include/shared"
    else
        echo "ERROR: Could not detect HDF5 include path"
    fi
    export MGLET_FFLAGS="${MGLET_FFLAGS} -I${H5_INCPATH}"
    export MGLET_CFLAGS="${MGLET_CFLAGS} -I${H5_INCPATH} -DHDF5_COMPATIBILITY"
    
    # Find name of HDF5 library
    if [ -e "$H5PATH/lib/libhdf5_fortran.so" ] ; then
        export MGLET_LFLAGS="${MGLET_LFLAGS} -L${H5PATH}/lib -lhdf5_fortran"
    elif [ -e "$H5PATH/lib/libhdf5_fortran-shared.so" ] ; then
        export MGLET_LFLAGS="${MGLET_LFLAGS} -L${H5PATH}/lib -lhdf5_fortran-shared"
    else
        echo "ERROR: Could not detect HDF5 library name"
    fi
}

if [ ! -z "$HDF5_DIR" ]; then
    mglet_set_hdf5 $HDF5_DIR
elif [ ! -z "$HDF5_BASE" ]; then
    mglet_set_hdf5 $HDF5_BASE
elif [ ! -z "$EBROOTHDF5" ]; then
    mglet_set_hdf5 $EBROOTHDF5
fi

# Which flags to use when preprocessing
export MGLET_CPPFLAGS="-P -traditional ${MGLET_MPI_INC}"

# Which preprocessor to use
export MGLET_CPP=`which cpp`

# Name of executable
export MGLET_EXEC_NAME="mg41.exe"

# Set this variable to "No" to explicitly NOT build & link against QPack
# (any other variable, even an empty string, means building with QPack)
export MGLET_BUILD_QPACK="Yes"

# QPack version
export MGLET_QPACK_VERSION=4

