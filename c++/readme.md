# cloudmrhub Reconstructors

## how to
- convert your data with siemens_to ismrmrd your data
- read signal ISMRMRD file to mha with  ismrmrdToMha -i /data/MYDATA/H5example/n.h5 -o /data/MYDATA/H5example/N.mha
- read noise ISMRMRD file to mha with  ismrmrdToMha -i /data/MYDATA/H5example/n.h5 -o /data/MYDATA/H5example/N.mha


## programs list

- [x] ismrmrdToMha
    - [ ] Oversapling
- [ ] Reconstructors
	- [ ] base class
            - [ ] CalculaltenoiseBW()
            - [ ] CalculatenNoiseCorrelationMatrix()
	- [ ] RSS
	- [ ] B1
	- [ ] SENSE
	- [ ] GRAPPA


## requirement
- ITK >= 4.9
- boost
- hdf5
- ismrmrd

## note
MHa file and don't support std::vector<complex<float>> tensors so we developed 2 method to interact with the classes, though a h5 or though the images

## Classes

### Recon
- take a kspace matrix mha of complex kspace data and a noise one
- Noise is a VNL matrix coils x number of pixels
- when the noise kspace is set the Noise is updated
- everytime the Noise is updated the Noiseconvaiace matrix is calculated and the signal should be prewhitened







### problems

// HDF5 C++ Wrapper compiler. Used only to detect HDF5 compile flags.
HDF5_CXX_COMPILER_EXECUTABLE:FILEPATH=/usr/bin/h5c++

// HDF5 Wrapper compiler. Used only to detect HDF5 compile flags.
HDF5_C_COMPILER_EXECUTABLE:FILEPATH=/usr/bin/h5cc

// Path to a file.
HDF5_C_INCLUDE_DIR:PATH=/usr/include

// HDF5 file differencing tool.
HDF5_DIFF_EXECUTABLE:FILEPATH=HDF5_DIFF_EXECUTABLE-NOTFOUND

// The directory containing a CMake configuration file for HDF5.
HDF5_DIR:PATH=HDF5_DIR-NOTFOUND

// HDF5 Fortran Wrapper compiler. Used only to detect HDF5 compile flags.
HDF5_Fortran_COMPILER_EXECUTABLE:FILEPATH=/usr/bin/h5fc

// HDF5 library compiled with parallel IO support
HDF5_IS_PARALLEL:BOOL=TRUE

// Path to a library.
HDF5_dl_LIBRARY_DEBUG:FILEPATH=HDF5_dl_LIBRARY_DEB UG-NOTFOUND

// Path to a library.
HDF5_dl_LIBRARY_RELEASE:FILEPATH=/usr/lib/x86_64-linux-gnu/libdl.so

// Path to a library.
HDF5_hdf5_LIBRARY_DEBUG:FILEPATH=HDF5_hdf5_LIBRARY _DEBUG-NOTFOUND

// Path to a library.
HDF5_hdf5_LIBRARY_RELEASE:FILEPATH=/usr/lib/x86_64-linux-gnu/libhdf5.so

// Path to a library.
HDF5_m_LIBRARY_DEBUG:FILEPATH=HDF5_m_LIBRARY_DEBUG-NOTFOUND

// Path to a library.
HDF5_m_LIBRARY_RELEASE:FILEPATH=/usr/lib/x86_64-linux-gnu/libm.so

// Path to a library.
HDF5_pthread_LIBRARY_DEBUG:FILEPATH=HDF5_pthread_L IBRARY_DEBUG-NOTFOUND

// Path to a library.
HDF5_pthread_LIBRARY_RELEASE:FILEPATH=/usr/lib/x86_64-linux-gnu/libpthread.so

// Path to a library.
HDF5_z_LIBRARY_DEBUG:FILEPATH=HDF5_z_LIBRARY_DEBUG-NOTFOUND

// Path to a library.
HDF5_z_LIBRARY_RELEASE:FILEPATH=/usr/lib/x86_64-linux-gnu/libz.so
for hdf5

Now how do I use the -D directive to override these so that program with compile with my local copy of hdf5 in $HOME/opt/hdf5? I am not very familiar with cmake. Thanks
