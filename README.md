# virtual-casing

## Installation
### C++ Executable
In addition to having C++ compiler (tested with GCC-9 and newer) that supports C++11 standard and OpenMP, install FFTW and LAPACK libraries. Intel MKL provides both of them.
Also install cmake package. The CMake setup is written in a monolithic way so we have to install python dependencies also,
even though they are strictly not required. On the python side, we require python header files and Numpy. All of them can be installed by using pip.

    pip install cmake numpy
  
Create a folder called ``build`` and switch to ``build`` directory. 

    mkdir build; cd build
    
Run cmake from within the build folder.

    cmake ..
  
If the above command succeeds, we can go to the next step. Build the ``vc_testing`` executable by running

    make vc_testing
  
If the compilation is successful, you will have ``vc_testing`` in the ``build`` folder. 


### Python package

We also need ninja package, which can be installed with pip

    pip install ninja
  
Then run 

    python setup.py install

