# virtual-casing

## Installation
### C++ Executable
In addition to having C++ compiler (tested with GCC-9 and newer) that supports C++17 standard and OpenMP, install FFTW and LAPACK libraries. Intel MKL provides both of them.
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

To install the python package use pip:

    python -m pip install .

To build the wheel use:

    python -m build --wheel

which can be installed by:

    pip install dist/virtual_casing*.whl

## Testing

Run the tests using

    python -m unittest test -v