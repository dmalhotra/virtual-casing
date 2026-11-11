# virtual-casing

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23619289.svg)](https://doi.org/10.5281/zenodo.23619289)

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

Optionally, set the environment variable ``LIBXSMM_DIR`` to a build of [LIBXSMM](https://github.com/libxsmm/libxsmm)
(version 2.0 or newer) before running cmake, or before installing the Python package, to use it for small matrix
products. The Makefile uses it in the same way.


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

## Citing virtual-casing

If you use virtual-casing, or code derived from it, please cite the software and the papers below (also listed in [CITATION.cff](CITATION.cff)):

- D. Malhotra, B. Medasani, M. Landreman, C. Smiet, A. Cerfon, *virtual-casing: the virtual-casing principle on toroidal surfaces*, https://github.com/hiddenSymmetries/virtual-casing, [doi:10.5281/zenodo.23619289](https://doi.org/10.5281/zenodo.23619289).
- D. Malhotra, *BIEST: Boundary Integral Equation Solver for Taylor states*, https://github.com/dmalhotra/BIEST, [doi:10.5281/zenodo.23511944](https://doi.org/10.5281/zenodo.23511944).
- D. Malhotra, T. Li, R. Blackwell, A. Rahimian, W. Yan, X. Gao, *SCTL: Scientific Computing Template Library*, https://github.com/dmalhotra/SCTL, [doi:10.5281/zenodo.23511629](https://doi.org/10.5281/zenodo.23511629).
- D. Malhotra, A. J. Cerfon, M. O'Neil, E. Toler, "Efficient high-order singular quadrature schemes in magnetic fusion," *Plasma Physics and Controlled Fusion* 62 (2020) 024004, [doi:10.1088/1361-6587/ab57f4](https://doi.org/10.1088/1361-6587/ab57f4).
- D. Malhotra, A. Cerfon, L.-M. Imbert-Gérard, M. O'Neil, "Taylor states in stellarators: A fast high-order boundary integral solver," *Journal of Computational Physics* 397 (2019) 108791, [doi:10.1016/j.jcp.2019.06.067](https://doi.org/10.1016/j.jcp.2019.06.067).

<details>
<summary>BibTeX</summary>

```bibtex
@software{virtualcasing,
  author    = {Malhotra, Dhairya and Medasani, Bharat and Landreman, Matt and Smiet, Chris and Cerfon, Antoine},
  title     = {virtual-casing: the virtual-casing principle on toroidal surfaces},
  url       = {https://github.com/hiddenSymmetries/virtual-casing},
  doi       = {10.5281/zenodo.23619289},
  version   = {0.3.1},
  year      = {2026}
}

@software{biest,
  author    = {Malhotra, Dhairya},
  title     = {{BIEST}: Boundary Integral Equation Solver for {Taylor} states},
  url       = {https://github.com/dmalhotra/BIEST},
  doi       = {10.5281/zenodo.23511944},
  version   = {2.0.0},
  year      = {2026}
}

@software{sctl,
  author    = {Malhotra, Dhairya and Li, Tianyue and Blackwell, Robert and Rahimian, Abtin and Yan, Wen and Gao, Xuanzhao},
  title     = {{SCTL}: Scientific Computing Template Library},
  url       = {https://github.com/dmalhotra/SCTL},
  doi       = {10.5281/zenodo.23511629},
  version   = {2.1.0},
  year      = {2026}
}

@article{malhotra2020quadrature,
  author  = {Malhotra, Dhairya and Cerfon, Antoine J. and O'Neil, Michael and Toler, Evan},
  title   = {Efficient high-order singular quadrature schemes in magnetic fusion},
  journal = {Plasma Physics and Controlled Fusion},
  volume  = {62},
  number  = {2},
  pages   = {024004},
  year    = {2020},
  doi     = {10.1088/1361-6587/ab57f4}
}

@article{malhotra2019taylor,
  author  = {Malhotra, Dhairya and Cerfon, Antoine and Imbert-G{\'e}rard, Lise-Marie and O'Neil, Michael},
  title   = {Taylor states in stellarators: A fast high-order boundary integral solver},
  journal = {Journal of Computational Physics},
  volume  = {397},
  pages   = {108791},
  year    = {2019},
  doi     = {10.1016/j.jcp.2019.06.067}
}
```

</details>

## License

virtual-casing is licensed under the [Apache License 2.0](LICENSE). The [NOTICE](NOTICE) file lists the authors and must be kept in any redistribution, including ports to other languages. It also has the authors' requests on crediting them and on AI use.
