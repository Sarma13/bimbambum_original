# BIMBAMBUM v1.0

## BIMBAMBUM

**BIMBAMBUM** stands for **B**oundary **I**ntegral **M**ethod for **B**ubble **A**nalysis and **M**odeling 
in **B**ounded and **U**nbounded **M**edia. It is a potential flow solver mainly written in C++ intended at modeling the behaviour 
 of a single cavitation bubble in the vicinity of an initially flat fluid-fluid interface. The software solves the
 Laplace equation using a boundary integral formulation in axisymmetric coordinates. The code has the vocation to allow 
 research teams, students and any other interested person to quickly get insights into the behaviour of cavitation bubbles
and hopefully to provide a basis for an accelerated and collaborative developments of this solver.


## What can I do with the solver?

With **BIMBAMBUM** you can solve:

* The dynamics of a cavitation bubble near a fluid-fluid interface of any density ratio.
* The pressure and velocity fields associated with the bubble dynamics.

## Repository Structure

**BIMBAMBUM** is a two-stage numerical flow solver. The first stage computes the time evolution of the surfaces of the bubble and 
nearby boundary. The associated code is available in the `Bubble_Dynamics` folder. The second stage computes the velocity and pressure 
fields associated with the bubble dynamics at any selected time point in the bubble lifetime. 
The associated code is available in the `Flow_Field_Quantities` folder. Both stages of the solver must be built and executed separately.
The first stage of the solver may be used independently while the second stage needs
the simualtion results of the first stage as inputs.

    BIMBAMBUM v1.0
    ├── Bubble_Dynamics         # compuation of the bubble and fluid-fluid interface dynamics
    │   ├── src                 # software c++ source files
    │   ├── include             # software c++ header files 
    │   ├── solver              # main function for solver execution, called at program startup.
    │   ├── config_examples     # examples configuration files (software inputs)
    │   ├── post_processing     # python scripts for rudimentary post-processing of the simulation
    │   ├── doc                 # doxygen folder (compilation optional)
    │   └── CMakeLists.txt 
    ├── Flow_Field_Quantities   # compuation of the velocity and pressure fields
    │   ├── solver_codes        # software scripts (c++ and Python)
    │   ├── simulation_examples # simualtion result examples of first stage (needed as inputs in second stage)
    │   ├── post_processing     # python scripts for rudimentary post-processing of the simulation
    │   └── CMakeLists.txt 
    ├── license.md
    └── README.md
    
## Dependencies

The solver relies on the following dependencies.

* c++:
    * a C++ compiller with preferably OpenMP support (tested: GCC 11.3 on Linux Ubuntu and Apple clang 14.0.3 and GCC 13.1.0 on macOS Ventura)
    * a CMake build system (https://cmake.org/, version 3.18 or higher)
    * the Armadillo library for linear algebra (https://arma.sourceforge.net/docs.html; version tested 10.8.2 and 12.4.0).
    * the GNU scientific library for numerical computing (https://www.gnu.org/software/gsl/; version tested 2.7).
    * the BOOST library for input parsing (https://www.boost.org/; version tested 1.74.0 and 1.81.0)

* Python:
    * Python 3 (version tested 3.11.3)
    * Numpy (version tested 1.24.3)
    * Matplotlib (version tested 3.7.1)
    
* Python and c++ binding:
    * pybind11 (https://github.com/pybind/pybind11; version tested 2.10.4)


## Building the software and running the examples

### Bubble dynamics

#### A) Building the software

This portion of the solver computes the temporal evolution of the bubble and fluid-fluid interface boundaries. It may be 
run as stand-alone and can be built from the terminal with the following command lines on Linux systems:

```
$ cd Bubble_Dynamics/
$ mkdir build
$ cd build
$ cmake -DBUILD_OPENMP=<OpenMP option> -DBUILD_DOXYGEN=<Doxygen option> ../
$ make install       # optional -> install the solver as a library (note that special permission may be needed for installation)
$ make solver
$ make doxydoc       # optional -> write the doxygen documentation
```
where `OpenMP option` and `Doxygen option` are either 0 (false) or 1 (true) to enable parallel processing
and the built of the Doxygen documentation, respectively.

On macOS systems the same commands lines may be used with an Apple clang or a GCC compiler. In the later case, the `cmake` command 
must fitted with the following additional flags:

```
$ cmake -DBUILD_OPENMP=<OpenMP option> -DBUILD_DOXYGEN=<Doxygen option> -DCMAKE_C_COMPILER=/your/path/to/gcc-XX -DCMAKE_CXX_COMPILER=/your/path/to/g++-XX ../
```

where `XX` stands for the version of the GCC compiler. Note that `OpenMP` may not be supported on every system in which case the software
must be run in series, i.e. with  `-DBUILD_OPENMP=0`.

#### B) Running an example

The executable is installed in a `bin` folder alongside a set of `.json` configuration files provided as examples to run
the simulations. These configurations files contain the input data defining the nature of the simulation. The set of 
input data is provided in the table below:

| Parameters         | Type    | Short description                                                                  |
| -------------------| ------- | ---------------------------------------------------------------------------------- |
| Nb                 | int     | The number of elements for the bubble discretization                               |
| Ns                 | int     | The number of elements for the fluid-fluid interface discretization              |
| bubble_dynamics    | string  | The physics governing the bubble (`Rayleigh_Bubble` or `Rayleigh_Plesset_Bubble`)  |
| zeta               | double  | Buoyancy parameter                                                                 |
| gamma              | double  | Stand-off distance                                                                 |
| alpha              | double  | Density ratio of the fluid-fluid interface                                         |
| epsilon            | double  | Strength parameter                                                                 |
| k                  | double  | Specific heat ratio                                                                |
| surface_elasticity | boolean | Whether or not the interface has elasticity (`true` or `false`)                    |
| sigma_s            | double  | Fluid-fluid interface surface tension                                              |
| boundary           | string  | Nature of the initial geometry of the boundary (`from_code` must be used in v1.0)  |
| temporal_solver    | string  | Order of the temporal discretization (`RK1` or `RK2`)                              |
| dumper_filename    | string  | Name of the dumper file (stores the solution results)                              |
| delta_phi          | double  | Constant for adaptive time stepping                                                |
| filtering_freq     | int     | Filtering frequency (i.e. filtering every nth time steps)                          |
| n_threads          | int     | Number of threads used for parallel processing                                     |


The solver may then be run with the following command lines:

```
$ cd bin/
$ ./main <name of config file>.json
```

#### C) Post-processing the results

The results of the simulation are written in a `.txt` file whose name is defined by `dumper_filename`. The 10 first lines 
of this file contain information regarding key input parameter of the simulation. For lines 11 and below, the first column 
indicates the time step and the second column the simulation time. The remainder of the columns store the position of the node
points of the discretized bubble surface and fluid-fluid interface as well as the values of the normal velocities 
and the potential values at those points. Further information regarding the structure of this file is documented in 
`Bubble_Dynamics/include/BIM_solver.hpp` (in the method `write_solution`) and in `Bubble_Dynamics/post_processing/quick_output_dispaly.py`. 

Simple post-processing of the simulation results can be done with the `Python` scripts available in the folder `Bubble_Dynamics/post_processing`. 
These short scripts are provided alongside a sample `.txt` simulation result file and may be run as such. 

### Flow field quantities

#### A) Building the software

This portion of the solver computes the flow field quantities (velocity and pressure) associated with the 
bubble dynamics at any selected time point in the bubble lifetime. To build this portion of the code, the 
user must first download pybind11 (https://github.com/pybind/pybind11) . This can be done as follow 
from the terminal with the following command lines:

```
$ cd Flow_Field_Quantities/
$ git clone https://github.com/pybind/pybind11.git
```
Once pybind11 is downloaded, the code may be built with the following command lines:
```
$ mkdir build
$ cd build
$ cmake ../
$ make
```
Note that the `source_dir` of the command `add_subdirectory(<source_dir>)` may have to be modified in the `CMakeLists.txt` file 
to match the name of the directory containing `pybind11`.  

#### B) Running an example

All files needed to run this portion of the solver are copied in the `build` folder alongside sample 
results of the first processing stage which are needed as inputs for the second processing stage. The remainder
of the inputs needed for the simulation may be changed by the user in the `inputs.py` script. The set of 
input data is provided in the table below: 

| Parameters         | Type         | Short description                                                                                                        |
| -------------------| ------------ | -------------------------------------------------------------------------------------------------------------------------|
| file_name          | string       | Input file name (i.e. results files from first processing stage)                                                         |
| dumper_name        | string       | Output file name (no need to provide an extension)                                                                       |
| Nb                 | int          | The number of elements for the bubble discretization (must be the same as the one written in `file_name`)                |
| Ns                 | int          | The number of elements for the fluid-fluid interface discretization (must be the same as the one written in `file_name`) |
| zeta               | double       | Buoyancy parameter (must be the same as the one written in `file_name`)                                                  |
| time_step          | list of int  | List of the time steps (from the first processing phase) for which the flow fields are computed                          |
| h_grid             | double       | Grid spacing                                                                                                             |
| s_grid             | double       | Offsetting constant for shifting the domain boundaries away from the original boundaries                                 |
| r_max              | double       | Extend of the computational domain in the r-direction                                                                    |
| z_min              | double       | Extend of the computational domain in the z-direction (min) (N.B. the bubble is located at z<0)                          |
| z_max              | double       | Extend of the computational domain in the z-direction (max)                                                              |


An example may then directly be run with the following command line:  

```
$ python3 main.py
```

#### C) Post-processing the results

The results of the simulation are written in a `<dumper_name><time_step>.txt` file.

We provide a simple `Python` script in `Flow_Field_Quantities/post_processing` to display the computed velocity and pressure
fields. This short script is provided alongside a sample `.txt` simulation result file and may be run as such. 

Alternatively, the `<dumper_name><time_step>.txt` file may also be read in the open source post-processing visualization 
engine `Paraview` (https://www.paraview.org/) with the following procedure: 
1) open the file `<dumper_name><time_step>.txt` in `Paraview`, 
2) apply the `Table To Points` filter and assign the appropriate columns to the `X`, `Y` and `Z` coordinates, 
3) apply the `Delaunay 2D` filter and,
4) select the flow field quantity to be displayed.

## Funding

We gratefully acknowledge the support of the Swiss National Science Foundation (grants No. 179018 and No. 186158) and the MSCA-ITN-ETN of the European Union's H2020 program (REA grant agreement no. 813766).

## License

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

**BIMBAMBUM**  is a free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
**BIMBAMBUM**  is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details http://www.gnu.org/licenses.

