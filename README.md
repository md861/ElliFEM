# ElliFEM
A Fortran code to solve a boundary-value elliptic problem in 2D. The solver uses lagrange polynomials based Finite Element Method (FEM) for space discretizations. The solver also provides an option to enrich the solution space for better convergence rates [[1]](#1). The code is written for the Linux API, however, the system calls could be modified to be run over other operating systems as well. 

<img src="https://github.com/md861/ElliFEM/blob/main/images/mesh_p4.png" width="400" height="450"> <img src="https://github.com/md861/ElliFEM/blob/main/images/30pi_p4.png" width="400" height="450">
## Features
* Directly import 2D meshes generated in [Gmsh](https://gmsh.info/) - an open source mesh generator.
* Dynamic allocation of all variables and arrays, depending on the imported mesh and order of FEM polynomials.
* Plot the mesh and numerical solutions in [Paraview](https://www.paraview.org/) - an open source alternative to Tecplot.
* Apply Neumann, Dirichlet and Robin boundary conditions on edges marked respectively in Gmsh. Automated calculation of edge normals for using Neumann boundaries.
* LAPACK libraries for matrix solutions.

A complete documentation for the description and usage of the package is under way. However, for the time-being, this read-me file should provide a succinct manual.

## Version system
The version number is based on the following format: *y*:*[m]m*:*n[s]*, for the *n*-th iteration of the code, in the *m*-th month of the *y*-th year in the current decade. The *s* denotes a stable version. Download the latest version available from the repository. The old(er) versions may not be archived for long periods in the future, depending on the number of iterations and the size of updates.

## Setting up the API
### Compiler 
Make sure that the compiler supports Fortran. For Linux based systems, generally this is provided by the GNU compiler. The latest version of ElliFEM was compiled using gcc 9.3.0
### LAPACK
Download the latest version from [netlib.org](http://www.netlib.org/lapack/#_release_history). The latest version of ElliFEM was compiled using LAPACK, version 3.9.0. Make sure the path of the library is resolved properly. In Linux, this can be easily checked via running `ld -lblas -llapack --verbose` in the terminal.
### [Gmsh](https://gmsh.info/)
This open-source meshing tool is required to produce the requisite format of data file, that is fed to the ElliFEM for computations. Read the instructions below on how to prepare the mesh and define the computational domain.
### [Paraview](https://www.paraview.org/)
This tool is not required to run the ElliFEM code, although, it does provide a way to view the output numerical results (if plotted) using this open-source software. 

## Directory tree and output files
The output files are stored in a folder named *case_default* created at runtime. The following files are created:
* *logfile.txt* - This file logs the information about the problem solved such as total degrees of freedom, mesh coordinates, node and edge mappings, integration points used, etc.
* *Plots* - Paraview plot files that contain the numerical (and analytical if available) solutions computed over the mesh.
* *error_data* - If analytical solutions are given, then the normed errors in numerical solution are stored in these files.

Towards the end of the run, the ElliFEM code collates all the output files generated and stores them inside the *case_default* directory. By default, at compilation, the code cleans any previous existing folder named *case_default*, including its contents, thus caution must be taken (by renaiming the *case_default* folder, or backing it up, from any previous runs)

## Mesh preparation
The ElliFEM code searches for a file named *dat*, inside the current working directory, during runtime to get the mesh information. If this file is not present in the folder, the code would compile correctly, but would generate an error during runtime (NB: the type of error generated is API dependent, and is generally a segmentation fault. However future development of the code would consist a proper handler for this, and other potential, errors).

To prepare the mesh, create a 2D geometry in [Gmsh](https://gmsh.info/) and export the mesh (only meshed using quadrilateral elements) as the *SU2* format and name it as *dat*. Some sample *dat* files are given in the [Example](https://github.com/md861/ElliFEM/tree/main/Example) folder. Then make sure that this *dat* file (or a copy of it) resides in the same folder as the compiled ElliFEM code. 

Make sure to mark the boundary curves/nodes for each type of boundary condition. The boundaries can be described in any sequence, however do make note of the sequence of their definition as this would provide the index (for *pellib_Hlmhltz.f90* file) to implement the corresponding conditions for each boundadry type. 

NB: At the moment, only non-homogenous Neumann and homogenous Dirichlet conditions may be applied. Although, the future versions would contain the provision for non-homogenous mixed boundaries (even with enriched solution basis).

## Implementing problem sources
### Boundary condition(s)
* Non-homogeneous Neumann boundary: 
    * In *pellib_Hlmhltz.f90* under the "!Integrate over edges" section set the `NBC_pos` index as the same as the corresponding index of the boundary type defined during meshing (see [Mesh preparation](#mesh-preparation) for details). In this sections, update the `GZ_Phi` variable that is used to implement the relevant condition for the given boundary type (chosen using `NBC_pos` index). NB: `n(1)` and `n(2)` store the x and y components of the normal vector to a given edge. 
    * For each different boundary condition (of non-homogeneous Neumann type), copy and paste the entire "!Integrate over edges" section (i.e. all the 4 subsections integrating over the four sides of a quadrilaterl) and set the `NBC_pos` index appropriately. 
* Homogeneous Neumann boundary:
    * Since these boundaries do not require integration, simply do not include any "!Integrate over edges" section for this `NBC_pos` index in the *pellib_Hlmhltz.f90* file.
* Homogeneous Dirichlet boundary: 
   * In *femSolver.f90* under the "!Build marked indices for DBC" section, set the `DBC_pos` index as the same as the corresponding index of the Dirichlet boundary type defined during meshing (see [Mesh preparation](#mesh-preparation) for details). 
### Source term(s)
In the *pellib_Hlmhltz.f90* file, under the "!Integrate inside domain" section, update the `FZ_Phi` variable that evaluates the source terms for the given problem. 

## Implementing numerical parameters
Numerical parameters, e.g. the quadrature points for integration, etc could be set in the *femSolver.f90* file. Some (most used) parameters are detailed below:
* `NGAUSS` - this variable defines the number of quadrature points to be used for line integrals, and the square of this for integrations inside the element.
* `NPOINTPLOT` - similar to `NGAUSS`, defines the number of quadrature points for plotting the numerical results. 
* `NANGL` - number of enrichment functions per node. Set `NANGL = 1` for purely p-FEM solutions.
* `K_W` - wavenumber of the enrichment functions. NB: these are not used if `NANGL = 1`.

 ## Description of files:
 A short summary of some files (which are usually modified the most) is presented here. 
 * *dat* - This is the "su2" format mesh file supplied by the user, generated from Gmsh. The file should be named as "dat" to be read by the solver.
 * *femSolver.f90* - The main code that coordinates the subroutines and functions. This is where you could change some numerical parameters e.g. 
    * the wavenumber of the problem (if solving for Helmholtz), 
    * number of integration points to be used, etc.
 * *pellib_Hlmhltz.f90* - This file allows to specify the boundary sources as well as any sources (`FZ_Phi`) inside the domain.
 * *ln_norm.f90* and *proslib.f90* - These two files are used to specify the analytical solution (if available) for the computation of normed errors and plotting of analytical values over mesh, respectively.
 
## Usage:
All the files should be in the same folder. Open a terminal in the code folder, and type the following for:
* Compilation: `./CleanNCompile`
* Run: `./femSolver`

## Example files:
An example *dat* file that has a 2D mesh with 4-th order quadrilateral elements (generated with [Gmsh](https://gmsh.info/)) is located in the [Example/30pi_p4](https://github.com/md861/ElliFEM/tree/main/Example/30pi_p4) folder. The Paraview plots of the mesh and an example numerical solution for a Helmholtz problem (for wavenumber = 30<img src="http://latex2png.com/pngs/ed7e6344a97275d7ed3b0ab52b9f9eb9.png" width="10" height="10">) with non-zero Neuman (inner) boundaries and Homogeneous Dirichlet (outer) boundaries solved over this mesh, are shown at the beginning of this read-me file. 
 
## References
<a id="1">[1]</a> 
M. Drolia, *et al*. Enriched finite elements for initial-value problem of transverse electromagnetic waves in time domain. *Computers & Structures*, 182, 354-367, 2017.


