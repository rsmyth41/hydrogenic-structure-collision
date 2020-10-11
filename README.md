# Hydrogenic structure and electron impact excitation

## Project overview
This is a project which models the atomic structure of an arbirary single electron atom or ion and describes the electron-impact excitation of that model for any incident electron energy. It works for any reasonable arbitrary central potential. However, a Colulomb potential is the default and is currently the only potential which is coded into the project. 

The `structure.f95` part of the project calculates the energy levels (in Rydbergs) and radial functions (using the Numerov numerical method for second order differential equations) for two arbitrary gross structure levels, for any given (`n`, `l`) pair. It also calculates electric dipole (E1) radiative transition data (radiative transition rates, line strengths, oscillator strengths) between those states (where the trapezoidal rule is used for single integrations when necessary) and prints the data for the user. The calculated radial atomic orbitals are also written to files for use in the collision part of this project. 

The `collision.f95` part of this project uses the radial parts of the atomic orbitals written from the atomic structure calculation to calculate the electric dipole collision strength for the transition. The code makes use of the First-Born Approximation and the trapezoidal rule for the double integrations. The collision strength for each energy point is written to a file for easy plotting.

## Compiling and executing
To compile the atomic structure part of the project:

    gfortran subroutines/subroutinesMod.f95 structure.f95 -o structure.x -mcmodel=medium -O3

To run: 

    ./structure.x

The user is presented with an interactive calculation where they must provide a principal quantum number `n` and an angular momentum quantum number `l`. These are integers which must satisfy `n > 0` and `l >= n`.

Similarly, to compile the electron impact excitation part of the project:

    gfortran collision.f95 -o collision.x -mcmodel=medium -O3

To run: 

    ./collision.x
