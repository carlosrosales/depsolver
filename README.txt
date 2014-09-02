=== depSolver v2.0===

Multimaterial electrostatic solver with focus on dielectrophoretic force calculations. Released under GPLv3 (see file COPYING.txt).


=== INSTALLATION ===

A full installation requires working versions of:

 * gcc
 * make
 * doxygen (Reference Manual Only)

To install depSolver Type:

$ ./install.sh

This will compile all the binary files and move them to the depSolver/bin/ directory. It will also generate the Reference Manual and save it in the depSolver/docs/ directory.


=== CHANGELOG ===

Revision 2.0: (2008-06-16) Second public release of depSolver.
              * Strings read from the input file now have limited size to 
                avoid the possibility of stack overflow.
              * Header style now matches most OSS projects.
              * Related functions are now grouped in files shorter than 1000 lines.
              * Post-processing options include the vtk file format for Paraview.

Revision 1.0: (2006-09-01) First public release of depSolver

2008-06-22 Carlos Rosales Fernandez
