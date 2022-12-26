# Optimized Sparse Davidsons Algorithm
The Optimized Davidson's algorithm discussed here is an optimized version of the original Davidson algorithm, which is used to find extreme eigenpairs of Hermitian matrices. 

The FORTRAN 90 subroutines of the ODA and OSDA for both real and complex sparse Hermitian matrices are provided and can be readily used. To include these subroutines in the main program file (i.e., the matrix for which the user wishes to find the extreme eigenpairs), and compile the main program, the following command must be executed:

```
$ ifort myprogram.f90 sparse_davidson.f90 -mkl
```
If a segmentation fault is shown (which may be due to stack memory issues), use the command:

```
$ ifort myprogram.f90 sparse_davidson.f90 -mkl -heap-arrays
```
Note that, myprogram reads the input from the file 'fort.124' (which contains the nonzero elements of the matrix) and stores the eigenvalues in output file 'fort.201'. For more details regarding paramters, read the instruction_for_compilation.pdf
