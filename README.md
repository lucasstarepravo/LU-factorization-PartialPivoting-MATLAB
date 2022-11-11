#Low-memory LU-factorization-PartialPivoting-MATLAB

This code performs low-memory LU factorization partial pivoting

Gaussian weights and upper triangular matrix are both stored on the output matrix, while the permutation matrix is given in vector format.

Python and MATLAB versions are available. For the Python version, please go to corresponding repository

# Example
```ruby
# given a positive definite matrix A and a vector b which A*x=b

[A,P] = lu_pp(A)
b = lt_solve(A,b)
b = ut_solve(A,b)
```
