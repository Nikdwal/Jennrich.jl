# Jennrich

Compute the canonical polyadic decomposition or (L,L,1) block term decomposition of a third order tensor using a pencil-based alorithm, also known as Jennrich's algorithm.

## Usage
Use `jennrich(A, r)` for the rank `r` CPD or `jennrich(A, L)` for the `(L[1],L[1],1), ..., (L[R], L[R], 1)`-BTD of `A`.
The result is returned in factor matrix format.
```
𝒯 = reshape(sum([kron(randn(9), randn(8), randn(7)) for i = 1:r]), (7,8,9))
A, B, C = jennrich(𝒯, r)
```

## References
```
@article{DeLathauwer2008,
  author = {{De Lathauwer}, Lieven},
  doi = {10.1137/070690729},
  journal = {SIAM Journal on Matrix Analysis and Applications},
  number = {3},
  title = {{Decompositions of a higher-order tensor in block terms - Part II: Definitions and uniqueness}},
  url = {http://epubs.siam.org/doi/10.1137/070690729},
  volume = {30},
  year = {2008}
}

@article{harshman1970foundations,
  title={Foundations of the PARAFAC procedure: Models and conditions for an" explanatory" multimodal factor analysis},
  author={Harshman, Richard A and others},
  year={1970},
  publisher={University of California at Los Angeles Los Angeles, CA}
}
```
 
