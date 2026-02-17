# Lifting Wavelet Transform (Forward)

Performs the Forward Wavelet Transform using the Lifting Scheme.
Optimized with 'C++' backend.

## Usage

``` r
lwt(signal, scheme, levels = 1, extension = "symmetric")
```

## Arguments

- signal:

  Numeric vector containing the input signal.

- scheme:

  A `lifting_scheme` object.

- levels:

  Integer. Number of decomposition levels.

- extension:

  Boundary extension mode: "symmetric" (default), "periodic", or "zero".

## Value

An object of class 'lwt'. It is a list containing 'coeffs' (list of
details d1..dn and approximation an) and 'scheme' (the scheme object
used).

## Examples

``` r
data = c(1, 2, 3, 4, 5, 6, 7, 8)
sch = lifting_scheme("haar")
res = lwt(data, sch, levels = 2)
#> Warning: Residual signal at level 2 has only 2.0 samples.
print(res)
#> --- LWT Decomposition (C++ Accelerated) ---
#> Levels: 2
#> Wavelet: haar
#> Coefficients:
#>   a2: length 2
#>   d1: length 4
#>   d2: length 2
```
