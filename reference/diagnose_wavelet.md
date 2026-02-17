# Complete Wavelet Diagnosis

Runs a battery of physical and mathematical tests on a wavelet.

## Usage

``` r
diagnose_wavelet(wavelet_name, config, verbose = TRUE, plot = TRUE)
```

## Arguments

- wavelet_name:

  Name string or a `lifting_scheme` object.

- config:

  Configuration list (is_ortho, vm_degrees, max_taps).

- verbose:

  Print results to console handling? (Defaults to TRUE).

- plot:

  Boolean. Visualize basis functions during diagnosis? (Defaults to
  TRUE).

## Value

An object of class `wavelet_diagnosis` (S3), which is a list containing
the results of each test. The object has a dedicated `print` method.
