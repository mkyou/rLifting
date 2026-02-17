# Validate Shift Sensitivity (Shift Variance)

Decimated wavelets are not translation invariant. This test quantifies
the variation in detail energy when shifting the input signal by 1
sample.

## Usage

``` r
validate_shift_sensitivity(scheme)
```

## Arguments

- scheme:

  Object of class `lifting_scheme`.

## Value

List with status and percentage variation.
