# Validate Compact Support (FIR Compliance)

Verifies if the impulse response is finite (FIR Filter).

## Usage

``` r
validate_compact_support(scheme, max_width)
```

## Arguments

- scheme:

  Object of class `lifting_scheme`.

- max_width:

  Maximum expected width (number of taps).

## Value

List with status and number of active taps.
