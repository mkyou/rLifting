# Create a custom wavelet

Wrapper to create a `lifting_scheme` object from manual steps.

## Usage

``` r
custom_wavelet(name, steps, norm = c(1, 1))
```

## Arguments

- name:

  Identifier name for the wavelet.

- steps:

  List of steps created via `lift_step`.

- norm:

  Normalization vector c(K, 1/K).

## Value

An object of class `lifting_scheme`.

## Examples

``` r
p1 = lift_step("predict", c(1), position = "center")
u1 = lift_step("update", c(0.5), position = "center")
w = custom_wavelet("HaarManual", list(p1, u1), c(1.41, 0.707))
```
