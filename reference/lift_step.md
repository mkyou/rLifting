# Create an individual Lifting Step

Helper function to create prediction (P) or update (U) steps,
abstracting the complexity of index management.

## Usage

``` r
lift_step(
  type = c("predict", "update"),
  coeffs,
  start_idx = NULL,
  position = "center"
)
```

## Arguments

- type:

  Step type: "predict" (P) or "update" (U).

- coeffs:

  Numeric vector containing the filter coefficients.

- start_idx:

  (Optional) Manual start index. If provided, ignores the `position`
  parameter. Use this for fine-grained control.

- position:

  Automatic index adjustment (used only if `start_idx` is NULL):

  - "center": Centers the filter (default).

  - "left": Causal filter (looks into the past).

  - "right": Anti-causal filter (looks into the future).

## Value

A list formatted for the internal lifting engine.
