# Create an individual Lifting Step

Helper function to create prediction (P) or update (U) steps,
abstracting the complexity of index management.

## Usage

``` r
lift_step(
  type = c("predict", "update"),
  coeffs,
  start_idx = NULL,
  position = "center",
  degree = NULL
)
```

## Arguments

- type:

  Step type: `"predict"` (P) or `"update"` (U).

- coeffs:

  Numeric vector containing the filter coefficients.

- start_idx:

  (Optional) Manual start index. If provided, ignores the `position`
  parameter. Use this for fine-grained control. The filter reads
  neighbours at offsets `start_idx + 0..(length(coeffs) - 1)` relative
  to the current index.

- position:

  Automatic index adjustment, used only when `start_idx` is `NULL`:

  - `"center"`: centres the filter (default).
    `start_idx = -floor((length(coeffs) - 1) / 2)`.

  - `"left"`: causal filter (looks into the past).
    `start_idx = -length(coeffs) + 1`.

  - `"right"`: anti-causal filter (looks into the future).
    `start_idx = 0`.

- degree:

  (Optional) Polynomial degree the predict step reproduces exactly.
  Drives the irregular-grid Lagrange interpolation (see
  [`vignette("v06-extensions")`](https://mkyou.github.io/rLifting/articles/v06-extensions.md)).
  If `NULL`, inferred as `length(coeffs) - 1` when `type == "predict"`
  and `sum(coeffs) == 1` (interpolating filter); otherwise `-1` (filter
  not interpretable as polynomial interpolation, e.g. CDF 9/7 or DB2
  predicts). Update steps always carry `degree = -1`.

## Value

A list `list(type, coeffs, start_idx, degree)` formatted for the
internal lifting engine.
