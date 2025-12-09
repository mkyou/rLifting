## Resubmission

This is a resubmission. In this version I have:

* Added single quotes around software names (e.g., 'C++') in DESCRIPTION
  and roxygen documentation as requested.

* Added `\value` tags to all exported methods that were missing them:
  - `print.lifting_scheme()`
  - `print.lwt()`
  - `threshold()`
  - `visualize_wavelet_basis()`

* Fixed `visualize_wavelet_basis()` to use `on.exit()` for restoring
  graphical parameters (`par()`) when the function exits.

## Test environments
* local Windows 11 install, R 4.4.1
* win-builder (devel and release)

## R CMD check results
0 errors | 0 warnings | 1 note



## Explanation of Notes
* Spelling: The word 'denoising' in the DESCRIPTION file is a standard
  technical term in Signal Processing and is correct.
