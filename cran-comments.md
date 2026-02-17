## Resubmission (v0.9.0)
This is a resubmission.

### Response to previous Human Review (Benjamin Altmann):
* Added single quotes around software names (e.g., 'C++') in DESCRIPTION.
* Added \value tags to all exported methods that were missing them.
* Fixed visualize_wavelet_basis() to use on.exit() for restoring graphical parameters.

### Recent improvements (v0.9.0):
* Added counterfactual leakage tests to demonstrate causal correctness.
* Expanded vignettes with real-time streaming examples and custom wavelet construction.
* Removed unused dependencies from Suggests.
* Cleaned up orphan files.

### Test environments
* local Windows 11 install, R 4.4.1
* win-builder (devel and release)
* R-hub (Linux, macOS, Windows)

### R CMD check results
0 errors | 0 warnings | 1 note

### Explanation of Notes
* Spelling: The word 'denoising' in the DESCRIPTION file is a standard
  technical term in Signal Processing and is correct.
