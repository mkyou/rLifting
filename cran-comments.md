## Resubmission (v0.9.0)
This is a resubmission addressing the installation failure on CRAN pre-tests.

### Changes since last submission:
* Fixed signed/unsigned integer comparison warnings in C++ source files that
  caused compilation errors on GCC 14.3 (Windows) and Clang 21.1 (Debian)
  with `-Wall -pedantic` flags.

### Response to previous Human Review (Benjamin Altmann):
* Added single quotes around software names (e.g., 'C++') in DESCRIPTION.
* Added \value tags to all exported methods that were missing them.
* Fixed visualize_wavelet_basis() to use on.exit() for restoring graphical parameters.

### Test environments
* local Windows 11 install, R 4.4.1
* win-builder (devel and release)
* R-hub (Linux, macOS, Windows)

### R CMD check results
0 errors | 0 warnings | 1 note

### Explanation of Notes
* Spelling: The word 'denoising' in the DESCRIPTION file is a standard
  technical term in Signal Processing and is correct.
