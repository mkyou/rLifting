# Changelog

## rLifting 0.9.1

- Updated Roadmap to incorporate dissertation findings (irregular grids,
  multivariate denoising, and boundary refinements).
- Corrected complexity notation: removed misleading $O(1)$ mentions,
  clarifying the efficient $O(T)$ nature of the ring-buffer
  architecture.
- Improved documentation order: vignettes are now correctly ordered in
  the documentation site.

## rLifting 0.9.0

CRAN release: 2026-02-20

- Initial release of the package.
- Implements high-performance Wavelet Lifting Transforms.
- Features:
  - Unified offline (batch) denoising.
  - Causal (real-time) filtering using a ring buffer engine.
  - Adaptive recursive thresholding.
  - Zero-allocation ‘C++’ core.
