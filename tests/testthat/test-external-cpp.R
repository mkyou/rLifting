
test_that("External C++ package can headers from rLifting", {
  skip_on_cran()

  if (!requireNamespace("Rcpp", quietly = TRUE)) {
    skip("Rcpp not available")
  }

  # This test requires rLifting to be installed so Rcpp can find headers
  # via [[Rcpp::depends(rLifting)]]
  if (!requireNamespace("rLifting", quietly = TRUE)) {
    skip("rLifting not installed")
  }

  cpp_src = "
  // [[Rcpp::depends(rLifting)]]
  #include <Rcpp.h>
  #include <rLifting/WaveletEngine.h>

  using namespace Rcpp;

  // [[Rcpp::export]]
  bool check_engine_compiles() {
    // Just verifying the type and class are visible and linkable
    // We instantiate a pointer to verify the type definition is complete
    WaveletEngine* engine = nullptr;
    return (engine == nullptr);
  }
  "

  # Attempt to compile and run
  # We use tryCatch to properly report compilation errors as test failures
  result = tryCatch(
    {
      Rcpp::sourceCpp(code = cpp_src, env = environment())
      check_engine_compiles()
    },
    error = function(e) {
      fail(paste("Failed to compile external C++ source:", e$message))
      FALSE
    })

  expect_true(result)
})
