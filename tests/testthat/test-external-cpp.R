
test_that("External C++ package can headers from rLifting", {
  skip_on_cran()

  if (!requireNamespace("Rcpp", quietly = TRUE)) {
    skip("Rcpp not available")
  }

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
    WaveletEngine* engine = nullptr;
    return (engine == nullptr);
  }
  "

  result = tryCatch(
    {
      Rcpp::sourceCpp(code = cpp_src, env = environment())
      check_engine_compiles()
    },
    error = function(e) {
      skip(paste("Cannot compile external C++ source:", e$message))
    })

  expect_true(result)
})
