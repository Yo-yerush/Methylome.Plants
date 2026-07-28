#include <Rcpp.h>
#include <vector>
#include <limits>

using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

namespace {

inline bool is_na_like(double x) {
  return ISNAN(x);
}

}  // namespace

// Compute mean(values) by 1-based group id with NA handling compatible with
// aggregate(..., mean, na.rm=TRUE):
// - no hits for a bin -> NA
// - hits exist but all values NA/NaN -> NaN
// [[Rcpp::export]]
NumericVector meta_cpp_mean_by_group(
    const IntegerVector& groups,
    const NumericVector& values,
    int n_groups) {
  if (groups.size() != values.size()) {
    stop("groups and values must have the same length");
  }
  if (n_groups < 0) {
    stop("n_groups must be >= 0");
  }

  NumericVector sums(n_groups, 0.0);
  IntegerVector counts(n_groups, 0);
  IntegerVector seen(n_groups, 0);

  for (R_xlen_t i = 0; i < groups.size(); ++i) {
    const int g = groups[i];
    if (g < 1 || g > n_groups) {
      continue;
    }

    const int idx = g - 1;
    seen[idx] += 1;

    const double v = values[i];
    if (!is_na_like(v)) {
      sums[idx] += v;
      counts[idx] += 1;
    }
  }

  NumericVector out(n_groups, NA_REAL);
  for (int i = 0; i < n_groups; ++i) {
    if (seen[i] == 0) {
      out[i] = NA_REAL;
    } else if (counts[i] == 0) {
      out[i] = std::numeric_limits<double>::quiet_NaN();
    } else {
      out[i] = sums[i] / static_cast<double>(counts[i]);
    }
  }

  return out;
}

// Row means with na.rm=TRUE semantics:
// if a row has only NA/NaN -> NaN
// [[Rcpp::export]]
NumericVector meta_cpp_row_means_ignore_na(const NumericMatrix& mat) {
  const int nr = mat.nrow();
  const int nc = mat.ncol();

  NumericVector out(nr, std::numeric_limits<double>::quiet_NaN());
  if (nc == 0) {
    return out;
  }

  for (int i = 0; i < nr; ++i) {
    double acc = 0.0;
    int cnt = 0;

    for (int j = 0; j < nc; ++j) {
      const double v = mat(i, j);
      if (!is_na_like(v)) {
        acc += v;
        cnt += 1;
      }
    }

    if (cnt > 0) {
      out[i] = acc / static_cast<double>(cnt);
    }
  }

  return out;
}
