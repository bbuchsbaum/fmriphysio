// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include <algorithm>
#include <vector>

using namespace Rcpp;

namespace {

// Median of the finite entries of r2[idx]; NA when none are finite.
// Matches stats::median (mean of the two middle values for even counts).
double finite_median(const arma::vec& r2, const arma::uvec& idx) {
  std::vector<double> vals;
  vals.reserve(idx.n_elem);
  for (arma::uword i = 0; i < idx.n_elem; ++i) {
    const double v = r2[idx[i]];
    if (std::isfinite(v)) vals.push_back(v);
  }
  const std::size_t n = vals.size();
  if (n == 0) return NA_REAL;
  const std::size_t mid = n / 2;
  std::nth_element(vals.begin(), vals.begin() + mid, vals.end());
  const double upper = vals[mid];
  if (n % 2 == 1) return upper;
  const double lower = *std::max_element(vals.begin(), vals.begin() + mid);
  return 0.5 * (lower + upper);
}

} // namespace

// Voxelwise squared correlation of each component with the (row-centred)
// data, summarised as median R^2 in non-neuronal (wNN < 0.5) and neuronal
// (wNN > 0.5) voxels. Expects X_nt rows and T_tl columns already centred.
// Zero-variance voxels and components are judged relative to the largest
// one (tolerance `eps`), matching r2_maps() in R.
// [[Rcpp::export]]
Rcpp::List score_components_nn_nt_cpp(
  const arma::mat& X_nt,
  const arma::mat& T_tl,
  const arma::vec& wNN,
  const double eps = 1e-12
) {
  if (X_nt.n_rows != wNN.n_elem) {
    stop("Length of wNN must match voxel dimension.");
  }
  if (X_nt.n_cols != T_tl.n_rows) {
    stop("X_nt columns must match T_tl rows.");
  }

  const arma::uword L = T_tl.n_cols;
  NumericVector ratio(L, NA_REAL);
  NumericVector med_nn(L, NA_REAL);
  NumericVector med_nt(L, NA_REAL);
  NumericVector r2_mean(L, NA_REAL);

  if (L == 0) {
    return List::create(
      _["ratio"] = ratio, _["med_nn"] = med_nn,
      _["med_nt"] = med_nt, _["r2_mean"] = r2_mean
    );
  }

  const arma::uvec nn_idx = arma::find(wNN < 0.5);
  const arma::uvec nt_idx = arma::find(wNN > 0.5);
  if (nn_idx.n_elem == 0 || nt_idx.n_elem == 0) {
    stop("Need both non-neuronal (wNN < 0.5) and neuronal (wNN > 0.5) voxels for ratio scoring.");
  }

  const arma::vec xnorm2 = arma::sum(arma::square(X_nt), 1);
  const arma::rowvec tnorm2_all = arma::sum(arma::square(T_tl), 0);
  const double x_tol = eps * std::max(xnorm2.max(), 0.0);
  const double t_tol = eps * std::max(tnorm2_all.max(), 0.0);
  const arma::uvec x_null = arma::find(xnorm2 <= x_tol);

  for (arma::uword k = 0; k < L; ++k) {
    const double tnorm2 = tnorm2_all[k];
    if (!std::isfinite(tnorm2) || tnorm2 <= t_tol) continue;

    const arma::vec cvec = X_nt * T_tl.col(k);
    arma::vec r2 = arma::square(cvec) / (xnorm2 * tnorm2);
    r2.elem(x_null).fill(arma::datum::nan);

    const double mn = finite_median(r2, nn_idx);
    const double mt = finite_median(r2, nt_idx);
    med_nn[k] = mn;
    med_nt[k] = mt;
    if (!NumericVector::is_na(mn) && !NumericVector::is_na(mt)) {
      ratio[k] = mn / std::max(mt, eps);
    }

    const arma::vec finite_r2 = r2.elem(arma::find_finite(r2));
    if (finite_r2.n_elem > 0) r2_mean[k] = arma::mean(finite_r2);
  }

  return List::create(
    _["ratio"] = ratio, _["med_nn"] = med_nn,
    _["med_nt"] = med_nt, _["r2_mean"] = r2_mean
  );
}
