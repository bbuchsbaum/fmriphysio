// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

namespace {

inline arma::vec zscore_vec(const arma::vec& x, const double eps = 1e-12) {
  if (x.n_elem == 0) return arma::vec();
  const double mu = arma::mean(x);
  arma::vec xc = x - mu;
  const double sd = std::sqrt(arma::dot(xc, xc) / std::max(1.0, static_cast<double>(x.n_elem - 1)));
  if (!std::isfinite(sd) || sd < eps) {
    return arma::zeros<arma::vec>(x.n_elem);
  }
  return xc / sd;
}

inline double safe_abs_cor(const arma::vec& a, const arma::vec& b, const double eps = 1e-12) {
  if (a.n_elem != b.n_elem || a.n_elem < 2) return 0.0;
  arma::vec az = zscore_vec(a, eps);
  arma::vec bz = zscore_vec(b, eps);
  double denom = static_cast<double>(a.n_elem - 1);
  if (denom <= 0) return 0.0;
  double r = arma::dot(az, bz) / denom;
  if (!std::isfinite(r)) return 0.0;
  return std::abs(r);
}

} // namespace

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
  arma::vec ratio(L, arma::fill::zeros);
  arma::vec med_nn(L, arma::fill::zeros);
  arma::vec med_nt(L, arma::fill::zeros);
  arma::vec r2_mean(L, arma::fill::zeros);

  if (L == 0) {
    return Rcpp::List::create(
      _["ratio"] = ratio,
      _["med_nn"] = med_nn,
      _["med_nt"] = med_nt,
      _["r2_mean"] = r2_mean
    );
  }

  arma::uvec nn_idx = arma::find(wNN < 0.5);
  arma::uvec nt_idx = arma::find(wNN > 0.5);
  if (nn_idx.n_elem == 0 || nt_idx.n_elem == 0) {
    stop("Need both NN and NT voxels for ratio scoring.");
  }

  arma::vec xnorm2 = arma::sum(arma::square(X_nt), 1);

  for (arma::uword k = 0; k < L; ++k) {
    arma::vec tk = T_tl.col(k);
    const double tnorm2 = arma::dot(tk, tk);
    if (!std::isfinite(tnorm2) || tnorm2 < eps) {
      continue;
    }

    arma::vec cvec = X_nt * tk;
    arma::vec denom = xnorm2 * tnorm2;
    denom.transform([eps](double d) { return std::max(d, eps); });
    arma::vec r2 = arma::square(cvec) / denom;

    const double mn = arma::median(r2.elem(nn_idx));
    const double mt = arma::median(r2.elem(nt_idx));

    med_nn[k] = mn;
    med_nt[k] = mt;
    ratio[k] = mn / std::max(mt, eps);
    r2_mean[k] = arma::mean(r2);
  }

  return Rcpp::List::create(
    _["ratio"] = ratio,
    _["med_nn"] = med_nn,
    _["med_nt"] = med_nt,
    _["r2_mean"] = r2_mean
  );
}

// [[Rcpp::export]]
Rcpp::LogicalVector stability_filter_cpp(
  const arma::mat& T_tl,
  const arma::mat& S1,
  const arma::mat& S2,
  const Rcpp::IntegerVector& i1,
  const Rcpp::IntegerVector& i2,
  const double thresh = 0.30
) {
  if (T_tl.n_cols == 0) return Rcpp::LogicalVector(0);

  const arma::uword L = T_tl.n_cols;
  Rcpp::LogicalVector keep(L, false);

  arma::uvec idx1(i1.size());
  arma::uvec idx2(i2.size());
  for (int j = 0; j < i1.size(); ++j) {
    int id = i1[j] - 1;
    if (id < 0 || id >= static_cast<int>(T_tl.n_rows)) stop("i1 contains out-of-range indices.");
    idx1[j] = static_cast<arma::uword>(id);
  }
  for (int j = 0; j < i2.size(); ++j) {
    int id = i2[j] - 1;
    if (id < 0 || id >= static_cast<int>(T_tl.n_rows)) stop("i2 contains out-of-range indices.");
    idx2[j] = static_cast<arma::uword>(id);
  }

  for (arma::uword k = 0; k < L; ++k) {
    arma::vec tk = T_tl.col(k);
    arma::vec t1 = zscore_vec(tk.elem(idx1));
    arma::vec t2 = zscore_vec(tk.elem(idx2));

    double m1 = 0.0;
    double m2 = 0.0;

    for (arma::uword j = 0; j < S1.n_cols; ++j) {
      double c = safe_abs_cor(t1, S1.col(j));
      if (c > m1) m1 = c;
    }
    for (arma::uword j = 0; j < S2.n_cols; ++j) {
      double c = safe_abs_cor(t2, S2.col(j));
      if (c > m2) m2 = c;
    }

    keep[k] = std::isfinite(m1) && std::isfinite(m2) && (std::min(m1, m2) >= thresh);
  }

  return keep;
}
