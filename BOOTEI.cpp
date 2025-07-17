// =====================================================================
//  BOOTEI – Bootstrap Ensemble Inference (core C++/Rcpp engine)
//  Implements χ², Spearman and Mann–Whitney with bootstrap‑averaged
//  statistic and permutation p‑value. Supports alternative =
//  "two.sided" | "greater" | "less".
//
// Setting B = 1 disables bootstrap averaging: the test reduces to standard permutation.
//
// Compile in R with Rcpp::sourceCpp("bootei.cpp")
// =====================================================================

#include <Rcpp.h>
#include <algorithm>
#include <unordered_map>
#include <numeric>
#include <cmath>
using namespace Rcpp;

/* ================================================================
 Helper: average ranks (handles ties)
 ================================================================ */
NumericVector rank_numeric(const NumericVector &x) {
  int n = x.size();
  std::vector<std::pair<double,int>> ord(n);
  for (int i = 0; i < n; ++i) ord[i] = { x[i], i };
  std::sort(ord.begin(), ord.end(),
            [](const std::pair<double,int>& a,
               const std::pair<double,int>& b){ return a.first < b.first; });
  NumericVector rk(n);
  int i = 0;
  while (i < n) {
    int j = i;
    while (j + 1 < n && ord[j+1].first == ord[i].first) ++j;
    double r = (i + j + 2) / 2.0;  // average rank, 1‑based
    for (int k = i; k <= j; ++k) rk[ ord[k].second ] = r;
    i = j + 1;
  }
  return rk;
}

/* ================================================================
 χ² for independence (fixed levels to avoid NA on bootstrap)
 ================================================================ */

// builds a contingency table using pre‑computed maps
inline double chi2_from_table(const std::vector<std::vector<int>>& tab,
                              const std::vector<int>& rs,
                              const std::vector<int>& cs,
                              int tot){
  int nr = tab.size(), nc = tab[0].size();
  if (tot == 0) return 0.0;                       // all empty
  double chi2 = 0.0;
  for (int r=0;r<nr;++r) for (int c=0;c<nc;++c){
    double expct = (double)rs[r] * cs[c] / tot;
    if (expct <= 0.0) continue;                  // skip degenerate cells
    double diff = tab[r][c] - expct;
    chi2 += diff*diff/expct;
  }
  return chi2;
}

// factory: returns a lambda that computes χ² using fixed level maps
std::function<double(const CharacterVector&, const CharacterVector&)>
  chisq_fixed_factory(const CharacterVector &x_ref,
                      const CharacterVector &y_ref){
    // discover full set of levels from ORIGINAL data
    std::unordered_map<std::string,int> rmap, cmap;
    int n = x_ref.size();
    for (int i=0;i<n;++i){
      if (CharacterVector::is_na(x_ref[i]) || CharacterVector::is_na(y_ref[i])) continue;
      rmap[ as<std::string>(x_ref[i]) ];
      cmap[ as<std::string>(y_ref[i]) ];
    }
    int nr = rmap.size(), nc = cmap.size();
    // assign sequential indices
    int id = 0; for (auto &kv : rmap) kv.second = id++;
    id = 0;     for (auto &kv : cmap) kv.second = id++;
    
    // capture maps by value inside the lambda
    return [=](const CharacterVector &x, const CharacterVector &y){
      std::vector<std::vector<int>> tab(nr, std::vector<int>(nc, 0));
      std::vector<int> rs(nr,0), cs(nc,0);
      int tot = 0;
      int nloc = x.size();
      for (int i=0;i<nloc;++i){
        if (CharacterVector::is_na(x[i]) || CharacterVector::is_na(y[i])) continue;
        auto it_r = rmap.find( as<std::string>(x[i]) );
        auto it_c = cmap.find( as<std::string>(y[i]) );
        if (it_r == rmap.end() || it_c == cmap.end()) continue; // unseen level
        int r = it_r->second;
        int c = it_c->second;
        tab[r][c]++; rs[r]++; cs[c]++; tot++;
      }
      return chi2_from_table(tab, rs, cs, tot);
    };
  }

/* ================================================================
 Spearman ρ and Mann–Whitney U (as before)
 ================================================================ */

double spearman_stat(const NumericVector &x, const NumericVector &y){
  int n = x.size();
  if (n < 3) return NA_REAL;
  NumericVector rx = rank_numeric(x);
  NumericVector ry = rank_numeric(y);
  double mx = mean(rx), my = mean(ry);
  double num=0.0, dx=0.0, dy=0.0;
  for(int i=0;i<n;++i){
    double cx = rx[i]-mx, cy = ry[i]-my;
    num += cx*cy; dx += cx*cx; dy += cy*cy;
  }
  if (dx==0.0 || dy==0.0) return NA_REAL;
  return num / std::sqrt(dx*dy);
}

double mannwhitney_U(const NumericVector &x, const NumericVector &y){
  int nx = x.size(), ny = y.size();
  int n = nx + ny;
  NumericVector z(n); IntegerVector grp(n);
  for(int i=0;i<nx;++i){ z[i]=x[i]; grp[i]=0; }
  for(int j=0;j<ny;++j){ z[nx+j]=y[j]; grp[nx+j]=1; }
  NumericVector rk = rank_numeric(z);
  double R1 = 0.0;
  for(int i=0;i<n;++i) if (grp[i]==0) R1 += rk[i];
  return R1 - nx*(nx+1)/2.0; // già centrato per two-sided
}

/* ================================================================
 Generic bootstrap helpers
 ================================================================ */

template<typename VECX, typename VECY, typename STATFN>
inline double bagged_paired(const VECX &x, const VECY &y,
                            STATFN statfn, int B){
  if (B == 1) return statfn(x,y);
  RNGScope scope;
  int n = x.size();
  double acc = 0.0; int ok = 0;
  for (int b=0;b<B;++b){
    IntegerVector idx = Rcpp::sample(n, n, true); // 1‑based indices
    VECX xb(n); VECY yb(n);
    for(int i=0;i<n;++i){ xb[i] = x[idx[i]-1]; yb[i] = y[idx[i]-1]; }
    double s = statfn(xb,yb);
    if (!std::isnan(s)) { acc += s; ++ok; }
  }
  if (ok==0) return statfn(x,y);   // fallback to original stat
  return acc / ok;
}

inline double bagged_unpaired(const NumericVector &x, const NumericVector &y,
                              std::function<double(const NumericVector&,const NumericVector&)> statfn,
                              int B){
  if (B == 1) return statfn(x,y);
  RNGScope scope;
  int nx=x.size(), ny=y.size();
  double acc=0.0; int ok=0;
  for(int b=0;b<B;++b){
    IntegerVector idxx = Rcpp::sample(nx, nx, true);
    IntegerVector idxy = Rcpp::sample(ny, ny, true);
    NumericVector xb(nx), yb(ny);
    for(int i=0;i<nx;++i) xb[i] = x[idxx[i]-1];
    for(int j=0;j<ny;++j) yb[j] = y[idxy[j]-1];
    double s = statfn(xb,yb);
    if (!std::isnan(s)) { acc+=s; ++ok; }
  }
  if (ok==0) return statfn(x,y);
  return acc/ok;
}

/* ================================================================
 Permutation helpers
 ================================================================ */

template<typename VECX, typename VECY, typename STATFN>
double perm_paired(const VECX &x, const VECY &y, STATFN statfn,
                   int B, int R, const std::string &alt, double ref=0.0){
  double obs = bagged_paired(x,y,statfn,B);
  if (R<=0 || std::isnan(obs)) return NA_REAL;
  RNGScope scope;
  int n = x.size();
  int geq=0, valid=0;
  for(int r=0;r<R;++r){
    IntegerVector perm = Rcpp::sample(n, n, false);
    VECY yp(n);
    for(int i=0;i<n;++i) yp[i] = y[perm[i]-1];
    double s = bagged_paired(x,yp,statfn,B);
    if (std::isnan(s)) continue;
    ++valid;
    if (alt=="greater")      { if (s >= obs) ++geq; }
    else if (alt=="less")    { if (s <= obs) ++geq; }
    else { // two‑sided
      if (std::fabs(s-ref) >= std::fabs(obs-ref)) ++geq;
    }
  }
  if (valid==0) return NA_REAL;
  return (double)(geq + 1) / (valid + 1);  // add‑one rule
}

double perm_unpaired(const NumericVector &x, const NumericVector &y,
                     std::function<double(const NumericVector&,const NumericVector&)> statfn,
                     int B, int R, const std::string &alt){
  double obs = bagged_unpaired(x,y,statfn,B);
  if (R<=0 || std::isnan(obs)) return NA_REAL;
  RNGScope scope;
  int nx=x.size(), ny=y.size(); int N = nx+ny;
  NumericVector pool(N);
  for(int i=0;i<nx;++i) pool[i]=x[i];
  for(int j=0;j<ny;++j) pool[nx+j]=y[j];
  double ref = nx*ny/2.0;                   // expected U under H0
  int geq=0, valid=0;
  for(int r=0;r<R;++r){
    IntegerVector perm = Rcpp::sample(N, N, false);
    NumericVector xg(nx), yg(ny);
    for(int i=0;i<nx;++i) xg[i] = pool[perm[i]-1];
    for(int j=0;j<ny;++j) yg[j] = pool[perm[nx+j]-1];
    double s = bagged_unpaired(xg,yg,statfn,B);
    if (std::isnan(s)) continue;
    ++valid;
    if (alt=="greater")      { if (s >= obs) ++geq; }
    else if (alt=="less")    { if (s <= obs) ++geq; }
    else { // two‑sided
      if (std::fabs(s-ref) >= std::fabs(obs-ref)) ++geq;
    }
  }
  if (valid==0) return NA_REAL;
  return (double)(geq + 1) / (valid + 1);
}

/* ================================================================
 Exported wrapper
 ================================================================ */
// [[Rcpp::export]]
List bootei(SEXP x, SEXP y,
            std::string test = "chisq",
            int B = 100, int R = 100,
            std::string alternative = "two.sided"){
  
  if (alternative != "two.sided" && alternative != "greater" && alternative != "less")
    stop("alternative must be 'two.sided', 'greater', or 'less'");
  
  if (test == "chisq") {
    CharacterVector xc(x), yc(y);
    if (xc.size() != yc.size())
      stop("For chisq test x and y must have equal length.");
    
    // build statfn with fixed levels
    auto statfn = chisq_fixed_factory(xc,yc);
    
    double stat = bagged_paired(xc,yc,statfn,B);
    double pval = perm_paired(xc,yc,statfn,B,R,"greater"); // χ² is one‑sided
    return List::create(_["statistic"] = stat,
                        _["p.value"]  = pval,
                        _["method"]   = "BOOTEI χ² test",
                        _["alternative"] = "greater");
  }
  
  else if (test == "spearman") {
    NumericVector xn(x), yn(y);
    if (xn.size() != yn.size())
      stop("For spearman test x and y must have equal length.");
    auto statfn = [&](const NumericVector &a, const NumericVector &b){ return spearman_stat(a,b); };
    
    double stat = bagged_paired(xn,yn,statfn,B);
    double pval = perm_paired(xn,yn,statfn,B,R,alternative,0.0);
    return List::create(_["statistic"] = stat,
                        _["p.value"]  = pval,
                        _["method"]   = "BOOTEI Spearman test",
                        _["alternative"] = alternative);
  }
  
  else if (test == "mannwhitney") {
    NumericVector x1(x), y1(y);
    auto statfn = [&](const NumericVector &a, const NumericVector &b){ return mannwhitney_U(a,b); };
    
    double stat = bagged_unpaired(x1,y1,statfn,B);
    double pval = perm_unpaired(x1,y1,statfn,B,R,alternative);
    return List::create(_["statistic"] = stat,
                        _["p.value"]  = pval,
                        _["method"]   = "BOOTEI Mann–Whitney test",
                        _["alternative"] = alternative);
  }
  
  else {
    stop("Unknown test type. Use 'chisq', 'spearman', or 'mannwhitney'.");
  }
}
