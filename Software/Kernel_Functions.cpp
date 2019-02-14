// load Rcpp
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

//================================================================================================
//(1) Listed below are wrapper functions to sample from different distributions

// [[Rcpp::export]]
//random multivariate normal sample generator using RcppArmadillo
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

// [[Rcpp::export]]
double rgammadouble(int a, double b, double c){
  Rcpp::NumericVector x = rgamma(a,b,1/c);
  return x(0);
}

// [[Rcpp::export]]
double rInvGamma(int n, double shape,double scale){
  double x = rgammadouble(n, shape, scale);
  return 1.0/x;
}

// [[Rcpp::export]]
double rScaledInvChiSq(int n, double nu, double tau2){
  double x = rInvGamma(n, nu/2,(nu * tau2)/2);
  return x;
}

// [[Rcpp::export]]
double median_rcpp(NumericVector x) {
  NumericVector y = clone(x);
  int n, half;
  double y1, y2;
  n = y.size();
  half = n / 2;
  if(n % 2 == 1) {
    // median for odd length vector
    std::nth_element(y.begin(), y.begin()+half, y.end());
    return y[half];
  } else {
    // median for even length vector
    std::nth_element(y.begin(), y.begin()+half, y.end());
    y1 = y[half];
    std::nth_element(y.begin(), y.begin()+half-1, y.begin()+half);
    y2 = y[half-1];
    return (y1 + y2) / 2.0;
  }
}

// [[Rcpp::export]]
double mad_rcpp(NumericVector x, double scale_factor = 1.4826) {
  // scale_factor = 1.4826; default for normal distribution consistent with R
  return median_rcpp(abs(x - median_rcpp(x))) * scale_factor;
}

// norm_rs(a, b)
// generates a sample from a N(0,1) RV restricted to be in the interval
// (a,b) via rejection sampling.
// ======================================================================

// [[Rcpp::export]]

double norm_rs(double a, double b)
{
    double  x;
    x = Rf_rnorm(0.0, 1.0);
    while( (x < a) || (x > b) ) x = norm_rand();
    return x;
}

// half_norm_rs(a, b)
// generates a sample from a N(0,1) RV restricted to the interval
// (a,b) (with a > 0) using half normal rejection sampling.
// ======================================================================

// [[Rcpp::export]]

double half_norm_rs(double a, double b)
{
    double   x;
    x = fabs(norm_rand());
    while( (x<a) || (x>b) ) x = fabs(norm_rand());
    return x;
}

// unif_rs(a, b)
// generates a sample from a N(0,1) RV restricted to the interval
// (a,b) using uniform rejection sampling.
// ======================================================================

// [[Rcpp::export]]

double unif_rs(double a, double b)
{
    double xstar, logphixstar, x, logu;
    
    // Find the argmax (b is always >= 0)
    // This works because we want to sample from N(0,1)
    if(a <= 0.0) xstar = 0.0;
    else xstar = a;
    logphixstar = R::dnorm(xstar, 0.0, 1.0, 1.0);
    
    x = R::runif(a, b);
    logu = log(R::runif(0.0, 1.0));
    while( logu > (R::dnorm(x, 0.0, 1.0,1.0) - logphixstar))
    {
        x = R::runif(a, b);
        logu = log(R::runif(0.0, 1.0));
    }
    return x;
}

// exp_rs(a, b)
// generates a sample from a N(0,1) RV restricted to the interval
// (a,b) using exponential rejection sampling.
// ======================================================================

// [[Rcpp::export]]

double exp_rs(double a, double b)
{
    double  z, u, rate;
    
    //  Rprintf("in exp_rs");
    rate = 1/a;
    //1/a
    
    // Generate a proposal on (0, b-a)
    z = R::rexp(rate);
    while(z > (b-a)) z = R::rexp(rate);
    u = R::runif(0.0, 1.0);
    
    while( log(u) > (-0.5*z*z))
    {
        z = R::rexp(rate);
        while(z > (b-a)) z = R::rexp(rate);
        u = R::runif(0.0,1.0);
    }
    return(z+a);
}

// rnorm_trunc( mu, sigma, lower, upper)
//
// generates one random normal RVs with mean 'mu' and standard
// deviation 'sigma', truncated to the interval (lower,upper), where
// lower can be -Inf and upper can be Inf.
//======================================================================

// [[Rcpp::export]]
double rnorm_trunc (double mu, double sigma, double lower, double upper)
{
    int change;
    double a, b;
    double logt1 = log(0.150), logt2 = log(2.18), t3 = 0.725;
    double z, tmp, lograt;
    
    change = 0;
    a = (lower - mu)/sigma;
    b = (upper - mu)/sigma;
    
    // First scenario
    if( (a == R_NegInf) || (b == R_PosInf))
    {
        if(a == R_NegInf)
        {
            change = 1;
            a = -b;
            b = R_PosInf;
        }
        
        // The two possibilities for this scenario
        if(a <= 0.45) z = norm_rs(a, b);
        else z = exp_rs(a, b);
        if(change) z = -z;
    }
    // Second scenario
    else if((a * b) <= 0.0)
    {
        // The two possibilities for this scenario
        if((R::dnorm(a, 0.0, 1.0,1.0) <= logt1) || (R::dnorm(b, 0.0, 1.0, 1.0) <= logt1))
        {
            z = norm_rs(a, b);
        }
        else z = unif_rs(a,b);
    }
    // Third scenario
    else
    {
        if(b < 0)
        {
            tmp = b; b = -a; a = -tmp; change = 1;
        }
        
        lograt = R::dnorm(a, 0.0, 1.0, 1.0) - R::dnorm(b, 0.0, 1.0, 1.0);
        if(lograt <= logt2) z = unif_rs(a,b);
        else if((lograt > logt1) && (a < t3)) z = half_norm_rs(a,b);
        else z = exp_rs(a,b);
        if(change) z = -z;
    }
    double output;
    output = sigma*z + mu;
    return (output);
}

//================================================================================================
//(2) Listed below are functions to create different kernel matrices (e.g. approximate kernel, Gaussian kernel, linear or additive kernel)

// [[Rcpp::export]]
arma::mat ApproxGaussKernel(arma::mat X, double iter, double h){
    int i;
    double ncov = X.n_rows, samp_size = X.n_cols;
    mat zeta(iter,samp_size);
    for(i = 0; i<iter; i++){
        vec omega = arma::randn(ncov)*sqrt(2/h);
        vec b = runif(1,0,2*datum::pi);
        zeta.row(i) = sqrt(2/iter)*cos(trans(omega)*X+as_scalar(b));
    }
    mat K_hat = zeta.t()*zeta;
    return K_hat;
}

// [[Rcpp::export]]
arma::mat GaussKernel(arma::mat X, double h = 1){
    int i,j;
    double n = X.n_cols;
    double p = X.n_rows;
    mat K = zeros<mat>(n,n);
    for (i = 0; i<n; i++){
        for(j = 0; j<n; j++){
            if(i==j){
                break;
            }
            else{
                K(i,j) = exp(-h/(p)*sum(pow(X.col(i)-X.col(j),2)));
            }
        }
    }
    return K + trans(K);
}

// [[Rcpp::export]]
arma::mat CauchyKernel(arma::mat X, double h = 1){
    int i,j;
    double n = X.n_cols;
    double p = X.n_rows;
    mat K = zeros<mat>(n,n);
    for (i = 0; i<n; i++){
        for(j = 0; j<n; j++){
            if(i==j){
                break;
            }
            else{
                K(i,j) = 1/(1+h/(p)*sum(pow(X.col(i)-X.col(j),2)));
            }
        }
    }
    return K + trans(K);
}

// [[Rcpp::export]]
arma::mat LogKernel(arma::mat X, double h = 1){
    int i,j;
    double n = X.n_cols;
    mat K = zeros<mat>(n,n);
    for (i = 0; i<n; i++){
        for(j = 0; j<n; j++){
            if(i==j){
                break;
            }
            else{
                K(i,j) = -log(sum(pow(X.col(i)-X.col(j),h))+1);
            }
        }
    }
    return K + trans(K);
}

// [[Rcpp::export]]
arma::mat LinearKernel(arma::mat X,double h = 1){
    return X.t()*X/h;
}

// [[Rcpp::export]]
arma::mat SigmoidKernel(arma::mat X, double alpha = 1,double h = 1){
    return tanh(X.t()*X/alpha+h);
}

// [[Rcpp::export]]
List EigDecomp(arma::mat X){
    mat U;
    vec lambda;
    mat V;
    svd(U,lambda,V,X);
    return Rcpp::List::create(Rcpp::Named("U") = U,Rcpp::Named("V") = V,Rcpp::Named("lambda") = lambda);
}

// [[Rcpp::export]]
arma::mat ComputePCs(arma::mat X,int top = 10){
    mat U;
    vec s;
    mat V;
    svd(U,s,V,X);
    
    mat PCs = U*diagmat(s);
    return PCs.cols(0,top-1);
}
