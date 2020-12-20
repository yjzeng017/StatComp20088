#include <Rcpp.h>
using namespace Rcpp;

//' @title Laplace density
//' @description Laplace density
//' @param x real value
//' @return the density at x
//' @export
// [[Rcpp::export]]
double denLa(double x){
  return exp(-abs(x));
}

#include <Rcpp.h>
using namespace Rcpp;

//' @title Random walk
//' @description  Random walk with Metropolis sampling method using Rcpp
//' @param sigma variance
//' @param x0 initial value
//' @param N sample size
//' @return a random sample of size N
//' @export

//[[Rcpp::export]]
List RWcpp(double sigma, double x0, int N) {
  NumericVector x(N);
  x[0]=x0;
  NumericVector u(N);
  u=runif(N);
  int k=0;
  for(int i=1; i<N; i++){
    double y;
    y=rnorm(1,x[i-1],sigma)[0];
    if(u[i]<=(denLa(y)/denLa(x[i-1])))
      x[i]=y;
    else{
      x[i]=x[i-1];
      k++;
    }
  }
  List out;
  out["x"]=x;
  out["k"]=k;
  return out;
  //return x;
}
