#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

inline double d(const double A, const double B, const double y, const double T, const double q) {
  return A * exp(B) * y * exp(T * log(q) / 10);
}

// [[Rcpp::export]]
Rcpp::List derivs(double t, NumericVector y, NumericVector parms) {
  
  const double A1 = parms[0];
  const double B1 = parms[1];
  const double q1 = parms[2];
  const double A2 = parms[3];
  const double B2 = parms[4];
  const double q2 = parms[5];
  const double A3 = parms[6];
  const double B3 = parms[7];
  const double q3 = parms[8];
  const double A4 = parms[9];
  const double B4 = parms[10];
  const double q4 = parms[11];
  const double A11 = parms[12];
  const double B11 = parms[13];
  const double q11 = parms[14];
  const double A21 = parms[15];
  const double B21 = parms[16];
  const double q21 = parms[17];
  const double A31 = parms[18];
  const double B31 = parms[19];
  const double q31 = parms[20];
  const double A41 = parms[21];
  const double B41 = parms[22];
  const double q41 = parms[23];
  const double A51 = parms[24];
  const double B51 = parms[25];
  const double q51 = parms[26];
  const double A52 = parms[27];
  const double B52 = parms[28];
  const double q52 = parms[29];
  const double A53 = parms[30];
  const double B53 = parms[31];
  const double q53 = parms[32];
  const double A61 = parms[33];
  const double B61 = parms[34];
  const double q61 = parms[35];
  const double A62 = parms[36];
  const double B62 = parms[37];
  const double q62 = parms[38];
  const double A63 = parms[39];
  const double B63 = parms[40];
  const double q63 = parms[41];
  const double Kmax = parms[42];
  const double Ku = parms[43];
  const double n = parms[44];
  const double halfmax = parms[45];
  const double Kt = parms[46];
  const double Vhalf = parms[47];
  const double T = parms[48];
  const double timeout = parms[49];
  const double starttime = parms[50];
  
  const double T20 = T - 20;
  
  const double dA21 = d(A21, B21 * y[10], y[1], T20, q21);
  const double dA11 = d(A11, B11 * y[10], y[0], T20, q11);
  
  const double dA2  = d(A2, B2 * y[10], y[3], T20, q2);
  const double dA1  = d(A1, B1 * y[10], y[2], T20, q1);
  
  const double dA51 = d(A51, B51 * y[10], y[2], T20, q51);
  const double dA61 = d(A61, B61 * y[10], y[0], T20, q61);
  
  const double dA3  = d(A3, B3 * y[10], y[1], T20, q3);
  const double dA4  = d(A4, B4 * y[10], y[5], T20, q4);
  
  const double dA52 = d(A52, B52 * y[10], y[3], T20, q52);
  const double dA62 = d(A62, B62 * y[10], y[1], T20, q62);
  
  const double dA31 = d(A31, B31 * y[10], y[3], T20, q31);
  const double dA41 = d(A41, B41 * y[10], y[4], T20, q41);
  
  const double dA53 = d(A53, B53 * y[10], y[4], T20, q53);
  const double dA63 = d(A63, B63 * y[10], y[5], T20, q63);
  
  const double dA53_2 = A53 * exp(B53 * y[10]) * exp(T20 * log(q53) / 10);
  const double dA63_2 = A63 * exp(B63 * y[10]) * exp(T20 * log(q63) / 10);
  
  const double dKmax = Kmax * Ku * exp(n * log(y[9])) / (exp(n * log(y[9])) + halfmax);
  const double ddKmax = dKmax * y[5]-Ku * dA53_2/(dA63_2)*y[6];
  
  const double dKt = Kt / (1 + exp(-(y[10] - Vhalf) / 6.789));
  
  NumericVector ydot(11);
  
  ydot[0] = -(dA11 - dA21) + (dA51 - dA61);
  ydot[1] = (dA11 - dA21) - (dA3 - dA4) + (dA52 - dA62);
  ydot[2] = -(dA1 - dA2) - (dA51 - dA61);
  ydot[3] = (dA1 - dA2) - (dA31 - dA41) - (dA52 - dA62);
  ydot[4] = (dA31 - dA41) - (dA53 - dA63) - (dKmax * y[4] - Ku * y[7]);
  ydot[5] = (dA3 - dA4) + (dA53 - dA63) - ddKmax;
  ydot[6] = ddKmax + (dKt * y[8] - Kt * y[6]);
  ydot[7] = (dKmax * y[4] - Ku * y[7]) + (dKt * y[8] - Kt * y[7]);
  ydot[8] = -(dKt * y[8] - Kt * y[7]) - (dKt * y[8] - Kt * y[6]);
  ydot[9] = 0;
  ydot[10] = 0;
  
  List result(1);
  result[0] = ydot;
  return result;
}
