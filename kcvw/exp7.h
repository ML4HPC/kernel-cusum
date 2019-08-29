const char * prefix = "change_in_std_gaussian";
const float mu1 = 0;
const float sigma1 = 1;  // <--- The std will be 1 and the variance will be 1.
const float p_change = 0.02;
const float mu2 = 0;
const float sigma2 = 1.5; // <-- The std will be 1.5 and the variance will be 1.5^2

const int famax = 1000;

const float h_delta = 0.4;
const int w1_h0 = -2;
const int w5_h0 = -2;
const int w10_h0 = -2;
const int w15_h0 = -2 ;
const int w20_h0 = -2;
const int w40_h0 = -2;

const double cu_delta = 1./8.;

double k(double x, double y, double alpha=2, double sigma2=2){
  return exp(-pow(fabs(x-y),alpha)/sigma2);
}
