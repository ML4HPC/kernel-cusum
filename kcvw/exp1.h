const char * prefix = "change_in_mean";
const float mu1 = 0;
const float sigma1 = 3; // <--- The std will be 3 and the variance will be 9
const float p_change = 0.02;
const float mu2 = 1;
const float sigma2 = 3; // <-- The std will be 3 and the variance will be 9
const int nref =  5; // <-- 5 is also used by M stat.

const int famax = 1000;

const float h_delta = 0.1;
const int w1_h0 = -2;
const int w5_h0 = -2 ;
const int w10_h0 = -2 ;
const int w15_h0 = -2 ;
const int w20_h0 = -2 ;
const int w40_h0 = -2;

const double cu_delta = 0.1;

double k(double x, double y, double alpha=1, double sigma2=2){
  return exp(-pow(fabs(x-y),alpha)/sigma2);
}
